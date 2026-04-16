#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Have to use `box` instead of `library` and `source` so that custom scripts can be loaded when executed by Nextflow
box::use(tidyverse[...])
box::use(data.table[...])
box::use(dada2[...])

# Custom function for tracking reads to their DADA2-generated ASVs (bin/read_asv_tracking.R)
box::use(./read_asv_tracking[...])
# Custom function for automatic truncation of reads based on quality scores (bin/trunc_len_automation.R)
box::use(./trunc_len_automation[...])

# Expects four arguments: prefix, forward fastq, reverse fastq (or "NA" for SE), merge_mode
# merge_mode options: "standard" (default), "gap", "separate"
args       = commandArgs(trailingOnly=TRUE)
prefix     = args[1] # Prefix
path_f     = args[2] # Forward fastq
path_r     = args[3] # Reverse fastq, or "NA" for single-end
merge_mode = if (length(args) >= 4) args[4] else "standard"
is_paired  = !is.na(path_r) && path_r != "NA"

# different tax ranks for silva/pr2
silva_tax_vec = c("Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
pr2_tax_vec = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

count_fastq_reads <- function(reads_path){
  decompressed_reads <- gzfile(reads_path, "rt")
  on.exit(close(decompressed_reads))
  total_lines <- 0L
  while(length(chunk <- readLines(decompressed_reads, n = 100000)) > 0){
    total_lines <- total_lines + length(chunk)
  }
  return(as.integer(total_lines/4))
}

# read counts for reads will only be based on the forward strand
# (at this stage we assume numbers should be similar for forward and reverse)
f_reads_count <- count_fastq_reads(path_f)

# Identify truncLen parameter for filterAndTrim function
final_where_to_cut_f = trunc_len_automation(path_f)
if (is_paired){
  final_where_to_cut_r = trunc_len_automation(path_r)
}

# Do some quality filtering
filt_f =  paste0("./", prefix, "_1", "_filt.fastq.gz")
tryCatch(
  {
    if (is_paired){
      filt_r =  paste0("./", prefix, "_2", "_filt.fastq.gz")
      print(paste0("The forward strand truncation point is: ", final_where_to_cut_f))
      print(paste0("The reverse strand truncation point is: ", final_where_to_cut_r))
      out = filterAndTrim(path_f, filt_f, path_r, filt_r, rm.phix=TRUE, maxEE=c(2,5), truncQ=2, truncLen=c(final_where_to_cut_f,final_where_to_cut_r), compress=TRUE, multithread=TRUE)
    } else{
      print(paste0("The forward strand truncation point is: ", final_where_to_cut_f))
      out = filterAndTrim(path_f, filt_f, rm.phix=TRUE, maxEE=2, truncQ=2, truncLen=final_where_to_cut_f, compress=TRUE, multithread=TRUE)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `filterAndTrim` stage:\n", msg))
    quit()
  }
)

f_trimmed_reads_count <- count_fastq_reads(filt_f)

tryCatch(
  {
    # Learn error model
    err_f = learnErrors(filt_f, multithread=TRUE)
    if (is_paired){
      err_r = learnErrors(filt_r, multithread=TRUE)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `learnErrors` stage:\n", msg))
    quit()
  }
)

tryCatch(
  {
    # Dereplicate sequences
    drp_f = derepFastq(filt_f)
    if (is_paired){
      drp_r = derepFastq(filt_r)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `derepFastq` stage:\n", msg))
    quit()
  }
)

drp_reads_count <- sum(drp_f$uniques)

tryCatch(
  {
    # Generate stranded ASVs
    dada_f = dada(drp_f, err=err_f, multithread=TRUE)
    if (is_paired){
      dada_r = dada(drp_r, err=err_r, multithread=TRUE)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `dada` stage:\n", msg))
    quit()
  }
)

# ---- Separate mode: run F and R as independent SE analyses then combine ----
if (merge_mode == "separate" && is_paired) {

  tryCatch(
    {
      # Forward
      seqtab_f        = makeSequenceTable(dada_f)
      seqtab_f.nochim = removeBimeraDenovo(seqtab_f, method="consensus", multithread=TRUE, verbose=TRUE)

      chimera_ids_f     = which(colnames(seqtab_f) %in% colnames(seqtab_f.nochim) == FALSE)
      duplicated_asvs_f = which(duplicated(colnames(seqtab_f)))
      ids_to_remove_f   = unique(c(chimera_ids_f, duplicated_asvs_f))

      final_f_map = read_asv_tracking(dada_f, drp_f, dada_f, "single", ids_to_remove_f)

      # Reverse
      seqtab_r        = makeSequenceTable(dada_r)
      seqtab_r.nochim = removeBimeraDenovo(seqtab_r, method="consensus", multithread=TRUE, verbose=TRUE)

      chimera_ids_r     = which(colnames(seqtab_r) %in% colnames(seqtab_r.nochim) == FALSE)
      duplicated_asvs_r = which(duplicated(colnames(seqtab_r)))
      ids_to_remove_r   = unique(c(chimera_ids_r, duplicated_asvs_r))

      final_r_map = read_asv_tracking(dada_r, drp_r, dada_r, "single", ids_to_remove_r)
    }, error = function(msg){
      message(paste("Caught an error at the `separate` mode processing stage:\n", msg))
      quit()
    }
  )

  # Collect valid ASV indices (excluding 0 = unassigned)
  asvs_left_f = sort(as.numeric(unique(unlist(lapply(final_f_map, `[[`, 1)))))
  asvs_left_f = asvs_left_f[asvs_left_f > 0]

  asvs_left_r = sort(as.numeric(unique(unlist(lapply(final_r_map, `[[`, 1)))))
  asvs_left_r = asvs_left_r[asvs_left_r > 0]

  if (length(asvs_left_f) == 0 && length(asvs_left_r) == 0) {
    message("Caught an error - No ASVs in either strand - stopping script early.")
    quit()
  }

  # Write forward map
  fwrite(lapply(final_f_map, `[[`, 1), file = paste0("./", prefix, "_1_map.txt"), sep="\n")

  # Write reverse map with index offset so IDs don't collide with forward IDs in the combined FASTA
  n_f_asvs = ncol(seqtab_f.nochim)
  final_r_map_offset = lapply(final_r_map, function(x) {
    v = x[[1]]
    if (!is.na(v) && v > 0) v + n_f_asvs else v
  })
  fwrite(final_r_map_offset, file = paste0("./", prefix, "_2_map.txt"), sep="\n")

  # Combined FASTA: forward ASVs labelled seq_f_N, reverse ASVs labelled seq_r_N
  unqs_f       = getUniques(seqtab_f)[asvs_left_f]
  ids_f        = paste("seq_f", asvs_left_f, sep="_")
  unqs_r       = getUniques(seqtab_r)[asvs_left_r]
  ids_r        = paste("seq_r", asvs_left_r, sep="_")
  combined_seqs = c(unqs_f, unqs_r)
  combined_ids  = c(ids_f, ids_r)
  uniquesToFasta(combined_seqs, paste0("./", prefix, "_asvs.fasta"), combined_ids)

  # Combined ASV count table
  combined_seqtab = cbind(seqtab_f.nochim, seqtab_r.nochim)
  write.table(combined_seqtab, file = paste0("./", prefix, "_asv_counts.tsv"), sep="\t", row.names=FALSE)

  # Stats
  total_dada2_reads_f = sum(seqtab_f.nochim)
  total_dada2_reads_r = sum(seqtab_r.nochim)
  prop_chim_f = 1 - (sum(seqtab_f.nochim) / sum(seqtab_f))
  prop_chim_r = 1 - (sum(seqtab_r.nochim) / sum(seqtab_r))
  merged_read_count = length(asvs_left_f) + length(asvs_left_r)

  output_report_df <- data.frame(
    names = c(
      "initial_read_count",
      "filtered_trimmed_read_count",
      "dereplicated_read_count",
      "merged_read_count",
      "with_asv_read_count",
      "proportion_reads_matched",
      "proportion_reads_chimeric",
      "final_read_count",
      "truncation_point_forward",
      "truncation_point_reverse"
    ),
    values = c(
      f_reads_count,
      f_trimmed_reads_count,
      drp_reads_count,
      merged_read_count,
      length(final_f_map),
      NA,
      (prop_chim_f + prop_chim_r) / 2,
      total_dada2_reads_f + total_dada2_reads_r,
      final_where_to_cut_f,
      final_where_to_cut_r
    )
  )
  write.table(output_report_df, file = paste0("./", prefix, "_dada2_stats.tsv"),
              sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  quit()
}

# ---- Standard and gap modes ----
tryCatch(
  {
    if (is_paired){
      if (merge_mode == "gap") {
        # Join non-overlapping read pairs with a NNNNNNNNNN spacer
        merged = mergePairs(dada_f, drp_f, dada_r, drp_r, verbose=TRUE, justConcatenate=TRUE)
      } else {
        # Standard: merge overlapping pairs, discard those that don't merge
        merged = mergePairs(dada_f, drp_f, dada_r, drp_r, verbose=TRUE)
      }
    } else{
      merged = dada_f
    }
  }, error = function(msg){
    message(paste("Caught an error at the `mergePairs` stage:\n", msg))
    quit()
  }
)

merged_read_count <- length(merged$sequence)
if (merged_read_count == 0){
  message("Caught an error - No ASVs - stopping script early.")
  quit()
}

tryCatch(
  {
    # Make ASV count table
    seqtab = makeSequenceTable(merged)
  }, error = function(msg){
    message(paste("Caught an error at the `makeSequenceTable` stage:\n", msg))
    quit()
  }
)

tryCatch(
  {
    # Remove chimeras
    seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  }, error = function(msg){
    message(paste("Caught an error at the `removeBimeraDenovo` stage:\n", msg))
    quit()
  }
)

chimera_ids = which(colnames(seqtab) %in% colnames(seqtab.nochim) == FALSE)
duplicated_asvs = which(duplicated(merged$sequence))
ids_to_remove = unique(c(rbind(chimera_ids, duplicated_asvs)))

tryCatch(
  {
    # Track reads to their ASVs
    if (is_paired){
      final_f_map = read_asv_tracking(dada_f, drp_f, merged, "forward", ids_to_remove)
      final_r_map = read_asv_tracking(dada_r, drp_r, merged, "reverse", ids_to_remove)
    } else{
      final_f_map = read_asv_tracking(dada_f, drp_f, merged, "single", ids_to_remove)
    }
  }, error = function(msg){
    message(paste("Caught an error at the `read_asv_tracking` stage:\n", msg))
    quit()
  }
)

final_f_output = list()
if (is_paired){
  final_r_output = list()
}

unmatched_asvs = character()

for (i in 1:length(final_f_map)){
  f_map_list = final_f_map[[i]]
  if (is_paired){
    r_map_list = final_r_map[[i]]
    overlap = intersect(f_map_list, r_map_list)

    if (length(overlap) == 0){
      unmatched_asvs = append(unmatched_asvs, i)
    } else if (overlap == 0){
      unmatched_asvs = append(unmatched_asvs, i)
    }
  }

  final_f_output[[i]] = f_map_list[1]
  if (is_paired){
    final_r_output[[i]] = f_map_list[1]
  }
}

# The extremely vast majority of forwards+reverse pairs should be assigned the same ASV. This checks it
traced_remainder = length(final_f_map) - length(unmatched_asvs)
total_dada2_reads = sum(seqtab.nochim)
final_matched_perc = traced_remainder / total_dada2_reads

# Save map to file, with each line representing an ASV's read
if (is_paired){
  fwrite(final_f_output, file = paste0("./", prefix, "_1_map.txt"), sep="\n")
  fwrite(final_r_output, file = paste0("./", prefix, "_2_map.txt"), sep="\n")
} else{
  fwrite(final_f_output, file = paste0("./", prefix, "_map.txt"), sep="\n")
}

# Save ASV count table
write.table(seqtab.nochim, file = paste0("./", prefix, "_asv_counts.tsv"), sep = "\t", row.names=FALSE)

# Write proportion of chimeric reads into a file
seqtab_read_count = sum(seqtab)
seqtab.nochim_read_count = sum(seqtab.nochim)
proportion_chimeric = 1 - (seqtab.nochim_read_count / seqtab_read_count)

# Get count of unique ASVs left after all types of filtering
asvs_left = sort(as.numeric(unique(unlist(final_f_output))))
asvs_left = asvs_left[2:length(asvs_left)]

# Save ASV sequences to FASTA file
seqtab.length = length(seqtab)
num_list = as.character(1:seqtab.length)
id_list = paste("seq", asvs_left, sep="_")
unqs = getUniques(seqtab)[asvs_left]
uniquesToFasta(unqs, paste0("./", prefix, "_asvs.fasta"), id_list)

output_report_df <- data.frame(
  names = c(
    "initial_read_count",
    "filtered_trimmed_read_count",
    "dereplicated_read_count",
    "merged_read_count",
    "with_asv_read_count",
    "proportion_reads_matched",
    "proportion_reads_chimeric",
    "final_read_count",
    "truncation_point_forward",
    "truncation_point_reverse"
  ),
  values = c(
    f_reads_count,
    f_trimmed_reads_count,
    drp_reads_count,
    merged_read_count,
    length(final_f_map),
    final_matched_perc,
    proportion_chimeric,
    total_dada2_reads,
    final_where_to_cut_f,
    if (is_paired) final_where_to_cut_r else NA
  )
)
write.table(output_report_df, file = paste0("./", prefix, "_dada2_stats.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
