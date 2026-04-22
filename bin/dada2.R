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
if (!is_paired) merge_mode = "standard"  # gap/separate only meaningful for PE

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
filt_f = paste0("./", prefix, "_1_filt.fastq.gz")
tryCatch(
  {
    if (is_paired){
      filt_r = paste0("./", prefix, "_2_filt.fastq.gz")
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

# ---- Mode-specific blocks: each sets a common set of output variables ----
#
# Common variables produced by both blocks:
#   merged_read_count  — unique ASVs before chimera removal (or combined strand total for separate)
#   seqtab_out         — final nochim sequence table
#   asv_ids            — ASV names for FASTA (seq_N, or seq_f_N / seq_r_N for separate)
#   asv_seqs           — getUniques() subsetted to kept ASVs, for FASTA
#   f_map_out          — per-read ASV index list, for _1_map.txt / _map.txt
#   r_map_out          — per-read ASV index list for _2_map.txt; NULL for SE
#   total_dada2_reads  — sum of reads in seqtab_out
#   proportion_chimeric
#   final_matched_perc — fraction of reads tracked to an ASV; NA for separate mode

if (merge_mode == "separate") {

  # Helper: run the single-end DADA2 post-processing for one strand
  make_strand_outputs = function(dada_obj, drp_obj) {
    tryCatch(
      { seqtab = makeSequenceTable(dada_obj) },
      error = function(msg){
        message(paste("Caught an error at the `makeSequenceTable` stage:\n", msg))
        quit()
      }
    )
    tryCatch(
      { seqtab_nc = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) },
      error = function(msg){
        message(paste("Caught an error at the `removeBimeraDenovo` stage:\n", msg))
        quit()
      }
    )
    tryCatch(
      {
        chimera_ids   = which(colnames(seqtab) %in% colnames(seqtab_nc) == FALSE)
        dup_ids       = which(duplicated(colnames(seqtab)))
        ids_to_remove = unique(c(chimera_ids, dup_ids))
        asv_map       = read_asv_tracking(dada_obj, drp_obj, dada_obj, "single", ids_to_remove)
      }, error = function(msg){
        message(paste("Caught an error at the `read_asv_tracking` stage:\n", msg))
        quit()
      }
    )
    list(seqtab=seqtab, seqtab_nc=seqtab_nc, asv_map=asv_map)
  }

  out_f = make_strand_outputs(dada_f, drp_f)
  out_r = make_strand_outputs(dada_r, drp_r)

  asvs_left_f = sort(as.numeric(unique(unlist(lapply(out_f$asv_map, `[`, 1)))))
  asvs_left_f = asvs_left_f[asvs_left_f > 0]
  asvs_left_r = sort(as.numeric(unique(unlist(lapply(out_r$asv_map, `[`, 1)))))
  asvs_left_r = asvs_left_r[asvs_left_r > 0]

  merged_read_count = length(asvs_left_f) + length(asvs_left_r)
  if (merged_read_count == 0){
    message("Caught an error - No ASVs in either strand - stopping script early.")
    quit()
  }

  seqtab_out = cbind(out_f$seqtab_nc, out_r$seqtab_nc)
  asv_seqs   = c(getUniques(out_f$seqtab)[asvs_left_f], getUniques(out_r$seqtab)[asvs_left_r])
  asv_ids    = c(paste("seq_f", asvs_left_f, sep="_"), paste("seq_r", asvs_left_r, sep="_"))

  f_map_out = lapply(out_f$asv_map, `[`, 1)
  n_f_asvs  = ncol(out_f$seqtab_nc)
  r_map_out = lapply(out_r$asv_map, function(x) { v = x[1]; if (!is.na(v) && v > 0) v + n_f_asvs else v })

  total_dada2_reads   = sum(out_f$seqtab_nc) + sum(out_r$seqtab_nc)
  prop_chim_f         = 1 - sum(out_f$seqtab_nc) / sum(out_f$seqtab)
  prop_chim_r         = 1 - sum(out_r$seqtab_nc) / sum(out_r$seqtab)
  proportion_chimeric = (prop_chim_f + prop_chim_r) / 2
  final_matched_perc  = NA

} else {

  # Standard and gap modes: merge PE pairs (with or without gap-filling), or use dada_f for SE
  tryCatch(
    {
      if (is_paired){
        merged = mergePairs(dada_f, drp_f, dada_r, drp_r, verbose=TRUE,
                            justConcatenate=(merge_mode == "gap"))
      } else{
        merged = dada_f
      }
    }, error = function(msg){
      message(paste("Caught an error at the `mergePairs` stage:\n", msg))
      quit()
    }
  )

  merged_read_count = length(merged$sequence)
  if (merged_read_count == 0){
    message("Caught an error - No ASVs - stopping script early.")
    quit()
  }

  tryCatch(
    { seqtab = makeSequenceTable(merged) },
    error = function(msg){
      message(paste("Caught an error at the `makeSequenceTable` stage:\n", msg))
      quit()
    }
  )

  tryCatch(
    { seqtab_out = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) },
    error = function(msg){
      message(paste("Caught an error at the `removeBimeraDenovo` stage:\n", msg))
      quit()
    }
  )

  chimera_ids   = which(colnames(seqtab) %in% colnames(seqtab_out) == FALSE)
  dup_ids       = which(duplicated(merged$sequence))
  ids_to_remove = unique(c(rbind(chimera_ids, dup_ids)))

  tryCatch(
    {
      tracking_strand = if (is_paired) "forward" else "single"
      final_f_map = read_asv_tracking(dada_f, drp_f, merged, tracking_strand, ids_to_remove)
      if (is_paired){
        final_r_map = read_asv_tracking(dada_r, drp_r, merged, "reverse", ids_to_remove)
      }
    }, error = function(msg){
      message(paste("Caught an error at the `read_asv_tracking` stage:\n", msg))
      quit()
    }
  )

  # Count reads where forward and reverse map to different ASVs (quality check; expected to be near zero)
  unmatched_count = if (is_paired){
    sum(mapply(function(f, r){
      ov = intersect(f, r)
      length(ov) == 0 || ov == 0
    }, final_f_map, final_r_map))
  } else 0

  f_map_out = lapply(final_f_map, `[`, 1)
  r_map_out = if (is_paired) lapply(final_f_map, `[`, 1) else NULL  # both strands carry the same merged ASV index

  asvs_left = sort(as.numeric(unique(unlist(f_map_out))))
  asvs_left = asvs_left[asvs_left > 0]
  asv_seqs  = getUniques(seqtab)[asvs_left]
  asv_ids   = paste("seq", asvs_left, sep="_")

  total_dada2_reads   = sum(seqtab_out)
  proportion_chimeric = 1 - total_dada2_reads / sum(seqtab)
  final_matched_perc  = (length(final_f_map) - unmatched_count) / total_dada2_reads

}

# ---- Shared output section ----

# Save map files (one entry per read → its ASV index)
if (is_paired){
  fwrite(f_map_out, file=paste0("./", prefix, "_1_map.txt"), sep="\n")
  fwrite(r_map_out, file=paste0("./", prefix, "_2_map.txt"), sep="\n")
} else{
  fwrite(f_map_out, file=paste0("./", prefix, "_map.txt"), sep="\n")
}

# Save ASV sequences to FASTA
uniquesToFasta(asv_seqs, paste0("./", prefix, "_asvs.fasta"), asv_ids)

# Save ASV count table
write.table(seqtab_out, file=paste0("./", prefix, "_asv_counts.tsv"), sep="\t", row.names=FALSE)

# Save stats report
report_where_to_cut_f = if (final_where_to_cut_f == 0) -1 else final_where_to_cut_f
report_where_to_cut_r = if (is_paired) { if (final_where_to_cut_r == 0) -1 else final_where_to_cut_r } else NA

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
    length(f_map_out),
    final_matched_perc,
    proportion_chimeric,
    total_dada2_reads,
    report_where_to_cut_f,
    report_where_to_cut_r
  )
)
write.table(output_report_df, file=paste0("./", prefix, "_dada2_stats.tsv"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
