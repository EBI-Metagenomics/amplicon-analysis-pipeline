#!/usr/bin/env python3
"""Render a pipeline metro map (nf-core workflow-schematic style) from a TOML spec.

The spec gives the lines, the lanes, the stages and the stations. The script
calculates all the positions: the columns from the label widths, the lanes from
the space that the labels need, and the curves between the lanes.

The style is the style of the MGnify pipeline maps (miassembler, metaViraVerse):

- one coloured line for each input or sequence type, in parallel tracks;
- grey boxes for the stages;
- the tool name (bold) above each station, a short description below;
- a black dot for a mandatory step, an empty station for an optional step;
- file icons, input triangles, output bars, a title and a legend.

By default the map is one long row. A ``[layout]`` table with ``down`` (the first
stages) and ``up`` (the last stages) folds the map into a U, so that it is more
square: the ``down`` stages are a column on the left (the lines go down), the other
stages are the bottom row, and the ``up`` stages are a column on the right (the lines
go up). The lines bend through 90 degrees between the columns and the row. The title
and the legend go in the space above the bottom row.

The script writes two themes: MAP.svg (light) and MAP_dark.svg (dark), and the
same two PNGs with ``--png``.

This script comes from the mi-metro-map skill of EBI-Metagenomics/mi-agentic
(``skills/mi-metro-map/scripts/metromap.py``). The spec format is in
``skills/mi-metro-map/references/spec.md`` of that repository. Change the spec, not
this copy. To change the script, change it in the mi-metro-map skill of
EBI-Metagenomics/mi-agentic, run the tests of that repository, and copy the script again.

The script needs Python 3.11 or later. It uses only the Python standard library.

Usage::

    make_metromap.py SPEC.toml [--out MAP.svg] [--font MavenPro.ttf] [--png MAP.png]
"""

import argparse
import base64
import math
import re
import shutil
import subprocess
import sys
import textwrap
import tomllib
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import Path

FONT = "'Maven Pro', 'Helvetica Neue', Arial, sans-serif"

# Lines and stations
GAP = 14  # distance between parallel tracks
LW = 9  # line width
R = 13  # station radius
DOT = 5.5  # radius of the black dot of a mandatory station

# Labels: font size, line height and colour of each kind of text
SIZE = {"name": 28, "desc": 26, "note": 22}
LH = {"name": 32, "desc": 30, "note": 28}
FILL = {"name": "#000", "desc": "#000", "note": "#555"}
# Other text: the legend and the map title (the pipeline name keeps its size)
NAME_SIZE = [86]  # the pipeline name (a list, so that the style can change it)
TSIZE = {"legend": 34, "item": 28, "lnote": 24, "org": 46, "subtitle": 28}
BASE_SIZES = [dict(d) for d in (SIZE, LH, TSIZE)]  # the sizes for font_delta = 0
NAME_DY = 30  # last tool-name baseline, above the highest track in use at the station
DESC_DY = 44  # first description baseline, below the lowest track in use at the station
FORK_DY = 45  # branch offset of a fork with two branches
FORK_PITCH = 80  # minimum branch distance of a fork with three or more branches

# Horizontal layout
COL_GAP = 24  # minimum space between the labels of two columns
MIN_PITCH = 150  # minimum distance between two columns (room for a curve)
CURVE = 90  # width of the curve from one lane to another
STAGE_PAD, STAGE_GAP = 20, 15
TERMINUS_GAP = 20
MARGIN = 16

# Vertical layout
LANE_GAP = 120  # minimum distance between the last track of a lane and the next lane
CLEAR = 16  # minimum space between the labels and tracks of different lanes
TITLE_LH = 32

# File icons: (width, height) of the map inputs and outputs, and of the station inputs
ICON = (62 * 1.25, 80 * 1.25)
SMALL_ICON = (62, 80)

# Folded layout (a [layout] table with 'down' and 'up' stages)
BEND = 90  # radius of the innermost line in a bend between a column and the bottom row
VGAP = 24  # minimum space between the labels of two stations in a column
MIN_VPITCH = 96  # minimum distance between two stations in a column
ROW_GAP = 40  # minimum space between the column boxes and the bottom row boxes
FILE_GAP = 30  # space between two file icons in a row
LEGEND_GROUP_GAP = 16  # space between the lines and the symbols in a one-column legend
COL_TITLE_GAP = 30  # extra space between the title of a column box and its first station
COL_EXTRA = 120  # extra width of the column boxes, on the side towards the middle of the map


class SpecError(Exception):
    """The spec has errors. The message lists all of them."""


@dataclass
class Station:
    names: list  # tool name lines (stacked above the station)
    branches: list  # for a fork: [(name lines, desc lines)], else []
    desc: list
    note: list
    lane: str
    stage: str
    col: int
    lines: list
    optional: bool
    input: tuple | None = None  # (label, caption lines, note lines)
    rewrap: bool = True  # False: a folded map keeps the line breaks of the description and note
    x: float = 0.0
    top: float = 0.0  # highest track in use at the station (set by the layout)
    bottom: float = 0.0  # lowest track in use at the station (set by the layout)
    icon_y: float = 0.0  # top of the station input icon (set by the layout)

    @property
    def fork(self):
        return bool(self.branches)


@dataclass
class Model:
    title: dict
    style: dict
    legend: dict
    lines: dict  # id -> {"label", "color", "start", "end"}, in legend order
    lanes: dict  # id -> list of tracks (each a list of line ids), lanes from top to bottom
    stages: list  # [(id, title lines)], from left to right
    stations: list
    inputs: list  # [(label, caption lines, lane or None, row)]
    outputs: list
    fold: tuple | None = None  # (down stage ids, up stage ids) of a folded map
    weight: int = field(init=False)

    def __post_init__(self):
        self.weight = self.style["tool_weight"]


def set_font_delta(delta):
    """Make all the text ``delta`` larger, except the pipeline name and the labels on the file icons.

    The line heights grow by the same amount. The delta is relative to the default sizes.
    """
    for d, base in zip((SIZE, LH, TSIZE), BASE_SIZES):
        for k in d:
            d[k] = base[k] + delta


def text_lines(value):
    """Split a label into lines. A label is a string with ``\\n`` or a list."""
    if not value:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    return str(value).split("\n")


def text_width(s, size, weight=400):
    """Estimate the width of a text in Maven Pro (no font metrics needed)."""
    return len(s) * size * (0.53 if weight < 700 else 0.6 if weight < 800 else 0.63)


def lane_lines(tracks):
    return [ln for t in tracks for ln in t]


def parse(spec):
    """Check the spec and make the model.

    :param spec: The spec, as a dict from ``tomllib``.
    :return: A :class:`Model`.
    :raises SpecError: If the spec has errors.
    """
    errors = []
    lines = {}
    for ln in spec.get("lines", []):
        if not ln.get("id") or not ln.get("color"):
            errors.append(f"line {ln}: needs 'id' and 'color'")
            continue
        if ln["id"] in lines:
            errors.append(f"line '{ln['id']}': duplicate id")
        start, end = ln.get("start", "input"), ln.get("end", "output")
        if start not in ("input", "station") or end not in ("output", "station"):
            errors.append(f"line '{ln['id']}': 'start' must be input or station, 'end' must be output or station")
        lines[ln["id"]] = {"label": ln.get("label", ln["id"]), "color": ln["color"], "start": start, "end": end}
    if not lines:
        errors.append("the spec has no [[lines]]")

    lanes = {}
    for lane in spec.get("lanes", []):
        lid = lane.get("id")
        if not lid or lid in lanes:
            errors.append(f"lane {lane}: needs a unique 'id'")
            continue
        tracks = [t if isinstance(t, list) else [t] for t in lane.get("lines", [])]
        flat = lane_lines(tracks)
        unknown = [x for x in flat if x not in lines]
        if unknown or not flat:
            errors.append(f"lane '{lid}': unknown or no lines {unknown}")
        elif len(flat) != len(set(flat)):
            errors.append(f"lane '{lid}': a line is on two tracks")
        lanes[lid] = tracks
    if not lanes:
        errors.append("the spec has no [[lanes]]")

    style = {"numbered_stages": False, "tool_weight": 700, "font_delta": 0, "name_size": 86}
    style.update(spec.get("style", {}))
    stages = []
    for st in spec.get("stages", []):
        if not st.get("id") or st["id"] in [s[0] for s in stages]:
            errors.append(f"stage {st}: needs a unique 'id'")
            continue
        title = text_lines(st.get("title", st["id"]))
        if style["numbered_stages"]:
            title[0] = f"{len(stages) + 1}. {title[0]}"
        stages.append((st["id"], title))
    stage_ids = [s[0] for s in stages]

    layout = spec.get("layout", {})
    down, up = list(layout.get("down", [])), list(layout.get("up", []))
    fold = None
    if down or up:
        n = len(stage_ids)
        if not down or not up:
            errors.append("layout: give both 'down' and 'up'")
        elif down != stage_ids[:len(down)] or up != stage_ids[n - len(up):] or len(down) + len(up) >= n:
            errors.append("layout: 'down' must be the first stages and 'up' the last stages, "
                          "with one or more stages between them")
        else:
            fold = (down, up)

    stations = []
    last_col = {}
    for i, raw in enumerate(spec.get("stations", []), 1):
        branches = [(text_lines(b.get("tool")), text_lines(b.get("desc"))) for b in raw.get("branches", [])]
        names = text_lines(raw.get("tool"))
        label = " / ".join(" ".join(n) for n in [names] + [b[0] for b in branches] if n) or "no tool"
        where = f"station {i} ({label})"
        if bool(names) == bool(branches):
            errors.append(f"{where}: needs 'tool' or 'branches', not both")
        if branches and (len(branches) < 2 or not all(b[0] for b in branches)):
            errors.append(f"{where}: 'branches' needs two or more tables, each with 'tool'")
        lane, stage = raw.get("lane"), raw.get("stage")
        if lane not in lanes:
            errors.append(f"{where}: unknown lane '{lane}'")
            continue
        if stage not in stage_ids:
            errors.append(f"{where}: unknown stage '{stage}'")
            continue
        # Default column: the next column of this lane in this stage.
        col = raw.get("col", last_col.get((stage, lane), 0) + 1)
        if not isinstance(col, int) or col < 1:
            errors.append(f"{where}: 'col' must be a positive integer")
            continue
        last_col[(stage, lane)] = col
        st_lines = raw.get("lines", lane_lines(lanes[lane]))
        outside = [x for x in st_lines if x not in lane_lines(lanes[lane])]
        if outside or not st_lines:
            errors.append(f"{where}: lines {outside or st_lines} are not in lane '{lane}'")
            continue
        icon = raw.get("input")
        if icon is not None and not icon.get("label"):
            errors.append(f"{where}: 'input' needs 'label'")
        icon = icon and (icon.get("label", ""), text_lines(icon.get("caption")), text_lines(icon.get("note")))
        stations.append(Station(names=names, branches=branches, desc=text_lines(raw.get("desc")),
                                note=text_lines(raw.get("note")), lane=lane, stage=stage, col=col,
                                lines=list(st_lines), optional=bool(raw.get("optional")), input=icon,
                                rewrap=bool(raw.get("rewrap", True))))

    seen = {}
    for s in stations:
        key = (s.stage, s.col, s.lane)
        if key in seen:
            errors.append(f"two stations share lane '{s.lane}', stage '{s.stage}' and col {s.col}")
        seen[key] = True
    for lid in lines:
        cols = [(stage_ids.index(s.stage), s.col) for s in stations if lid in s.lines]
        if not cols:
            errors.append(f"line '{lid}' has no stations")
        elif len(cols) != len(set(cols)):
            errors.append(f"line '{lid}' has two stations in the same column")
    for sid, _ in stages:
        if not any(s.stage == sid for s in stations):
            errors.append(f"stage '{sid}' has no stations")
    if fold:
        vertical = set(fold[0]) | set(fold[1])
        cols = [(s.stage, s.col) for s in stations if s.stage in vertical]
        if len(cols) != len(set(cols)):
            errors.append("two stations of a 'down' or 'up' stage share a col")
        if any(s.fork for s in stations if s.stage in vertical):
            errors.append("a station of a 'down' or 'up' stage cannot have 'branches'")
        for lid, info in lines.items():
            where = {"down" if s.stage in fold[0] else "up" if s.stage in fold[1] else "mid"
                     for s in stations if lid in s.lines}
            left = "down" in where or info["start"] == "input"
            right = "up" in where or info["end"] == "output"
            if left and right and "mid" not in where:
                errors.append(f"line '{lid}' goes from the left column to the right column: "
                              "it needs a station in a stage between them")

    def files(key):
        out = []
        for f in spec.get(key, []):
            if not f.get("label"):
                errors.append(f"{key} {f}: needs 'label'")
            if f.get("lane") is not None and f["lane"] not in lanes:
                errors.append(f"{key} {f}: unknown lane '{f['lane']}'")
            row = f.get("row", 1)
            if not isinstance(row, int) or row < 1:
                errors.append(f"{key} {f}: 'row' must be a positive integer")
            out.append((f.get("label", ""), text_lines(f.get("caption")), f.get("lane"), row))
        return out

    inputs, outputs = files("inputs"), files("outputs")
    if errors:
        raise SpecError("\n".join(f"- {e}" for e in errors))
    legend = {"mandatory": "Mandatory", "optional": "Optional", "input": "Input", "output": "Output", "note": ""}
    legend.update(spec.get("legend", {}))
    return Model(title=spec.get("title", {}), style=style, legend=legend, lines=lines, lanes=lanes,
                 stages=stages, stations=stations, inputs=inputs, outputs=outputs, fold=fold)


# ---------------------------------------------------------------------------
# Geometry relative to a lane (y = 0 is the first track of the lane)
# ---------------------------------------------------------------------------


def track(m, lane, line):
    return next(i for i, t in enumerate(m.lanes[lane]) if line in t) * GAP


def station_span(m, s):
    """y of the first and the last track of the station lines."""
    ys = [track(m, s.lane, ln) for ln in s.lines]
    return min(ys), max(ys)


def block_height(items):
    """Height from the first to the last baseline of a block of (kind, text) items."""
    return sum(LH[k] for k, _ in items[1:])


def fork_offsets(m, s):
    """y offset of each branch of a station. A normal station has one branch."""
    n = len(s.branches)
    if n < 2:
        return [0.0]
    y0, y1 = station_span(m, s)
    if n == 2:  # the station description and note go between the two branches
        middle = [("desc", d) for d in s.desc] + [("note", t) for t in s.note]
        need = block_height(middle) + (LH[middle[0][0]] if middle else 0) + 2 * R + 24
        pitch = max(2 * FORK_DY, need)
    else:  # the name and description of each branch go above it
        above = max(len(nm) * LH["name"] + len(d) * LH["desc"] for nm, d in s.branches[1:])
        pitch = max(FORK_PITCH, above + 2 * R + 16)
    pitch += y1 - y0
    return [(i - (n - 1) / 2) * pitch for i in range(n)]


def stack(x, base, items, up):
    """Baselines for a block of (kind, text) items.

    :param up: True: the last item is on ``base`` and the block goes up. False: the first
        item is on ``base`` and the block goes down.
    """
    out, y = [], base
    seq = reversed(items) if up else items
    for i, (kind, t) in enumerate(seq):
        if i:
            y += -LH[kind] if up else LH[kind]
        out.append((x, y, t, kind))
    return out


def station_labels(m, s):
    """Labels of a station: (x, baseline y, text, kind). kind is name, desc or note."""
    below = [("desc", d) for d in s.desc] + [("note", t) for t in s.note]
    if not s.fork:
        return (stack(s.x, s.top - NAME_DY, [("name", n) for n in s.names], up=True)
                + stack(s.x, s.bottom + DESC_DY, below, up=False))
    y0, y1 = station_span(m, s)
    offs = fork_offsets(m, s)

    def block(b):
        names, desc = s.branches[b]
        return [("name", n) for n in names] + [("desc", d) for d in desc]

    if len(s.branches) == 2:  # labels outside the branches, station description between them
        out = stack(s.x, offs[0] + y0 - NAME_DY, block(0), up=True)
        out += stack(s.x, offs[1] + y1 + DESC_DY + 2, block(1), up=False)
        if below:
            first = (y0 + y1) / 2 - block_height(below) / 2 + 0.3 * SIZE[below[0][0]]
            out += stack(s.x, first, below, up=False)
        return out
    out = []
    for b, off in enumerate(offs):  # labels above each branch, station description below the last one
        out += stack(s.x, off + y0 - NAME_DY, block(b), up=True)
    return out + stack(s.x, offs[-1] + y1 + DESC_DY, below, up=False)


def label_rect(m, x, y, t, kind):
    w = text_width(t, SIZE[kind], m.weight if kind == "name" else 400)
    return (x - w / 2, x + w / 2, y - 0.8 * SIZE[kind], y + 0.25 * SIZE[kind])


def station_width(m, s):
    return max([2 * R] + [r[1] - r[0] for r in (label_rect(m, *lb) for lb in station_labels(m, s))])


def icon_block(m, s):
    """x range and height of the input icon of a station, with its caption on the right."""
    w, h = SMALL_ICON
    texts = [text_width(c, SIZE["desc"]) for c in s.input[1]] + [text_width(c, SIZE["note"]) for c in s.input[2]]
    return s.x - w / 2 - 8, s.x + w / 2 + 18 + max(texts, default=0), h


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


class Layout:
    """All the positions of a map."""

    def __init__(self, m):
        self.m = m
        # Pass 1 finds the gaps between columns that have curves. Pass 2 makes these gaps wide
        # enough for a curve between the labels. The order of the columns does not change.
        self.place_columns({})
        self.make_routes()
        need = {i: (self.col_w[i] + self.col_w[i + 1]) / 2 + CURVE + 2 * CLEAR for i in self.curve_gaps}
        self.place_columns(need)
        self.make_routes()
        self.find_tracks_in_use()
        self.place_station_icons()
        self.place_lanes()

    # -- x ------------------------------------------------------------------
    def place_columns(self, need):
        """Set the x of each column and stage box.

        :param need: Minimum distance between column i and column i + 1, for the gaps with curves.
        """
        m = self.m
        widths = {}
        for s in m.stations:
            key = (s.stage, s.col)
            widths[key] = max(widths.get(key, 0), station_width(m, s))

        self.in_w = self.file_column_width(m.inputs)
        self.x_start = MARGIN + (self.in_w + TERMINUS_GAP if m.inputs else 0)
        cursor = self.x_start + TERMINUS_GAP
        col_x, boxes = {}, []
        prev_x = prev_w = None
        gap = 1  # gap 0 is between the input terminus and the first column
        for sid, title in m.stages:
            x0 = cursor
            for i, c in enumerate(sorted(c for (st, c) in widths if st == sid)):
                w = widths[(sid, c)]
                if prev_x is None:
                    x = x0 + STAGE_PAD + w / 2
                elif i == 0:
                    x = max(x0 + STAGE_PAD + w / 2, prev_x + MIN_PITCH, prev_x + need.get(gap, 0))
                    gap += 1
                else:
                    x = prev_x + max((prev_w + w) / 2 + COL_GAP, MIN_PITCH, need.get(gap, 0))
                    gap += 1
                col_x[(sid, c)] = x
                prev_x, prev_w = x, w
            title_w = max(text_width(t, SIZE["name"], 700) for t in title) + 48
            x1 = max(prev_x + prev_w / 2 + STAGE_PAD, x0 + title_w)
            boxes.append([x0, x1, title])
            cursor = x1 + STAGE_GAP
        for s in m.stations:
            s.x = col_x[(s.stage, s.col)]
        self.boxes = boxes
        self.x_end = boxes[-1][1] + 12
        keys = sorted(col_x, key=lambda k: col_x[k])
        self.columns = [self.x_start] + [col_x[k] for k in keys] + [self.x_end]
        self.col_w = [0] + [widths[k] for k in keys] + [0]
        self.out_w = self.file_column_width(m.outputs)
        self.width = self.x_end + 11 + (TERMINUS_GAP + self.out_w if m.outputs else 0) + MARGIN

    @staticmethod
    def file_column_width(files):
        if not files:
            return 0
        return max([ICON[0] + 16] + [text_width(c, SIZE["desc"]) for _, caps, *_ in files for c in caps])

    def curve(self, gap):
        """x range of a curve in the gap between column gap and column gap + 1.

        The curve is in the middle of the free space between the labels of the two columns.
        """
        self.curve_gaps.add(gap)
        a = self.columns[gap] + self.col_w[gap] / 2
        b = self.columns[gap + 1] - self.col_w[gap + 1] / 2
        mid, w = (a + b) / 2, min(CURVE, self.columns[gap + 1] - self.columns[gap] - 2 * R - 20)
        return mid - w / 2, mid + w / 2

    # -- routes -------------------------------------------------------------
    def make_routes(self):
        """For each line, one route per fork branch: a list of (x, lane, y in lane)."""
        m = self.m
        self.routes = {}
        self.curve_gaps = set()
        for ln, info in m.lines.items():
            stations = sorted((s for s in m.stations if ln in s.lines), key=lambda s: s.x)
            n_branches = max(len(fork_offsets(m, s)) for s in stations)
            routes = []
            for b in range(n_branches):
                pts = []
                prev = None
                for s in stations:
                    offs = fork_offsets(m, s)
                    y = track(m, s.lane, ln) + offs[min(b, len(offs) - 1)]
                    if prev is None:
                        if info["start"] == "input":
                            pts.append((self.x_start, s.lane, y))
                    elif (prev[1], prev[2]) != (s.lane, y):
                        # After a fork, the curve is just before the next station; else just after the last one.
                        gap = self.columns.index(s.x) - 1 if prev[3].fork else self.columns.index(prev[0])
                        t0, t1 = self.curve(gap)
                        pts += [(t0, prev[1], prev[2]), (t1, s.lane, y)]
                    pts.append((s.x, s.lane, y))
                    prev = (s.x, s.lane, y, s)
                if info["end"] == "output":
                    pts.append((self.x_end, prev[1], prev[2]))
                routes.append(pts)
            self.routes[ln] = routes

    def segments(self, lane):
        """Straight parts of the routes in a lane: (x0, x1, y)."""
        for routes in self.routes.values():
            for pts in routes:
                for (xa, la, ya), (xb, lb, yb) in pairwise(pts):
                    if la == lb == lane and ya == yb:
                        yield xa, xb, ya

    def find_tracks_in_use(self):
        """Set the highest and lowest track in use at each station: the labels go outside them."""
        m = self.m
        for s in self.m.stations:
            tracks = {i * GAP for i in range(len(m.lanes[s.lane]))}
            used = {y for x0, x1, y in self.segments(s.lane) if x0 <= s.x <= x1 and y in tracks}
            used.update(station_span(m, s))
            s.top, s.bottom = min(used), max(used)

    def place_station_icons(self):
        """Put each station input icon above the labels of the stations under it."""
        m = self.m
        for s in m.stations:
            if not s.input:
                continue
            x0, x1, h = icon_block(m, s)
            tops = [r[2] for o in m.stations if o.lane == s.lane
                    for r in (label_rect(m, *lb) for lb in station_labels(m, o)) if r[0] < x1 and r[1] > x0]
            tops.append(s.top - R)
            s.icon_y = min(tops) - 14 - h

    # -- y ------------------------------------------------------------------
    def lane_rects(self, lane):
        """Rectangles (x0, x1, y0, y1) of labels, stations, icons and tracks, relative to the lane."""
        m = self.m
        rects = []
        for s in m.stations:
            if s.lane != lane:
                continue
            rects += [label_rect(m, *lb) for lb in station_labels(m, s)]
            y0, y1 = station_span(m, s)
            for off in fork_offsets(m, s):
                rects.append((s.x - R, s.x + R, y0 + off - R, y1 + off + R))
            if s.input:
                x0, x1, h = icon_block(m, s)
                rects.append((x0, x1, s.icon_y, s.icon_y + h))
        rects += [(x0, x1, y - LW / 2, y + LW / 2) for x0, x1, y in self.segments(lane)]
        return rects

    def place_lanes(self):
        m = self.m
        rects = {lane: self.lane_rects(lane) for lane in m.lanes}
        self.lane_y = {}
        prev = None
        for lane in m.lanes:
            y = 0.0 if prev is None else self.lane_y[prev] + (len(m.lanes[prev]) - 1) * GAP + LANE_GAP
            for done, yd in self.lane_y.items():
                for a in rects[done]:
                    for b in rects[lane]:
                        if a[0] < b[1] + CLEAR and b[0] < a[1] + CLEAR:
                            y = max(y, yd + a[3] - b[2] + CLEAR)
            self.lane_y[lane] = y
            prev = lane

        # Stage boxes: the same top and bottom for all; room for the titles at the top.
        absolute = [(r[0], r[1], r[2] + self.lane_y[ln], r[3] + self.lane_y[ln]) for ln in m.lanes for r in rects[ln]]
        top = math.inf
        for x0, x1, title in self.boxes:
            inside = [r[2] for r in absolute if r[0] < x1 and r[1] > x0]
            top = min(top, min(inside) - (46 + (len(title) - 1) * TITLE_LH + 20))
        bottom = max(r[3] for r in absolute) + 30
        shift = 20 - top
        for ln in self.lane_y:
            self.lane_y[ln] += shift
        self.box_top, self.box_bottom = 20, bottom + shift

    def y(self, lane, rel):
        return self.lane_y[lane] + rel


# ---------------------------------------------------------------------------
# SVG
# ---------------------------------------------------------------------------


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def fmt(v):
    return f"{v:.1f}".rstrip("0").rstrip(".")


def text(x, y, s, size=SIZE["desc"], weight=400, anchor="middle", fill="#000"):
    return (f'<text x="{fmt(x)}" y="{fmt(y)}" font-family="{FONT}" font-size="{fmt(size)}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>')


def path_through(points):
    """Polyline through points; a change of y is a smooth S-curve."""
    d = f"M{fmt(points[0][0])},{fmt(points[0][1])}"
    for (x0, y0), (x1, y1) in pairwise(points):
        if y0 == y1:
            d += f" L{fmt(x1)},{fmt(y1)}"
        else:
            xm = (x0 + x1) / 2
            d += f" C{fmt(xm)},{fmt(y0)} {fmt(xm)},{fmt(y1)} {fmt(x1)},{fmt(y1)}"
    return d


def pill(x, y0, y1, r, **attrs):
    extra = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in attrs.items())
    return (f'<rect x="{fmt(x - r)}" y="{fmt(y0 - r)}" width="{fmt(2 * r)}" height="{fmt(y1 - y0 + 2 * r)}" '
            f'rx="{fmt(r)}" {extra}/>')


def draw_station(lay, s):
    m = lay.m
    y0, y1 = station_span(m, s)
    out = []
    for off in fork_offsets(m, s):
        a, b = lay.y(s.lane, y0 + off), lay.y(s.lane, y1 + off)
        out.append(pill(s.x, a, b, R, fill="#fff", stroke="#000", stroke_width=4))
        if not s.optional:
            out.append(pill(s.x, a, b, DOT, fill="#000"))
    for x, y, t, kind in station_labels(m, s):
        out.append(text(x, lay.y(s.lane, y), t, SIZE[kind], m.weight if kind == "name" else 400, fill=FILL[kind]))
    if s.input:
        label, caption, note = s.input
        out += file_icon(s.x - SMALL_ICON[0] / 2, lay.y(s.lane, s.icon_y), label, caption, note, SMALL_ICON, "right")
    return out


def file_icon(x, y, label, caption, note=(), size=ICON, caption_side="below"):
    """nf-core style document icon; (x, y) is the top-left corner.

    The caption goes below the icon, or on its right or left (for a station input). The
    note lines (grey) follow the caption.
    """
    w, h = size
    sc = w / 62
    fold = 18 * sc
    page = f"M{fmt(x)},{fmt(y)} H{fmt(x + w - fold)} L{fmt(x + w)},{fmt(y + fold)} V{fmt(y + h)} H{fmt(x)} Z"
    corner = f"M{fmt(x + w - fold)},{fmt(y)} V{fmt(y + fold)} H{fmt(x + w)}"
    out = [f'<path d="{page}" fill="#fff" stroke="#000" stroke-width="5" stroke-linejoin="round"/>',
           f'<path d="{corner}" fill="none" stroke="#000" stroke-width="4" stroke-linejoin="round"/>']
    bw = max(w + 16, text_width(label, 18 * sc, 700) + 16)
    out.append(f'<rect x="{fmt(x + w / 2 - bw / 2)}" y="{fmt(y + h * 0.42)}" width="{fmt(bw)}" '
               f'height="{fmt(26 * sc)}" rx="4" fill="#000"/>')
    out.append(text(x + w / 2, y + h * 0.42 + 19 * sc, label, 18 * sc, 700, fill="#fff"))
    items = [("desc", c) for c in caption] + [("note", n) for n in note]
    if caption_side in ("right", "left"):
        first = y + h / 2 - block_height(items) / 2 + 8
        cx, anchor = (x + w + 18, "start") if caption_side == "right" else (x - 18, "end")
        for _, yy, t, kind in stack(0, first, items, up=False):
            out.append(text(cx, yy, t, SIZE[kind], anchor=anchor, fill=FILL[kind]))
    else:
        for _, yy, t, kind in stack(0, y + h + 32, items, up=False):
            out.append(text(x + w / 2, yy, t, SIZE[kind], fill=FILL[kind]))
    return out


def file_column(lay, files, centre_x):
    """Icons in a column. An icon with a lane is level with the lane; the others are stacked.

    :return: A tuple ``(svg parts, bottom y)``.
    """
    gap = 40
    heights = [ICON[1] + 40 + LH["desc"] * len(c) for _, c, *_ in files]
    y = lay.box_top
    if not any(f[2] for f in files):  # centre the stack on the boxes if it fits
        y = max(y, (lay.box_top + lay.box_bottom - sum(heights) - gap * (len(files) - 1)) / 2)
    out = []
    for (label, caption, lane, _), h in zip(files, heights):
        if lane:
            centre = lay.y(lane, (len(lay.m.lanes[lane]) - 1) * GAP / 2)
            y = max(y, centre - h / 2 + 20)
        out += file_icon(centre_x - ICON[0] / 2, y, label, caption)
        y += h + gap
    return out, y - gap


def termini(lay, end):
    """Input triangles or output bars: one for each lane where lines start or end."""
    groups = {}
    for ln, routes in lay.routes.items():
        if lay.m.lines[ln]["end" if end else "start"] == "station":
            continue
        _, lane, y = routes[0][-1 if end else 0]
        groups.setdefault(lane, []).append(lay.y(lane, y))
    out = []
    for ys in groups.values():
        y0, y1 = min(ys) - 12, max(ys) + 12
        if end:
            out.append(f'<rect x="{fmt(lay.x_end)}" y="{fmt(y0)}" width="11" height="{fmt(y1 - y0)}" fill="#000"/>')
        else:
            xs = lay.x_start
            out.append(f'<path d="M{fmt(xs)},{fmt(y0)} L{fmt(xs + 26)},{fmt((y0 + y1) / 2)} '
                       f'L{fmt(xs)},{fmt(y1)} Z" fill="#000"/>')
    return out


def legend_symbols(m):
    symbols = [("mandatory", m.legend["mandatory"])]
    if any(s.optional for s in m.stations):
        symbols.append(("optional", m.legend["optional"]))
    return symbols + [("in", m.legend["input"]), ("out", m.legend["output"])]


def legend_note(m, folded):
    """Lines of the legend note. In a folded map the note wraps to the width of the legend items."""
    note = m.legend["note"]
    if not folded or not note:
        return [note] if note else []
    items_w = 120 + max(text_width(lab, TSIZE["item"]) for lab in
                        [v["label"] for v in m.lines.values()] + [lab for _, lab in legend_symbols(m)])
    return textwrap.wrap(note, max(20, int(items_w / text_width("x", TSIZE["lnote"]))))


def legend_size(m, folded=False):
    """Widths of the title and the legend box, and height of the legend box.

    :param folded: True: the lines and the symbols are in one column (a taller, narrower box).
    """
    t = m.title
    # In a folded map the width of the pipeline name is measured better (0.51 em per character in the
    # 800 weight of Maven Pro), because the legend box is as wide as the name.
    name_w = (len(t.get("name", "")) * NAME_SIZE[0] * 0.51 if folded
              else text_width(t.get("name", ""), NAME_SIZE[0], 800))
    title_w = max(text_width(t.get("org", ""), TSIZE["org"], 700), name_w,
                  text_width(t.get("subtitle", ""), TSIZE["subtitle"]))
    symbols = legend_symbols(m)
    lines_w = 120 + max(text_width(v["label"], TSIZE["item"]) for v in m.lines.values())
    sym_w = 52 + max(text_width(lab, TSIZE["item"]) for _, lab in symbols)
    note = legend_note(m, folded)
    if folded:  # the box is at least as wide as the pipeline name (it goes below it)
        bw = max([30 + max(lines_w, 120 + sym_w - 52) + 30, name_w]
                 + [text_width(n, TSIZE["lnote"]) + 60 for n in note])
        rows = len(m.lines) + len(symbols)
        bh = 100 + (rows - 1) * 42 + LEGEND_GROUP_GAP + 40 + LH["note"] * len(note) + (12 if note else 0)
    else:
        bw = max([30 + lines_w + 40 + sym_w + 30] + [text_width(n, TSIZE["lnote"]) + 60 for n in note])
        rows = max(len(m.lines), len(symbols))
        bh = 100 + (rows - 1) * 42 + 40 + (40 if note else 0)
    return title_w, bw, bh, lines_w


def legend(m, left, right, top, folded=False):
    """Title (on the left) and legend box (on the right, or below the title if there is no room).

    :param folded: True: the input and output symbols point down and up.
    :return: A tuple ``(svg parts, bottom y)``.
    """
    out = []
    t = m.title
    title_w, bw, bh, lines_w = legend_size(m, folded)
    if t.get("org"):
        out.append(text(left, top + 70, t["org"], TSIZE["org"], 700, "start", "#555"))
    if t.get("name"):
        out.append(text(left, top + 140, t["name"], NAME_SIZE[0], 800, "start"))
    if t.get("subtitle"):
        out.append(text(left, top + 190, t["subtitle"], TSIZE["subtitle"], 400, "start", "#444"))

    symbols = legend_symbols(m)
    note = legend_note(m, folded)
    bx, by = right - bw, top
    if bx < left + title_w + 60:  # no room on the right of the title: put the legend below it
        bx, by = left, top + 230
    out.append(f'<rect x="{fmt(bx)}" y="{fmt(by)}" width="{fmt(bw)}" height="{fmt(bh)}" rx="24" fill="#E9E9E9"/>')
    out.append(text(bx + 30, by + 50, "Legend", TSIZE["legend"], 700, "start"))
    for i, v in enumerate(m.lines.values()):
        y = by + 100 + i * 42
        out.append(f'<line x1="{fmt(bx + 30)}" y1="{fmt(y)}" x2="{fmt(bx + 100)}" y2="{fmt(y)}" '
                   f'stroke="{v["color"]}" stroke-width="{LW + 3}"/>')
        out.append(text(bx + 120, y + 8, v["label"], TSIZE["item"], 400, "start"))
    sx = bx + 30 + lines_w + 40
    for i, (kind, lab) in enumerate(symbols):
        y, cx = by + 100 + i * 42, sx + 20
        if folded:  # one column: the symbols below the lines, centred on the line samples
            y = by + 100 + (len(m.lines) + i) * 42 + LEGEND_GROUP_GAP
            cx, sx = bx + 65, bx + 120 - 52
        if kind in ("mandatory", "optional"):
            out.append(f'<circle cx="{fmt(cx)}" cy="{fmt(y)}" r="{R}" fill="#fff" stroke="#000" stroke-width="4"/>')
            if kind == "mandatory":
                out.append(f'<circle cx="{fmt(cx)}" cy="{fmt(y)}" r="{DOT}" fill="#000"/>')
        elif kind == "in" and folded:
            out.append(f'<path d="M{fmt(cx - 14)},{fmt(y - 10)} L{fmt(cx)},{fmt(y + 12)} '
                       f'L{fmt(cx + 14)},{fmt(y - 10)} Z" fill="#000"/>')
        elif kind == "in":
            out.append(f'<path d="M{fmt(cx - 10)},{fmt(y - 14)} L{fmt(cx + 12)},{fmt(y)} '
                       f'L{fmt(cx - 10)},{fmt(y + 14)} Z" fill="#000"/>')
        elif folded:
            out.append(f'<rect x="{fmt(cx - 15)}" y="{fmt(y - 5)}" width="30" height="10" fill="#000"/>')
        else:
            out.append(f'<rect x="{fmt(cx - 5)}" y="{fmt(y - 15)}" width="10" height="30" fill="#000"/>')
        out.append(text(sx + 52, y + 9, lab, TSIZE["item"], 400, "start"))
    for i, n in enumerate(note):
        out.append(text(bx + 30, by + bh - 25 - (len(note) - 1 - i) * LH["note"], n, TSIZE["lnote"], 400, "start", "#555"))
    return out, max(by + bh, top + 200)


def render(m, font_b64=None):
    """Make the SVG.

    :param m: The :class:`Model`.
    :param font_b64: A TTF font as base64, to embed in the SVG, or ``None``.
    :return: A tuple ``(svg text, width, height)``.
    """
    set_font_delta(m.style["font_delta"])
    NAME_SIZE[0] = m.style["name_size"]
    if m.fold:
        return render_folded(m, font_b64)
    lay = Layout(m)
    boxes, lines, stations = draw_row(lay)
    body = boxes + lines + termini(lay, end=False) + termini(lay, end=True) + stations
    # The legend is below the stage boxes. The file icons are on the left and right of the boxes,
    # so they can go lower than the boxes.
    parts, bottom = legend(m, lay.boxes[0][0], lay.boxes[-1][1], lay.box_bottom + 40)
    body += parts
    for files, x in [(m.inputs, MARGIN + lay.in_w / 2), (m.outputs, lay.x_end + 11 + TERMINUS_GAP + lay.out_w / 2)]:
        if files:
            parts, icons_bottom = file_column(lay, files, x)
            body += parts
            bottom = max(bottom, icons_bottom)

    return svg_document(body, math.ceil(lay.width), math.ceil(bottom + 40), font_b64)


def svg_document(body, w, h, font_b64):
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">']
    if font_b64:
        head.append("<defs><style>@font-face{font-family:'Maven Pro';font-weight:400 900;"
                    f"src:url(data:font/ttf;base64,{font_b64}) format('truetype');}}</style></defs>")
    head.append(f'<rect width="{w}" height="{h}" fill="#fff"/>')
    return "\n".join(head + body + ["</svg>"]), w, h


def draw_row(lay):
    """Stage boxes, lines and stations of a row. Returns three lists of svg parts."""
    m = lay.m
    boxes, lines, stations = [], [], []
    for x0, x1, title in lay.boxes:
        boxes.append(f'<rect x="{fmt(x0)}" y="{fmt(lay.box_top)}" width="{fmt(x1 - x0)}" '
                     f'height="{fmt(lay.box_bottom - lay.box_top)}" rx="24" fill="#E9E9E9"/>')
        boxes += [text(x0 + 24, lay.box_top + 46 + i * TITLE_LH, t, SIZE["name"], 700, "start", "#8A8A8A")
                  for i, t in enumerate(title)]
    for ln, routes in lay.routes.items():
        for pts in routes:
            d = path_through([(x, lay.y(lane, y)) for x, lane, y in pts])
            lines.append(line_path(d, m.lines[ln]["color"]))
    for s in m.stations:
        stations += draw_station(lay, s)
    return boxes, lines, stations


def line_path(d, color):
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{LW}" stroke-linejoin="round"/>'


def group(parts, dx, dy):
    return [f'<g transform="translate({fmt(dx)},{fmt(dy)})">'] + parts + ["</g>"]


# ---------------------------------------------------------------------------
# Folded layout: a column down on the left, the bottom row, a column up on the right
# ---------------------------------------------------------------------------


def hpill(x0, x1, y, r, **attrs):
    extra = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in attrs.items())
    return (f'<rect x="{fmt(x0 - r)}" y="{fmt(y - r)}" width="{fmt(x1 - x0 + 2 * r)}" height="{fmt(2 * r)}" '
            f'rx="{fmt(r)}" {extra}/>')


class Column:
    """The stations of the 'down' stages (left column) or of the 'up' stages (right column).

    The lines are vertical, at the x given by ``xs``. In the left column the lines go down
    and the labels are on the right of the lines; in the right column the lines go up and
    the labels are on the left. The positions are relative until :meth:`move`.
    """

    def __init__(self, m, stage_ids, xs, down, head_extra=0.0):
        """:param head_extra: Extra space between each box title and the first station of the box."""
        self.m, self.xs, self.down = m, xs, down
        self.lo, self.hi = min(xs.values()), max(xs.values())
        self.anchor = "start" if down else "end"
        self.lx = self.hi + R + 14 if down else self.lo - R - 14
        titles = dict(m.stages)
        self.y, self.stations, self.boxes = {}, [], []
        cursor = 0.0
        for sid in stage_ids:
            title = [" ".join(titles[sid])]  # one line: the box is as wide as the title needs
            head = 46 + (len(title) - 1) * TITLE_LH + 26 + COL_TITLE_GAP + head_extra
            prev = None
            for s in sorted((s for s in m.stations if s.stage == sid), key=lambda s: s.col):
                rects = self.rects(s, 0)
                t, b = min(r[2] for r in rects), max(r[3] for r in rects)
                if down:  # (y, bottom of the station) of the previous station
                    y = cursor + head - t if prev is None else max(prev[0] + MIN_VPITCH, prev[1] + VGAP - t)
                    prev = (y, y + b)
                else:  # (y, top of the station) of the previous station
                    y = cursor - 24 - b if prev is None else min(prev[0] - MIN_VPITCH, prev[1] - VGAP - b)
                    prev = (y, y + t)
                self.y[id(s)] = y
                self.stations.append(s)
            if down:
                self.boxes.append([cursor, prev[1] + 24, title])
                cursor = prev[1] + 24 + STAGE_GAP
            else:
                self.boxes.append([prev[1] - head, cursor, title])
                cursor = prev[1] - head - STAGE_GAP
        rects = [r for s in self.stations for r in self.rects(s, self.y[id(s)])]
        for _, _, title in self.boxes:
            for t in title:
                w = text_width(t, SIZE["name"], 700)
                rects.append((self.lx if down else self.lx - w, self.lx + w if down else self.lx, 0, 0))
        self.x0 = min([r[0] for r in rects] + [self.lo - R]) - STAGE_PAD
        self.x1 = max([r[1] for r in rects] + [self.hi + R]) + STAGE_PAD
        if down:
            self.x1 += COL_EXTRA
        else:
            self.x0 -= COL_EXTRA

    def move(self, dy):
        for k in self.y:
            self.y[k] += dy
        for box in self.boxes:
            box[0] += dy
            box[1] += dy

    @property
    def top(self):
        return min(b[0] for b in self.boxes)

    @property
    def bottom(self):
        return max(b[1] for b in self.boxes)

    def labels(self, s, y):
        """[((x, baseline y, text, kind), rect)]: the tool names next to the station, the description below them."""
        items = [("name", n) for n in s.names] + [("desc", d) for d in s.desc] + [("note", t) for t in s.note]
        first = y - (len(s.names) - 1) * LH["name"] / 2 + 0.35 * SIZE["name"]
        out = []
        for lb in stack(self.lx, first, items, up=False):
            x, yy, t, kind = lb
            w = text_width(t, SIZE[kind], self.m.weight if kind == "name" else 400)
            x0 = x if self.down else x - w
            out.append((lb, (x0, x0 + w, yy - 0.8 * SIZE[kind], yy + 0.25 * SIZE[kind])))
        return out

    def icon(self, s, y):
        """Top-left corner of the station input icon, on the other side of the lines than the labels."""
        w, h = SMALL_ICON
        return (self.lo - R - 18 - w if self.down else self.hi + R + 18), y - h / 2

    def rects(self, s, y):
        xs = [self.xs[ln] for ln in s.lines]
        rects = [r for _, r in self.labels(s, y)] + [(min(xs) - R, max(xs) + R, y - R, y + R)]
        if s.input:
            (ix, iy), (w, h) = self.icon(s, y), SMALL_ICON
            items = [("desc", c) for c in s.input[1]] + [("note", n) for n in s.input[2]]
            cw = max([text_width(t, SIZE[k]) for k, t in items], default=0)
            ch = max(h, block_height(items) + 40)
            cx0 = ix - 18 - cw if self.down else ix + w + 18
            rects += [(ix, ix + w, iy, iy + h), (cx0, cx0 + cw, y - ch / 2, y + ch / 2)]
        return rects

    def draw_boxes(self):
        out = []
        for y0, y1, title in self.boxes:
            out.append(f'<rect x="{fmt(self.x0)}" y="{fmt(y0)}" width="{fmt(self.x1 - self.x0)}" '
                       f'height="{fmt(y1 - y0)}" rx="24" fill="#E9E9E9"/>')
            out += [text(self.lx, y0 + 46 + i * TITLE_LH, t, SIZE["name"], 700, self.anchor, "#8A8A8A")
                    for i, t in enumerate(title)]
        return out

    def draw_stations(self):
        out = []
        for s in self.stations:
            y = self.y[id(s)]
            xs = [self.xs[ln] for ln in s.lines]
            out.append(hpill(min(xs), max(xs), y, R, fill="#fff", stroke="#000", stroke_width=4))
            if not s.optional:
                out.append(hpill(min(xs), max(xs), y, DOT, fill="#000"))
            for (x, yy, t, kind), _ in self.labels(s, y):
                out.append(text(x, yy, t, SIZE[kind], self.m.weight if kind == "name" else 400,
                                self.anchor, FILL[kind]))
            if s.input:
                label, caption, note = s.input
                ix, iy = self.icon(s, y)
                out += file_icon(ix, iy, label, caption, note, SMALL_ICON, "left" if self.down else "right")
        return out


def rewrap(lines, max_w, size):
    """Wrap the words of the lines again, in as few lines as fit ``max_w``. Keeps the lines if that is not fewer."""
    out = []
    for word in " ".join(lines).split():
        if out and text_width(out[-1] + " " + word, size) <= max_w:
            out[-1] += " " + word
        else:
            out.append(word)
    return out if len(out) < len(lines) else lines


def file_width(f):
    return max([ICON[0] + 16] + [text_width(c, SIZE["desc"]) for c in f[1]])


def file_row(files, centre, y, slot=None):
    """File icons side by side, centred on ``centre``.

    :param slot: The width of each icon with its caption (default: its own width).
    :return: A tuple ``(svg parts, (x0, x1))``.
    """
    widths = [slot or file_width(f) for f in files]
    total = sum(widths) + FILE_GAP * (len(files) - 1)
    x = x0 = centre - total / 2
    out = []
    for (label, caption, *_), w in zip(files, widths):
        out += file_icon(x + w / 2 - ICON[0] / 2, y, label, caption)
        x += w + FILE_GAP
    return out, (x0, x0 + total)


def file_rows(files, centre, y0=0.0):
    """File icons in rows (the ``row`` of each file, from the top), each row centred on ``centre``.

    :param y0: Top of the first row.
    :return: A tuple ``(svg parts, (x0, x1), height)``.
    """
    rows = [[f for f in files if f[3] == row] for row in sorted({f[3] for f in files})]
    # All the icons have the same width, so that the icons of the rows are in columns.
    slot = max((file_width(f) for f in files), default=0)
    out, x0, x1, y = [], math.inf, -math.inf, 0
    for i, in_row in enumerate(rows):
        parts, (a2, z2) = file_row(in_row, centre, y0 + y, slot=slot)
        out += parts
        x0, x1, y = min(x0, a2), max(x1, z2), y + file_row_height(in_row) + (FILE_GAP if i < len(rows) - 1 else 0)
    return out, (x0, x1), y


def file_row_height(files):
    return max((ICON[1] + 32 + LH["desc"] * (len(c) - 1) + 10 for _, c, *_ in files), default=0)


def render_folded(m, font_b64=None):
    """Make the SVG of a folded map (see :func:`render`)."""
    down, up = m.fold
    vertical = set(down) | set(up)
    mid = [s for s in m.stations if s.stage not in vertical]
    left_st = [s for s in m.stations if s.stage in down]
    right_st = [s for s in m.stations if s.stage in up]

    def on(ln, stations):
        return [s for s in stations if ln in s.lines]

    in_left = [ln for ln, i in m.lines.items() if on(ln, left_st) or i["start"] == "input"]
    in_right = [ln for ln, i in m.lines.items() if on(ln, right_st) or i["end"] == "output"]

    # The bottom row: the usual layout. The lines from the left column start at its left end,
    # the lines to the right column end at its right end.
    row_lines = {ln: dict(i, start="input" if ln in in_left else "station",
                          end="output" if ln in in_right else "station")
                 for ln, i in m.lines.items() if on(ln, mid)}
    row = Model(title=m.title, style=m.style, legend=m.legend, lines=row_lines,
                lanes={k: v for k, v in m.lanes.items() if any(s.lane == k for s in mid)},
                stages=[st for st in m.stages if st[0] not in vertical], stations=mid, inputs=[], outputs=[])
    # Fewer line breaks in the station descriptions: one line in the columns (the labels are next to
    # the lines, so there is room), and as few lines as fit the column width in the bottom row.
    for s in [s for s in left_st + right_st if s.rewrap]:
        s.desc = [" ".join(s.desc)] if s.desc else []
        s.note = [" ".join(s.note)] if s.note else []
    col_w = {}
    for s in mid:
        col_w[(s.stage, s.col)] = max(col_w.get((s.stage, s.col), 0), station_width(row, s))
    for s in mid:
        if s.rewrap and not s.fork:
            s.desc = rewrap(s.desc, col_w[(s.stage, s.col)], SIZE["desc"])
            s.note = rewrap(s.note, col_w[(s.stage, s.col)], SIZE["note"])
    lay = Layout(row)
    entry = {ln: lay.y(*lay.routes[ln][0][0][1:]) for ln in row_lines if ln in in_left}
    exit_ = {ln: lay.y(*lay.routes[ln][0][-1][1:]) for ln in row_lines if ln in in_right}

    # The order of the tracks in a column keeps the lines parallel in the bends: the line that is
    # highest in the bottom row is on the inside of the bend.
    order = list(m.lines)

    def ranked(lines, ends):
        return sorted(lines, key=lambda ln: (0, ends[ln]) if ln in ends else (1, order.index(ln)))

    left_rank, right_rank = ranked(in_left, entry), ranked(in_right, exit_)
    dx = BEND - lay.x_start  # the left column is at x = 0
    x_end = lay.x_end + dx
    left_xs = {ln: -k * GAP for k, ln in enumerate(left_rank)}
    right_xs = {ln: x_end + BEND + k * GAP for k, ln in enumerate(right_rank)}
    left = Column(m, down, left_xs, down=True)
    right = Column(m, up, right_xs, down=False)

    # y: the file icons at the top, then the input and output termini, then the columns.
    xin = [left.xs[ln] for ln in in_left if m.lines[ln]["start"] == "input"]
    xout = [right.xs[ln] for ln in in_right if m.lines[ln]["end"] == "output"]
    # The file icons are centred on the input triangle and on the output bar.
    out_centre = (min(xout or [right.lo]) + max(xout or [right.hi])) / 2
    in_svg, in_x, in_h = file_rows(m.inputs, (min(xin or [0]) + max(xin or [0])) / 2)
    _, out_x, out_h = file_rows(m.outputs, out_centre)
    y_in, y_out = in_h + 14, out_h + 14
    col_top = y_in + 26 + TERMINUS_GAP
    title_w, bw, bh, _ = legend_size(m, folded=True)

    def radius(rank, ln):
        return BEND + rank.index(ln) * GAP

    def place(left, right):
        """b (the top of the bottom row boxes), c (the bottom of the right column is at b + c) and
        (x, width, height) of the title and legend between the columns (None if they go below the map)."""
        # The row is below the left column, and each bend starts below the last station of its line.
        b = left.bottom + ROW_GAP
        for ln, ye in entry.items():
            last = max(left.y[id(s)] for s in on(ln, left_st)) if on(ln, left_st) else col_top
            b = max(b, last + R + 10 + radius(left_rank, ln) - (ye - lay.box_top))
        # The right column is above the bottom row, with room for the output bar and icons above it.
        lowest = max(right.y[id(s)] for s in right_st)
        c = min([-ROW_GAP] + [exit_[ln] - lay.box_top - radius(right_rank, ln) - 10 - R - lowest for ln in exit_])
        b = max(b, y_out + 11 + TERMINUS_GAP - right.top - c)
        # The title and the legend go in the middle of the white space between the columns, above the
        # bottom row, if they are not too wide. The bottom row goes down if they are too tall.
        x0 = left.x1 + 40
        x1 = min(right.x0, out_x[0] if m.outputs else right.x0) - 40
        for w, h in [(title_w + 60 + bw, max(bh, 200)), (max(title_w, bw), 230 + bh)]:
            if w <= x1 - x0:
                b = max(b, h + 2 * ROW_GAP)
                return b, c, ((x0 + x1 - w) / 2, w, h)
        return b, c, None

    # If the bottom row is lower than the left column needs, all the column boxes get taller (the same
    # extra space below each title), so that the left column stays close to the bottom row.
    extra = 0.0
    for _ in range(50):
        left = Column(m, down, left_xs, down=True, head_extra=extra)
        left.move(col_top)
        right = Column(m, up, right_xs, down=False, head_extra=extra)
        b, c, corner = place(left, right)
        free = b - ROW_GAP - left.bottom
        if free <= 1:
            break
        extra += free / len(down)
    right.move(b + c)
    dy = b - lay.box_top
    # The output bar and icons are just above the right column.
    y_out = right.top - TERMINUS_GAP - 11
    out_top = y_out - 14 - out_h
    out_svg, out_x, _ = file_rows(m.outputs, out_centre, out_top)

    # Lines: the left column with its bends, the bottom row, the right column with its bends.
    lines = []
    for ln in in_left:
        info, x = m.lines[ln], left.xs[ln]
        ys = [left.y[id(s)] for s in on(ln, left_st)]
        ya = y_in if info["start"] == "input" else min(ys)
        if ln in entry:
            r, ye = radius(left_rank, ln), entry[ln] + dy
            d = f"M{fmt(x)},{fmt(ya)} V{fmt(ye - r)} A{fmt(r)},{fmt(r)} 0 0 0 {fmt(x + r)},{fmt(ye)}"
        else:
            d = f"M{fmt(x)},{fmt(ya)} V{fmt(max(ys))}"
        lines.append(line_path(d, info["color"]))
    for ln in in_right:
        info, x = m.lines[ln], right.xs[ln]
        ys = [right.y[id(s)] for s in on(ln, right_st)]
        yt = y_out + 11 if info["end"] == "output" else min(ys)
        if ln in exit_:
            r, ye = radius(right_rank, ln), exit_[ln] + dy
            d = f"M{fmt(x_end)},{fmt(ye)} A{fmt(r)},{fmt(r)} 0 0 0 {fmt(x_end + r)},{fmt(ye - r)} V{fmt(yt)}"
        else:
            d = f"M{fmt(x)},{fmt(max(ys))} V{fmt(yt)}"
        lines.append(line_path(d, info["color"]))

    # Termini: input triangles above the left column, output bar above the right column.
    termini_svg = []
    if xin:
        a, z = min(xin) - 12, max(xin) + 12
        termini_svg.append(f'<path d="M{fmt(a)},{fmt(y_in)} L{fmt(z)},{fmt(y_in)} '
                           f'L{fmt((a + z) / 2)},{fmt(y_in + 26)} Z" fill="#000"/>')
    if xout:
        termini_svg.append(f'<rect x="{fmt(min(xout) - 12)}" y="{fmt(y_out)}" '
                           f'width="{fmt(max(xout) - min(xout) + 24)}" height="11" fill="#000"/>')

    row_boxes, row_lines_svg, row_stations = draw_row(lay)
    body = (left.draw_boxes() + group(row_boxes, dx, dy) + right.draw_boxes()
            + lines + group(row_lines_svg, dx, dy) + termini_svg
            + left.draw_stations() + group(row_stations, dx, dy) + right.draw_stations())

    xs, ys = [left.x0, right.x1, lay.boxes[0][0] + dx, lay.boxes[-1][1] + dx], [0, lay.box_bottom + dy]
    for parts, (a, z) in [(in_svg, in_x), (out_svg, out_x)]:
        if parts:
            body += parts
            xs += [a, z]

    # Title and legend: in the middle of the white space (see above), else below the bottom row.
    if corner:
        lx, w, h = corner
        top = (b - h) / 2 - 30  # the title text starts about 30 below ``top``
        parts, bottom = legend(m, lx, lx + w, top, folded=True)
    else:
        lx, w = left.x0, max(right.x1 - left.x0, title_w + 60 + bw)
        parts, bottom = legend(m, lx, lx + w, lay.box_bottom + dy + 40, folded=True)
    body += parts
    xs += [lx, lx + w]
    ys.append(bottom)

    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    w, h = math.ceil(x1 - x0 + 2 * MARGIN), math.ceil(y1 - y0 + MARGIN + 40)
    return svg_document(group(body, MARGIN - x0, MARGIN - y0), w, h, font_b64)


# The dark theme swaps the neutral colours of the light theme. The line colours do not change.
# The roles swap together: for example, a station is a dark circle with a light ring.
DARK = {
    "#fff": "#0D1117",  # background, station and icon fill (the GitHub dark background)
    "#000": "#E6EDF3",  # text, rings, dots, icon outlines, input triangles, output bars
    "#444": "#B1BAC4",  # subtitle
    "#555": "#9DA7B3",  # notes, org
    "#8A8A8A": "#8B949E",  # stage titles
    "#E9E9E9": "#21262D",  # stage boxes, legend box
}


def dark(svg):
    """The dark theme of a light SVG."""
    return re.sub(r'(fill|stroke)="(#[0-9A-Fa-f]+)"', lambda m: f'{m[1]}="{DARK.get(m[2], m[2])}"', svg)


def dark_path(path):
    """``map.svg`` -> ``map_dark.svg``."""
    return path.with_name(f"{path.stem}_dark{path.suffix}")


CHROME = ["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome", "google-chrome", "chromium",
          "chromium-browser", "chrome"]


def render_png(svg, png, w, h):
    """Take a screenshot of the SVG with headless Chrome (it shows the embedded font)."""
    chrome = next((c if Path(c).exists() else shutil.which(c) for c in CHROME
                   if Path(c).exists() or shutil.which(c)), None)
    if not chrome:
        sys.exit("no Chrome or Chromium found: open the SVG in a browser and export the PNG by hand")
    subprocess.run([chrome, "--headless", f"--screenshot={png.resolve()}", f"--window-size={w},{h}",
                    "--hide-scrollbars", svg.resolve().as_uri()], check=True, capture_output=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("spec", type=Path, help="TOML spec of the map")
    ap.add_argument("--out", type=Path, help="SVG file (default: the spec path with .svg). "
                    "The dark theme goes to the same name with _dark.")
    ap.add_argument("--font", type=Path, help="TTF to embed in the SVG (Maven Pro)")
    ap.add_argument("--png", type=Path, help="also write PNGs, light and _dark (needs Chrome or Chromium)")
    a = ap.parse_args(argv)
    try:
        with open(a.spec, "rb") as f:
            model = parse(tomllib.load(f))
    except (SpecError, tomllib.TOMLDecodeError) as e:
        sys.exit(f"{a.spec}: errors in the spec:\n{e}")
    font = base64.b64encode(a.font.read_bytes()).decode() if a.font else None
    svg, w, h = render(model, font)
    out = a.out or a.spec.with_suffix(".svg")
    for svg_path, png_path, content in [(out, a.png, svg), (dark_path(out), a.png and dark_path(a.png), dark(svg))]:
        svg_path.write_text(content)
        print(f"wrote {svg_path} ({w} x {h})")
        if png_path:
            render_png(svg_path, png_path, w, h)
            print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
