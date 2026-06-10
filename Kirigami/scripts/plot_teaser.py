#!/usr/bin/env python3
"""Render the poster teaser panels (closed / half-open / 0.9*Theta_max) directly
from the kiri_export SVGs for hero2 (export/hero2/*.svg, "cut" layer = one
closed straight-line polygon path per face).

Each state is rendered into its own square canvas, filled faces in mid-blue,
dark thin edges, gaps between faces left white (the paper background), no
axes, tight per-state framing (equal aspect, centered square crop so the
shape fills as much of the panel as its own aspect ratio allows).

Run under the arm64 interpreter:
    arch -arm64 /usr/local/bin/python3 code/scripts/plot_teaser.py
"""
import math
import re
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

NS = "{http://www.w3.org/2000/svg}"
NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")

FACE_COLOR = "#5F82A8"
EDGE_COLOR = "#1a1a1a"
EDGE_LW = 0.6
SIZE_PX = 1400
DPI = 300
MARGIN_FRAC = 0.04  # padding around the shape's own bounding box

# Each "cut_face" path traces the face's macro triangle plus small in-and-out
# detours where a living-hinge neck notch is cut into the boundary (a few mm
# of perpendicular deviation). Left alone, those detours make the path
# self-intersecting, which under the nonzero fill rule hollows out most of
# the face instead of just nicking the corner. Simplify with a small
# Ramer-Douglas-Peucker tolerance (in mm, same units as the SVG) to drop the
# neck detail and recover the true macro polygon before filling.
RDP_EPS_MM = 3.0


def _perp_dist(pt, a, b):
    (x, y), (x1, y1), (x2, y2) = pt, a, b
    dx, dy = x2 - x1, y2 - y1
    if dx == 0 and dy == 0:
        return math.hypot(x - x1, y - y1)
    t = ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy)
    px, py = x1 + t * dx, y1 + t * dy
    return math.hypot(x - px, y - py)


def rdp(points, eps):
    if len(points) < 3:
        return points
    dmax, idx = 0.0, 0
    for i in range(1, len(points) - 1):
        d = _perp_dist(points[i], points[0], points[-1])
        if d > dmax:
            dmax, idx = d, i
    if dmax > eps:
        left = rdp(points[:idx + 1], eps)
        right = rdp(points[idx:], eps)
        return left[:-1] + right
    return [points[0], points[-1]]

STATES = [
    ("closed", "hero2_130_sigma_mc_closed.svg", "teaser_closed.png"),
    ("half", "hero2_130_sigma_mc_open_half.svg", "teaser_half.png"),
    ("open90", "hero2_130_sigma_mc_open_0.9tm.svg", "teaser_open90.png"),
]

SRC_DIR = "export/hero2"
OUT_DIR = "docs/poster/figs"


def face_polygons(svg_path):
    root = ET.parse(svg_path).getroot()
    cut = next(g for g in root.findall(f"{NS}g") if g.get("id") == "cut")
    polys = []
    for el in cut.findall(f"{NS}path"):
        d = el.get("d", "")
        tokens = re.findall(r"[MLZmlz]|" + NUM.pattern, d)
        pts, i = [], 0
        while i < len(tokens):
            t = tokens[i]
            if t in "MLml":
                pts.append((float(tokens[i + 1]), float(tokens[i + 2])))
                i += 3
            else:
                i += 1
        if len(pts) >= 3:
            simplified = rdp(pts + [pts[0]], RDP_EPS_MM)[:-1]
            if len(simplified) >= 3:
                polys.append(simplified)
            else:
                polys.append(pts)
    return polys


def render(name, svg_name, out_name):
    svg_path = f"{SRC_DIR}/{svg_name}"
    polys = face_polygons(svg_path)
    xs = [x for poly in polys for x, y in poly]
    ys = [y for poly in polys for x, y in poly]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    w, h = xmax - xmin, ymax - ymin
    half = max(w, h) * (1.0 + MARGIN_FRAC) / 2.0
    cx, cy = (xmin + xmax) / 2.0, (ymin + ymax) / 2.0

    fig = plt.figure(figsize=(SIZE_PX / DPI, SIZE_PX / DPI), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("white")
    for poly in polys:
        ax.add_patch(Polygon(poly, closed=True, facecolor=FACE_COLOR,
                              edgecolor=EDGE_COLOR, linewidth=EDGE_LW))
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy + half, cy - half)  # SVG y points down
    ax.set_aspect("equal")
    ax.axis("off")
    out_path = f"{OUT_DIR}/{out_name}"
    fig.savefig(out_path, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"{out_path}: {len(polys)} faces, bbox {w:.1f} x {h:.1f} mm")


if __name__ == "__main__":
    for name, svg_name, out_name in STATES:
        render(name, svg_name, out_name)
