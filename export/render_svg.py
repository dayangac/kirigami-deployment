#!/usr/bin/env python3
"""Render a kiri_export SVG to a PNG preview.

Only the subset of SVG that kiri_export writes is understood: absolute
"M x y L x y ... Z" paths, circular "A" arcs, text, and circles, grouped
into the layers cut / score / engrave. Run under the arm64 interpreter:

    arch -arm64 /usr/local/bin/python3 export/render_svg.py in.svg -o out.png
"""
import argparse
import math
import re
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle

NS = "{http://www.w3.org/2000/svg}"
NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")

LAYER_STYLE = {
    "cut": dict(edgecolor="#111111", facecolor="#cfd8e3", lw=0.6, zorder=2),
    "score": dict(edgecolor="#1f6feb", facecolor="none", lw=0.9, zorder=3),
    "engrave": dict(edgecolor="#d1242f", facecolor="none", lw=0.5, zorder=4),
}


def _arc_points(p0, rx, ry, large, sweep, p1, n=48):
    """Sample a circular SVG arc (rx == ry) using the centre parameterization."""
    (x0, y0), (x1, y1) = p0, p1
    if abs(x1 - x0) < 1e-12 and abs(y1 - y0) < 1e-12:
        return [p1]
    r = max(rx, math.hypot(x1 - x0, y1 - y0) / 2.0)
    mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    dx, dy = (x1 - x0) / 2.0, (y1 - y0) / 2.0
    h2 = max(0.0, r * r / (dx * dx + dy * dy) - 1.0)
    sign = 1.0 if large != sweep else -1.0
    cx = mx + sign * math.sqrt(h2) * dy
    cy = my - sign * math.sqrt(h2) * dx
    a0 = math.atan2(y0 - cy, x0 - cx)
    a1 = math.atan2(y1 - cy, x1 - cx)
    da = a1 - a0
    if sweep and da < 0:
        da += 2 * math.pi
    if not sweep and da > 0:
        da -= 2 * math.pi
    return [(cx + r * math.cos(a0 + da * k / n), cy + r * math.sin(a0 + da * k / n))
            for k in range(1, n + 1)]


def path_points(d):
    """Vertices of one path; circular arcs are sampled."""
    pts, i = [], 0
    tokens = re.findall(r"[MLAZmlaz]|" + NUM.pattern, d)
    while i < len(tokens):
        t = tokens[i]
        if t in "MLml":
            pts.append((float(tokens[i + 1]), float(tokens[i + 2])))
            i += 3
        elif t in "Aa":  # rx ry rot large sweep x y
            args = [float(v) for v in tokens[i + 1:i + 8]]
            start = pts[-1] if pts else (args[5], args[6])
            pts.extend(_arc_points(start, args[0], args[1], int(args[3]), int(args[4]),
                                   (args[5], args[6])))
            i += 8
        elif t in "Zz":
            i += 1
        else:
            i += 1
    return pts


def render(svg_path, png_path, dpi=200):
    root = ET.parse(svg_path).getroot()
    vb = [float(v) for v in root.get("viewBox").split()]
    fig_w = max(3.0, min(14.0, vb[2] / 25.4))
    fig, ax = plt.subplots(figsize=(fig_w, fig_w * vb[3] / vb[2]))
    counts = {}
    for g in root.findall(f"{NS}g"):
        layer = g.get("id", "cut")
        style = LAYER_STYLE.get(layer, LAYER_STYLE["cut"])
        n = 0
        for el in g:
            if el.tag == f"{NS}path":
                pts = path_points(el.get("d", ""))
                if len(pts) < 2:
                    continue
                closed = el.get("d", "").rstrip().endswith("Z")
                if closed and len(pts) >= 3:
                    ax.add_patch(Polygon(pts, closed=True, **style))
                else:
                    xs, ys = zip(*pts)
                    ax.plot(xs, ys, color=style["edgecolor"], lw=style["lw"],
                            zorder=style["zorder"])
                n += 1
            elif el.tag == f"{NS}circle":
                ax.add_patch(Circle((float(el.get("cx")), float(el.get("cy"))),
                                    float(el.get("r")),
                                    edgecolor=style["edgecolor"], facecolor="white",
                                    lw=style["lw"], zorder=style["zorder"] + 1))
                n += 1
            elif el.tag == f"{NS}text" and el.text:
                ax.text(float(el.get("x")), float(el.get("y")), el.text,
                        fontsize=max(2.0, float(el.get("font-size", "3")) * 1.6),
                        color=style["edgecolor"], ha="center",
                        zorder=style["zorder"] + 1)
                n += 1
        counts[layer] = n
    ax.set_xlim(vb[0], vb[0] + vb[2])
    ax.set_ylim(vb[1] + vb[3], vb[1])  # SVG y points down
    ax.set_aspect("equal")
    ax.set_xlabel("mm")
    ax.set_ylabel("mm")
    desc = root.find(f"{NS}desc")
    ax.set_title(f"{svg_path.split('/')[-1]}\n{desc.text if desc is not None else ''}",
                 fontsize=7)
    fig.tight_layout()
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)
    print(f"{png_path}  " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("svg", nargs="+")
    ap.add_argument("-o", "--out", default=None, help="output PNG (single input only)")
    a = ap.parse_args()
    if a.out and len(a.svg) != 1:
        ap.error("--out takes a single input")
    for s in a.svg:
        render(s, a.out or s[:-4] + ".png")
