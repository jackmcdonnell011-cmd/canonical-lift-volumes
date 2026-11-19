"""The three figures of Section 3 from the data file of the release.

usage:  python make_figures.py modular_geodesics_L25.csv.gz <out_dir>
writes data_vol_vs_wordlength.png, data_vol_vs_length.png, data_twist_families.png
"""
import os
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
df = pd.read_csv(csv, usecols=["word", "word_length", "hyperbolic_length", "volume"])
df = df.rename(columns={"word_length": "L", "hyperbolic_length": "geodesic_length"})

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})

# 1. volume vs word length
rng = np.random.default_rng(0)
g = df.groupby("L").volume.agg(["min", "max"])
fig, ax = plt.subplots(figsize=(6, 3.6))
ax.scatter(df.L + rng.uniform(-0.3, 0.3, len(df)), df.volume, s=0.3, color=BLUE, alpha=0.05,
           linewidths=0, rasterized=True)
ax.plot(g.index, g["max"], color=ORANGE, lw=2, marker="o", ms=4, label="largest volume")
ax.plot(g.index, g["min"], color=AQUA, lw=2, marker="o", ms=4, label="smallest volume")
ax.set_xlabel(r"word length $|w|$")
ax.set_ylabel("volume")
ax.set_xticks(range(2, 26, 2))
ax.legend(loc="upper left")
fig.savefig(f"{outdir}/data_vol_vs_wordlength.png")
plt.close(fig)

# 2. volume vs hyperbolic length (density)
fig, ax = plt.subplots(figsize=(6, 3.6))
hb = ax.hexbin(df.geodesic_length, df.volume, gridsize=160, bins="log", cmap="Blues", mincnt=1,
               linewidths=0)
cb = fig.colorbar(hb, ax=ax, pad=0.02)
cb.set_label("number of knots", color=INK)
cb.outline.set_edgecolor(MUTED)
ax.set_xlabel(r"hyperbolic length $\ell(\gamma_w)$")
ax.set_ylabel("volume")
fig.savefig(f"{outdir}/data_vol_vs_length.png")
plt.close(fig)

# 3. twist families in the low-volume region
VMAX, LMAX = 11.5, 10
sub = df[(df.volume < VMAX) & (df.geodesic_length < LMAX)]
idx = df.set_index("word")
def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))
fams = [
    (r"$X^kY$", [canon("X" * k + "Y") for k in range(1, 25)], BLUE),
    (r"$X^kY^2$", [canon("X" * k + "YY") for k in range(1, 24)], ORANGE),
    (r"$X^kY^3$", [canon("X" * k + "YYY") for k in range(1, 23)], AQUA),
]
fig, ax = plt.subplots(figsize=(6, 3.6))
ax.scatter(sub.geodesic_length, sub.volume, s=4, color="#b9b8b2", linewidths=0, rasterized=True,
           label="all knots in region")
for name, words, col in fams:
    pts = idx.loc[[w for w in words if w in idx.index]]
    pts = pts[(pts.volume < VMAX) & (pts.geodesic_length < LMAX)]
    ax.plot(pts.geodesic_length, pts.volume, color=col, lw=2, marker="o", ms=4, label=name,
            markeredgecolor="white", markeredgewidth=0.8)
ax.set_xlabel(r"hyperbolic length $\ell(\gamma_w)$")
ax.set_ylabel("volume")
ax.set_xlim(right=LMAX); ax.set_ylim(top=VMAX)
ax.legend(loc="lower right")
fig.savefig(f"{outdir}/data_twist_families.png")
plt.close(fig)
print(len(sub), "knots in zoom region")
