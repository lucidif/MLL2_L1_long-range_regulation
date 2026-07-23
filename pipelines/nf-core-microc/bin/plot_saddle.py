#!/usr/bin/env python3
"""
plot_saddle.py - Generate saddle plots from cooltools .saddledump.npz files.

Usage:
    python3 plot_saddle.py \
        --npz   path/to/sample.saddledump.npz \
        --vecs  path/to/sample.cis.vecs.tsv \
        --out   path/to/output.pdf \
        --title "My sample - saddle plot (cis)" \
        [--fasta path/to/genome.fa] \
        [--qlo 0.05] [--qhi 0.95] [--vmin 0.5] [--vmax 2.0] \
        [--cmap coolwarm] [--color steelblue]

n_bins is inferred automatically from the saddledump matrix shape and must
not be specified -- it is determined by --n-bins used in cooltools saddle.

If --fasta is provided, GC content per bin is computed from the fasta and used
as phasing track to orient E1 (positive E1 = high GC = compartment A).
Otherwise the raw E1 values are used as margin track without orientation.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec


def compute_gc_per_bin(vecs_df, fasta_path):
    """
    Compute GC content for each bin defined in vecs_df using pyfaidx.

    Parameters
    ----------
    vecs_df : pd.DataFrame
        DataFrame with columns chrom, start, end (from .cis.vecs.tsv).
    fasta_path : str
        Path to genome fasta (must be indexed or indexable by pyfaidx).

    Returns
    -------
    gc : np.ndarray
        GC content (0-1) per bin, NaN for bins with no valid sequence.
    """
    try:
        from pyfaidx import Fasta
    except ImportError:
        raise ImportError("pyfaidx is required for --fasta. Install with: pip install pyfaidx")

    print(f"Loading fasta: {fasta_path}")
    fa = Fasta(fasta_path, as_raw=True)

    gc_values = []
    for _, row in vecs_df.iterrows():
        chrom, start, end = row["chrom"], int(row["start"]), int(row["end"])
        if chrom not in fa:
            gc_values.append(np.nan)
            continue
        seq = fa[chrom][start:end].upper()
        total = len(seq)
        if total == 0:
            gc_values.append(np.nan)
            continue
        gc = (seq.count("G") + seq.count("C")) / total
        gc_values.append(gc)

    fa.close()
    return np.array(gc_values)


def orient_e1_by_gc(e1, gc):
    """
    Orient E1 so that positive values correlate with high GC (compartment A).
    If the correlation is negative, flip the sign of E1.

    Parameters
    ----------
    e1 : np.ndarray
        Raw E1 values (may contain NaN).
    gc : np.ndarray
        GC content per bin (may contain NaN).

    Returns
    -------
    e1_oriented : np.ndarray
    flipped : bool
        True if the sign was flipped.
    """
    mask = np.isfinite(e1) & np.isfinite(gc)
    corr = np.corrcoef(e1[mask], gc[mask])[0, 1]
    print(f"Pearson correlation E1 vs GC: {corr:.4f}")
    if corr < 0:
        print("Flipping E1 sign to align with GC content (A = positive E1).")
        return -e1, True
    return e1.copy(), False


def compute_groupmean(e1, q_lo, q_hi, n_bins):
    """
    Compute the mean E1 value per quantile bin (used for margin bar plots).

    Parameters
    ----------
    e1 : np.ndarray
        E1 values (NaN will be dropped).
    q_lo, q_hi : float
        Quantile range.
    n_bins : int
        Number of bins (inferred from saddledata matrix shape).

    Returns
    -------
    groupmean : np.ndarray of shape (n_bins,)
    binedges_plot : np.ndarray of shape (n_bins+1,)
    """
    e1_valid = e1[np.isfinite(e1)]
    lo_val = np.quantile(e1_valid, q_lo)
    hi_val = np.quantile(e1_valid, q_hi)
    binedges_plot = np.linspace(q_lo, q_hi, n_bins + 1)

    e1_clipped = e1_valid[(e1_valid >= lo_val) & (e1_valid <= hi_val)]
    bin_idx = np.digitize(e1_clipped, np.quantile(e1_clipped, binedges_plot))
    bin_idx = np.clip(bin_idx, 1, n_bins)
    groupmean = np.array(
        [e1_clipped[bin_idx == i].mean() for i in range(1, n_bins + 1)]
    )
    return groupmean, binedges_plot


def plot_saddle(
    npz_path,
    vecs_path,
    out_path,
    fasta_path=None,
    title=None,
    q_lo=0.05,
    q_hi=0.95,
    vmin=0.5,
    vmax=2.0,
    cmap="coolwarm",
    color="steelblue",
):
    """
    Plot a saddle plot from a cooltools saddledump .npz file.

    Parameters
    ----------
    npz_path : str
        Path to the .saddledump.npz file produced by cooltools saddle.
    vecs_path : str
        Path to the .cis.vecs.tsv file produced by cooltools eigs-cis.
    out_path : str
        Output file path (pdf, png, svg, etc.).
    fasta_path : str, optional
        Path to genome fasta. If provided, GC content is computed per bin
        and used to orient E1 (E1 > 0 = compartment A = high GC).
    title : str, optional
        Plot title.
    q_lo, q_hi : float
        Quantile range used when running cooltools saddle (must match).
    vmin, vmax : float
        Color scale limits (log scale).
    cmap : str
        Matplotlib colormap name.
    color : str
        Face color for margin bar plots.
    """

    # --- load data ---
    d = np.load(npz_path, allow_pickle=True)
    saddledata = d["saddledata"]

    track = pd.read_csv(vecs_path, sep="\t")

    # remove boundary bins: n_bins is inferred from the matrix
    C = saddledata[1:-1, 1:-1]
    n_bins = C.shape[0]
    print(f"n_bins inferred from saddledata: {n_bins}")

    # --- E1 orientation ---
    e1_raw = track["E1"].values.astype(float)

    if fasta_path is not None:
        gc = compute_gc_per_bin(track, fasta_path)
        e1_oriented, flipped = orient_e1_by_gc(e1_raw, gc)
        if flipped:
            C = C[::-1, ::-1]
    else:
        e1_oriented = e1_raw
        print("No fasta provided: E1 used as-is (orientation not guaranteed).")

    # --- compute groupmean for margin bars ---
    groupmean, binedges_plot = compute_groupmean(e1_oriented, q_lo, q_hi, n_bins)

    # --- layout ---
    fig = plt.figure(figsize=(5, 5))
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    X, Y = np.meshgrid(binedges_plot, binedges_plot)
    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)

    # heatmap
    ax_heatmap = plt.subplot(gs[4])
    img = ax_heatmap.pcolormesh(X, Y, C, norm=norm, cmap=cmap, rasterized=True)
    ax_heatmap.yaxis.set_visible(False)
    ax_heatmap.set_xlim(q_lo, q_hi)
    ax_heatmap.set_ylim(q_hi, q_lo)
    ax_heatmap.grid(False)
    ax_heatmap.set_xlabel("E1 quantiles")

    # left margin
    ax_margin_y = plt.subplot(gs[3], sharey=ax_heatmap)
    ax_margin_y.barh(
        binedges_plot[:-1],
        width=groupmean,
        height=1 / n_bins,
        align="edge",
        edgecolor="k",
        facecolor=color,
        linewidth=0.5,
    )
    ax_margin_y.set_xlim(ax_margin_y.get_xlim()[1], ax_margin_y.get_xlim()[0])
    ax_margin_y.set_ylim(q_hi, q_lo)
    ax_margin_y.spines[["top", "bottom", "left"]].set_visible(False)
    ax_margin_y.xaxis.set_visible(False)
    ax_margin_y.set_ylabel("E1 quantiles")

    # top margin
    ax_margin_x = plt.subplot(gs[1], sharex=ax_heatmap)
    ax_margin_x.bar(
        binedges_plot[:-1],
        height=groupmean,
        width=1 / n_bins,
        align="edge",
        edgecolor="k",
        facecolor=color,
        linewidth=0.5,
    )
    ax_margin_x.set_xlim(q_lo, q_hi)
    ax_margin_x.spines[["top", "right", "left"]].set_visible(False)
    ax_margin_x.xaxis.set_visible(False)
    ax_margin_x.yaxis.set_visible(False)
    if title:
        ax_margin_x.set_title(title)

    # colorbar
    ax_cbar = plt.subplot(gs[5])
    plt.colorbar(img, cax=ax_cbar, label="obs/exp contact frequency")

    plt.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot saddle plot from cooltools npz.")
    parser.add_argument("--npz",   required=True, help="Path to .saddledump.npz")
    parser.add_argument("--vecs",  required=True, help="Path to .cis.vecs.tsv")
    parser.add_argument("--out",   required=True, help="Output file path")
    parser.add_argument("--fasta", default=None,  help="Genome fasta for GC-based E1 orientation")
    parser.add_argument("--title", default=None,  help="Plot title")
    parser.add_argument("--qlo",   type=float, default=0.05)
    parser.add_argument("--qhi",   type=float, default=0.95)
    parser.add_argument("--vmin",  type=float, default=0.5)
    parser.add_argument("--vmax",  type=float, default=2.0)
    parser.add_argument("--cmap",  default="coolwarm")
    parser.add_argument("--color", default="steelblue")
    args = parser.parse_args()

    plot_saddle(
        npz_path=args.npz,
        vecs_path=args.vecs,
        out_path=args.out,
        fasta_path=args.fasta,
        title=args.title,
        q_lo=args.qlo,
        q_hi=args.qhi,
        vmin=args.vmin,
        vmax=args.vmax,
        cmap=args.cmap,
        color=args.color,
    )