#!/usr/bin/env python

import argparse
from pathlib import Path

import colorcet as cc
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from celluloid import Camera

SCATTER_KWARGS = {"lw": 0, "s": 5}


def main():
    parser = argparse.ArgumentParser(description="Generate gif from batches")
    parser.add_argument("input", help="Path to the parquet file.", type=Path)
    parser.add_argument("output", help="Path to gif file", type=Path)
    parser.add_argument(
        "--batch",
        help="Name of the column identifying the batch variable",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--width", help="Figure width", required=False, type=float, default=8
    )
    parser.add_argument(
        "--height", help="Figure height", required=False, type=float, default=8
    )
    parser.add_argument(
        "--fps", help="Frames per second", required=False, type=int, default=2
    )

    args = parser.parse_args()

    labels = pd.read_parquet(args.input).astype({"spatialleiden": "category"})

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    camera = Camera(fig)
    for _, sdf in labels.groupby(args.batch, observed=True):
        scatter_labels(sdf, "spatialleiden", ax)
        fig.tight_layout()
        camera.snap()

    animation = camera.animate()
    animation.save(args.output, fps=args.fps)


def remove_tick_and_label(ax):
    ax.set(xticklabels=[], yticklabels=[], xlabel=None, ylabel=None)
    ax.tick_params(left=False, bottom=False)


def scatter_labels(df, name, ax):
    _ = sns.scatterplot(
        df,
        x="x",
        y="y",
        hue=name,
        ax=ax,
        legend=False,
        palette=sns.color_palette(cc.glasbey, n_colors=df[name].cat.categories.size),
        **SCATTER_KWARGS,
    )
    ax.set(aspect=1)
    remove_tick_and_label(ax)


if __name__ == "__main__":
    main()
