import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D


def plot_conservation(csv_path: str, output_path: str = None, title: str = None):

    df = pd.read_csv(csv_path)
    df = df.sort_values("position").reset_index(drop=True)

    # --- Figure setup ---
    fig, ax1 = plt.subplots(figsize=(14, 4.5), dpi=300)
    fig.patch.set_facecolor("white")

    # --- Font sizes ---
    FS_LABEL    = 15   # axis labels
    FS_TICK     = 13   # tick labels
    FS_LEGEND   = 13   # legend
    FS_TITLE    = 17   # suptitle

    # --- Color palette ---
    COLOR_CONS  = "#2E6DA4"   # steel blue for conservation
    COLOR_COV   = "#C0392B"   # muted red for coverage
    FILL_CONS   = "#AEC6E0"
    FILL_COV    = "#F1A89B"

    # --- Left axis: conservation ---
    ax1.set_xlabel("Residue Position", fontsize=FS_LABEL, labelpad=8)
    ax1.set_ylabel("Conservation Score", color=COLOR_CONS, fontsize=FS_LABEL, labelpad=8)
    ax1.tick_params(axis="y", labelcolor=COLOR_CONS, labelsize=FS_TICK)
    ax1.tick_params(axis="x", labelsize=FS_TICK)
    ax1.set_ylim(0, 1.05)
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))

    line1, = ax1.plot(
        df["position"], df["conservation_score"],
        color=COLOR_CONS, linewidth=1.5, alpha=0.95, zorder=3
    )
    ax1.fill_between(
        df["position"], df["conservation_score"],
        color=FILL_CONS, alpha=0.35, zorder=2
    )

    # --- Right axis: coverage ---
    ax2 = ax1.twinx()
    ax2.set_ylim(0, df["n_seqs"].max() * 1.05)
    ax2.set_ylabel("Coverage (n seqs)", color=COLOR_COV, fontsize=FS_LABEL, labelpad=8)
    ax2.tick_params(axis="y", labelcolor=COLOR_COV, labelsize=FS_TICK)
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    line2, = ax2.plot(
        df["position"], df["n_seqs"],
        color=COLOR_COV, linewidth=1.5, alpha=0.8, zorder=3
    )
    ax2.fill_between(
        df["position"], df["n_seqs"],
        color=FILL_COV, alpha=0.15, zorder=1
    )

    # --- Grid and spines ---
    ax1.set_xlim(df["position"].min(), df["position"].max())
    ax1.grid(axis="y", linestyle=":", linewidth=0.6, color="#CCCCCC", zorder=0)
    ax1.grid(axis="x", linestyle=":", linewidth=0.4, color="#DDDDDD", zorder=0)
    ax1.set_axisbelow(True)

    for spine in ["top"]:
        ax1.spines[spine].set_visible(False)
        ax2.spines[spine].set_visible(False)

    ax1.spines["left"].set_color(COLOR_CONS)
    ax2.spines["right"].set_color(COLOR_COV)
    ax1.spines["left"].set_linewidth(1.2)
    ax2.spines["right"].set_linewidth(1.2)
    """
    # --- Legend ---
    legend_elements = [
        Line2D([0], [0], color=COLOR_CONS, linewidth=1.5, label="Conservation Score"),
        Line2D([0], [0], color=COLOR_COV,  linewidth=1.5, linestyle="--", label="Coverage"),
    ]
    ax1.legend(
        handles=legend_elements,
        loc="upper left",
        fontsize=FS_LEGEND,
        framealpha=0.85,
        edgecolor="#CCCCCC",
        fancybox=False,
    )
    """
    # --- Title ---
    plot_title = title or os.path.splitext(os.path.basename(csv_path))[0]
    fig.suptitle(plot_title, fontsize=FS_TITLE, fontweight="bold", y=1.01)

    plt.tight_layout()

    # --- Output ---
    if output_path is None:
        output_path = os.path.splitext(csv_path)[0] + "_plot.png"

    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Plot saved to: {output_path}")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot conservation score and MSA coverage from a conservation CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_conservation.py conservation.csv
  python plot_conservation.py conservation.csv --output my_plot.png --title "EXAMPLE Conservation"
        """
    )
    parser.add_argument("csv", type=str, help="Path to conservation CSV file")
    parser.add_argument("--output", type=str, default=None, help="Output image path (default: <csv>_plot.png)")
    parser.add_argument("--title",  type=str, default=None, help="Plot title (default: CSV filename)")
    args = parser.parse_args()

    plot_conservation(args.csv, args.output, args.title)