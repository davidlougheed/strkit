import re
from collections import deque, defaultdict
from logging import Logger
from pathlib import Path
from strkit.call.types import LocusResult
from strkit.utils import cat_strs

__all__ = ["plot_locus"]


def normalize_motif(normalized_motifs_cache: set, motif: str):
    motif_deque = deque(motif)

    if len(motif) == 3 and "C" in motif:
        if motif.count("C") == 2:
            # special case: if len(motif) == 3 and we have two C, rotate until we start with two Cs
            while not (motif_deque[0] == "C" and motif_deque[1] == "C"):
                motif_deque.rotate(1)
            return cat_strs(motif_deque)
        else:
            # special case: if len(motif) == 3 and we have a C, rotate until we start with C
            while motif_deque[0] != "C":
                motif_deque.rotate(1)
            return cat_strs(motif_deque)
    else:
        for _ in range(len(motif)):
            if (c := cat_strs(motif_deque)) in normalized_motifs_cache:
                return c
            motif_deque.rotate(1)
        normalized_motifs_cache.add(motif)
        return motif


def plot_locus(
    locus_data: LocusResult,
    allele: str,
    allele_label: str,
    n_bins: int | None,
    x_label: str,
    y_label: str,
    annotations: list[str] | None,
    out_file: str | Path,
    logger: Logger,
) -> None:
    """
    TODO
    :param locus_data: TODO
    :param allele: TODO
    :param allele_label: TODO
    :param n_bins: TODO
    :param x_label: TODO
    :param y_label: TODO
    :param annotations: TODO
    :param out_file: TODO
    :param logger: TODO
    :return: TODO
    """

    try:
        import seaborn.objects as so
    except ImportError:
        so = None

    if so is None:
        logger.error("Could not import seaborn. Make sure to install 'strkit[plot]'!")
        exit(1)

    import matplotlib as mpl
    import pandas as pd
    import seaborn as sns

    from matplotlib import figure as mpl_figure

    font_size = 10
    font_rc = {
        "font.family": "Arial",
        "font.size": font_size,
        "axes.titlesize": font_size,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.title_fontsize": font_size,
        "legend.fontsize": font_size,
    }
    mpl.rcParams.update(font_rc)

    def add_text(plot_, x: int, y: int, text: str):
        return plot_.add(
            so.Text(halign="right"),
            data=pd.DataFrame.from_records([{"x": x, "y": y, "text": text}]),
            x="x", y="y", text="text"
        )

    def themed(plot_):
        return plot_.theme({
            **sns.axes_style("white"),
            "axes.spines.top": False,
            "axes.spines.right": False,
            **font_rc,
            "patch.linewidth": 0,
        })

    # ------------------------------------------------------------------------------------------------------------------

    if "reads" not in locus_data or "peaks" not in locus_data:
        logger.error("No reads and/or peaks in specified locus. This locus may not have been called.")
        exit(1)

    records = []

    for read in locus_data["reads"].values():
        records.append({
            allele_label: str(read.get("p", -1)),
            "Copy Number": read["cn"],
            # "Sequence Length": # TODO,  in future STRKit, include read sequence length!
        })

    df = pd.DataFrame.from_records(records)

    if df[df[allele_label] == -1].shape[0] == df.shape[0]:
        logger.warning("No reads with peaks in specified locus; ignoring specified allele and plotting everything.")

    # TODO: add k-mer plot

    if re.match(r"\d+", allele):
        allele_idx = int(allele)
        df = df[df[allele_label] == allele_idx]

    ppi = 300

    n_bins = n_bins or 90  # TODO: normal default

    f = mpl_figure.Figure(figsize=(6.5, 8))  # tight_layout=True
    subfigs = f.subfigures(2, 1)

    plot = (
        themed(so.Plot(df, x="Copy Number"))
        .scale(color=sns.color_palette("muted", 2))
        .add(so.Bars(edgewidth=0), so.Hist(bins=n_bins), color=allele_label)
        .on(subfigs[0])
    )

    def bad_annot(a: str):
        logger.error("Bad annotation (format: x y message): %s", a)

    for annot in (annotations or []):
        ad = annot.split(" ", maxsplit=2)
        if len(ad) != 3:
            bad_annot(annot)
            exit(1)

        try:
            ad_x = int(ad[0])
            ad_y = int(ad[1])
        except ValueError:
            bad_annot(annot)
            exit(1)

        ad_msg = ad[2]

        plot = plot.add(so.Line(color="#666666"), data=pd.DataFrame({"y": [0, ad_y]}).assign(x=ad_x), x="x", y="y")
        plot = add_text(plot, ad_x, ad_y, ad_msg)

    plot = plot.label(x=x_label, y=y_label or "# reads")
    ptr = plot.plot()
    ptr._figure.legends[0].set_bbox_to_anchor((0.8, 0.95))

    normalize: bool = True

    if "kmers" in locus_data["peaks"]:
        normalized_motifs_cache = set()
        kmers_dict: dict[tuple[str, str], float] = defaultdict(lambda: 0.0)
        for pi, p in enumerate(locus_data["peaks"]["kmers"]):
            for k, v in p.items():
                vnorm = v / locus_data["peaks"]["n_reads"][pi] / (len(k) if normalize else 1)
                kmers_dict[str(pi), normalize_motif(normalized_motifs_cache, k) if normalize else k] += vnorm
        kmers_df = pd.DataFrame.from_records([
            {allele_label: pi, "k-mer": kmer, "Count": count} for (pi, kmer), count in kmers_dict.items()
        ])
        kmers_df_top5 = kmers_df.nlargest(5, "Count", keep="all")
        kmers_plot = (
            themed(so.Plot(kmers_df_top5, x="k-mer", y="Count"))
            .facet(col=allele_label)
            .add(so.Bar(), so.Dodge())
            .on(subfigs[1])
        )
        kmers_ptr = kmers_plot.plot()
        kmers_ptr._figure.axes[1].xaxis.set_tick_params(rotation=90)
        kmers_ptr._figure.axes[2].xaxis.set_tick_params(rotation=90)

    # plot = add_text(plot, 107, 80, "Main expansion peak (De Luca et al.)")
    #
    # plot = plot.add(so.Line(color="#666666"), data=pd.DataFrame({"y": [0, 95]}).assign(x=134), x="x", y="y")
    # plot = add_text(plot, 134, 95, "First mosaic peak (De Luca et al.)")
    #
    # plot = plot.add(so.Line(color="#666666"), data=pd.DataFrame({"y": [0, 110]}).assign(x=175), x="x", y="y")
    # plot = add_text(plot, 175, 110, "Second mosaic peak (De Luca et al.)")

    # ptr = plot.plot()
    # # noinspection PyProtectedMember
    # ptr._figure.legends[0].set_bbox_to_anchor((0.8, 0.95))
    #
    # if kmers_plot:
    #     kmers_plot.plot()

    f.savefig(out_file, dpi=ppi)

    # ptr.save(out_file, dpi=ppi)

    # TODO: panel two: overall k-mer distribution
    # TODO: panel three?: comparison with other tools?
