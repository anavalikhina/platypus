import numpy as np
import pandas as pd
import plotly.graph_objects as go


def plot_volcano(
    df: pd.DataFrame,
    logfc_thr: float = 1.5,
    pvalue_thr: float = 0.05,
    logfc: str = 'logFC',
    pvalue: str = 'adj.P.Val',
    name: str = 'ID',
    title: str = "",
    filename: str = "volcano",
    file_format: str = "html",
):
    """
    Makes Volcano plot for DE analysis results and saves it under filename.

    Args:
        df: DataFrame with DE analysis results.
        logfc_thr: Threshold for logFC value used as decision boundary on plot. Default to 1.5.
        pvalue_thr: Threshold for p-value used as decision boundary on plot. Default to 0.05.
        logfc: Name of column containing logFC values. Default to 'logFC'.
        pvalue: Name of column containing p-values. Default to 'adj.P.Val'.
        name: Name of column containing gene names. Default to 'ID'.
        title: Title of the figure. Default to ''.
        filename: Name of file to store figure.
        file_format: Format of the file in which plot will be saved. Either 'json' or 'html'. Default to 'html'.
    """

    # determine DE genes (1 - upregulated, -1 - downregulated, 0 - no changes) with pvalue_thr and logfc_thr
    df['DE'] = df.apply(
            lambda row: 1
            if row[pvalue] < pvalue_thr and row[logfc] >= logfc_thr
            else -1
            if row[pvalue] < pvalue_thr and row[logfc] <= -logfc_thr
            else 0,
            axis=1,
        )

    # calculate -log10(p-val)
    df["logp"] = -(np.log10(df[pvalue]))

    # plot volcano
    trace_up = go.Scatter(
        x=df.loc[df['DE'] == 1, logfc],
        y=df.loc[df['DE'] == 1, "logp"],
        text=df.loc[df['DE'] == 1, name],
        mode="markers",
        hoverinfo="text+x+y",
        marker=dict(color="red"),
        textposition="top right",
        textfont=dict(size=7),
    )
    trace_down = go.Scatter(
        x=df.loc[df['DE'] == -1, logfc],
        y=df.loc[df['DE'] == -1, "logp"],
        text=df.loc[df['DE'] == -1, name],
        mode="markers",
        hoverinfo="text+x+y",
        marker=dict(color="blue"),
        textposition="top right",
        textfont=dict(size=7),
    )
    trace_neut = go.Scatter(
        x=df.loc[df['DE'] == 0, logfc],
        y=df.loc[df['DE'] == 0, "logp"],
        text=df.loc[df['DE'] == 0, name],
        hoverinfo="text+x+y",
        mode="markers",
        marker=dict(color="grey"),
    )

    fig = go.Figure(
        data=[trace_up, trace_down, trace_neut],
        layout=go.Layout(
            title=f"<br>{title}",
            xaxis_title="log\u2082(Fold Change)",
            yaxis_title="-log\u2081\u2080(P-value)",
            font=dict(size=16),
            showlegend=False,
        ),
    )
    fig.add_vline(x=logfc_thr, line_width=1, line_dash="dash", line_color="#7d7d7d")
    fig.add_vline(x=-logfc_thr, line_width=1, line_dash="dash", line_color="#7d7d7d")
    fig.add_hline(
        y=-np.log10(pvalue_thr), line_width=1, line_dash="dash", line_color="#7d7d7d"
    )

    if file_format == "json":
        fig.write_json(f"{filename}.json")
    else:
        fig.write_html(f"{filename}.html")
