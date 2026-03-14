"""
Plotly visualization functions for the MORPHIC Portal.
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .config import (
    CONDITION_COLORS,
    KNOCKDOWN_COLORS,
    LINEAGE_COLORS,
    LINEAGES,
)


def plot_lineage_heatmap(df: pd.DataFrame, gene: str, metric: str = "observed_e_score_std") -> go.Figure:
    """
    Create heatmap of lineage effects for a gene across conditions.

    Args:
        df: Lineage analysis DataFrame filtered to gene
        gene: Gene symbol for title
        metric: Column to use for color (e.g., observed_e_score_std, z_score_e_score_std)

    Returns:
        Plotly figure
    """
    # Pivot to get condition x lineage matrix
    pivot = df.pivot_table(
        index="condition",
        columns="lineage",
        values=metric,
        aggfunc="mean"
    )

    # Reorder columns to match LINEAGES order
    cols = [c for c in LINEAGES if c in pivot.columns]
    pivot = pivot[cols]

    # Create heatmap
    fig = px.imshow(
        pivot,
        color_continuous_scale="RdBu_r",
        color_continuous_midpoint=0,
        labels=dict(x="Lineage", y="Condition", color="Score"),
        aspect="auto",
    )

    fig.update_layout(
        title=f"Lineage Effects: {gene}",
        height=250,
        margin=dict(l=10, r=10, t=40, b=10),
    )

    # Add text annotations
    for i, condition in enumerate(pivot.index):
        for j, lineage in enumerate(pivot.columns):
            val = pivot.loc[condition, lineage]
            if not np.isnan(val):
                fig.add_annotation(
                    x=j, y=i,
                    text=f"{val:.2f}",
                    showarrow=False,
                    font=dict(size=10, color="white" if abs(val) > 2 else "black"),
                )

    return fig


def plot_knockdown_comparison(df: pd.DataFrame, gene: str) -> go.Figure:
    """
    Bar chart comparing knockdown efficiency across conditions.

    Args:
        df: Knockdown efficiency DataFrame filtered to gene
        gene: Gene symbol for title

    Returns:
        Plotly figure
    """
    fig = go.Figure()

    for condition in ["iPSC", "EBs"]:
        cond_data = df[df["condition"] == condition]
        if len(cond_data) == 0:
            continue

        # Take first row if multiple perturbations
        row = cond_data.iloc[0]

        fig.add_trace(go.Bar(
            x=[condition],
            y=[row["knockdown_pct"]],
            name=condition,
            marker_color=CONDITION_COLORS.get(condition, "#888"),
            text=[f"{row['knockdown_pct']:.1f}%"],
            textposition="outside",
        ))

    fig.update_layout(
        title=f"Knockdown Efficiency: {gene}",
        yaxis_title="Knockdown %",
        yaxis_range=[0, 100],
        showlegend=False,
        height=300,
        margin=dict(l=10, r=10, t=40, b=10),
    )

    return fig


def plot_viability_rank(df: pd.DataFrame, gene: str, condition: str = "EBs") -> go.Figure:
    """
    Show gene's position in the LFC distribution.

    Args:
        df: Viability DataFrame
        gene: Gene symbol to highlight
        condition: Which condition to show

    Returns:
        Plotly figure
    """
    cond_df = df[df["condition"] == condition].copy()

    if len(cond_df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No viability data for {condition}", showarrow=False)
        return fig

    # Sort by LFC
    cond_df = cond_df.sort_values("median_lfc")
    cond_df["rank"] = range(1, len(cond_df) + 1)

    # Find the gene
    gene_row = cond_df[cond_df["gene"] == gene]

    fig = go.Figure()

    # All genes as scatter
    fig.add_trace(go.Scatter(
        x=cond_df["rank"],
        y=cond_df["median_lfc"],
        mode="markers",
        marker=dict(size=4, color="#cccccc"),
        name="All genes",
        hovertemplate="%{customdata}<br>LFC: %{y:.2f}<extra></extra>",
        customdata=cond_df["gene"],
    ))

    # Highlight the gene
    if len(gene_row) > 0:
        fig.add_trace(go.Scatter(
            x=gene_row["rank"],
            y=gene_row["median_lfc"],
            mode="markers",
            marker=dict(size=12, color="#e74c3c", symbol="diamond"),
            name=gene,
            hovertemplate=f"{gene}<br>Rank: %{{x}}<br>LFC: %{{y:.2f}}<extra></extra>",
        ))

        rank = gene_row["rank"].values[0]
        lfc = gene_row["median_lfc"].values[0]
        total = len(cond_df)

        fig.add_annotation(
            x=rank, y=lfc,
            text=f"{gene}: Rank {rank}/{total}",
            showarrow=True,
            arrowhead=2,
            yshift=20,
        )

    # Add reference line at 0
    fig.add_hline(y=0, line_dash="dash", line_color="gray")

    fig.update_layout(
        title=f"Viability Rank ({condition})",
        xaxis_title="Rank (depleted → enriched)",
        yaxis_title="Median LFC",
        height=350,
        showlegend=False,
        margin=dict(l=10, r=10, t=40, b=10),
    )

    return fig


def plot_deg_volcano(df: pd.DataFrame, perturbation: str, padj_threshold: float = 0.05, lfc_threshold: float = 1.0) -> go.Figure:
    """
    Volcano plot of DEGs for a perturbation.

    Args:
        df: DEG table DataFrame
        perturbation: Perturbation ID for title
        padj_threshold: Adjusted p-value threshold for significance
        lfc_threshold: Log2 fold change threshold for labeling

    Returns:
        Plotly figure
    """
    df = df.copy()

    # Calculate -log10 padj
    df["neg_log_padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    # Categorize points
    df["category"] = "Not significant"
    df.loc[(df["padj"] < padj_threshold) & (df["log2FoldChange"] > lfc_threshold), "category"] = "Up"
    df.loc[(df["padj"] < padj_threshold) & (df["log2FoldChange"] < -lfc_threshold), "category"] = "Down"

    color_map = {
        "Not significant": "#b0aba5",
        "Up": "#e74c3c",
        "Down": "#4a8eb5",
    }

    fig = px.scatter(
        df,
        x="log2FoldChange",
        y="neg_log_padj",
        color="category",
        color_discrete_map=color_map,
        hover_data=["gene", "baseMean", "padj"],
        labels={
            "log2FoldChange": "Log2 Fold Change",
            "neg_log_padj": "-Log10 Adjusted P-value",
        },
    )

    # Add threshold lines
    fig.add_hline(y=-np.log10(padj_threshold), line_dash="dash", line_color="gray")
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")

    # Label top genes
    top_genes = df.nlargest(5, "neg_log_padj")
    for _, row in top_genes.iterrows():
        if pd.notna(row.get("gene")):
            gene_name = row["gene"]
        else:
            gene_name = row.name if isinstance(row.name, str) else ""

        fig.add_annotation(
            x=row["log2FoldChange"],
            y=row["neg_log_padj"],
            text=gene_name[:10],
            showarrow=True,
            arrowhead=0,
            font=dict(size=9),
        )

    n_up = (df["category"] == "Up").sum()
    n_down = (df["category"] == "Down").sum()

    fig.update_layout(
        title=f"DEGs: {perturbation} (↑{n_up} ↓{n_down})",
        height=400,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        margin=dict(l=10, r=10, t=60, b=10),
    )

    return fig


def plot_timecourse_expression(df: pd.DataFrame, gene: str) -> go.Figure:
    """
    Line plot of gene expression across timepoints.

    Args:
        df: Timecourse expression DataFrame
        gene: Gene symbol to plot

    Returns:
        Plotly figure
    """
    gene_df = df[df["gene"] == gene].copy()

    if len(gene_df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No timecourse data for {gene}", showarrow=False)
        return fig

    # Sort by day with iPSC first (before day 0)
    day_order = ["iPSC", "0", "2", "4", "6", "8", "10", "12", "14"]
    gene_df["day_str"] = gene_df["day"].astype(str)
    gene_df["day_order"] = gene_df["day_str"].apply(
        lambda x: day_order.index(x) if x in day_order else 999
    )
    gene_df = gene_df.sort_values("day_order")

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Mean Expression", "% Cells Expressing"),
        horizontal_spacing=0.12,
    )

    # Use day_str for x-axis so iPSC displays correctly
    x_values = gene_df["day_str"].tolist()

    # Mean expression
    fig.add_trace(
        go.Scatter(
            x=x_values,
            y=gene_df["mean_expression"],
            mode="lines+markers",
            name="Expression",
            line=dict(color="#5a9e8f"),
        ),
        row=1, col=1
    )

    # Percent expressing at different thresholds
    # Check if threshold columns exist (new format) or fall back to old format
    if "pct_expr_gt1" in gene_df.columns:
        # New format with multiple thresholds - show >1 TPM and >5 TPM
        fig.add_trace(
            go.Scatter(
                x=x_values,
                y=gene_df["pct_expr_gt1"],
                mode="lines+markers",
                name=">1 TPM",
                line=dict(color="#2ecc71"),
            ),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=x_values,
                y=gene_df["pct_expr_gt5"],
                mode="lines+markers",
                name=">5 TPM",
                line=dict(color="#e74c3c"),
            ),
            row=1, col=2
        )
        show_legend = True
    else:
        # Old format - just show pct_expressing
        fig.add_trace(
            go.Scatter(
                x=x_values,
                y=gene_df["pct_expressing"],
                mode="lines+markers",
                name="% Expressing",
                line=dict(color="#2ecc71"),
            ),
            row=1, col=2
        )
        show_legend = False

    fig.update_layout(
        title=f"Timecourse Expression: {gene}",
        height=300,
        showlegend=show_legend,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=10, r=10, t=60, b=10),
    )

    # Force category order so iPSC appears first
    fig.update_xaxes(title_text="Day", categoryorder="array", categoryarray=x_values, row=1, col=1)
    fig.update_xaxes(title_text="Day", categoryorder="array", categoryarray=x_values, row=1, col=2)
    fig.update_yaxes(title_text="Mean Expression", row=1, col=1)
    fig.update_yaxes(title_text="% Cells", range=[0, 100], row=1, col=2)

    return fig


def plot_compositional_bars(df: pd.DataFrame, gene: str) -> go.Figure:
    """
    Bar chart showing compositional changes (chi-square residuals) per lineage.

    Args:
        df: Compositional analysis DataFrame filtered to gene
        gene: Gene symbol for title

    Returns:
        Plotly figure
    """
    if len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No compositional data for {gene}", showarrow=False)
        return fig

    # Get residual columns
    residual_cols = [c for c in df.columns if c.startswith("residual_")]
    lineages = [c.replace("residual_", "") for c in residual_cols]

    fig = go.Figure()

    for condition in ["iPSC", "EBs"]:
        cond_df = df[df["condition"] == condition]
        if len(cond_df) == 0:
            continue

        # Take first perturbation
        row = cond_df.iloc[0]
        values = [row[col] for col in residual_cols]

        fig.add_trace(go.Bar(
            x=lineages,
            y=values,
            name=condition,
            marker_color=CONDITION_COLORS.get(condition, "#888"),
        ))

    fig.update_layout(
        title=f"Compositional Changes: {gene}",
        xaxis_title="Lineage",
        yaxis_title="Chi-Square Residual",
        barmode="group",
        height=300,
        margin=dict(l=10, r=10, t=40, b=10),
    )

    # Add reference line at 0
    fig.add_hline(y=0, line_dash="dash", line_color="gray")

    return fig


def plot_marker_barplot(df: pd.DataFrame, perturbation: str) -> go.Figure:
    """
    Bar plot showing marker gene expression for a perturbation.

    Args:
        df: Marker counts DataFrame for the perturbation
        perturbation: Perturbation ID for title

    Returns:
        Plotly figure
    """
    if df is None or len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No marker data for {perturbation}", showarrow=False)
        return fig

    # Define marker gene groups with colors
    marker_groups = {
        "Pluripotency": {"genes": ["POU5F1", "NANOG", "SOX2", "DNMT3B", "DPPA3"], "color": "#5a9e8f"},
        "Ectoderm": {"genes": ["PAX6", "SOX1", "NES", "OTX2", "TFAP2A"], "color": "#9b59b6"},
        "Endoderm": {"genes": ["SOX17", "FOXA2", "GATA4", "HNF4A", "CXCR4"], "color": "#f1c40f"},
        "Mesoderm": {"genes": ["T", "MIXL1", "MESP1", "TBX6", "HAND1"], "color": "#e74c3c"},
        "Trophoblast": {"genes": ["CDX2", "GATA3", "KRT7", "TFAP2C", "TP63"], "color": "#2ecc71"},
    }

    fig = go.Figure()

    # Add bars for each group
    for group_name, group_info in marker_groups.items():
        group_genes = [g for g in group_info["genes"] if g in df["gene"].values]
        if not group_genes:
            continue

        group_df = df[df["gene"].isin(group_genes)].copy()
        group_df = group_df.set_index("gene").loc[group_genes].reset_index()

        fig.add_trace(go.Bar(
            x=group_df["gene"],
            y=group_df["cpm"],
            name=group_name,
            marker_color=group_info["color"],
            hovertemplate="<b>%{x}</b><br>CPM: %{y:.1f}<extra></extra>",
        ))

    fig.update_layout(
        title=f"Marker Expression: {perturbation}",
        xaxis_title="Marker Gene",
        yaxis_title="CPM",
        height=350,
        margin=dict(l=10, r=10, t=40, b=80),
        xaxis=dict(tickangle=45),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        barmode="group",
    )

    return fig


def plot_marker_dotplot(df: pd.DataFrame, perturbation: str, ntc_df: pd.DataFrame = None) -> go.Figure:
    """
    Dotplot showing marker gene expression for a perturbation vs NTC.

    Args:
        df: Marker counts DataFrame for the perturbation
        perturbation: Perturbation ID for title
        ntc_df: Optional NTC marker counts for comparison

    Returns:
        Plotly figure
    """
    if df is None or len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No marker data for {perturbation}", showarrow=False)
        return fig

    # Define marker gene groups
    marker_genes = {
        "Pluripotency": ["POU5F1", "NANOG", "SOX2", "DNMT3B", "DPPA3"],
        "Ectoderm": ["PAX6", "SOX1", "NES", "OTX2", "TFAP2A"],
        "Endoderm": ["SOX17", "FOXA2", "GATA4", "HNF4A", "CXCR4"],
        "Mesoderm": ["T", "MIXL1", "MESP1", "TBX6", "HAND1"],
        "Trophoblast": ["CDX2", "GATA3", "KRT7", "TFAP2C", "TP63"],
    }

    # Flatten marker genes in order
    ordered_genes = []
    gene_to_group = {}
    for group, genes in marker_genes.items():
        for g in genes:
            if g in df["gene"].values:
                ordered_genes.append(g)
                gene_to_group[g] = group

    if not ordered_genes:
        # Just use all genes if no markers found
        ordered_genes = df["gene"].tolist()[:20]

    # Filter to ordered genes
    df_filt = df[df["gene"].isin(ordered_genes)].copy()

    # Set gene order
    df_filt["gene"] = pd.Categorical(df_filt["gene"], categories=ordered_genes, ordered=True)
    df_filt = df_filt.sort_values("gene")

    # Calculate log CPM
    df_filt["log_cpm"] = np.log1p(df_filt["cpm"])

    fig = go.Figure()

    # Add perturbation as dots
    fig.add_trace(go.Scatter(
        x=df_filt["gene"],
        y=[perturbation] * len(df_filt),
        mode="markers",
        marker=dict(
            size=df_filt["log_cpm"] * 3,
            color=df_filt["log_cpm"],
            colorscale="Reds",
            showscale=True,
            colorbar=dict(title="Log CPM"),
        ),
        name=perturbation,
        hovertemplate="<b>%{x}</b><br>CPM: %{customdata:.1f}<extra></extra>",
        customdata=df_filt["cpm"],
    ))

    # Add NTC comparison if available
    if ntc_df is not None and len(ntc_df) > 0:
        ntc_filt = ntc_df[ntc_df["gene"].isin(ordered_genes)].copy()
        ntc_filt["gene"] = pd.Categorical(ntc_filt["gene"], categories=ordered_genes, ordered=True)
        ntc_filt = ntc_filt.sort_values("gene")
        ntc_filt["log_cpm"] = np.log1p(ntc_filt["cpm"])

        fig.add_trace(go.Scatter(
            x=ntc_filt["gene"],
            y=["NTC"] * len(ntc_filt),
            mode="markers",
            marker=dict(
                size=ntc_filt["log_cpm"] * 3,
                color=ntc_filt["log_cpm"],
                colorscale="Blues",
                showscale=False,
            ),
            name="NTC",
            hovertemplate="<b>%{x}</b><br>CPM: %{customdata:.1f}<extra></extra>",
            customdata=ntc_filt["cpm"],
        ))

    # Add group dividers
    x_pos = 0
    for group in marker_genes.keys():
        group_genes = [g for g in ordered_genes if gene_to_group.get(g) == group]
        if group_genes:
            mid_pos = x_pos + len(group_genes) / 2 - 0.5
            fig.add_annotation(
                x=mid_pos,
                y=1.15,
                yref="paper",
                text=group,
                showarrow=False,
                font=dict(size=10, color="gray"),
            )
            if x_pos > 0:
                fig.add_vline(x=x_pos - 0.5, line_dash="dot", line_color="lightgray")
            x_pos += len(group_genes)

    fig.update_layout(
        title=f"Marker Expression: {perturbation}",
        xaxis_title="Marker Gene",
        height=300,
        margin=dict(l=10, r=10, t=60, b=80),
        xaxis=dict(tickangle=45),
    )

    return fig


def plot_umap_highlight(
    umap_df: pd.DataFrame,
    perturbation: str,
    color_by: str = "cell_type"
) -> go.Figure:
    """
    UMAP plot with contour density for background and highlighted perturbation cells.

    Args:
        umap_df: DataFrame with umap_1, umap_2, perturbation, cell_type columns
        perturbation: Perturbation to highlight
        color_by: Column to color background cells by

    Returns:
        Plotly figure
    """
    if umap_df is None or len(umap_df) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="UMAP data not available. Run: python scripts/extract_umap_coordinates.py",
            showarrow=False
        )
        return fig

    # Find perturbation column
    pert_col = None
    for col in ['perturbation', 'gene', 'target_gene', 'guide']:
        if col in umap_df.columns:
            pert_col = col
            break

    if pert_col is None:
        fig = go.Figure()
        fig.add_annotation(text="No perturbation info in UMAP data", showarrow=False)
        return fig

    fig = go.Figure()

    # Background cells as density contour
    other_cells = umap_df[umap_df[pert_col] != perturbation]

    # Add 2D histogram contour for background density
    fig.add_trace(go.Histogram2dContour(
        x=other_cells["umap_1"],
        y=other_cells["umap_2"],
        colorscale="Greys",
        reversescale=True,
        showscale=False,
        contours=dict(
            showlabels=False,
            coloring="fill",
        ),
        line=dict(width=0.5, color="lightgray"),
        hoverinfo="skip",
        name="All cells",
        opacity=0.6,
    ))

    # Highlighted cells as scatter points
    pert_cells = umap_df[umap_df[pert_col] == perturbation]

    if len(pert_cells) > 0:
        fig.add_trace(go.Scattergl(
            x=pert_cells["umap_1"],
            y=pert_cells["umap_2"],
            mode="markers",
            marker=dict(size=6, color="#e74c3c", line=dict(width=0.5, color="white")),
            name=perturbation,
            hovertemplate=f"{perturbation}<br>UMAP1: %{{x:.2f}}<br>UMAP2: %{{y:.2f}}<extra></extra>",
        ))

    n_pert = len(pert_cells)
    n_total = len(umap_df)

    fig.update_layout(
        title=f"UMAP: {perturbation} ({n_pert:,} / {n_total:,} cells)",
        xaxis_title="UMAP 1",
        yaxis_title="UMAP 2",
        height=500,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            font=dict(size=9),
        ),
        margin=dict(l=40, r=10, t=60, b=40),
    )

    # Force square aspect ratio
    fig.update_xaxes(showgrid=False, zeroline=False, scaleanchor="y", scaleratio=1)
    fig.update_yaxes(showgrid=False, zeroline=False)

    return fig


def plot_edist_rank(df: pd.DataFrame, gene: str, condition: str = "EBs") -> go.Figure:
    """
    E-distance rank plot — gene highlighted in the all-gene distribution.

    Args:
        df: Transcriptome E-distance DataFrame (all perturbations for a condition)
        gene: Gene symbol to highlight
        condition: Condition label for title
    """
    if df is None or len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No E-distance data for {condition}", showarrow=False)
        return fig

    df = df.copy()
    # Extract gene from perturbation
    df["gene"] = df["perturbation"].str.split("_").str[0]
    df = df.sort_values("edist_observed")
    df["rank"] = range(1, len(df) + 1)

    gene_rows = df[df["gene"] == gene]

    fig = go.Figure()

    # All perturbations
    fig.add_trace(go.Scatter(
        x=df["rank"],
        y=df["edist_observed"],
        mode="markers",
        marker=dict(size=4, color="#cccccc"),
        name="All perturbations",
        hovertemplate="%{customdata}<br>E-dist: %{y:.3f}<extra></extra>",
        customdata=df["perturbation"],
    ))

    # Highlight gene
    if len(gene_rows) > 0:
        fig.add_trace(go.Scatter(
            x=gene_rows["rank"],
            y=gene_rows["edist_observed"],
            mode="markers",
            marker=dict(size=12, color="#e74c3c", symbol="diamond"),
            name=gene,
            hovertemplate=f"{gene}<br>Rank: %{{x}}<br>E-dist: %{{y:.3f}}<extra></extra>",
        ))

        for _, row in gene_rows.iterrows():
            fig.add_annotation(
                x=row["rank"], y=row["edist_observed"],
                text=f"{row['perturbation']}",
                showarrow=True, arrowhead=2, yshift=15,
                font=dict(size=9),
            )

    fig.update_layout(
        title=f"Transcriptome E-distance ({condition})",
        xaxis_title="Rank",
        yaxis_title="E-distance (observed)",
        height=350,
        showlegend=False,
        margin=dict(l=10, r=10, t=40, b=10),
    )

    return fig


def plot_dose_response_scatter(df: pd.DataFrame, gene: str, condition: str = "EBs") -> go.Figure:
    """
    Scatter plot: knockdown_pct vs lineage_effect, gene's perturbations highlighted.

    Args:
        df: Dose response DataFrame
        gene: Gene symbol to highlight
        condition: Condition label for title
    """
    if df is None or len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No dose response data for {condition}", showarrow=False)
        return fig

    df = df.copy()
    df["gene"] = df["perturbation"].str.split("_").str[0]
    df["is_gene"] = df["gene"] == gene

    fig = go.Figure()

    # Background
    bg = df[~df["is_gene"]]
    fig.add_trace(go.Scatter(
        x=bg["knockdown_pct"],
        y=bg["lineage_effect"],
        mode="markers",
        marker=dict(size=4, color="#cccccc", opacity=0.5),
        name="All perturbations",
        hovertemplate="%{customdata}<br>KD: %{x:.1f}%<br>Effect: %{y:.3f}<extra></extra>",
        customdata=bg["perturbation"],
    ))

    # Highlighted gene
    fg = df[df["is_gene"]]
    if len(fg) > 0:
        fig.add_trace(go.Scatter(
            x=fg["knockdown_pct"],
            y=fg["lineage_effect"],
            mode="markers+text",
            marker=dict(size=10, color="#e74c3c"),
            name=gene,
            text=fg["perturbation"],
            textposition="top center",
            textfont=dict(size=9),
            hovertemplate=f"{gene}<br>KD: %{{x:.1f}}%<br>Effect: %{{y:.3f}}<extra></extra>",
        ))

    fig.update_layout(
        title=f"Dose Response ({condition})",
        xaxis_title="Knockdown %",
        yaxis_title="Lineage Effect",
        height=350,
        showlegend=False,
        margin=dict(l=10, r=10, t=40, b=10),
    )

    return fig


def plot_pathway_enrichment_bars(df: pd.DataFrame, perturbation: str, n_top: int = 15) -> go.Figure:
    """
    Horizontal bar chart of top enriched terms by -log10(p-value).

    Args:
        df: Pathway enrichment DataFrame filtered to perturbation
        perturbation: Perturbation ID for title
        n_top: Number of top terms to show
    """
    if df is None or len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No pathway enrichment data for {perturbation}", showarrow=False)
        return fig

    df = df.copy()

    # Use adjusted_p_value if available, else p_value
    p_col = "adjusted_p_value" if "adjusted_p_value" in df.columns else "p_value"
    df["neg_log_p"] = -np.log10(df[p_col].clip(lower=1e-300))

    # Top terms
    top = df.nlargest(n_top, "neg_log_p")
    top = top.sort_values("neg_log_p", ascending=True)

    # Truncate long term names
    top["term_short"] = top["term"].str[:60]

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=top["neg_log_p"],
        y=top["term_short"],
        orientation="h",
        marker_color="#c8a57b",
        hovertemplate="<b>%{y}</b><br>-log10(p): %{x:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title=f"Top Enriched Pathways: {perturbation}",
        xaxis_title=f"-log10({p_col})",
        height=max(300, n_top * 25),
        margin=dict(l=10, r=10, t=40, b=10),
        yaxis=dict(tickfont=dict(size=10)),
    )

    return fig


def plot_lineage_de_volcano(
    df: pd.DataFrame,
    perturbation: str,
    lineage: str,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> go.Figure:
    """
    Volcano plot for lineage-specific DE results.

    Uses log2fc/pval/padj column names (lineage DE format).
    Gene labels show ENSG IDs.
    """
    if df is None or len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(text=f"No lineage DE data", showarrow=False)
        return fig

    df = df.copy()

    # Determine column names — lineage DE uses log2fc/padj or log2FoldChange/padj
    lfc_col = "log2fc" if "log2fc" in df.columns else "log2FoldChange"
    padj_col = "padj" if "padj" in df.columns else "pvalue"
    gene_col = "gene" if "gene" in df.columns else df.columns[0]

    if lfc_col not in df.columns or padj_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="Missing required columns (log2fc, padj)", showarrow=False)
        return fig

    df["neg_log_padj"] = -np.log10(df[padj_col].clip(lower=1e-300))

    df["category"] = "Not significant"
    df.loc[(df[padj_col] < padj_threshold) & (df[lfc_col] > lfc_threshold), "category"] = "Up"
    df.loc[(df[padj_col] < padj_threshold) & (df[lfc_col] < -lfc_threshold), "category"] = "Down"

    color_map = {
        "Not significant": "#b0aba5",
        "Up": "#e74c3c",
        "Down": "#4a8eb5",
    }

    fig = px.scatter(
        df,
        x=lfc_col,
        y="neg_log_padj",
        color="category",
        color_discrete_map=color_map,
        hover_data=[gene_col],
        labels={
            lfc_col: "Log2 Fold Change",
            "neg_log_padj": "-Log10 Adjusted P-value",
        },
    )

    fig.add_hline(y=-np.log10(padj_threshold), line_dash="dash", line_color="gray")
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")

    # Label top genes (ENSG IDs)
    top_genes = df.nlargest(5, "neg_log_padj")
    for _, row in top_genes.iterrows():
        gene_label = str(row.get(gene_col, ""))[:15]
        fig.add_annotation(
            x=row[lfc_col],
            y=row["neg_log_padj"],
            text=gene_label,
            showarrow=True,
            arrowhead=0,
            font=dict(size=8),
        )

    n_up = (df["category"] == "Up").sum()
    n_down = (df["category"] == "Down").sum()

    fig.update_layout(
        title=f"Lineage DE: {perturbation} / {lineage} ({n_up} up, {n_down} down)",
        height=400,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        margin=dict(l=10, r=10, t=60, b=10),
    )

    return fig
