"""
Publication-quality plotting functions for COVID-19 single-cell analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
from typing import Optional, List, Dict, Union, Tuple
from pathlib import Path

# Set publication style defaults
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.transparent': False,
})

# Color palettes
COVID_COLORS = {
    'COVID': '#E74C3C',
    'Control': '#3498DB',
}

CELL_TYPE_COLORS = {
    # Major types
    'Epithelial': '#2ECC71',
    'Myeloid': '#E74C3C',
    'Fibroblast': '#9B59B6',
    'Endothelial': '#F39C12',
    'T_cell': '#3498DB',
    'B_cell': '#1ABC9C',
    'NK_cell': '#E91E63',
    'Plasma_cell': '#00BCD4',
    'Mast_cell': '#795548',

    # Epithelial subtypes
    'AT1': '#27AE60',
    'AT2': '#2ECC71',
    'DATP': '#F1C40F',
    'Basal': '#8BC34A',
    'Club': '#CDDC39',
    'Ciliated': '#009688',

    # Myeloid subtypes
    'AM': '#C0392B',
    'MDM': '#E74C3C',
    'Classical_mono': '#D32F2F',
    'DC': '#FF5722',
    'pDC': '#FF7043',
    'Neutrophil': '#FF8A65',

    # Fibroblast subtypes
    'Alveolar_FB': '#7B1FA2',
    'Adventitial_FB': '#9C27B0',
    'pFB': '#E91E63',
    'Myofibroblast': '#F48FB1',

    # T cell subtypes
    'CD4_T': '#1976D2',
    'CD8_T': '#2196F3',
    'Treg': '#03A9F4',
    'Exhausted_T': '#00BCD4',
    'Cytotoxic_T': '#0097A7',

    # B cell subtypes
    'Naive_B': '#00897B',
    'Memory_B': '#009688',
    'Plasma': '#4DB6AC',
}


# =============================================================================
# UMAP Plots
# =============================================================================

def plot_umap_celltype(
    adata,
    color_col: str = 'cell_type',
    title: str = 'Cell Types',
    palette: Optional[Dict] = None,
    figsize: Tuple[int, int] = (8, 6),
    legend_loc: str = 'right margin',
    save_path: Optional[str] = None,
    **kwargs
) -> None:
    """
    Plot UMAP colored by cell type.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with UMAP coordinates
    color_col : str
        Column to color by
    title : str
        Plot title
    palette : dict, optional
        Color palette
    figsize : tuple
        Figure size
    legend_loc : str
        Legend location
    save_path : str, optional
        Path to save figure
    """
    if palette is None:
        palette = CELL_TYPE_COLORS

    fig, ax = plt.subplots(figsize=figsize)

    sc.pl.umap(
        adata,
        color=color_col,
        palette=palette,
        legend_loc=legend_loc,
        title=title,
        frameon=True,
        ax=ax,
        show=False,
        **kwargs
    )

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


def plot_umap_condition(
    adata,
    condition_col: str = 'condition',
    title: str = 'COVID vs Control',
    figsize: Tuple[int, int] = (10, 5),
    save_path: Optional[str] = None
) -> None:
    """Plot UMAP split by condition."""
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    for ax, cond in zip(axes, ['Control', 'COVID']):
        mask = adata.obs[condition_col] == cond
        sc.pl.umap(
            adata[mask],
            color=condition_col,
            palette=COVID_COLORS,
            title=f'{cond} (n={mask.sum():,})',
            ax=ax,
            show=False,
            legend_loc=None
        )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


def plot_umap_expression(
    adata,
    genes: List[str],
    ncols: int = 4,
    figsize: Tuple[int, int] = (4, 4),
    cmap: str = 'viridis',
    save_path: Optional[str] = None,
    **kwargs
) -> None:
    """Plot UMAP colored by gene expression."""
    # Filter to available genes
    available = [g for g in genes if g in adata.var_names]
    if not available:
        print("No genes found in data")
        return

    nrows = int(np.ceil(len(available) / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize[0] * ncols, figsize[1] * nrows)
    )
    axes = np.atleast_2d(axes).flatten()

    for i, gene in enumerate(available):
        sc.pl.umap(
            adata,
            color=gene,
            cmap=cmap,
            ax=axes[i],
            show=False,
            frameon=True,
            **kwargs
        )

    # Hide unused axes
    for j in range(len(available), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


# =============================================================================
# Expression Plots
# =============================================================================

def plot_dotplot(
    adata,
    marker_genes: Dict[str, List[str]],
    groupby: str = 'cell_type',
    figsize: Optional[Tuple[int, int]] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> None:
    """
    Create a dotplot of marker gene expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    marker_genes : dict
        Dictionary of {cell_type: [genes]}
    groupby : str
        Column to group by
    figsize : tuple, optional
        Figure size
    save_path : str, optional
        Path to save figure
    """
    sc.pl.dotplot(
        adata,
        marker_genes,
        groupby=groupby,
        dendrogram=True,
        standard_scale='var',
        show=False,
        **kwargs
    )

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


def plot_violin_comparison(
    adata,
    genes: List[str],
    groupby: str = 'cell_type',
    split_by: str = 'condition',
    figsize: Tuple[int, int] = (3, 3),
    save_path: Optional[str] = None,
    **kwargs
) -> None:
    """
    Plot violin plots comparing expression across conditions.
    """
    available = [g for g in genes if g in adata.var_names]
    if not available:
        print("No genes found")
        return

    ncols = min(4, len(available))
    nrows = int(np.ceil(len(available) / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize[0] * ncols, figsize[1] * nrows)
    )
    axes = np.atleast_2d(axes).flatten()

    for i, gene in enumerate(available):
        sc.pl.violin(
            adata,
            gene,
            groupby=groupby,
            ax=axes[i],
            show=False,
            rotation=90,
            **kwargs
        )

    for j in range(len(available), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


def plot_heatmap_top_genes(
    adata,
    groupby: str = 'cell_type',
    n_genes: int = 10,
    figsize: Tuple[int, int] = (12, 10),
    save_path: Optional[str] = None,
    **kwargs
) -> None:
    """Plot heatmap of top marker genes per group."""
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=n_genes,
        groupby=groupby,
        show_gene_labels=True,
        show=False,
        **kwargs
    )

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


# =============================================================================
# Cell Proportion Plots
# =============================================================================

def plot_cell_proportions(
    adata,
    cell_type_col: str = 'cell_type',
    condition_col: str = 'condition',
    sample_col: str = 'sample_id',
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None
) -> None:
    """
    Plot cell type proportions per sample.
    """
    # Calculate proportions
    props = adata.obs.groupby([sample_col, cell_type_col]).size().unstack(fill_value=0)
    props = props.div(props.sum(axis=1), axis=0)

    # Add condition info
    sample_condition = adata.obs.groupby(sample_col)[condition_col].first()

    # Sort by condition
    props['condition'] = sample_condition
    props = props.sort_values('condition')
    condition = props.pop('condition')

    # Plot stacked bar
    fig, ax = plt.subplots(figsize=figsize)

    props.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=[CELL_TYPE_COLORS.get(c, '#888888') for c in props.columns],
        edgecolor='white',
        linewidth=0.5
    )

    # Add condition labels
    colors = [COVID_COLORS[c] for c in condition]
    for i, (sample, color) in enumerate(zip(props.index, colors)):
        ax.get_xticklabels()[i].set_color(color)

    ax.set_xlabel('Sample')
    ax.set_ylabel('Proportion')
    ax.set_title('Cell Type Proportions per Sample')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

    # Add condition legend
    patches = [mpatches.Patch(color=c, label=l) for l, c in COVID_COLORS.items()]
    ax.legend(
        handles=ax.get_legend_handles_labels()[0] + patches,
        bbox_to_anchor=(1.02, 1),
        loc='upper left'
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


def plot_proportion_comparison(
    adata,
    cell_type_col: str = 'cell_type',
    condition_col: str = 'condition',
    sample_col: str = 'sample_id',
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None
) -> None:
    """
    Box plot comparing cell type proportions between conditions.
    """
    # Calculate proportions per sample
    props = adata.obs.groupby([sample_col, cell_type_col]).size().unstack(fill_value=0)
    props = props.div(props.sum(axis=1), axis=0)

    # Add condition
    sample_condition = adata.obs.groupby(sample_col)[condition_col].first()
    props['condition'] = sample_condition

    # Melt for plotting
    props_melt = props.melt(
        id_vars='condition',
        var_name='cell_type',
        value_name='proportion'
    )

    fig, ax = plt.subplots(figsize=figsize)

    sns.boxplot(
        data=props_melt,
        x='cell_type',
        y='proportion',
        hue='condition',
        palette=COVID_COLORS,
        ax=ax
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Proportion')
    ax.set_title('Cell Type Proportions: COVID vs Control')
    ax.legend(title='Condition')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


# =============================================================================
# Differential Expression Plots
# =============================================================================

def plot_volcano(
    de_results: pd.DataFrame,
    logfc_col: str = 'log2FoldChange',
    pval_col: str = 'padj',
    gene_col: Optional[str] = None,
    logfc_thresh: float = 1.0,
    pval_thresh: float = 0.05,
    top_n_labels: int = 10,
    figsize: Tuple[int, int] = (8, 6),
    title: str = 'Differential Expression',
    save_path: Optional[str] = None
) -> None:
    """
    Create a volcano plot of differential expression results.
    """
    df = de_results.copy()

    if gene_col is None:
        df['gene'] = df.index
    else:
        df['gene'] = df[gene_col]

    # Calculate -log10 p-value
    df['-log10_pval'] = -np.log10(df[pval_col].clip(lower=1e-300))

    # Classify significance
    df['significance'] = 'NS'
    df.loc[
        (df[logfc_col] > logfc_thresh) & (df[pval_col] < pval_thresh),
        'significance'
    ] = 'Up'
    df.loc[
        (df[logfc_col] < -logfc_thresh) & (df[pval_col] < pval_thresh),
        'significance'
    ] = 'Down'

    colors = {'Up': '#E74C3C', 'Down': '#3498DB', 'NS': '#AAAAAA'}

    fig, ax = plt.subplots(figsize=figsize)

    for sig, color in colors.items():
        mask = df['significance'] == sig
        ax.scatter(
            df.loc[mask, logfc_col],
            df.loc[mask, '-log10_pval'],
            c=color,
            alpha=0.5,
            s=10,
            label=f'{sig} ({mask.sum()})'
        )

    # Add threshold lines
    ax.axhline(-np.log10(pval_thresh), ls='--', c='gray', alpha=0.5)
    ax.axvline(logfc_thresh, ls='--', c='gray', alpha=0.5)
    ax.axvline(-logfc_thresh, ls='--', c='gray', alpha=0.5)

    # Label top genes
    top_up = df[df['significance'] == 'Up'].nlargest(top_n_labels, '-log10_pval')
    top_down = df[df['significance'] == 'Down'].nlargest(top_n_labels, '-log10_pval')

    for _, row in pd.concat([top_up, top_down]).iterrows():
        ax.annotate(
            row['gene'],
            (row[logfc_col], row['-log10_pval']),
            fontsize=8,
            alpha=0.8
        )

    ax.set_xlabel(f'Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value')
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


# =============================================================================
# Trajectory Plots
# =============================================================================

def plot_trajectory_genes(
    adata,
    genes: List[str],
    pseudotime_col: str = 'dpt_pseudotime',
    color_by: str = 'cell_type',
    ncols: int = 3,
    figsize: Tuple[int, int] = (4, 3),
    save_path: Optional[str] = None
) -> None:
    """
    Plot gene expression along pseudotime.
    """
    available = [g for g in genes if g in adata.var_names]
    if not available:
        print("No genes found")
        return

    nrows = int(np.ceil(len(available) / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize[0] * ncols, figsize[1] * nrows)
    )
    axes = np.atleast_2d(axes).flatten()

    pseudotime = adata.obs[pseudotime_col]

    for i, gene in enumerate(available):
        ax = axes[i]

        if gene in adata.var_names:
            expr = adata[:, gene].X.toarray().flatten()
        else:
            continue

        # Scatter plot
        scatter = ax.scatter(
            pseudotime,
            expr,
            c=adata.obs[color_by].astype('category').cat.codes,
            alpha=0.3,
            s=5,
            cmap='tab20'
        )

        # Add loess smoothing (approximate with rolling mean)
        sorted_idx = np.argsort(pseudotime)
        window = max(len(pseudotime) // 50, 10)
        smoothed = pd.Series(expr[sorted_idx]).rolling(window, center=True).mean()
        ax.plot(pseudotime[sorted_idx], smoothed, 'r-', linewidth=2)

        ax.set_xlabel('Pseudotime')
        ax.set_ylabel('Expression')
        ax.set_title(gene)

    for j in range(len(available), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


# =============================================================================
# QC Plots
# =============================================================================

def plot_qc_metrics(
    adata,
    groupby: str = 'sample_id',
    figsize: Tuple[int, int] = (15, 4),
    save_path: Optional[str] = None
) -> None:
    """
    Plot QC metric distributions.
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    titles = ['Genes per Cell', 'UMIs per Cell', 'Mitochondrial %']

    for ax, metric, title in zip(axes, metrics, titles):
        if metric not in adata.obs:
            continue

        sc.pl.violin(
            adata,
            metric,
            groupby=groupby,
            ax=ax,
            show=False,
            rotation=90
        )
        ax.set_title(title)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"Saved: {save_path}")

    plt.show()


# =============================================================================
# Save Utilities
# =============================================================================

def save_figure(
    fig,
    path: str,
    formats: List[str] = ['png', 'pdf'],
    dpi: int = 300
) -> None:
    """
    Save figure in multiple formats.
    """
    path = Path(path)

    for fmt in formats:
        save_path = path.with_suffix(f'.{fmt}')
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {save_path}")
