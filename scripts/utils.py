"""
Utility functions for COVID-19 single-cell analysis.
"""

import os
import gzip
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
from typing import Optional, List, Dict, Union
from scipy import sparse


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_csv_to_anndata(
    file_path: str,
    sample_id: str,
    transpose: bool = True
) -> ad.AnnData:
    """
    Load a gzipped CSV file into AnnData format.

    Parameters
    ----------
    file_path : str
        Path to the CSV.gz file
    sample_id : str
        Sample identifier to add to obs
    transpose : bool
        Whether to transpose (genes as rows -> cells as rows)

    Returns
    -------
    AnnData object with counts
    """
    # Read CSV (genes as rows, cells as columns)
    df = pd.read_csv(file_path, index_col=0)

    if transpose:
        df = df.T

    # Create AnnData
    adata = ad.AnnData(
        X=sparse.csr_matrix(df.values.astype(np.float32)),
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns)
    )

    # Add sample info
    adata.obs['sample_id'] = sample_id

    return adata


def load_all_samples(
    data_dir: str,
    metadata_path: str,
    verbose: bool = True
) -> ad.AnnData:
    """
    Load all samples and concatenate into a single AnnData.

    Parameters
    ----------
    data_dir : str
        Directory containing raw CSV.gz files
    metadata_path : str
        Path to sample metadata CSV
    verbose : bool
        Print progress

    Returns
    -------
    Concatenated AnnData with all samples
    """
    # Load metadata
    metadata = pd.read_csv(metadata_path)

    adatas = []
    for _, row in metadata.iterrows():
        file_path = os.path.join(data_dir, row['file_name'])
        sample_id = row['sample_id']

        if verbose:
            print(f"Loading {sample_id}...")

        adata = load_csv_to_anndata(file_path, sample_id)
        adata.obs['condition'] = row['condition']
        adata.obs['patient_id'] = row['patient_id']

        adatas.append(adata)

    # Concatenate
    if verbose:
        print("Concatenating samples...")

    adata = ad.concat(
        adatas,
        join='outer',
        merge='same',
        uns_merge='same',
        index_unique='_'
    )

    # Make var_names unique
    adata.var_names_make_unique()

    if verbose:
        print(f"Total cells: {adata.n_obs:,}")
        print(f"Total genes: {adata.n_vars:,}")

    return adata


# =============================================================================
# Quality Control Functions
# =============================================================================

def calculate_qc_metrics(
    adata: ad.AnnData,
    mito_prefix: str = 'MT-',
    ribo_prefix: str = 'RPS|RPL'
) -> ad.AnnData:
    """
    Calculate standard QC metrics for snRNA-seq data.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    mito_prefix : str
        Prefix for mitochondrial genes
    ribo_prefix : str
        Pattern for ribosomal genes

    Returns
    -------
    AnnData with QC metrics in obs
    """
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith(mito_prefix)

    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.contains(ribo_prefix)

    # Calculate metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    return adata


def mad_filter(
    values: np.ndarray,
    n_mads: float = 5.0,
    direction: str = 'both'
) -> np.ndarray:
    """
    Calculate MAD-based outlier thresholds.

    Parameters
    ----------
    values : np.ndarray
        Values to analyze
    n_mads : float
        Number of MADs for threshold
    direction : str
        'upper', 'lower', or 'both'

    Returns
    -------
    Boolean mask of cells to keep
    """
    median = np.median(values)
    mad = np.median(np.abs(values - median))

    if mad == 0:
        return np.ones(len(values), dtype=bool)

    if direction == 'upper':
        return values <= median + n_mads * mad
    elif direction == 'lower':
        return values >= median - n_mads * mad
    else:
        return (values >= median - n_mads * mad) & (values <= median + n_mads * mad)


def filter_cells_mad(
    adata: ad.AnnData,
    n_genes_mads: float = 5.0,
    n_counts_mads: float = 5.0,
    mito_threshold: float = 20.0,
    min_genes: int = 200,
    min_counts: int = 500,
    verbose: bool = True
) -> ad.AnnData:
    """
    Filter cells using MAD-based thresholds (appropriate for snRNA-seq).

    Parameters
    ----------
    adata : AnnData
        Input AnnData with QC metrics
    n_genes_mads : float
        Number of MADs for gene count filtering
    n_counts_mads : float
        Number of MADs for UMI count filtering
    mito_threshold : float
        Maximum mitochondrial percentage
    min_genes : int
        Minimum genes per cell
    min_counts : int
        Minimum UMIs per cell
    verbose : bool
        Print filtering statistics

    Returns
    -------
    Filtered AnnData
    """
    n_before = adata.n_obs

    # MAD-based upper bounds
    genes_keep = mad_filter(adata.obs['n_genes_by_counts'], n_genes_mads, 'upper')
    counts_keep = mad_filter(adata.obs['total_counts'], n_counts_mads, 'upper')

    # Fixed lower bounds and mito threshold
    min_genes_keep = adata.obs['n_genes_by_counts'] >= min_genes
    min_counts_keep = adata.obs['total_counts'] >= min_counts
    mito_keep = adata.obs['pct_counts_mt'] < mito_threshold

    # Combine filters
    keep = genes_keep & counts_keep & min_genes_keep & min_counts_keep & mito_keep

    adata = adata[keep].copy()

    if verbose:
        print(f"Cells before filtering: {n_before:,}")
        print(f"Cells after filtering: {adata.n_obs:,}")
        print(f"Cells removed: {n_before - adata.n_obs:,} ({(n_before - adata.n_obs)/n_before*100:.1f}%)")

    return adata


def run_scrublet(
    adata: ad.AnnData,
    batch_key: str = 'sample_id',
    threshold: Optional[float] = None,
    verbose: bool = True
) -> ad.AnnData:
    """
    Run Scrublet doublet detection per sample.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    batch_key : str
        Column for sample batches
    threshold : float, optional
        Manual doublet score threshold
    verbose : bool
        Print progress

    Returns
    -------
    AnnData with doublet scores and predictions
    """
    import scrublet as scr

    doublet_scores = np.zeros(adata.n_obs)
    predicted_doublets = np.zeros(adata.n_obs, dtype=bool)

    for sample in adata.obs[batch_key].unique():
        if verbose:
            print(f"Running Scrublet on {sample}...")

        mask = adata.obs[batch_key] == sample
        sample_adata = adata[mask].copy()

        # Run Scrublet
        scrub = scr.Scrublet(sample_adata.X)
        scores, predictions = scrub.scrub_doublets(verbose=False)

        if threshold is not None:
            predictions = scores > threshold

        doublet_scores[mask.values] = scores
        predicted_doublets[mask.values] = predictions

    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets

    if verbose:
        n_doublets = predicted_doublets.sum()
        print(f"Total predicted doublets: {n_doublets:,} ({n_doublets/adata.n_obs*100:.1f}%)")

    return adata


# =============================================================================
# Gene Scoring Functions
# =============================================================================

def score_gene_signature(
    adata: ad.AnnData,
    gene_list: List[str],
    score_name: str,
    ctrl_size: int = 50,
    use_raw: bool = False
) -> ad.AnnData:
    """
    Score cells for a gene signature.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    gene_list : list
        List of genes in signature
    score_name : str
        Name for the score column
    ctrl_size : int
        Number of control genes
    use_raw : bool
        Use raw counts

    Returns
    -------
    AnnData with signature score
    """
    # Filter to available genes
    available_genes = [g for g in gene_list if g in adata.var_names]

    if len(available_genes) == 0:
        print(f"Warning: No genes from {score_name} found in data")
        adata.obs[score_name] = 0
        return adata

    sc.tl.score_genes(
        adata,
        available_genes,
        score_name=score_name,
        ctrl_size=ctrl_size,
        use_raw=use_raw
    )

    return adata


# =============================================================================
# Pseudobulk Functions
# =============================================================================

def create_pseudobulk(
    adata: ad.AnnData,
    group_by: List[str],
    layer: Optional[str] = None,
    use_raw: bool = True
) -> pd.DataFrame:
    """
    Create pseudobulk expression matrix.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    group_by : list
        Columns to group by (e.g., ['sample_id', 'cell_type'])
    layer : str, optional
        Layer to use
    use_raw : bool
        Use raw counts

    Returns
    -------
    DataFrame with pseudobulk counts (samples x genes)
    """
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    elif layer is not None:
        X = adata.layers[layer]
        var_names = adata.var_names
    else:
        X = adata.X
        var_names = adata.var_names

    # Convert to dense if sparse
    if sparse.issparse(X):
        X = X.toarray()

    # Create group labels
    groups = adata.obs[group_by].astype(str).agg('_'.join, axis=1)

    # Aggregate
    df = pd.DataFrame(X, index=adata.obs_names, columns=var_names)
    df['group'] = groups.values

    pseudobulk = df.groupby('group').sum()

    return pseudobulk


# =============================================================================
# File I/O
# =============================================================================

def save_adata(
    adata: ad.AnnData,
    path: str,
    compression: str = 'gzip'
) -> None:
    """Save AnnData with compression."""
    adata.write_h5ad(path, compression=compression)
    print(f"Saved: {path}")


def load_adata(path: str) -> ad.AnnData:
    """Load AnnData from h5ad file."""
    return sc.read_h5ad(path)


# =============================================================================
# Statistical Helpers
# =============================================================================

def calculate_cell_proportions(
    adata: ad.AnnData,
    cell_type_col: str,
    group_col: str
) -> pd.DataFrame:
    """
    Calculate cell type proportions per group.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    cell_type_col : str
        Column with cell type labels
    group_col : str
        Column to group by (e.g., 'sample_id')

    Returns
    -------
    DataFrame with proportions (groups x cell types)
    """
    counts = adata.obs.groupby([group_col, cell_type_col]).size().unstack(fill_value=0)
    proportions = counts.div(counts.sum(axis=1), axis=0)
    return proportions
