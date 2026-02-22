"""
Cell type marker genes and gene signatures for COVID-19 lung analysis.

Based on Melms et al. 2021 "A molecular single-cell lung atlas of lethal COVID-19"
"""

from typing import Dict, List

# =============================================================================
# Major Cell Type Markers
# =============================================================================

MAJOR_CELL_TYPES: Dict[str, List[str]] = {
    # Epithelial
    'Epithelial': ['EPCAM', 'KRT18', 'KRT19', 'CDH1'],

    # Myeloid
    'Myeloid': ['CD68', 'CD14', 'FCGR3A', 'LYZ', 'CSF1R'],

    # Fibroblasts
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA'],

    # Endothelial
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'CLDN5'],

    # T cells
    'T_cell': ['CD3D', 'CD3E', 'CD3G', 'TRAC'],

    # B cells
    'B_cell': ['CD79A', 'MS4A1', 'CD19', 'PAX5'],

    # NK cells
    'NK_cell': ['NCAM1', 'NKG7', 'GNLY', 'KLRD1'],

    # Plasma cells
    'Plasma_cell': ['JCHAIN', 'MZB1', 'XBP1', 'IGHG1'],

    # Mast cells
    'Mast_cell': ['TPSAB1', 'CPA3', 'KIT', 'HPGDS'],
}

# =============================================================================
# Fine Cell Type Markers
# =============================================================================

# Epithelial subtypes
EPITHELIAL_SUBTYPES: Dict[str, List[str]] = {
    # Alveolar Type 1
    'AT1': ['AGER', 'PDPN', 'CLIC5', 'HOPX', 'CAV1', 'AQP5'],

    # Alveolar Type 2
    'AT2': ['SFTPC', 'SFTPA1', 'SFTPA2', 'SFTPB', 'ABCA3', 'SLC34A2', 'LAMP3'],

    # Damage-Associated Transient Progenitors (DATP/KRT8+)
    'DATP': ['KRT8', 'CLDN4', 'CDKN1A', 'KRT17', 'SOX4', 'TP63'],

    # Basal cells
    'Basal': ['KRT5', 'KRT14', 'TP63', 'NGFR'],

    # Club cells
    'Club': ['SCGB1A1', 'SCGB3A1', 'CYP2F1'],

    # Ciliated cells
    'Ciliated': ['FOXJ1', 'PIFO', 'CAPS', 'TPPP3', 'DNAH5'],

    # Secretory cells
    'Secretory': ['MUC5B', 'MUC5AC', 'SPDEF'],
}

# Myeloid subtypes
MYELOID_SUBTYPES: Dict[str, List[str]] = {
    # Alveolar macrophages
    'AM': ['MARCO', 'FABP4', 'MCEMP1', 'PPARG', 'SIGLEC1'],

    # Monocyte-derived macrophages
    'MDM': ['FCN1', 'S100A8', 'S100A9', 'S100A12', 'VCAN'],

    # Classical monocytes
    'Classical_mono': ['CD14', 'LYZ', 'S100A8', 'S100A9'],

    # Non-classical monocytes
    'Nonclassical_mono': ['FCGR3A', 'MS4A7', 'CX3CR1'],

    # Dendritic cells
    'DC': ['CD1C', 'CLEC10A', 'FCER1A'],

    # Plasmacytoid DC
    'pDC': ['LILRA4', 'IL3RA', 'CLEC4C', 'IRF7'],

    # Neutrophils
    'Neutrophil': ['FCGR3B', 'CSF3R', 'CXCR2', 'NAMPT'],
}

# Fibroblast subtypes
FIBROBLAST_SUBTYPES: Dict[str, List[str]] = {
    # Alveolar fibroblasts
    'Alveolar_FB': ['NPNT', 'WNT2', 'LIMCH1'],

    # Adventitial fibroblasts
    'Adventitial_FB': ['PI16', 'MFAP5', 'SFRP2'],

    # Pathological fibroblasts (COVID-associated)
    'pFB': ['CTHRC1', 'COL1A1', 'POSTN', 'TNC', 'COL3A1'],

    # Myofibroblasts
    'Myofibroblast': ['ACTA2', 'MYH11', 'TAGLN', 'CNN1'],

    # Lipofibroblasts
    'Lipofibroblast': ['PLIN2', 'APOE', 'TCF21'],
}

# T cell subtypes
TCELL_SUBTYPES: Dict[str, List[str]] = {
    # CD4 T cells
    'CD4_T': ['CD4', 'IL7R', 'TCF7', 'LEF1'],

    # CD8 T cells
    'CD8_T': ['CD8A', 'CD8B', 'GZMK', 'GZMB'],

    # Regulatory T cells
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],

    # Exhausted T cells
    'Exhausted_T': ['PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'TOX'],

    # Cytotoxic T cells
    'Cytotoxic_T': ['GZMB', 'PRF1', 'GNLY', 'NKG7', 'IFNG'],

    # Proliferating T cells
    'Proliferating_T': ['MKI67', 'TOP2A', 'PCNA'],

    # Naive T cells
    'Naive_T': ['CCR7', 'SELL', 'TCF7', 'LEF1'],
}

# B cell subtypes
BCELL_SUBTYPES: Dict[str, List[str]] = {
    # Naive B cells
    'Naive_B': ['MS4A1', 'IGHD', 'TCL1A'],

    # Memory B cells
    'Memory_B': ['CD27', 'AIM2', 'TNFRSF13B'],

    # Germinal center B cells
    'GC_B': ['BCL6', 'AICDA', 'RGS13'],

    # Plasma cells
    'Plasma': ['JCHAIN', 'MZB1', 'XBP1', 'IGHG1', 'IGHA1'],

    # Plasmablasts
    'Plasmablast': ['PRDM1', 'IRF4', 'TNFRSF17'],
}

# =============================================================================
# COVID-19 Specific Gene Signatures
# =============================================================================

# Myeloid dysregulation signature (COVID)
MYELOID_DYSFUNCTION: Dict[str, List[str]] = {
    'lncRNA_high': ['NEAT1', 'MALAT1'],
    'efferocytosis_low': ['AXL', 'MERTK', 'GAS6'],
    'pro_inflammatory': ['IL1B', 'IL6', 'TNF', 'CCL2', 'CCL3', 'CXCL8'],
}

# DATP signature (transitional AT2)
DATP_SIGNATURE: List[str] = [
    'KRT8', 'CLDN4', 'CDKN1A', 'KRT17', 'SOX4',
    'AREG', 'HBEGF', 'MMP7', 'CXCL8', 'IL33'
]

# Fibrosis signature
FIBROSIS_SIGNATURE: List[str] = [
    'COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL6A1',
    'FN1', 'POSTN', 'TNC', 'CTHRC1', 'ACTA2'
]

# T cell exhaustion signature
EXHAUSTION_SIGNATURE: List[str] = [
    'PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'CTLA4',
    'TOX', 'ENTPD1', 'CXCL13'
]

# T cell cytotoxicity signature
CYTOTOXICITY_SIGNATURE: List[str] = [
    'GZMB', 'GZMA', 'GZMH', 'GZMK', 'PRF1',
    'GNLY', 'NKG7', 'IFNG', 'FASLG'
]

# Interferon response signature
IFN_RESPONSE_SIGNATURE: List[str] = [
    'ISG15', 'IFI6', 'IFI27', 'IFI44L', 'IFIT1',
    'IFIT2', 'IFIT3', 'MX1', 'MX2', 'OAS1',
    'OAS2', 'OAS3', 'STAT1', 'IRF7'
]

# Inflammatory cytokines
INFLAMMATORY_CYTOKINES: List[str] = [
    'IL1B', 'IL6', 'TNF', 'CXCL8', 'CCL2',
    'CCL3', 'CCL4', 'CXCL10', 'IL1A'
]

# =============================================================================
# Transcription Factor Lists
# =============================================================================

# Key TFs in COVID fibrosis
FIBROSIS_TFS: List[str] = [
    'JUNB', 'JUND', 'FOS', 'JUN', 'FOSB',
    'ATF3', 'ATF4', 'CEBPB', 'NFKB1', 'RELA'
]

# Key TFs in DATP
DATP_TFS: List[str] = [
    'TP63', 'SOX4', 'SOX9', 'FOXA2', 'NKX2-1',
    'HIF1A', 'CEBPD', 'YAP1'
]

# =============================================================================
# Cell Communication Genes
# =============================================================================

# Key ligands in COVID
COVID_LIGANDS: List[str] = [
    'IL1B', 'IL6', 'TNF', 'TGFB1', 'TGFB2',
    'PDGFA', 'PDGFB', 'FGF2', 'FGF7', 'WNT5A',
    'CCL2', 'CCL3', 'CXCL8', 'CXCL10', 'CXCL12'
]

# Key receptors in COVID
COVID_RECEPTORS: List[str] = [
    'IL1R1', 'IL6R', 'IL6ST', 'TNFRSF1A', 'TGFBR1',
    'TGFBR2', 'PDGFRA', 'PDGFRB', 'FGFR1', 'FGFR2',
    'CCR2', 'CXCR4', 'ACE2'
]

# =============================================================================
# Helper Functions
# =============================================================================

def get_all_markers() -> Dict[str, List[str]]:
    """Return all marker genes as a single dictionary."""
    all_markers = {}
    all_markers.update(MAJOR_CELL_TYPES)
    all_markers.update(EPITHELIAL_SUBTYPES)
    all_markers.update(MYELOID_SUBTYPES)
    all_markers.update(FIBROBLAST_SUBTYPES)
    all_markers.update(TCELL_SUBTYPES)
    all_markers.update(BCELL_SUBTYPES)
    return all_markers


def get_markers_for_celltype(cell_type: str) -> List[str]:
    """Get marker genes for a specific cell type."""
    all_markers = get_all_markers()
    return all_markers.get(cell_type, [])


def get_signature_genes(signature_name: str) -> List[str]:
    """Get genes for a named signature."""
    signatures = {
        'DATP': DATP_SIGNATURE,
        'fibrosis': FIBROSIS_SIGNATURE,
        'exhaustion': EXHAUSTION_SIGNATURE,
        'cytotoxicity': CYTOTOXICITY_SIGNATURE,
        'interferon': IFN_RESPONSE_SIGNATURE,
        'inflammatory': INFLAMMATORY_CYTOKINES,
    }
    return signatures.get(signature_name, [])


# Combine all markers for easy iteration
ALL_CELL_TYPE_MARKERS = get_all_markers()
