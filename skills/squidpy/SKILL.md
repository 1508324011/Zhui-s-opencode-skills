---
name: squidpy
description: Spatial transcriptomics analysis toolkit. Use for analyzing spatial gene expression data, cell-cell interactions in tissue context, spatial statistics, and integrative analysis of single-cell and spatial data. Best for understanding RBP spatial distribution, tissue architecture, and neighborhood analysis. Part of the scverse ecosystem.
license: BSD-3
compatibility: opencode
metadata:
  category: spatial-transcriptomics
  language: python
  skill-author: K-Dense Inc.
---

# Squidpy: Spatial Transcriptomics Analysis

## Overview

Squidpy is a Python tool for the analysis and visualization of spatial transcriptomics data. It is part of the scverse ecosystem and integrates seamlessly with Scanpy and AnnData. Squidpy enables the analysis of gene expression in its spatial context, allowing researchers to study tissue architecture, cell-cell interactions, and spatial patterns of gene expression.

**Current Version**: Squidpy 1.4.0
**Key Features**: Spatial statistics, neighborhood analysis, image analysis, spatial graphs

## When to Use This Skill

Use Squidpy when you need to:

- **Analyze Visium/Slide-seq/MERFISH data** and other spatial transcriptomics platforms
- **Study spatial distribution of RBP expression** in tissues
- **Identify spatially variable genes** (SVGs)
- **Analyze cell-cell interactions** in spatial context
- **Perform neighborhood enrichment analysis**
- **Integrate single-cell and spatial data**
- **Extract features from histological images**
- **Analyze spatial patterns** and tissue organization

## Installation and Setup

```bash
# Basic installation
pip install squidpy

# With all optional dependencies (recommended)
pip install squidpy[all]

# Core dependencies
pip install scanpy anndata squidpy
```

## Core Capabilities

### 1. Spatial Data Import

**Supported Technologies**:
- 10x Genomics Visium
- Slide-seq / Slide-seqV2
- MERFISH / seqFISH
- CODEX / PhenoCycler
- STARmap / osmFISH
- Custom spatial data

### 2. Spatial Statistics

- **Ripley's K**: Point pattern analysis
- **Moran's I**: Spatial autocorrelation
- **Geary's c**: Spatial heterogeneity
- **Variogram**: Spatial continuity
- **Co-occurrence**: Cell type spatial relationships

### 3. Neighborhood Analysis

- **Cell type deconvolution**: Estimate cell proportions
- **Ligand-receptor analysis**: Predict cell-cell communication
- **Spatially variable genes**: Find genes with spatial patterns
- **Niche analysis**: Identify cellular neighborhoods

### 4. Image Analysis

- **Feature extraction**: Texture, intensity, morphology
- **Segmentation**: Cell/nucleus segmentation
- **Tissue classification**: Identify tissue regions
- **Registration**: Align images with expression data

## Using This Skill

### Example 1: Load and Visualize Visium Data

```python
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt

# Load Visium data
adata = sq.datasets.visium_fluo_adata()
img = sq.datasets.visium_fluo_image()

# Visualize spatial expression
sq.pl.spatial_scatter(
    adata, 
    color='cluster',
    library_id='V1_Adult_Mouse_Brain',
    img=img
)

# Show RBP expression spatially
sq.pl.spatial_scatter(
    adata,
    color=['FUS', 'TARDBP', 'HNRNPA1'],  # RBP genes
    cmap='magma'
)
```

### Example 2: Identify Spatially Variable Genes

```python
# Compute spatial statistics
sq.gr.spatial_neighbors(adata)

# Find spatially variable genes
sq.gr.spatial_autocorr(adata, mode='moran')

# Top spatially variable RBPs
spatial_rbps = adata.uns['moranI'].head(20)
print("Top spatially variable RBPs:")
print(spatial_rbps[spatial_rbps.index.str.contains('RBP|FUS|TARDBP')])

# Visualize
sc.pl.spatial(adata, color=['FUS', 'moranI'])
```

### Example 3: Cell-Cell Interaction Analysis

```python
# Compute spatial graph
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)

# Calculate connectivity
sq.gr.centrality_scores(adata, cluster_key='cluster')

# Co-occurrence analysis
sq.gr.co_occurrence(adata, cluster_key='cluster')
sq.pl.co_occurrence(adata, cluster_key='cluster', figsize=(10, 5))

# Neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key='cluster')
sq.pl.nhood_enrichment(adata, cluster_key='cluster')
```

### Example 4: Ligand-Receptor Analysis

```python
# Import ligand-receptor database
import liana as li

# Run ligand-receptor inference
li.mt.rank_aggregate(adata, groupby='cluster', expr_prop=0.1)

# Visualize interactions
li.pl.dotplot(adata, colour='magnitude', size='specificity')

# Spatially restricted interactions
sq.gr.ligrec(adata, n_perms=1000, cluster_key='cluster')
sq.pl.ligrec(adata, source_groups=['Cluster_A'], target_groups=['Cluster_B'])
```

### Example 5: Integrate with Single-Cell Data

```python
# Load single-cell reference
sc_adata = sc.datasets.pbmc3k_processed()

# Transfer labels to spatial data
sq.tl.transfer_labels(adata, sc_adata, label_key='cell_type')

# Visualize deconvolution
sq.pl.spatial_scatter(adata, color=['cell_type', 'cell_type_probability'])
```

### Example 6: Image Feature Extraction

```python
# Calculate image features
sq.im.calculate_image_features(
    adata,
    img,
    features='summary',  # summary, histogram, texture
    key_added='image_features',
    n_jobs=4
)

# Use features for analysis
sc.pp.pca(adata, use_highly_variable=False)
sc.pl.pca(adata, color='image_features_summary_mean')
```

## Integration with Other Skills

- **scanpy**: Single-cell analysis backbone
- **anndata**: Data structure for spatial data
- **scvi-tools**: Deep learning for spatial deconvolution
- **cellxgene-census**: Reference atlases
- **matplotlib/scientific-visualization**: Publication figures
- **pandas**: Data manipulation

## Best Practices

### Data Preprocessing
1. **Quality control**: Filter low-quality spots
2. **Normalization**: Use scanpy's normalization functions
3. **Feature selection**: Select highly variable genes
4. **Dimensionality reduction**: PCA before spatial analysis

### Spatial Analysis Workflow
1. **Visual inspection**: Always plot raw data first
2. **Spatial graphs**: Define neighborhood structure
3. **Statistical testing**: Use permutation tests for significance
4. **Multiple corrections**: Adjust p-values for multiple testing

### Visualization Tips
1. **Color schemes**: Use perceptually uniform colormaps
2. **Scale bars**: Always include for spatial context
3. **Annotations**: Label tissue regions and landmarks
4. **Resolution**: Export high-resolution figures (300+ DPI)

## Advanced Topics

### Custom Spatial Graphs

```python
# Create custom Delaunay triangulation
sq.gr.spatial_neighbors(
    adata,
    coord_type='generic',
    spatial_key='spatial',
    delaunay=True
)

# Radius-based neighbors
sq.gr.spatial_neighbors(
    adata,
    coord_type='generic',
    radius=50  # microns
)
```

### Spatial Trajectory Inference

```python
# Pseudotime analysis in spatial context
import scvelo as scv

# Compute velocity
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# Spatial projection
scv.pl.velocity_embedding_grid(adata, basis='spatial', color='cluster')
```

## Troubleshooting

**Issue**: Image alignment problems
- **Solution**: Check coordinate system and scale factors

**Issue**: Memory errors with large datasets
- **Solution**: Use `inplace=False` or process in chunks

**Issue**: Slow computation
- **Solution**: Use `n_jobs` parameter for parallelization

## Resources

- **Documentation**: https://squidpy.readthedocs.io/
- **Tutorials**: https://squidpy.readthedocs.io/en/stable/tutorials.html
- **GitHub**: https://github.com/scverse/squidpy
- **Citation**: Palla et al. (2022) Nature Methods

## Summary

Squidpy provides a comprehensive toolkit for spatial transcriptomics analysis, enabling researchers to study gene expression in its spatial context. Its integration with the scverse ecosystem makes it a powerful tool for understanding tissue architecture and cellular interactions.
