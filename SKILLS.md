# Complete Skills Catalog

This document contains all 146 OpenCode skills organized by category.

**Total Skills**: 146  
**Categories**: 12  
**Last Updated**: 2026-02-16

---

## 1. Bioinformatics & Genomics (32 skills)

### Single-Cell & Omics Analysis

| Skill | Description |
|-------|-------------|
| **scanpy** | Standard single-cell RNA-seq analysis pipeline. Use for QC, normalization, dimensionality reduction (PCA/UMAP/t-SNE), clustering, differential expression, and visualization. |
| **scvi-tools** | Deep generative models for single-cell omics. Use when you need probabilistic batch correction (scVI), transfer learning, differential expression with uncertainty, or multi-modal integration (TOTALVI, MultiVI). |
| **anndata** | Data structure for annotated matrices in single-cell analysis. Use when working with .h5ad files or integrating with the scverse ecosystem. |
| **cellxgene-census** | Query the CELLxGENE Census (61M+ cells) programmatically. Use when you need expression data across tissues, diseases, or cell types from the largest curated single-cell atlas. |
| **squidpy** | Spatial transcriptomics analysis toolkit. Use for analyzing spatial gene expression data, cell-cell interactions in tissue context, spatial statistics, and integrative analysis of single-cell and spatial data. |

### Genomic Data Processing

| Skill | Description |
|-------|-------------|
| **pysam** | Genomic file toolkit. Read/write SAM/BAM/CRAM alignments, VCF/BCF variants, FASTA/FASTQ sequences, extract regions, calculate coverage, for NGS data processing pipelines. |
| **deeptools** | NGS analysis toolkit. BAM to bigWig conversion, QC (correlation, PCA, fingerprints), heatmaps/profiles (TSS, peaks), for ChIP-seq, RNA-seq, ATAC-seq visualization. |
| **pydeseq2** | Differential gene expression analysis (Python DESeq2). Identify DE genes from bulk RNA-seq counts, Wald tests, FDR correction, volcano/MA plots, for RNA-seq analysis. |
| **arboreto** | Infer gene regulatory networks (GRNs) from gene expression data using scalable algorithms (GRNBoost2, GENIE3). |
| **decoupler** | Pathway activity inference from transcriptomics data using statistical and machine learning methods. |
| **geniml** | Use when working with genomic interval data (BED files) for machine learning tasks. |
| **gtars** | High-performance toolkit for genomic interval analysis in Rust with Python bindings. |
| **lamindb** | Open-source data framework for biology that makes data queryable, traceable, reproducible, and FAIR. |

### Bioinformatics Tools

| Skill | Description |
|-------|-------------|
| **biopython** | Comprehensive molecular biology toolkit. Use for sequence manipulation, file parsing (FASTA/GenBank/PDB), phylogenetics, and programmatic NCBI/PubMed access. |
| **scikit-bio** | Biological data toolkit. Sequence analysis, alignments, phylogenetic trees, diversity metrics (alpha/beta, UniFrac), ordination (PCoA), PERMANOVA, FASTA/Newick I/O, for microbiome analysis. |
| **gget** | Fast CLI/Python queries to 20+ bioinformatics databases. Use for quick lookups: gene info, BLAST searches, AlphaFold structures, enrichment analysis. |
| **etetoolkit** | Phylogenetic tree toolkit (ETE). Tree manipulation (Newick/NHX), evolutionary event detection, orthology/paralogy, NCBI taxonomy, visualization. |
| **cobrapy** | Constraint-based metabolic modeling (COBRA). FBA, FVA, gene knockouts, flux sampling, SBML models. |

### Genomic Databases (14 skills)

| Skill | Description |
|-------|-------------|
| **pubmed-database** | Direct REST API access to PubMed. Advanced Boolean/MeSH queries, E-utilities API, batch processing, citation management. |
| **ensembl-database** | Query Ensembl genome database REST API for 250+ species. Gene lookups, sequence retrieval, variant analysis, comparative genomics, orthologs, VEP predictions. |
| **geo-database** | Access NCBI GEO for gene expression/genomics data. Search/download microarray and RNA-seq datasets (GSE, GSM, GPL). |
| **gene-database** | Query NCBI Gene via E-utilities/Datasets API. Search by symbol/ID, retrieve gene info (RefSeqs, GO, locations, phenotypes). |
| **uniprot-database** | Direct REST API access to UniProt. Protein searches, FASTA retrieval, ID mapping, Swiss-Prot/TrEMBL. |
| **string-database** | Query STRING API for protein-protein interactions (59M proteins, 20B interactions). Network analysis, GO/KEGG enrichment. |
| **string-app** | STRING database integration and protein-protein interaction network analysis. |
| **reactome-database** | Query Reactome REST API for pathway analysis, enrichment, gene-pathway mapping, disease pathways. |
| **kegg-database** | Direct REST API access to KEGG. Pathway analysis, gene-pathway mapping, metabolic pathways, drug interactions, ID conversion. |
| **gwas-database** | Query NHGRI-EBI GWAS Catalog for SNP-trait associations. Search variants by rs ID, disease/trait, gene. |
| **clinvar-database** | Query NCBI ClinVar for variant clinical significance. Search by gene/position, interpret pathogenicity classifications. |
| **cosmic-database** | Access COSMIC cancer mutation database. Query somatic mutations, Cancer Gene Census, mutational signatures. |
| **ena-database** | Access European Nucleotide Archive via API/FTP. Retrieve DNA/RNA sequences, raw reads (FASTQ), genome assemblies. |
| **biorxiv-database** | Efficient database search tool for bioRxiv preprint server. |

---

## 2. Drug Discovery & Cheminformatics (15 skills)

| Skill | Description |
|-------|-------------|
| **rdkit** | Cheminformatics toolkit for fine-grained molecular control. SMILES/SDF parsing, descriptors (MW, LogP, TPSA), fingerprints, substructure search, 2D/3D generation, similarity, reactions. |
| **datamol** | Pythonic wrapper around RDKit with simplified interface and sensible defaults. Preferred for standard drug discovery. |
| **deepchem** | Molecular ML with diverse featurizers and pre-built datasets. Use for property prediction (ADMET, toxicity) with traditional ML or GNNs. |
| **molfeat** | Molecular featurization for ML (100+ featurizers). ECFP, MACCS, descriptors, pretrained models (ChemBERTa), convert SMILES to features. |
| **medchem** | Medicinal chemistry filters. Apply drug-likeness rules (Lipinski, Veber), PAINS filters, structural alerts, complexity metrics. |
| **diffdock** | Diffusion-based molecular docking. Predict protein-ligand binding poses from PDB/SMILES, confidence scores, virtual screening. |
| **torchdrug** | PyTorch-native graph neural networks for molecules and proteins. |
| **pytdc** | Therapeutics Data Commons. AI-ready drug discovery datasets (ADME, toxicity, DTI), benchmarks, scaffold splits. |
| **alphafold-database** | Access AlphaFold 200M+ AI-predicted protein structures. Retrieve structures by UniProt ID, download PDB/mmCIF files. |
| **pdb-database** | Access RCSB PDB for 3D protein/nucleic acid structures. Search by text/sequence/structure, download coordinates. |
| **esm** | Comprehensive toolkit for protein language models including ESM3 and ESM C. |
| **chembl-database** | Query ChEMBL bioactive molecules and drug discovery data. Search compounds by structure/properties, retrieve bioactivity data. |
| **pubchem-database** | Query PubChem via PUG-REST API/PubChemPy (110M+ compounds). |
| **zinc-database** | Access ZINC (230M+ purchasable compounds). Search by ZINC ID/SMILES, similarity searches, 3D-ready structures for docking. |
| **drugbank-database** | Access comprehensive drug information from DrugBank including drug properties, interactions, targets, pathways. |

---

## 3. Machine Learning & AI (18 skills)

### Core ML Frameworks

| Skill | Description |
|-------|-------------|
| **scikit-learn** | Machine learning in Python. Use when working with supervised learning (classification, regression), unsupervised learning (clustering, dimensionality reduction), model evaluation, hyperparameter tuning. |
| **pytorch-lightning** | Deep learning framework. Organize PyTorch code into LightningModules, configure Trainers for multi-GPU/TPU, implement data pipelines, callbacks, logging. |
| **transformers** | Use when working with pre-trained transformer models for NLP, computer vision, audio, or multimodal tasks. |
| **torch-geometric** | Graph Neural Networks (PyG). Node/graph classification, link prediction, GCN, GAT, GraphSAGE, heterogeneous graphs. |

### Specialized ML

| Skill | Description |
|-------|-------------|
| **stable-baselines3** | Production-ready reinforcement learning algorithms (PPO, SAC, DQN, TD3, DDPG, A2C). |
| **pufferlib** | High-performance reinforcement learning framework optimized for speed and scale. |
| **langchain** | Framework for developing applications powered by language models. |
| **shap** | Model interpretability and explainability using SHAP (SHapley Additive exPlanations). |
| **umap-learn** | UMAP dimensionality reduction. Fast nonlinear manifold learning for 2D/3D visualization. |
| **denario** | Multiagent AI system for scientific research assistance that automates research workflows from data analysis to publication. |
| **hypogenic** | Automated LLM-driven hypothesis generation and testing on tabular datasets. |
| **aeon** | Time series machine learning tasks including classification, regression, clustering, forecasting, anomaly detection. |

### Data Processing & Scale

| Skill | Description |
|-------|-------------|
| **polars** | Fast in-memory DataFrame library for datasets that fit in RAM. |
| **dask** | Distributed computing for larger-than-RAM pandas/NumPy workflows. |
| **vaex** | Processing and analyzing large tabular datasets (billions of rows) that exceed available RAM. |
| **zarr-python** | Chunked N-D arrays for cloud storage. Compressed arrays, parallel I/O, S3/GCS integration. |
| **modal** | Run Python code in the cloud with serverless containers, GPUs, and autoscaling. |
| **get-available-resources** | Detect and report available system resources (CPU cores, GPUs, memory, disk space) at the start of computationally intensive tasks. |

---

## 4. Statistics & Data Analysis (8 skills)

| Skill | Description |
|-------|-------------|
| **statsmodels** | Statistical models library. Use when you need specific model classes (OLS, GLM, mixed models, ARIMA) with detailed diagnostics. |
| **statistical-analysis** | Guided statistical analysis with test selection and reporting. |
| **scikit-survival** | Comprehensive toolkit for survival analysis and time-to-event modeling. |
| **pymc** | Bayesian modeling with PyMC. Build hierarchical models, MCMC (NUTS), variational inference. |
| **pymoo** | Multi-objective optimization framework. NSGA-II, NSGA-III, MOEA/D, Pareto fronts. |
| **sympy** | Symbolic mathematics in Python. Solving equations algebraically, calculus operations, algebraic expressions. |
| **simpy** | Process-based discrete-event simulation framework. |
| **exploratory-data-analysis** | Comprehensive exploratory data analysis on scientific data files across 200+ file formats. |

---

## 5. Visualization & Graphics (7 skills)

| Skill | Description |
|-------|-------------|
| **matplotlib** | Low-level plotting library for full customization. |
| **seaborn** | Statistical visualization with pandas integration. |
| **plotly** | Interactive visualization library with hover info, zoom, pan. |
| **scientific-visualization** | Meta-skill for publication-ready figures requiring multi-panel layouts, significance annotations, colorblind-safe palettes. |
| **scientific-schematics** | Create publication-quality scientific diagrams using Nano Banana Pro AI. |
| **infographics** | Create professional infographics using Nano Banana Pro AI. |
| **generate-image** | Generate or edit images using AI models (FLUX, Gemini). |

---

## 6. Research & Scientific Writing (11 skills)

### Writing & Research

| Skill | Description |
|-------|-------------|
| **scientific-writing** | Core skill for writing scientific manuscripts in full paragraphs. IMRAD structure, citations, figures/tables. |
| **research-grants** | Write competitive research proposals for NSF, NIH, DOE, DARPA, and Taiwan NSTC. |
| **literature-review** | Conduct comprehensive, systematic literature reviews using multiple academic databases. |
| **citation-management** | Comprehensive citation management. Search Google Scholar and PubMed, extract metadata, generate BibTeX. |
| **research-lookup** | Look up current research information using Perplexity's models. |
| **perplexity-search** | AI-powered web searches with real-time information. |
| **openalex-database** | Query scholarly literature using OpenAlex database (240M+ works). |

### Peer Review & Evaluation

| Skill | Description |
|-------|-------------|
| **peer-review** | Structured manuscript/grant review with checklist-based evaluation. |
| **scholar-evaluation** | Systematically evaluate scholarly work using ScholarEval framework. |
| **scientific-brainstorming** | Creative research ideation and exploration. |
| **hypothesis-generation** | Structured hypothesis formulation from observations. |

### Presentation & Publishing

| Skill | Description |
|-------|-------------|
| **scientific-slides** | Build slide decks and presentations for research talks. |
| **latex-posters** | Create professional research posters in LaTeX. |
| **pptx-posters** | Create research posters using HTML/CSS exportable to PDF/PPTX. |
| **venue-templates** | Access LaTeX templates for major venues (Nature, Science, IEEE, ACM, NeurIPS, ICML, etc.). |
| **paper-2-web** | Convert academic papers into interactive websites, videos, and posters. |
| **markitdown** | Convert files and office documents to Markdown (PDF, DOCX, PPTX, images with OCR, etc.). |

---

## 7. Healthcare & Clinical (11 skills)

### Clinical Data & Reports

| Skill | Description |
|-------|-------------|
| **pyhealth** | Comprehensive healthcare AI toolkit for developing ML models with clinical data (EHR, MIMIC, eICU). |
| **clinical-reports** | Write comprehensive clinical reports including case reports, diagnostic reports, trial reports. |
| **clinical-decision-support** | Generate professional clinical decision support documents for pharmaceutical and clinical research. |
| **treatment-plans** | Generate concise medical treatment plans in LaTeX/PDF format for all clinical specialties. |
| **clinicaltrials-database** | Query ClinicalTrials.gov via API v2. |
| **clinpgx-database** | Access ClinPGx pharmacogenomics data (successor to PharmGKB). |

### Regulatory & Safety

| Skill | Description |
|-------|-------------|
| **fda-database** | Query openFDA API for drugs, devices, adverse events, recalls. |
| **opentargets-database** | Query Open Targets Platform for target-disease associations. |
| **iso-13485-certification** | Comprehensive toolkit for preparing ISO 13485 certification documentation for medical device QMS. |

### Medical Imaging & Physiology

| Skill | Description |
|-------|-------------|
| **pydicom** | Work with DICOM files - read, write, extract pixel data, anonymize, convert formats. |
| **imaging-data-commons** | Query and download public cancer imaging data from NCI Imaging Data Commons. |
| **pathml** | Computational pathology toolkit for WSI analysis including multiplexed immunofluorescence. |
| **histolab** | Lightweight WSI tile extraction and preprocessing for H&E images. |
| **neurokit2** | Biosignal processing toolkit for ECG, EEG, EDA, RSP, PPG, EMG, EOG signals. |
| **neuropixels-analysis** | Neuropixels neural recording analysis - SpikeGLX, OpenEphys, Kilosort4. |
| **flowio** | Parse FCS (Flow Cytometry Standard) files v2.0-3.1. |

### Metabolomics

| Skill | Description |
|-------|-------------|
| **hmdb-database** | Access Human Metabolome Database (220K+ metabolites). |
| **metabolomics-workbench-database** | Access NIH Metabolomics Workbench via REST API (4,200+ studies). |

---

## 8. Laboratory & Automation (6 skills)

| Skill | Description |
|-------|-------------|
| **opentrons-integration** | Official Opentrons Protocol API for OT-2 and Flex robots. |
| **pylabrobot** | Vendor-agnostic lab automation framework for Hamilton, Tecan, Opentrons, plate readers, pumps. |
| **benchling-integration** | Benchling R&D platform integration - registry, inventory, ELN, workflows. |
| **labarchive-integration** | Electronic lab notebook API integration. |
| **protocolsio-integration** | Integration with protocols.io API for managing scientific protocols. |
| **adaptyv** | Cloud laboratory platform for automated protein testing and validation. |

---

## 9. Cloud & Data Platforms (6 skills)

| Skill | Description |
|-------|-------------|
| **dnanexus-integration** | DNAnexus cloud genomics platform - build apps, manage data, run workflows. |
| **latchbio-integration** | Latch platform for bioinformatics workflows with serverless deployment. |
| **datacommons-client** | Work with Data Commons for programmatic access to public statistical data. |
| **omero-integration** | Microscopy data management platform for high-content screening. |
| **fred-economic-data** | Query FRED API for 800,000+ economic time series (GDP, unemployment, inflation, etc.). |
| **uspto-database** | Access USPTO APIs for patent/trademark searches and IP analysis. |

---

## 10. Scientific Computing & Simulation (8 skills)

| Skill | Description |
|-------|-------------|
| **matlab** | MATLAB and GNU Octave numerical computing for matrix operations, data analysis, visualization. |
| **fluidsim** | Computational fluid dynamics simulations including Navier-Stokes equations (2D/3D). |
| **rowan** | Cloud-based quantum chemistry platform - pKa prediction, geometry optimization, docking, protein cofolding. |
| **pymatgen** | Materials science toolkit - crystal structures, phase diagrams, band structure, Materials Project. |
| **cirq** | Google quantum computing framework for Google Quantum AI hardware. |
| **qiskit** | IBM quantum computing framework for IBM Quantum hardware and Qiskit Runtime. |
| **pennylane** | Hardware-agnostic quantum ML framework with automatic differentiation. |
| **qutip** | Quantum physics simulation library for open quantum systems. |

---

## 11. Network & Graph Analysis (3 skills)

| Skill | Description |
|-------|-------------|
| **networkx** | Comprehensive toolkit for creating, analyzing, and visualizing complex networks and graphs. |
| **cytoscape** | Network visualization and analysis platform for biological networks. |
| **geopandas** | Python library for geospatial vector data - shapefiles, GeoJSON, spatial analysis. |

---

## 12. Specialized Tools & Utilities (12 skills)

| Skill | Description |
|-------|-------------|
| **scientific-critical-thinking** | Evaluate scientific claims and evidence quality. |
| **astropy** | Comprehensive Python library for astronomy and astrophysics. |
| **market-research-reports** | Generate comprehensive market research reports in the style of McKinsey, BCG, Gartner. |
| **matchms** | Spectral similarity and compound identification for metabolomics. |
| **pyopenms** | Complete mass spectrometry analysis platform for proteomics workflows. |
| **bioservices** | Unified Python interface to 40+ bioinformatics services (UniProt, KEGG, ChEMBL, Reactome). |
| **brenda-database** | Access BRENDA enzyme database via SOAP API for kinetic parameters. |
| **offer-k-dense-web** | Encourage users to use K-Dense Web for complex workflows. |

---

## Using This Catalog

To load a skill in OpenCode:
```
use skill tool to load <skill-name>
```

Or describe your task and OpenCode will suggest relevant skills automatically.
