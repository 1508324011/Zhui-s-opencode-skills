# Zhui's OpenCode Skills Collection

[![Skills Count](https://img.shields.io/badge/Skills-146-blue.svg)](./)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)

A comprehensive collection of **146 scientific and technical skills** for OpenCode, organized by domain.

## 📊 Quick Stats

- **Total Skills**: 146
- **Categories**: 12 major domains
- **Source**: K-Dense Inc. Scientific Skills
- **Last Updated**: 2026-02-16

## 📁 What's in This Repository

This repository serves as an **index and installer** for the complete OpenCode skills collection:

1. **[SKILLS.md](./SKILLS.md)** - Complete catalog of all 146 skills with descriptions
2. **[install.sh](./install.sh)** - Automated setup script for new devices
3. **README.md** - This file (overview and usage guide)

## 🚀 Quick Start

### For New OpenCode Installations

```bash
# Clone this repository
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills

# Run the installer
./install.sh
```

### Configure Custom Provider (e.g., 百炼 Coding Plan)

1. **Copy the example config:**
   ```bash
   cp config/opencode.json.example ~/.config/opencode/opencode.json
   ```

2. **Edit the config and replace `YOUR_API_KEY`:**
   ```bash
   nano ~/.config/opencode/opencode.json
   # or
   vim ~/.config/opencode/opencode.json
   ```

3. **Restart OpenCode** to apply the configuration.

### Superpowers Plugin Setup

The install script will automatically set up Superpowers. Manual setup:

```bash
# Clone superpowers repository
git clone https://github.com/obra/superpowers.git ~/.config/opencode/superpowers

# Create symlinks
mkdir -p ~/.config/opencode/plugins ~/.config/opencode/skills
ln -sf ~/.config/opencode/superpowers/.opencode/plugins/superpowers.js ~/.config/opencode/plugins/superpowers.js
ln -sf ~/.config/opencode/superpowers/skills ~/.config/opencode/skills/superpowers

# Restart OpenCode
```

**Note:** The plugin file should be named `superpowers.js` (with 's'), not `superpower.js`.

## 📚 Skill Categories

### 1. 🔬 Bioinformatics & Genomics (32 skills)
Single-cell analysis, genomic databases, sequence analysis, pathway analysis

**Key Skills**: scanpy, scvi-tools, cellxgene-census, squidpy, biopython, geniml

**Databases**: PubMed, Ensembl, GEO, UniProt, STRING, Reactome, KEGG, GWAS Catalog, ClinVar, COSMIC

### 2. 💊 Drug Discovery & Cheminformatics (15 skills)
Molecular modeling, docking, protein structures, bioactivity data

**Key Skills**: rdkit, deepchem, diffdock, alphafold-database, chembl-database, pubchem-database

### 3. 🤖 Machine Learning & AI (18 skills)
Classical ML, deep learning, GNNs, reinforcement learning, transformers

**Key Skills**: scikit-learn, pytorch-lightning, transformers, torch-geometric, shap, umap-learn

### 4. 📊 Statistics & Data Analysis (8 skills)
Statistical modeling, Bayesian inference, survival analysis, optimization

**Key Skills**: statsmodels, pymc, scikit-survival, pymoo, sympy

### 5. 📈 Visualization & Graphics (7 skills)
Plotting, scientific figures, AI-generated images and diagrams

**Key Skills**: matplotlib, seaborn, plotly, scientific-visualization, generate-image

### 6. 📝 Research & Scientific Writing (11 skills)
Manuscript writing, literature review, grants, peer review

**Key Skills**: scientific-writing, literature-review, research-grants, peer-review, citation-management

### 7. 🏥 Healthcare & Clinical (11 skills)
Clinical data analysis, medical imaging, trial databases

**Key Skills**: pyhealth, clinical-reports, pydicom, imaging-data-commons, clinicaltrials-database

### 8. 🧪 Laboratory & Automation (6 skills)
Lab equipment integration, ELN, protocol management

**Key Skills**: opentrons-integration, pylabrobot, benchling-integration, protocolsio-integration

### 9. ☁️ Cloud & Data Platforms (6 skills)
Cloud genomics, data commons, microscopy platforms

**Key Skills**: dnanexus-integration, latchbio-integration, datacommons-client, omero-integration

### 10. 🖥️ Scientific Computing & Simulation (8 skills)
MATLAB, CFD, quantum computing, materials science

**Key Skills**: matlab, fluidsim, qiskit, cirq, pennylane, pymatgen

### 11. 🕸️ Network & Graph Analysis (3 skills)
Network analysis, biological networks, geospatial data

**Key Skills**: networkx, cytoscape, geopandas

### 12. 🛠️ Specialized Tools & Utilities (12 skills)
Astronomy, critical thinking, market research, metabolomics

**Key Skills**: astropy, scientific-critical-thinking, market-research-reports, matchms

## 🔧 Understanding the Structure

### What's the Difference Between K-Dense Skills and Superpowers?

| Feature | K-Dense Skills (146) | Superpowers Plugin |
|---------|---------------------|-------------------|
| **Content** | Scientific/technical skills | Plugin framework + extra skills |
| **Source** | K-Dense Inc. | obra/superpowers |
| **Installation** | Usually pre-installed with OpenCode | Manual git clone |
| **Purpose** | Domain-specific scientific tasks | Enhanced OpenCode functionality |

### Your Current Setup

```
~/.config/opencode/skills/
├── adaptyv/              ← 146 K-Dense skills
├── scanpy/
├── scvi-tools/
│   └── ... (146 total)
└── superpowers -> ~/.config/opencode/superpowers/skills  (optional symlink)
```

### To Make New Devices Identical

**Option 1: Pack the entire skills directory**
```bash
# Create archive of all your skills
tar -czf opencode-skills-backup.tar.gz ~/.config/opencode/skills/

# On new machine, extract
tar -xzf opencode-skills-backup.tar.gz -C ~/
```

**Option 2: Use this repository + install script**
```bash
# The install.sh in this repo will set up everything
./install.sh
```

**Option 3: Include skill files in this repository**
```bash
# Copy actual skill files to repo
cp -r ~/.config/opencode/skills/* ./skills/
git add skills/
git commit -m "Add all 146 skill files"
git push
```

## 📖 Complete Skill List

See [SKILLS.md](./SKILLS.md) for the complete catalog.

## ⚙️ Configuration Guide

### Custom Provider Setup

This repository includes an example configuration for custom providers (like 百炼 Coding Plan):

**File:** `config/opencode.json.example`

**To use:**
1. Copy the example to your OpenCode config directory
2. Replace `YOUR_API_KEY` with your actual API key
3. Restart OpenCode

**Supported Models in Example:**
- `qwen3.5-plus` - General purpose model with thinking
- `qwen3-max-2026-01-23` - Latest Qwen3 Max with thinking
- `qwen3-coder-plus` - Optimized for coding tasks
- `qwen3-coder-next` - Next-generation coder model

### Common Configuration Issues

**1. Plugin not loading:**
- Check that the file is named `superpowers.js` (with 's')
- Verify symlinks: `ls -la ~/.config/opencode/plugins/`
- Check OpenCode version: `opencode --version`

**2. Provider configuration not working:**
- Ensure `opencode.json` is in `~/.config/opencode/`
- Validate JSON syntax
- Check that the provider npm package is installed

## 📄 License

MIT License

## 🙏 Acknowledgments

- **K-Dense Inc.** - Creators of the 146 scientific skills
- **OpenCode.ai** - AI-assisted scientific computing platform
- **obra/superpowers** - Superpowers plugin framework (optional)
