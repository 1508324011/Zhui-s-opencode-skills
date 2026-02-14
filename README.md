# Zhui's OpenCode Skills Collection

Personal collection of 146 scientific skills for OpenCode AI assistant, migrated from Claude Scientific Skills.

## 📊 Skills Overview

This repository contains **146 specialized skills** covering:

### 🧬 Bioinformatics & Genomics (16+)
- biopython, scanpy, scvi-tools, anndata, arboreto
- pysam, deeptools, gget, scikit-bio
- Single-cell analysis, genomics tools, phylogenetics

### 🧪 Cheminformatics & Drug Discovery (11+)
- rdkit, datamol, deepchem, torchdrug
- diffdock, molfeat, medchem, pytdc
- Molecular property prediction, virtual screening

### 🔬 Proteomics & Mass Spectrometry (2)
- pyopenms, matchms

### 🏥 Clinical Research & Precision Medicine (12+)
- clinicaltrials-database, clinvar-database
- pyhealth, clinical-decision-support
- Healthcare AI, clinical documentation

### 🤖 Machine Learning & AI (15+)
- pytorch-lightning, transformers, scikit-learn
- stable-baselines3, pufferlib, torch-geometric
- Deep learning, reinforcement learning

### 📊 Data Analysis & Visualization (14+)
- matplotlib, seaborn, plotly
- networkx, umap-learn
- Scientific visualization

### 🔮 Materials Science & Chemistry (7)
- pymatgen, cobrapy, cirq, pennylane, qiskit

### 📚 Scientific Databases (28+)
- pubmed-database, chembl-database, uniprot-database
- openalex-database, string-database, kegg-database
- Reactome, Open Targets, COSMIC, ClinPGx

### 🔬 Multi-omics & Systems Biology (5+)
- decoupler (NEW), denario, hypogenic

### 🎓 Research Methodology & Planning (8+)
- literature-review, scientific-writing
- hypothesis-generation, research-grants

### 📚 Scientific Communication (20+)
- latex-posters, scientific-slides
- citation-management, infographics

### 🔧 Additional Tools
- langchain (NEW), cytoscape (NEW), squidpy (NEW)
- string-app (NEW)

## 🚀 Installation

### Method 1: Direct Copy
```bash
# Clone this repository
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git

# Copy to OpenCode skills directory
cp -r Zhui-s-opencode-skills/* ~/.config/opencode/skills/
```

### Method 2: Symlink (Recommended for updates)
```bash
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git ~/zhui-opencode-skills
ln -s ~/zhui-opencode-skills/* ~/.config/opencode/skills/
```

## 📖 Usage

After installation, OpenCode will automatically detect and use these skills when relevant.

Example usage in OpenCode:
```
# The assistant will automatically load appropriate skills
"Analyze this RNA-seq data with scanpy"
"Search for RBP interactions in STRING database"
"Create a literature review about phase separation"
```

## 🔄 Updates

To update skills from this repository:
```bash
cd ~/zhui-opencode-skills
git pull origin main
```

## 📝 Skill Format

Each skill follows the OpenCode skill specification:
- `SKILL.md` with YAML frontmatter
- `name`, `description`, `license`, `compatibility` fields
- Comprehensive documentation and examples

Example structure:
```
skill-name/
├── SKILL.md          # Main documentation
└── references/       # Optional reference materials
```

## 🏆 Original Source

These skills were originally from [K-Dense-AI/claude-scientific-skills](https://github.com/K-Dense-AI/claude-scientific-skills) and adapted for OpenCode compatibility.

## 📄 License

Each skill has its own license specified in the `license` field of its SKILL.md file. The repository structure follows the original MIT license where applicable.

## 🤝 Contributing

This is a personal collection. Feel free to fork and customize for your own use!

## 📊 Statistics

- **Total Skills**: 146
- **Categories**: 15+
- **Last Updated**: 2025-02-14
- **Compatibility**: OpenCode

---

*Made with ❤️ for scientific research automation*
