# Zhui's OpenCode Skills Collection

Personal collection of 146 scientific skills + Superpowers plugin for OpenCode AI assistant.

## 📊 What's Included

### 🎯 Scientific Skills (146 total)

#### 🧬 Bioinformatics & Genomics (16+)
biopython, scanpy, scvi-tools, anndata, arboreto, pysam, deeptools, gget, scikit-bio

#### 🧪 Cheminformatics & Drug Discovery (11+)
rdkit, datamol, deepchem, torchdrug, diffdock, molfeat, medchem, pytdc

#### 🔬 Proteomics & Mass Spectrometry (2)
pyopenms, matchms

#### 🏥 Clinical Research & Precision Medicine (12+)
clinicaltrials-database, clinvar-database, pyhealth, clinical-decision-support

#### 🤖 Machine Learning & AI (15+)
pytorch-lightning, transformers, scikit-learn, stable-baselines3, pufferlib, torch-geometric

#### 📊 Data Analysis & Visualization (14+)
matplotlib, seaborn, plotly, networkx, umap-learn

#### 🔮 Materials Science & Chemistry (7)
pymatgen, cobrapy, cirq, pennylane, qiskit

#### 📚 Scientific Databases (28+)
pubmed-database, chembl-database, uniprot-database, openalex-database, string-database, kegg-database

#### 🔬 Multi-omics & Systems Biology (5+)
decoupler, denario, hypogenic, cytoscape, squidpy, string-app

#### 🎓 Research Methodology & Planning (8+)
literature-review, scientific-writing, hypothesis-generation, research-grants

#### 📚 Scientific Communication (20+)
latex-posters, scientific-slides, citation-management, infographics

#### 🔧 Additional Tools
langchain

### ⚡ Superpowers Plugin

Advanced capabilities from [obra/superpowers](https://github.com/obra/superpowers):
- **brainstorming** - Creative ideation and exploration
- **dispatching-parallel-agents** - Execute multiple subagents concurrently
- **executing-plans** - Systematic plan execution
- **finishing-a-development-branch** - Branch completion workflow
- **receiving-code-review** - Handle code review feedback
- **requesting-code-review** - Request and manage code reviews
- **subagent-driven-development** - AI-driven development workflow
- **systematic-debugging** - Methodical debugging approach
- **test-driven-development** - TDD workflow
- **using-git-worktrees** - Advanced git worktree management
- **using-superpowers** - Guide to using superpowers
- **verification-before-completion** - Pre-completion verification
- **writing-plans** - Plan writing methodology
- **writing-skills** - Skill creation guide

## 🚀 Quick Setup (New Device)

### One-Command Installation

```bash
bash -c "$(curl -fsSL https://raw.githubusercontent.com/1508324011/Zhui-s-opencode-skills/main/install.sh)"
```

Or manually:

```bash
# 1. Clone this repository
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git ~/zhui-opencode-skills

# 2. Copy skills to OpenCode directory
mkdir -p ~/.config/opencode/skills
cp -r ~/zhui-opencode-skills/skills/* ~/.config/opencode/skills/

# 3. Install Superpowers plugin
git clone https://github.com/obra/superpowers.git ~/.config/opencode/superpowers

# 4. Create plugin symlink
mkdir -p ~/.config/opencode/plugins
ln -s ~/.config/opencode/superpowers/.opencode/plugins/superpowers.js ~/.config/opencode/plugins/superpowers.js

# 5. Create skills symlink for superpowers
ln -s ~/.config/opencode/superpowers/skills ~/.config/opencode/skills/superpowers
```

### Verification

Restart OpenCode and verify installation:

```
do you have superpowers?
```

The assistant should respond confirming it has superpowers capabilities.

## 📁 Repository Structure

```
Zhui-s-opencode-skills/
├── README.md              # This file
├── install.sh             # Automated installation script
├── skills/                # 146 scientific skills
│   ├── biopython/
│   ├── scanpy/
│   ├── rdkit/
│   └── ... (143 more)
└── superpowers/           # Reference to obra/superpowers
    └── (cloned separately)
```

## 🔄 Keeping Updated

### Update Skills
```bash
cd ~/zhui-opencode-skills
git pull origin main
# Re-copy if needed
cp -r skills/* ~/.config/opencode/skills/
```

### Update Superpowers
```bash
cd ~/.config/opencode/superpowers
git pull
```

## 📖 Usage Examples

After installation, OpenCode will automatically detect and use these skills:

```
# Scientific analysis
"Analyze this RNA-seq data with scanpy"
"Search for RBP interactions in STRING database"

# Research workflow
"Create a literature review about phase separation"
"Help me write a hypothesis for my experiment"

# Development with superpowers
"Use systematic debugging to find the bug in this code"
"Create a development plan for this feature"
"Request a code review for these changes"
```

## 🛠️ Troubleshooting

### Skills not detected
```bash
# Check skill files exist
ls ~/.config/opencode/skills/

# Verify symlinks (if using symlink method)
ls -la ~/.config/opencode/skills/
```

### Superpowers not loading
```bash
# Check plugin symlink
ls -la ~/.config/opencode/plugins/superpowers.js

# Verify source exists
ls ~/.config/opencode/superpowers/.opencode/plugins/superpowers.js
```

### Permission issues
```bash
# Fix permissions
chmod -R 755 ~/.config/opencode/skills/
chmod -R 755 ~/.config/opencode/plugins/
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

## 🏆 Credits

- **Scientific Skills**: Originally from [K-Dense-AI/claude-scientific-skills](https://github.com/K-Dense-AI/claude-scientific-skills)
- **Superpowers Plugin**: [obra/superpowers](https://github.com/obra/superpowers)

## 📄 License

Each skill has its own license specified in its SKILL.md. Repository structure follows MIT license where applicable.

## 📊 Statistics

- **Total Skills**: 146 scientific + 12 superpowers
- **Categories**: 15+ scientific domains
- **Last Updated**: 2025-02-16
- **Compatibility**: OpenCode

---

*Made with ❤️ for scientific research automation*
