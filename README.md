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

## 📖 Usage Strategies

### Strategy 1: Automatic Skill Detection

**Principle**: OpenCode automatically detects conversation content and loads matching skills.

**Best Practice**: Be specific about your domain and tools

```
❌ Less effective: "Help me analyze this data"
✅ More effective: "I have an RNA-seq count matrix. I need QC, 
    normalization, clustering, and differential expression analysis 
    to identify marker genes for each cluster"
    
→ Auto-loads: scanpy, pydeseq2, scientific-visualization
```

### Strategy 2: Superpowers for Complex Workflows

| Superpower | When to Use | Example Scenario |
|------------|-------------|------------------|
| `brainstorming` | Need creativity/exploration | "Brainstorm hypotheses for my experiment" |
| `writing-plans` | Large task breakdown | "Create a 3-day analysis plan" |
| `systematic-debugging` | Code issues | "Debug this error systematically" |
| `dispatching-parallel-agents` | Parallel tasks | "Analyze these 5 datasets simultaneously" |
| `subagent-driven-development` | Complex development | "Develop a literature crawler tool" |
| `requesting-code-review` | Code review | "Review this code quality" |
| `test-driven-development` | Need tests | "Write TDD approach for this project" |
| `verification-before-completion` | Quality check | "Verify all requirements before finishing" |

### Strategy 3: Combining Skills + Superpowers (Most Powerful)

**Standard Workflow for Complex Research Tasks**:

```
1. brainstorming → Explore ideas and hypotheses
2. writing-plans → Create detailed execution plan
3. 146 scientific skills → Execute domain-specific analysis
4. verification-before-completion → Validate results
5. scientific-writing → Document findings
```

### Strategy 4: Explicit Skill Requests

Ask OpenCode to use specific skills explicitly:

```
"Use the hypothesis-generation skill to build hypotheses for this experiment"

"Use literature-review skill to search recent 5-year literature on phase separation"

"Use scientific-critical-thinking to evaluate the credibility of this conclusion"
```

### Strategy 5: Skill Chaining

Chain multiple skills for complex workflows:

```
"First use gget to quickly query P53 gene info,
 then use string-database to find its interaction proteins,
 finally use cytoscape to visualize the network"
```

### Strategy 6: Let AI Guide You (Most Important!)

**In real-world problem solving, communicate with OpenCode and let AI tell you how to leverage skills and superpowers!**

```
"I want to analyze this metabolomics dataset but don't know where to start. 
 What skills do you have that can help? What's the best workflow?"

"I need to debug this complex bioinformatics pipeline. 
 Which superpowers should I use?"

"I'm planning a multi-omics study. Can you suggest a complete 
 analysis strategy using available skills?"

"I have this single-cell RNA-seq data with batch effects across 3 conditions.
 Walk me through the best approach step by step."
```

OpenCode will analyze your specific problem and recommend:
- Which skills to use
- The optimal execution order
- Which superpowers can enhance the workflow
- Step-by-step implementation plan



## 🔬 Real-World Usage Examples

### Example 1: Writing a Review Paper from Scratch

```
User: "I want to write a review about 'phase-separation proteins in RNA splicing'"

OpenCode automatically invokes:
├── literature-review → Systematic literature search
├── openalex-database → Retrieve latest papers from OpenAlex
├── string-database → Find protein-protein interactions
├── scientific-brainstorming → Build paper framework
├── scientific-writing → Write each section
├── scientific-visualization → Create publication-ready figures
└── citation-management → Manage references

Superpowers enhance:
├── writing-plans → Create writing schedule
├── executing-plans → Execute according to plan
└── verification-before-completion → Final quality check
```

### Example 2: Bioinformatics Analysis Pipeline

```
User: "I have single-cell RNA-seq data and want to:
1. QC and filtering
2. Clustering analysis
3. Find marker genes
4. Pathway enrichment"

Automatic skill chain:
scanpy (QC/clustering)
  → decoupler (pathway activity inference)
  → reactome-database (pathway database query)
  → scientific-visualization (publication-ready figures)
  → umap-learn (dimensionality reduction visualization)
```

### Example 3: Drug Discovery Project

```
User: "Design small molecule inhibitors targeting MDM2-p53 interaction"

Skill combination:
├── rdkit/datamol → Molecular manipulation
├── chembl-database → Find existing inhibitors
├── zinc-database → Virtual screening library
├── diffdock → Molecular docking
├── deepchem → Property prediction
├── medchem → Medicinal chemistry filtering
└── scientific-schematics → Create mechanism diagrams
```

## 💡 Pro Tips

### 1. Use Domain Keywords

Include professional terminology in questions to trigger relevant skills:

```
✅ "Use scanpy for batch effect correction"
✅ "Analyze crystal structure with pymatgen"
✅ "Flux balance analysis with cobrapy"
✅ "Batch integration using scvi-tools"
```

### 2. Let AI Guide You

**Most Important Strategy**: In actual problem-solving, communicate with OpenCode and let AI tell you how to leverage skills and superpowers!

```
"I want to analyze this metabolomics dataset but don't know where to start. 
 What skills do you have that can help? What's the best workflow?"

"I need to debug this complex bioinformatics pipeline. 
 Which superpowers should I use?"

"I'm planning a multi-omics study. Can you suggest a complete 
 analysis strategy using available skills?"
```

OpenCode will analyze your specific problem and recommend:
- Which skills to use
- The optimal execution order
- Which superpowers can enhance the workflow
- Step-by-step implementation plan

### 3. Daily Research Workflow

```
Morning:
├── "Check today's literature updates" → biorxiv-database
└── "Summarize these 3 relevant papers" → markitdown + scientific-writing

Lab work:
├── "Design this PCR experiment" → hypothesis-generation
└── "Record experimental steps" → protocolsio-integration

Analysis:
├── "Analyze the generated data" → [relevant analysis skills]
└── "Create result figures" → scientific-visualization

Writing:
├── "Summarize today's progress" → scientific-writing
└── "Plan tomorrow's work" → writing-plans
```

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
