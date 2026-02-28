# Zhui's OpenCode Configuration

[![Oh My OpenCode](https://img.shields.io/badge/Oh%20My%20OpenCode-v3.0-blue.svg)](https://github.com/code-yeongyu/oh-my-opencode)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)

OpenCode 配置和 Skills 集合 - 用于快速迁移和设置新设备

## 📋 目录

- [快速开始](#-快速开始)
- [配置说明](#-配置说明)
- [模型分配](#-模型分配)
- [Skills](#skills)
- [Superpowers](#superpowers)
- [设备迁移](#-设备迁移)
- [安全提交规范](#-安全提交规范)

---

## 🚀 快速开始

### 新设备一键安装

```bash
# 1. 克隆此仓库
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills

# 2. 运行安装脚本
./install.sh

# 3. 复制配置文件并填入 API Key
cp config/opencode.json.example ~/.config/opencode/opencode.json
nano ~/.config/opencode/opencode.json  # 替换 YOUR_API_KEY_HERE

# 4. 重启 OpenCode
```

### 手动配置

```bash
# 1. 复制配置文件
cp config/opencode.json.example ~/.config/opencode/opencode.json
cp config/oh-my-opencode.json ~/.config/opencode/oh-my-opencode.json

# 2. 编辑 opencode.json 填入你的 API Key
nano ~/.config/opencode/opencode.json

# 3. 重启 OpenCode
```

---

## ⚙️ 配置说明

### AI 提供商配置

#### 1. 百炼平台 Coding Plan（Qwen 系列）

编辑 `~/.config/opencode/opencode.json`：

```json
{
  "provider": {
    "bailian-coding-plan": {
      "npm": "@ai-sdk/anthropic",
      "name": "Model Studio Coding Plan",
      "options": {
        "baseURL": "https://coding.dashscope.aliyuncs.com/apps/anthropic/v1",
        "apiKey": "你的 API Key"
      }
    }
  }
}
```

**可用模型：**

| 模型 | 用途 |
|------|------|
| `qwen3.5-plus` | 主力模型（带思考模式） |
| `qwen3-max-2026-01-23` | 最强推理 |
| `qwen3-coder-plus` | 代码专用 |
| `qwen3-coder-next` | 新一代代码模型 |
| `MiniMax-M2.5` | 快速任务 |
| `glm-5` | 深度分析 |
| `glm-4.7` | 通用任务 |
| `kimi-k2.5` | 推理/分析 |

---

## 🎯 模型分配

### Agents 配置

| Agent | 模型 | 用途 |
|-------|------|------|
| **Sisyphus** | Qwen3.5 Plus | 主代理/协调员 |
| **Hephaestus** | Qwen3.5 Plus | 代码生成 |
| **Oracle** | Kimi K2.5 | 调试/架构分析 |
| **Librarian** | Qwen3.5 Plus | 文档/代码搜索 |
| **Explore** | MiniMax-M2.5 | 代码探索 |
| **Multimodal-looker** | Qwen3.5 Plus | 图像分析 |
| **Prometheus** | Qwen3.5 Plus | 规划师 |
| **Metis** | GLM-5 | 规划顾问 |
| **Momus** | Qwen3.5 Plus | 审查员 |
| **Atlas** | Qwen3.5 Plus | 执行者 |

### Categories 配置

| 类别 | 模型 | 用途 |
|------|------|------|
| `visual-engineering` | Kimi K2.5 | 前端/UI 工程 |
| `ultrabrain` | GLM-5 | 深度脑力 |
| `deep` | GLM-5 | 深度分析 |
| `artistry` | Qwen3.5 Plus | 创意任务 |
| `quick` | MiniMax-M2.5 | 快速任务 |
| `unspecified-low` | MiniMax-M2.5 | 低优先级任务 |
| `unspecified-high` | Qwen3.5 Plus | 高优先级任务 |
| `writing` | Kimi K2.5 | 写作任务 |

---

## 📚 Skills

此仓库包含完整的 Skills 集合，安装后自动可用。

### 主要技能类别

1. **生物信息学** (32 skills) - scanpy, scvi-tools, cellxgene-census, squidpy
2. **药物发现** (15 skills) - rdkit, deepchem, diffdock, alphafold-database
3. **机器学习** (18 skills) - scikit-learn, transformers, pytorch-lightning, torch-geometric
4. **统计分析** (8 skills) - statsmodels, pymc, scikit-survival, pymoo
5. **可视化** (7 skills) - matplotlib, seaborn, plotly, scientific-visualization
6. **科研写作** (11 skills) - scientific-writing, literature-review, research-grants
7. **医疗健康** (11 skills) - pyhealth, clinical-reports, pydicom, imaging-data-commons
8. **实验自动化** (6 skills) - opentrons-integration, benchling-integration, pylabrobot
9. **云计算** (6 skills) - dnanexus-integration, latchbio-integration, datacommons-client
10. **科学计算** (8 skills) - matlab, fluidsim, qiskit, cirq, pennylane, pymatgen
11. **网络分析** (3 skills) - networkx, cytoscape, geopandas
12. **专用工具** (12 skills) - astropy, scientific-critical-thinking, market-research-reports

完整列表见 [SKILLS.md](./SKILLS.md)

### github-safe-push 技能

**重要**：提交代码前必须运行安全检查，防止 API Key 泄露。

```bash
# 提交前检查
./skills/github-safe-push/pre-commit-check.sh

# 检查特定文件
./skills/github-safe-push/check-file.sh <file>

# 扫描仓库
./skills/github-safe-push/scan-repo.sh

# 安装 Git 钩子（自动检查）
./skills/github-safe-push/install-hook.sh
```

---

## 🦸 Superpowers

Superpowers 插件提供额外的开发增强功能。

### 核心技能

| 技能 | 用途 |
|------|------|
| `brainstorming` | 创意/功能设计 |
| `writing-plans` | 编写实现计划 |
| `test-driven-development` | TDD 开发 |
| `systematic-debugging` | 系统调试 |
| `verification-before-completion` | 完成前验证 |
| `requesting-code-review` | 请求代码审查 |
| `receiving-code-review` | 接收代码审查 |
| `using-git-worktrees` | Git 工作树管理 |
| `finishing-a-development-branch` | 开发分支完成 |

### 使用方式

- **自动加载** - 每次会话自动激活
- **创造性工作前** - 使用 brainstorming
- **复杂任务前** - 使用 writing-plans
- **完成前** - 使用 verification-before-completion

---

## 💡 Oh My OpenCode 使用指南

### 两种工作模式

#### 1. Ultrawork 模式（快速）

```bash
ulw [任务描述]
```

或

```bash
ultrawork [任务描述]
```

自动执行：探索→研究→实现→验证→完成

#### 2. Prometheus 模式（精确）

1. 按 **Tab** 进入 Prometheus 规划模式
2. 描述需求，回答澄清问题
3. 审查生成的计划（`.sisyphus/plans/*.md`）
4. 运行 `/start-work` 执行

### 常用命令

```bash
# 启动 OpenCode
opencode

# 查看可用模型
opencode models

# 连接 AI 提供商
opencode auth login

# 或在 TUI 中使用
/connect
/models
```

---

## 🔄 设备迁移

### 方法 1: 使用安装脚本（推荐）

```bash
# 新设备上运行
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills
./install.sh

# 然后手动填入 API Keys
cp config/opencode.json.example ~/.config/opencode/opencode.json
nano ~/.config/opencode/opencode.json
```

### 方法 2: 手动同步

```bash
# 旧设备导出（不含敏感文件）
rsync -av --exclude='auth.json' --exclude='opencode.json' \
  ~/.config/opencode/ /path/to/backup/

# 新设备导入
rsync -av /path/to/backup/ ~/.config/opencode/

# 手动配置 API Keys
```

---

## 📁 目录结构

```
~/.config/opencode/
├── opencode.json              # 提供商配置（含 API Key，不要提交）
├── oh-my-opencode.json        # Oh My OpenCode 配置（可提交）
├── auth.json                  # 认证信息（不要提交）
├── .gitignore                 # Git 忽略规则
├── plugins/
│   └── superpowers.js         # Superpowers 插件
└── skills/
    ├── scanpy/
    ├── scvi-tools/
    ├── superpowers/           # Superpowers 技能
    ├── github-safe-push/      # Git 安全提交技能
    └── ... (146 个技能)
```

---

## 🔒 安全提交规范

### ⚠️ 禁止提交的文件

以下文件**绝对不能**提交到 GitHub：

| 文件 | 原因 |
|------|------|
| `opencode.json` | 包含 API Key |
| `auth.json` | 包含认证信息 |
| `*.env` | 包含环境变量/密钥 |
| `*.local` | 本地配置文件 |

### ✅ 可以提交的文件

| 文件 | 说明 |
|------|------|
| `config/opencode.json.example` | 配置示例（不含真实 Key） |
| `config/oh-my-opencode.json` | 模型分配配置 |
| `README.md` | 配置说明文档 |
| `SKILLS.md` | Skills 目录 |
| `install.sh` | 安装脚本 |
| `skills/` | Skills 文件 |

### 提交流程

```bash
# 1. 先获取远程更新
git fetch origin

# 2. 合并远程更改
git pull --rebase origin main

# 3. 运行安全检查
./skills/github-safe-push/pre-commit-check.sh

# 4. 提交和推送
git add .
git commit -m "your message"
git push origin main
```

---

## 🛠️ 故障排除

### 插件未加载

```bash
# 检查插件文件
ls -la ~/.config/opencode/plugins/

# 确认文件名为 superpowers.js（带 s）
# 如果是 superpower.js，重命名：
mv ~/.config/opencode/plugins/superpower.js \
   ~/.config/opencode/plugins/superpowers.js
```

### 模型不可用

```bash
# 查看可用模型
opencode models

# 重新认证
opencode auth login
```

### 配置不生效

```bash
# 验证 JSON 语法
jq . ~/.config/opencode/opencode.json
jq . ~/.config/opencode/oh-my-opencode.json

# 重启 OpenCode
```

---

## 📄 License

MIT License

## 🙏 致谢

- **[Oh My OpenCode](https://github.com/code-yeongyu/oh-my-opencode)** - Oh My OpenCode 框架和命令系统
- **[OpenCode](https://github.com/anomalyco/opencode)** - AI 辅助科研计算平台
- **[Superpowers](https://github.com/obra/superpowers)** - Superpowers 插件框架和核心技能
- **K-Dense Inc. Scientific Skills** - 146 个科学和技术技能
