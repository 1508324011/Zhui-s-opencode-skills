# Zhui's OpenCode Configuration

[![Oh My OpenCode](https://img.shields.io/badge/Oh%20My%20OpenCode-v3.0-blue.svg)](https://github.com/code-yeongyu/oh-my-opencode)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)

OpenCode 配置和 Skills 集合 - 用于快速迁移和设置新设备

## 📋 目录

- [快速开始](#-快速开始)
- [配置说明](#-配置说明)
- [模型分配](#-模型分配)
- [OMO + Trellis 协同工作流](#-omo--trellis-协同工作流)
- [Skills](#skills)
- [已安装插件](#-已安装插件)
- [设备迁移](#-设备迁移)
- [安全提交规范](#-安全提交规范)

---

## 🚀 快速开始

### 新设备一键安装

```bash
# 1. 克隆此仓库
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills

# 2. 运行安装脚本（自动安装 npm 插件和 Skills）
./install.sh

# 3. 复制配置文件并填入 API Key
cp config/opencode.json.example ~/.config/opencode/opencode.json
nano ~/.config/opencode/opencode.json  # 替换 YOUR_API_KEY_HERE

# 4. 重启 OpenCode
opencode
```

**install.sh 自动安装以下 npm 插件**:
- `oh-my-opencode@latest` - 多 Agent 编排框架
- `@tarquinen/opencode-dcp@latest` - 动态上下文修剪

> `opencode-mem` 和 `Superpowers` 都已从这个仓库的默认迁移方案里移除。当前默认组合是 **OMO + DCP**，结构化多阶段工作流统一走 **Trellis**。

### 手动配置

```bash
# 1. 复制配置文件
cp config/opencode.json.example ~/.config/opencode/opencode.json
cp config/oh-my-opencode.json ~/.config/opencode/oh-my-opencode.json
cp config/oh-my-openagent.jsonc.example ~/.config/opencode/oh-my-openagent.jsonc
cp dcp.jsonc ~/.config/opencode/dcp.jsonc

# 2. 编辑 opencode.json 填入你的 API Key
nano ~/.config/opencode/opencode.json

# 3. 如需调整 DCP 行为，编辑 dcp.jsonc
nano ~/.config/opencode/dcp.jsonc

# 4. 如果使用 Trellis，编辑 categories.trellis-*
nano ~/.config/opencode/oh-my-openagent.jsonc

# 5. 重启 OpenCode
opencode
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

## 🌉 OMO + Trellis 协同工作流

### 直接从 OMO 输入 vs 从 Trellis 输入

| 入口 | 更适合的场景 | 优势 | 代价 |
|------|--------------|------|------|
| **直接和 OMO 对话** | 快速任务、探索式排障、需求还没想清楚时 | 启动快、交互自然、适合边想边做 | 过程结构化程度较低 |
| **`/trellis:start` / `/trellis:brainstorm`** | 多模块任务、长链路实现、需要 PRD / research / check 分阶段沉淀 | 流程稳定、便于复现、适合长期项目和复杂交付 | 前置准备更多，启动成本更高 |

### 推荐用法

1. **先用 OMO 定义问题**：当需求模糊、需要探索代码或先做小范围验证时，直接在 OpenCode 中对 OMO 下达目标，或者使用 `ultrawork` / `ulw`。
2. **再用 Trellis 稳定执行**：当任务已经清晰、涉及多个模块或希望走 research → implement → check → finish 流程时，切换到 `/trellis:start` 或 `/trellis:brainstorm`。
3. **遇到难点再回 OMO**：如果 Trellis 流程中的某个子问题特别棘手，可以把局部问题交回 OMO 深挖，再把结论写回 Trellis 任务目录继续推进。

### 当前仓库的工作流取舍

- `Superpowers` 已不再作为这个仓库的默认插件或默认技能来源
- 需要结构化、多阶段执行时，统一改走 Trellis 工作流
- 因此这里保留的是：**OMO 负责主代理能力，DCP 负责上下文管理，Trellis 负责复杂任务编排**

### Trellis 模型同步：把 `categories.trellis-*` 当成唯一来源

如果你想让 Trellis 子代理自动继承这套模型分配，使用仓库里的模板：

```bash
# 1. 复制 Trellis 模型同步模板
cp config/oh-my-openagent.jsonc.example ~/.config/opencode/oh-my-openagent.jsonc

# 2. 编辑五个 Trellis 类别
nano ~/.config/opencode/oh-my-openagent.jsonc
```

默认模板包含以下五个分类：

- `trellis-research`
- `trellis-implement`
- `trellis-check`
- `trellis-debug`
- `trellis-finish`

推荐把它们视为 **Trellis 的唯一模型源**。也就是说，平时只改 `~/.config/opencode/oh-my-openagent.jsonc` 里的这五项，然后让 Trellis 运行时自动同步，而不是手工去改活动中的 `.opencode/agents/*.md`。

### 如何执行同步

如果你的 Trellis 命令集已经包含 `/trellis:sync-models`，优先直接运行：

```text
/trellis:sync-models
```

如果还没有这个命令，就直接运行 Trellis 仓库里的同步脚本：

```bash
python3 /path/to/Trellis/.trellis/scripts/sync_trellis_models_from_omo.py --dry-run
python3 /path/to/Trellis/.trellis/scripts/sync_trellis_models_from_omo.py
```

同步完成后，再用一次 `--dry-run` 做回归检查，确认活动运行时已经和 `categories.trellis-*` 一致。

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

### 新增迁移技能

这次仓库额外收录了 3 个已经过本地适配和验证的技能，运行 `install.sh` 时会随 `skills/` 一起复制到新设备：

| 技能 | 仓库路径 | 用途 | 备注 |
|------|----------|------|------|
| `frontend-design` | `skills/frontend-design/` | 引入 Anthropic 的前端设计审美约束，生成更有辨识度的 UI / 页面实现 | 适合作为 OpenCode 原生 skill 直接使用 |
| `skill-creator` | `skills/skill-creator/` | 把 Anthropic 的技能设计/评审/打包工作流迁移到 OpenCode | 保留上游脚本；其中部分 trigger-eval 脚本仍是 Claude Code 参考实现 |
| `web-access` | `skills/web-access/` | 提供联网检索、登录态网页访问、浏览器 CDP 交互和站点经验积累能力 | 依赖 Node.js 22+ 与 Chrome remote debugging，见下方迁移说明 |

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

## 💾 已安装插件

此配置使用以下 OpenCode 插件：

### Oh My OpenCode

**仓库**: https://github.com/code-yeongyu/oh-my-opencode

**核心功能**:
- **多 Agent 编排** — Sisyphus（协调员）、Hephaestus（深度执行）、Prometheus（战略规划）、Oracle（调试）、Librarian（文档搜索）、Explore（代码探索）并行工作
- **`ultrawork`/`ulw` 命令** — 一键激活所有 Agent 直到任务完成
- **Hash-anchored 编辑工具** — LINE#ID 内容哈希验证防止编辑错误
- **LSP + AST-Grep 集成** — IDE 级精度工具（重命名、跳转定义、查找引用、诊断检查）
- **46 个生命周期 Hooks** — 可扩展的 Hook 系统自定义 Agent 行为
- **26 个内置工具** — 后台 Agent、tmux 集成、会话管理、技能嵌入 MCPs
- **内置 MCPs** — Exa（网络搜索）、Context7（官方文档）、Grep.app（GitHub 代码搜索）
- **Claude Code 兼容** — 完全支持现有 Hooks、Commands、Skills、MCPs 和插件
- **Ralph Loop** (`/ulw-loop`) — 自引用执行循环直到 100% 任务完成
- **Todo Enforcer** — 防止 Agent 在任务中途停止

### Dynamic Context Pruning (DCP)

**仓库**: https://github.com/Opencode-DCP/opencode-dynamic-context-pruning

**核心功能**:
- **三个核心工具** — `distill`（总结后删除）、`compress`（压缩对话段落）、`prune`（移除已完成内容）
- **自动去重** — 识别并删除重复的工具调用（如同一文件多次读取）
- **覆盖写入** — 删除后续被读取的写入工具调用，消除冗余
- **清除错误** — 在可配置轮次后删除错误工具输入（默认 4 轮）
- **零成本策略** — 每次请求自动运行 pruning 策略，无额外 token 消耗
- **可配置上下文限制** — Token 阈值（默认 100,000）或基于百分比的限制
- **斜杠命令** — `/dcp context`（token 使用分析）、`/dcp stats`（pruning 统计）
- **保护关键工具** — 保护重要工具（task、todowrite、distill、compress、prune 等）不被删除

**效果**:
- 减少 30-50% token 消耗
- 防止上下文窗口膨胀
- 更快的响应速度
- 更专注的 Agent 注意力

**当前仓库默认配置**:
- 默认安装 `@tarquinen/opencode-dcp@latest`
- 默认提供根目录 `dcp.jsonc`，可直接复制到 `~/.config/opencode/dcp.jsonc`
- 配置已切换到最新 schema，主要围绕 `compress`、`commands`、`manualMode`、`strategies` 四个配置块
- 常用命令是 `/dcp context`、`/dcp stats`、`/dcp manual on|off`

### 已移除的默认插件

- `opencode-mem` 已从仓库默认安装流程中移除
- 仓库不再提供 `opencode-mem.jsonc.example`
- 如果旧设备里还残留 `opencode-mem.jsonc`，请保留在 `.gitignore` 中，但不要再把它当成新设备的默认必装项
- `Superpowers` 已从仓库默认安装流程中移除
- 仓库不再分发 `plugins/superpowers.js`，也不再把 `./superpowers` 作为默认插件入口
- 原先依赖 Superpowers 的结构化开发流程，现改由 Trellis 负责
- 如果旧设备里还保留 `plugins/superpowers.js`、`skills/superpowers` 或 `opencode.json` 中的 `./superpowers`，请手动清理；重新运行 `install.sh` 时也会尝试删除这些遗留项

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

# 如需调整 DCP，复制当前模板
cp dcp.jsonc ~/.config/opencode/dcp.jsonc

# 如需 Trellis，同步 Trellis 模型模板
cp config/oh-my-openagent.jsonc.example ~/.config/opencode/oh-my-openagent.jsonc
```

`install.sh` 会自动把仓库中的 `skills/frontend-design`、`skills/skill-creator`、`skills/web-access` 一并复制到 `~/.config/opencode/skills/`，不需要额外登记插件。

其中 `web-access` 迁移后还需要额外确认：

1. `web-access` 默认要求新设备使用 **Node.js 22+**；低于 22 时，只有运行环境本身已经能解析 `ws` 模块时才可能工作，`install.sh` 不会自动补这个依赖。
2. Chrome 中已开启 `chrome://inspect/#remote-debugging` 并勾选 **Allow remote debugging for this browser instance**。
3. 如需手动检查环境，可在 `~/.config/opencode/skills/web-access/` 下运行 `node "./scripts/check-deps.mjs"`。

### 方法 2: 手动同步

```bash
# 旧设备导出（不含敏感文件）
rsync -av --exclude='auth.json' --exclude='opencode.json' --exclude='opencode-mem.jsonc' \
  ~/.config/opencode/ /path/to/backup/

# 新设备导入
rsync -av /path/to/backup/ ~/.config/opencode/

# 手动配置 API Keys，并根据需要清理旧的 opencode-mem.jsonc、plugins/superpowers.js，以及 opencode.json 中残留的 ./superpowers 插件入口
```

---

## 📁 目录结构

```
~/.config/opencode/
├── opencode.json              # 提供商配置（含 API Key，不要提交）
├── oh-my-opencode.json        # Oh My OpenCode 配置（可提交）
├── dcp.jsonc                  # DCP 插件配置（可提交）
├── oh-my-openagent.jsonc      # Trellis 模型同步源（可提交）
├── auth.json                  # 认证信息（不要提交）
├── .gitignore                 # Git 忽略规则
├── plugins/                   # OpenCode 插件目录
└── skills/
    ├── scanpy/
    ├── scvi-tools/
    ├── github-safe-push/      # Git 安全提交技能
    └── ... (包含 `frontend-design`、`skill-creator`、`web-access` 等技能)
```

---

## 🔒 安全提交规范

### ⚠️ 禁止提交的文件

以下文件**绝对不能**提交到 GitHub：

| 文件 | 原因 |
|------|------|
| `opencode.json` | 包含 API Key |
| `opencode-mem.jsonc` | 历史遗留 Memory 配置，若仍存在可能含敏感信息 |
| `auth.json` | 包含认证信息 |
| `*.env` | 包含环境变量/密钥 |
| `*.local` | 本地配置文件 |

### ✅ 可以提交的文件

| 文件 | 说明 |
|------|------|
| `config/opencode.json.example` | 配置示例（不含真实 Key） |
| `config/oh-my-opencode.json` | 模型分配配置 |
| `config/oh-my-openagent.jsonc.example` | Trellis 模型同步模板 |
| `dcp.jsonc` | DCP 插件配置模板 |
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

### Trellis / DCP 未生效

```bash
# 检查关键配置文件是否存在
ls ~/.config/opencode/dcp.jsonc
ls ~/.config/opencode/oh-my-openagent.jsonc
```

进入 OpenCode 后，先执行 `/dcp stats`；如果你使用 Trellis，再执行 `/trellis:sync-models` 或对应的 Trellis 同步脚本。

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

# 确认 DCP 配置位置
ls ~/.config/opencode/dcp.jsonc

# 快速检查 JSONC 关键项
grep -n 'model' ~/.config/opencode/oh-my-opencode.json
grep -n 'trellis-' ~/.config/opencode/oh-my-openagent.jsonc

# 重启 OpenCode
opencode
```

进入 OpenCode 后，再执行 `/dcp context` 或 `/dcp stats` 检查 DCP 是否正常加载。

---

## 📄 License

MIT License

## 🙏 致谢

本配置整合了以下优秀项目和贡献者：

### 框架和平台

- **[Oh My OpenCode](https://github.com/code-yeongyu/oh-my-opencode)** — 多 Agent 编排框架，提供 Sisyphus、Hephaestus、Prometheus 等纪律 Agent 和 `ultrawork` 自动化系统
- **[OpenCode](https://github.com/anomalyco/opencode)** — AI 辅助科研计算平台，提供基础 Agent 架构和插件系统
- **[Trellis](https://github.com/1508324011/Trellis)** — 结构化任务编排工作流，负责复杂任务的 research / implement / check / finish 流程

### 核心插件

- **[opencode-dynamic-context-pruning](https://github.com/Opencode-DCP/opencode-dynamic-context-pruning)** — 智能上下文管理插件，通过 distill/compress/prune 工具和自动去重减少 30-50% token 消耗

### 技能集合

- **K-Dense Inc. Scientific Skills** — 146 个科学和技术技能，覆盖生物信息学、药物发现、机器学习、统计分析、可视化、科研写作等领域

---
