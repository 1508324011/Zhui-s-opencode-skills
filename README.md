# Zhui's OpenCode Configuration

[![Oh My OpenCode](https://img.shields.io/badge/Oh%20My%20OpenCode-v3.0-blue.svg)](https://github.com/code-yeongyu/oh-my-opencode)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)

OpenCode 配置和 Skills 集合 - 用于快速迁移和设置新设备

## 📋 目录

- [快速开始](#快速开始)
- [配置说明](#配置说明)
- [模型分配](#模型分配)
- [Skills](#skills)
- [Superpowers](#superpowers)
- [设备迁移](#设备迁移)

---

## 🚀 快速开始

### 新设备一键安装

```bash
# 1. 克隆此仓库
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills

# 2. 运行安装脚本
./install.sh
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

#### 1. 百炼平台 Coding Plan（Qwen 3.5）

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
- `qwen3.5-plus` - 主力模型（带思考模式）
- `qwen3-max-2026-01-23` - 最强推理
- `qwen3-coder-plus` - 代码专用
- `qwen3-coder-next` - 新一代代码模型

#### 2. Moonshot AI（Kimi）

通过 OpenCode 内置命令配置：

```bash
opencode auth login
# 选择 Moonshot AI → 输入 API Key
```

或使用 `/connect` 命令在 TUI 中配置。

**可用模型：**
- `moonshotai/kimi-k2.5` - 最新 K2.5
- `moonshotai/kimi-k2-thinking` - 思考模式
- `moonshotai/kimi-k2-turbo-preview` - 快速版

---

## 🎯 模型分配

### Agents 配置

| Agent | 模型 | 用途 |
|-------|------|------|
| **Sisyphus** | Qwen3.5 Plus | 主代理/协调员 |
| **Oracle** | Kimi K2.5 | 调试/架构分析 |
| **Librarian** | Qwen3.5 Plus | 文档/代码搜索 |
| **Explore** | Qwen3.5 Plus | 代码探索 |
| **Prometheus** | Qwen3.5 Plus | 规划师 |
| **Metis** | Kimi K2.5 | 规划顾问 |
| **Atlas** | Qwen3.5 Plus | 执行者 |

### Categories 配置

| 类别 | 模型 | 用途 |
|------|------|------|
| **visual-engineering** | Kimi K2.5 | 前端/UI 工程 |
| **ultrabrain** | Kimi K2.5 | 深度脑力 |
| **deep** | Kimi K2.5 | 深度分析 |
| **writing** | Kimi K2.5 | 写作任务 |
| **其他** | Qwen3.5 Plus | 通用任务 |

---

## 📚 Skills

此仓库包含完整的 Skills 集合，安装后自动可用。

### 主要技能类别

1. **生物信息学** (32 skills) - scanpy, scvi-tools, cellxgene-census
2. **药物发现** (15 skills) - rdkit, deepchem, diffdock
3. **机器学习** (18 skills) - scikit-learn, transformers, pytorch-lightning
4. **统计分析** (8 skills) - statsmodels, pymc, scikit-survival
5. **可视化** (7 skills) - matplotlib, seaborn, plotly
6. **科研写作** (11 skills) - scientific-writing, literature-review
7. **医疗健康** (11 skills) - pyhealth, clinical-reports
8. **实验自动化** (6 skills) - opentrons, benchling
9. **云计算** (6 skills) - dnanexus, latchbio
10. **科学计算** (8 skills) - matlab, qiskit, fluidsim
11. **网络分析** (3 skills) - networkx, cytoscape
12. **专用工具** (12 skills) - astropy, market-research-reports

完整列表见 [SKILLS.md](./SKILLS.md)

---

## 🦸 Superpowers

Superpowers 插件提供额外的开发增强功能：

### 核心技能

| 技能 | 用途 |
|------|------|
| **brainstorming** | 创意/功能设计 |
| **writing-plans** | 编写实现计划 |
| **test-driven-development** | TDD 开发 |
| **systematic-debugging** | 系统调试 |
| **verification-before-completion** | 完成前验证 |
| **requesting-code-review** | 请求代码审查 |

### 使用方式

- **自动加载** - 每次会话自动激活
- **创造性工作前** - 使用 brainstorming
- **复杂任务前** - 使用 writing-plans
- **完成前** - 使用 verification-before-completion

---

## 💡 Oh My OpenCode 使用指南

### 两种工作模式

#### 1. Ultrawork 模式（快速）

```
ulw [任务描述]
```

或

```
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

## 🔒 安全提醒

### ⚠️ 不要提交敏感信息

此仓库已配置 `.gitignore`，以下文件**绝不会**被提交：

- `opencode.json` - 包含 API Key
- `auth.json` - 认证信息
- `oh-my-opencode.json` - 可提交（不含敏感信息）
- `.env*` - 环境变量

### 安全实践

1. **API Key 存储在本地** - 使用 `opencode auth login` 安全存储
2. **使用环境变量** - `apiKey: "{env:YOUR_API_KEY}"`
3. **定期轮换密钥** - 建议每 90 天更换

---

## 🔄 设备迁移

### 方法 1: 使用安装脚本（推荐）

```bash
# 新设备上运行
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills
./install.sh

# 然后手动填入 API Keys
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

### 方法 3: 使用 Git 同步配置

```bash
# 在 ~/.config/opencode/ 初始化 git
cd ~/.config/opencode
git init
git add oh-my-opencode.json skills/ plugins/
git commit -m "Backup config"

# 推送到私有仓库
git remote add origin <your-private-repo>
git push -u origin main
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
    └── ... (146 个技能)
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

- [Oh My OpenCode](https://github.com/code-yeongyu/oh-my-opencode)
- [OpenCode](https://github.com/anomalyco/opencode)
- [Superpowers](https://github.com/obra/superpowers)
- K-Dense Inc. Scientific Skills
