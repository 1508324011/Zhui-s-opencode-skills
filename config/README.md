# OpenCode 插件配置指南

本目录包含 OpenCode 核心插件的配置示例。

## 📦 核心插件

### 1. Oh My OpenCode

**npm 包**: `oh-my-opencode@latest`  
**GitHub**: https://github.com/code-yeongyu/oh-my-opencode

**功能**:
- 多 Agent 编排（Sisyphus、Hephaestus、Prometheus 等）
- `ultrawork`/`ulw` 一键执行命令
- Hash-anchored 编辑工具
- LSP + AST-Grep 集成
- 46 个生命周期 Hooks

**配置**: `oh-my-opencode.json` - 定义 OMO/OpenCode 的 Agent 和 Category 模型分配  
**Trellis 补充配置**: `oh-my-openagent.jsonc.example` - 定义 `categories.trellis-*` 作为 Trellis 同步源

### 2. Dynamic Context Pruning (DCP)

**npm 包**: `@tarquinen/opencode-dcp@latest`  
**GitHub**: https://github.com/Opencode-DCP/opencode-dynamic-context-pruning

**功能**:
- `distill` - 总结后删除冗余内容
- `compress` - 压缩对话段落
- `prune` - 移除已完成内容
- 自动去重、覆盖写入、清除错误
- 减少 30-50% token 消耗

**配置**: `dcp.jsonc` - 定义上下文限制和 pruning 策略

### 3. 仓库默认插件说明

此仓库当前默认只保留：

- `oh-my-opencode@latest`
- `@tarquinen/opencode-dcp@latest`

`opencode-mem` 已从默认迁移方案中移除，不再自动安装，也不再提供配置模板。

`Superpowers` 也已从默认迁移方案中移除；结构化、多阶段工作流统一改用 Trellis。

## 🔧 安装方法

### 方法 1: 使用 install.sh（推荐）

```bash
./install.sh
```

自动安装 OMO 与 DCP 两个核心插件。

### 方法 2: 手动安装

```bash
# 全局安装
npm install -g oh-my-opencode
npm install -g @tarquinen/opencode-dcp

# 或本地安装到 OpenCode 配置目录
cd ~/.config/opencode
npm install oh-my-opencode
npm install @tarquinen/opencode-dcp
```

## 📋 配置步骤

1. **编辑 opencode.json**
   ```bash
   nano ~/.config/opencode/opencode.json
   ```
   填入 API Key，并确保 plugin 数组至少包含 `oh-my-opencode@latest` 与 `@tarquinen/opencode-dcp@latest`

2. **编辑 oh-my-opencode.json**
    ```bash
    nano ~/.config/opencode/oh-my-opencode.json
    ```
    定义 Agent 和 Category 的模型分配

3. **如果使用 Trellis，复制并编辑 oh-my-openagent.jsonc**
   ```bash
   cp config/oh-my-openagent.jsonc.example ~/.config/opencode/oh-my-openagent.jsonc
   nano ~/.config/opencode/oh-my-openagent.jsonc
   ```
   把 `categories.trellis-research`、`trellis-implement`、`trellis-check`、`trellis-debug`、`trellis-finish` 作为 Trellis 的唯一模型源。

4. **复制并编辑 dcp.jsonc**
    ```bash
    cp dcp.jsonc ~/.config/opencode/dcp.jsonc
    nano ~/.config/opencode/dcp.jsonc
    ```
    当前仓库提供的是基于 DCP 最新 schema 的配置模板，重点使用 `compress`、`commands`、`manualMode` 和 `strategies` 四组设置。

5. **重启 OpenCode**
    ```bash
    opencode
    ```

6. **如果使用 Trellis，同步活动运行时模型**
   ```bash
   /trellis:sync-models
   ```
   如果你的 Trellis 环境尚未安装该命令，则直接运行 Trellis 仓库中的同步脚本。

## ✅ 验证安装

在 OpenCode 中：

```
/models           # 查看可用模型
/plugins          # 查看已安装插件
/dcp context      # 查看当前上下文使用情况
/dcp stats        # 查看 DCP 统计（如果已安装 DCP）
/trellis:sync-models  # 同步 Trellis 子代理模型（如果已安装 Trellis 命令）
```
