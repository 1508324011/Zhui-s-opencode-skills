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

**配置**: `oh-my-opencode.json` - 定义 Agent 和 Category 的模型分配

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

### 3. OpenCode Memory

**npm 包**: `opencode-mem@latest`  
**GitHub**: https://github.com/tickernelz/opencode-mem

**功能**:
- 本地向量数据库（SQLite + HNSW）
- 跨会话持久记忆
- 自动用户画像学习
- Web UI (http://127.0.0.1:4747)

**配置**: `opencode-mem.jsonc` - 定义记忆系统参数和 API Key

## 🔧 安装方法

### 方法 1: 使用 install.sh（推荐）

```bash
./install.sh
```

自动安装所有三个插件。

### 方法 2: 手动安装

```bash
# 全局安装
npm install -g oh-my-opencode
npm install -g @tarquinen/opencode-dcp
npm install -g opencode-mem

# 或本地安装到 OpenCode 配置目录
cd ~/.config/opencode
npm install oh-my-opencode
npm install @tarquinen/opencode-dcp
npm install opencode-mem
```

## 📋 配置步骤

1. **编辑 opencode.json**
   ```bash
   nano ~/.config/opencode/opencode.json
   ```
   填入 API Key 并确保 plugin 数组包含三个插件

2. **编辑 oh-my-opencode.json**
   ```bash
   nano ~/.config/opencode/oh-my-opencode.json
   ```
   定义 Agent 和 Category 的模型分配

3. **编辑 dcp.jsonc**（可选）
   ```bash
   nano ~/.config/opencode/dcp.jsonc
   ```
   配置上下文限制和 pruning 策略

4. **编辑 opencode-mem.jsonc**
   ```bash
   nano ~/.config/opencode/opencode-mem.jsonc
   ```
   配置记忆系统和 API Key

5. **重启 OpenCode**
   ```bash
   opencode
   ```

## ✅ 验证安装

在 OpenCode 中：

```
/models           # 查看可用模型
/plugins          # 查看已安装插件
/dcp stats        # 查看 DCP 统计（如果已安装 DCP）
```
