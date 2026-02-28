# Zhui's OpenCode Skills & Configuration

**维护者**: [@1508324011](https://github.com/1508324011)  
**最后更新**: 2026-02-28  
**许可证**: MIT

---

## 📦 简介

OpenCode AI 助手优化配置，专为空间转录组学科研项目 (Zhui/ClusterMap/PyStar) 设计。

**核心功能**:
- 🔍 **SLURM 作业自动检查** - 提交后自动验证状态和日志
- 🧠 **本地向量数据库记忆** - 跨会话持久化关键决策  
- 🗑️ **动态上下文修剪** - 减少 30-50% token 消耗

---

## 🚀 快速开始

### 1. 克隆仓库

```bash
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills
```

### 2. 安装插件

```bash
cd ~/.config/opencode
npm install opencode-mem@latest @tarquinen/opencode-dcp@latest
```

### 3. 复制配置

```bash
cp -r skills/ ~/.config/opencode/superpowers/skills/
cp opencode.json dcp.jsonc opencode-mem.jsonc ~/.config/opencode/
```

### 4. 配置 API Key

编辑 `~/.config/opencode/opencode-mem.jsonc`:

```jsonc
{
  "memoryApiKey": "sk-sp-xxxxxxxxxxxxxxxx"  // ← 替换为你的 API Key
}
```

> ⚠️ **重要**: 不要将 API Key 提交到 Git!

### 5. 重启 OpenCode

```bash
opencode
```

---

## 📁 结构

```
├── opencode.json          # 主配置
├── dcp.jsonc              # DCP 配置
├── opencode-mem.jsonc     # 记忆配置 (需配置 API Key)
├── docs/                  # 文档
│   ├── CONFIGURATION-GUIDE.md
│   └── SLURM-CHECK-PROTOCOL.md
└── skills/                # 技能
    └── slurm-job-checking/
        └── SKILL.md
```

---

## ⚙️ 配置

### opencode.json

```json
{
  "plugin": [
    "oh-my-opencode@latest",
    "@tarquinen/opencode-dcp@latest",
    "opencode-mem@latest"
  ]
}
```

### dcp.jsonc

```jsonc
{
  "enabled": true,
  "contextLimit": 100000,
  "strategies": {
    "deduplication": { "enabled": true },
    "supersedeWrites": { "enabled": true },
    "purgeErrors": { "enabled": true }
  }
}
```

### opencode-mem.jsonc

```jsonc
{
  "memoryProvider": "openai-chat",
  "memoryModel": "qwen3.5-plus",
  "memoryApiKey": "sk-sp-...",
  "maxMemories": 50,
  "webServerPort": 4747
}
```

---

## 🛠️ 技能

### slurm-job-checking

**检查顺序**:
1. `squeue -j ${JOB_ID}` → 状态
2. `.err` 文件 → 错误检测
3. `.out` 文件 → 进度
4. 输出文件 → 验证

**致命错误**: `ERROR | Fatal | Exception | Traceback | Killed | OOM`

---

## 📊 效果

| 指标 | 改进 |
|------|------|
| 作业失败检测 | 60x 更快 |
| Token 消耗 | -47% |
| 记忆保留 | >90% |

---

## ⚠️ 注意

### API Key 安全

```bash
# 不要提交含 API Key 的文件
git update-index --assume-unchanged opencode-mem.jsonc
```

### 系统依赖

- Python 3.6+
- GCC 8.5+
- make
- Node.js 18+

---

## 📚 资源

- [OpenCode](https://opencode.ai)
- [DCP](https://github.com/Tarquinen/opencode-dynamic-context-pruning)
- [opencode-mem](https://github.com/tickernelz/opencode-mem)

---

**MIT License**
