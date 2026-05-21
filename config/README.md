# Superpowers / OMO Line Templates

本目录现在只承载 **superpowers / OMO** 线路使用的模板文件。

## 这条线会安装什么

- `oh-my-opencode@latest`
- `@tarquinen/opencode-dcp@latest`
- 共享 `skills/`
- 官方 `superpowers` 仓库 clone
- 两个 symlink：
  - `~/.config/opencode/plugins/superpowers.js`
  - `~/.config/opencode/skills/superpowers`

## 本目录中的模板

- `oh-my-opencode.json`
- `opencode.json.example`
- `README.md`（本文件）

仓库根目录另外还会为这条线提供：

- `dcp.jsonc`

## 安装

```bash
./install.sh superpowers
```

不带模式参数时：

```bash
./install.sh
```

默认仍然安装 superpowers 线，以保持兼容旧行为。

## 手动安装（仅在需要时）

```bash
npm install -g oh-my-opencode
npm install -g @tarquinen/opencode-dcp
```

## 与 Trellis 线的区别

- Trellis 线模板位于 `../trellis/config/`
- Trellis 线安装 `oh-my-openagent`
- Trellis 线额外安装 `~/.opencode/` skeleton 和 `.trellis` project template

## 注意

- `opencode.json.example` 现在与根目录示例保持一致，只保留 `zhui` provider 示例，并展示 `baseURL`、`apiKey` 占位和 `setCacheKey`
- superpowers 插件通过 `~/.config/opencode/plugins/superpowers.js` symlink 接入
