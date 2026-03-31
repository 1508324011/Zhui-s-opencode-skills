# Skills 索引

这个仓库把可迁移的 OpenCode skills 放在根目录 `skills/` 下；运行 `./install.sh` 时，会整体复制到 `~/.config/opencode/skills/`。

## 本次新增并重点维护的迁移技能

### `frontend-design`

- 仓库路径：`skills/frontend-design/`
- 来源：Anthropic `claude-code/plugins/frontend-design`
- 作用：为前端页面、组件和产品界面提供更强的设计方向约束，避免 generic AI aesthetics。

### `skill-creator`

- 仓库路径：`skills/skill-creator/`
- 来源：Anthropic `skills/skill-creator`
- 作用：在 OpenCode 中创建、改进、验证和打包 skills。
- 说明：仓库保留了 upstream 的 `agents/`、`assets/`、`eval-viewer/`、`scripts/` 等资源；其中部分 trigger-eval 脚本仍然是 Claude Code 参考实现，而不是默认 OpenCode 路径。

### `web-access`

- 仓库路径：`skills/web-access/`
- 来源：<https://github.com/eze-is/web-access>
- 作用：联网检索、网页登录态访问、动态渲染页面探索、浏览器 CDP 交互、站点经验沉淀。
- 迁移前提：
  - 默认要求 Node.js **22+**；低于 22 时，只有运行环境本身已可解析 `ws` 模块时才可能工作
  - 需要在 Chrome 的 `chrome://inspect/#remote-debugging` 中允许当前实例 remote debugging
  - 安装后目录会位于 `~/.config/opencode/skills/web-access/`

## 使用说明

- 新设备安装时，直接运行仓库根目录的 `./install.sh`。
- 不要把 `opencode.json`、`auth.json`、`*.env` 等敏感文件提交到仓库。
- 对于 `web-access` 这类依赖本机运行环境的 skill，请在新设备上再次执行其自检脚本确认前置条件是否满足。
