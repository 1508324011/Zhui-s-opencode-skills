# Trellis Model Sync (Deprecated)

This repository no longer rewrites vendored Trellis runtime model values from `~/.config/opencode/oh-my-openagent.jsonc`.

The old sync layer depended on fixed file paths and Markdown `Task(... model: ...)` structure, which made future Trellis upstream updates unnecessarily fragile.

---

## What Changed

- `categories.trellis-*` in `~/.config/opencode/oh-my-openagent.jsonc` are no longer treated as the source of truth for vendored Trellis runtime files.
- `/trellis:sync-models` is kept only as a compatibility notice.
- Trellis runtime defaults now live directly in these files:
  - `.opencode/agents/dispatch.md`
  - `.opencode/agents/trellis-plan.md`
  - `.opencode/commands/trellis/start.md`

---

## If You Need to Change Trellis Model Defaults

1. Check the current upstream Trellis templates first.
2. Edit the vendored runtime files intentionally.
3. Verify the active `model:` lines you changed.

Example verification:

```bash
grep -n 'model:' /home/zhurui/.opencode/agents/dispatch.md /home/zhurui/.opencode/agents/trellis-plan.md /home/zhurui/.opencode/commands/trellis/start.md
```

---

## Compatibility Behavior

If you still run:

```bash
python3 /home/zhurui/.trellis/scripts/sync_trellis_models_from_omo.py
```

the script now prints a deprecation notice and exits without rewriting files.
