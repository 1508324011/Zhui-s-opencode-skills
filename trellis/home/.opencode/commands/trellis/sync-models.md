# Sync Trellis Models

Sync the active Trellis subagent model selections from `~/.config/opencode/oh-my-openagent.jsonc` into the active OpenCode runtime files.

[!] Use this command after editing `categories.trellis-*` in `~/.config/opencode/oh-my-openagent.jsonc`, or after reinstalling/updating Trellis.

---

## Sync Procedure

### Step 1: Preview the resolved model mapping

```bash
python3 /home/zhurui/.trellis/scripts/sync_trellis_models_from_omo.py --dry-run
```

Confirm the resolved models for:
- `research`
- `implement`
- `check`
- `debug`
- `finish`

If the script fails, stop and report which `categories.trellis-*` entries are missing or invalid.

### Step 2: Apply the sync

```bash
python3 /home/zhurui/.trellis/scripts/sync_trellis_models_from_omo.py
```

This sync updates the active runtime files:
- `.opencode/agents/dispatch.md`
- `.opencode/agents/trellis-plan.md`
- `.opencode/commands/trellis/start.md`

### Step 3: Verify the active runtime

```bash
python3 /home/zhurui/.trellis/scripts/sync_trellis_models_from_omo.py --dry-run
grep -n 'model:' /home/zhurui/.opencode/agents/dispatch.md /home/zhurui/.opencode/agents/trellis-plan.md /home/zhurui/.opencode/commands/trellis/start.md
```

Report:
- the resolved model mapping
- which files changed or remained unchanged
- whether the active Trellis files now match the configuration in `~/.config/opencode/oh-my-openagent.jsonc`

---

## Source of Truth

Do **not** manually maintain model selections in the active `.opencode` Trellis files.

The source of truth is:
- `~/.config/opencode/oh-my-openagent.jsonc`
- `categories.trellis-research`
- `categories.trellis-implement`
- `categories.trellis-check`
- `categories.trellis-debug`
- `categories.trellis-finish`

Update those values first, then run `/trellis:sync-models` again.
