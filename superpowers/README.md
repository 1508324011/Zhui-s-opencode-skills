# Superpowers Line

This line keeps the repository's shared skills and OpenAgent configuration, but the actual Superpowers plugin and skill tree come from the official upstream repository.

## What gets installed

- Shared skills from this repository
- `~/.config/opencode/oh-my-openagent.jsonc`
- `~/.config/opencode/dcp.jsonc`
- `~/.config/opencode/opencode.jsonc` (preserved if already present)
- `~/.config/opencode/superpowers` cloned from `https://github.com/obra/superpowers.git`
- Symlinks:
  - `~/.config/opencode/plugins/superpowers.js`
  - `~/.config/opencode/skills/superpowers`

No repository-local `./superpowers` plugin entry is required. The installer follows the upstream OpenCode layout and wires the plugin through the symlink in `~/.config/opencode/plugins/`.

## Install

```bash
./install.sh superpowers
```

If this server does not need bioinformatics / clinical shared skills:

```bash
./install.sh superpowers --exclude-bio-skills
```

Use `--skip-npm-install` when validating the file layout without mutating global npm state.

## Why this line does not vendor `skills/superpowers`

The upstream plugin expects to resolve `using-superpowers/SKILL.md` from the official Superpowers repository layout. Copying the plugin file without the matching upstream repo breaks that assumption, so this installer always uses the official clone + symlink flow.
