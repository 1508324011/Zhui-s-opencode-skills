# Trellis Line

This line packages the Trellis workflow skeleton plus the `oh-my-openagent` and DCP configuration used on this machine.

## What gets installed

- Shared skills from the repository `skills/` tree
- `~/.config/opencode/oh-my-openagent.jsonc`
- `~/.config/opencode/dcp.jsonc`
- `~/.config/opencode/opencode.json`
- `~/.opencode/` plugin and command skeleton
- `~/.config/opencode/trellis-project-template/.trellis/` project template

## Install

```bash
./install.sh trellis
```

Use `--skip-npm-install` when validating the file layout without mutating global npm state.

## Notes

- The generated `opencode.json` template uses only the `zhui` provider example and keeps `baseURL`, an `apiKey` placeholder, and `setCacheKey` without embedding a real key.
- The `.trellis` project template is stored under `~/.config/opencode/trellis-project-template/.trellis`; copy it into a project root when bootstrapping a new Trellis project.
