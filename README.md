# Zhui's OpenCode Configuration

This repository now maintains two incompatible OpenCode install lines:

- `superpowers`: OMO + DCP + upstream Superpowers repo
- `trellis`: OpenAgent + DCP + Trellis workflow skeleton

Both lines install the shared `skills/` tree from this repository. The difference is in the installer, config templates, and extra skeleton directories.

## Quick Start

```bash
git clone https://github.com/1508324011/Zhui-s-opencode-skills.git
cd Zhui-s-opencode-skills
```

### Superpowers Line

```bash
./install.sh
# equivalent to: ./install.sh superpowers
```

This installs:

- `oh-my-opencode@latest`
- `@tarquinen/opencode-dcp@latest`
- shared `skills/`
- upstream `superpowers` repository clone under `~/.config/opencode/superpowers`
- symlinks:
  - `~/.config/opencode/plugins/superpowers.js`
  - `~/.config/opencode/skills/superpowers`

Template sources:

- `config/oh-my-opencode.json`
- `config/opencode.json.example`
- `dcp.jsonc`

### Trellis Line

```bash
./install.sh trellis
```

This installs:

- `oh-my-openagent@latest`
- `@tarquinen/opencode-dcp@latest`
- shared `skills/`
- `~/.config/opencode/oh-my-openagent.jsonc`
- `~/.config/opencode/dcp.jsonc`
- `~/.config/opencode/opencode.json`
- `~/.opencode/` skeleton
- `~/.config/opencode/trellis-project-template/.trellis/` project template

Template sources:

- `trellis/config/oh-my-openagent.jsonc`
- `trellis/config/dcp.jsonc`
- `trellis/config/opencode.json.example`

### Layout-Only Validation

```bash
./install.sh superpowers --skip-npm-install
./install.sh trellis --skip-npm-install
```

Use this when you want to validate file layout without mutating global npm state.

## Repository Layout

```text
.
|-- skills/                     # shared skills for both lines
|-- config/                     # superpowers / OMO templates
|-- trellis/                    # trellis templates and skeletons
|   |-- config/
|   |-- home/.opencode/
|   +-- project/.trellis/
|-- superpowers/                # superpowers line notes
|-- install.sh                  # dual-line installer
|-- dcp.jsonc
`-- README.md
```

## Line Differences

| Item | superpowers | trellis |
|---|---|---|
| Primary plugin | `oh-my-opencode` | `oh-my-openagent` |
| DCP | install | install |
| Superpowers repo | upstream clone + symlinks | not installed |
| `~/.opencode` skeleton | not managed | installed |
| `.trellis` project template | not managed | installed |
| Provider template | `zhui` | `zhui` |

## Config Notes

### Superpowers / OMO

The superpowers line uses:

- `config/oh-my-opencode.json`
- `config/opencode.json.example`
- `dcp.jsonc`

This repository no longer vendors `skills/superpowers`. The installer follows the upstream OpenCode layout instead:

- clone upstream repo to `~/.config/opencode/superpowers`
- create `plugins/superpowers.js` symlink
- create `skills/superpowers` symlink

The config templates now use a single `zhui` provider example. Each template keeps `baseURL`, `apiKey` placeholder, and `setCacheKey` as structure examples, but never stores a real key.

### Trellis

The trellis line uses:

- `trellis/config/oh-my-openagent.jsonc`
- `trellis/config/dcp.jsonc`
- `trellis/config/opencode.json.example`

Its config templates are aligned to the `zhui/...` model names referenced by `oh-my-openagent.jsonc` and the vendored `.opencode` skeleton. The example `opencode.json` keeps `baseURL`, an `apiKey` placeholder, and `setCacheKey`, but no real credential.

## Shared Skills

`skills/` is the shared layer for both lines.

Notes:

- repo-only skills remain in this repository
- synchronized local skills are installed to both lines
- `skills/superpowers` is not treated as shared content; it comes from the upstream repo during superpowers installation

## Security

- Fill API keys manually after install.
- Example configs intentionally keep only the `zhui` provider, plus `baseURL`, an `apiKey` placeholder, and `setCacheKey`; replace the placeholder in your local copy.
- Do not commit real keys back into this repository.
- The installer creates a target-side `.gitignore` for `opencode.json`, `auth.json`, `*.env`, `*.local.*`, secrets, logs, caches, and `.claude_resources.json`.

## Line-Specific Notes

- `superpowers/README.md`
- `trellis/README.md`
