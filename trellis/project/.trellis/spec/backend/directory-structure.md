# Directory Structure

> How backend code is organized in this project.

---

## Overview

This project follows a **modular Python package structure** with clear separation between:

1. **Task Management System** (`.trellis/scripts/`) - Workflow automation scripts
2. **Scientific Package** (`ClusterMap/`) - Spatial transcriptomics analysis library

---

## Directory Layout

```
.
├── .trellis/
│   ├── scripts/
│   │   ├── common/          # Shared utilities (types, paths, logging, etc.)
│   │   ├── multi_agent/     # Multi-agent orchestration scripts
│   │   ├── hooks/           # Lifecycle hooks
│   │   └── *.py            # Top-level CLI scripts
│   ├── spec/               # Development guidelines
│   └── tasks/              # Task storage
│
├── ClusterMap/             # Main scientific package
│   ├── ClusterMap/         # Package source
│   │   ├── __init__.py
│   │   ├── clustermap.py   # Main analysis class
│   │   ├── preprocessing.py
│   │   ├── postprocessing.py
│   │   ├── utils.py        # Utility functions
│   │   └── ...
│   └── setup.py            # Package setup
│
└── workspace/              # Developer workspaces
    └── <developer>/
```

---

## Module Organization

### Common Utilities (`common/`)

All shared functionality lives in `.trellis/scripts/common/`:

| Module | Purpose |
|--------|---------|
| `types.py` | Type definitions (TypedDict, dataclasses) |
| `paths.py` | Path constants and resolution |
| `log.py` | Structured logging utilities |
| `io.py` | JSON file I/O |
| `git.py` | Git command execution |
| `config.py` | Configuration management |

**Rule**: New shared utilities go here, not in feature modules.

### Feature Scripts (`multi_agent/`, `hooks/`)

- Each major feature gets its own directory
- Imports from `common/` for shared functionality
- Scripts use `from __future__ import annotations` for forward references

### Scientific Package (`ClusterMap/`)

- Pure Python package for spatial transcriptomics
- Self-contained with its own `setup.py`
- Heavy use of scientific libraries (numpy, pandas, scanpy, etc.)

---

## Naming Conventions

### Files

- **Modules**: `snake_case.py` (e.g., `task_context.py`, `preprocessing.py`)
- **Scripts**: Descriptive verbs/nouns (e.g., `start.py`, `cleanup.py`)
- **Private**: Prefix with underscore for internal use (e.g., `_bootstrap.py`)

### Directories

- **Packages**: snake_case (e.g., `multi_agent/`, `common/`)
- **Constants**: UPPER_SNAKE_CASE in `paths.py`

---

## Examples

### Well-Organized Module

**`common/io.py`**:
```python
"""
JSON file I/O utilities.

Provides read_json and write_json as the single source of truth
for JSON file operations across all Trellis scripts.
"""

from __future__ import annotations

import json
from pathlib import Path


def read_json(path: Path) -> dict | None:
    """Read and parse a JSON file."""
    ...


def write_json(path: Path, data: dict) -> bool:
    """Write dict to JSON file with pretty formatting."""
    ...
```

### Package Structure

**`common/paths.py`**:
```python
# =============================================================================
# Path Constants (change here to rename directories)
# =============================================================================

# Directory names
DIR_WORKFLOW = ".trellis"
DIR_WORKSPACE = "workspace"
DIR_TASKS = "tasks"
```

---

## Common Mistakes

1. **Putting shared code in feature modules** - Use `common/` for reusable utilities
2. **Hardcoding paths** - Use `paths.py` constants
3. **Mixing CLI and library code** - Keep CLI entry points separate from business logic
