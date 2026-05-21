# Error Handling

> How errors are handled in this project.

---

## Overview

This project uses **Python exceptions** with a focus on:

1. **Specific exception types** - Always catch specific exceptions, never bare `except:`
2. **Graceful degradation** - Return `None` or sentinel values on recoverable errors
3. **Clear error messages** - Use `ValueError` with descriptive messages for invalid inputs
4. **Silent failures for I/O** - File operations return `None` or `False` instead of crashing

---

## Error Types

### Standard Exceptions

| Exception | When to Use |
|-----------|-------------|
| `ValueError` | Invalid argument values (e.g., wrong type, out of range) |
| `FileNotFoundError` | File doesn't exist |
| `json.JSONDecodeError` | Invalid JSON format |
| `OSError` / `IOError` | File system operations fail |
| `TypeError` | Wrong type passed (rarely raised explicitly) |

### Custom Patterns

No custom exception classes. Use standard library exceptions with descriptive messages.

---

## Error Handling Patterns

### Pattern 1: Try-Return-None

For I/O operations that may fail, return `None` on error:

```python
def read_json(path: Path) -> dict | None:
    """Read and parse a JSON file.
    
    Returns None if the file doesn't exist, is invalid JSON, or can't be read.
    """
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None
```

**Location**: `common/io.py:19-22`

### Pattern 2: Try-Return-Bool

For write operations, return `bool` indicating success:

```python
def write_json(path: Path, data: dict) -> bool:
    """Write dict to JSON file. Returns True on success, False on error."""
    try:
        path.write_text(
            json.dumps(data, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        return True
    except (OSError, IOError):
        return False
```

**Location**: `common/io.py:30-37`

### Pattern 3: Validate and Raise

For invalid inputs, raise `ValueError` with context:

```python
def _is_true_config_value(value: object) -> bool:
    """Return True when a config value represents an enabled flag."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return False


# In CLI adapter:
if len(args) < 2:
    raise ValueError(f"Usage: {command} <args>")
```

**Location**: `common/config.py:24-30`, `common/cli_adapter.py:300-358`

### Pattern 4: Try-Convert with Fallback

For type conversions with defaults:

```python
def get_max_journal_lines(repo_root: Path | None = None) -> int:
    """Get the maximum lines per journal file."""
    config = _load_config(repo_root)
    value = config.get("max_journal_lines", DEFAULT_MAX_JOURNAL_LINES)
    try:
        return int(value)
    except (ValueError, TypeError):
        return DEFAULT_MAX_JOURNAL_LINES
```

**Location**: `common/config.py:55-62`

---

## API Error Responses

### CLI Scripts

CLI scripts return exit codes:

```python
import sys

def main():
    try:
        result = do_work()
        return 0 if result else 1
    except Exception as e:
        log_error(str(e))
        return 1

sys.exit(main())
```

### Library Functions

Library functions use the patterns above (return `None`/`False` or raise).

---

## Common Mistakes

1. **Bare `except:` clauses** - Always catch specific exceptions
   ```python
   # BAD
   try: ...
   except: ...
   
   # GOOD
   try: ...
   except (OSError, IOError): ...
   ```

2. **Raising generic `Exception`** - Use specific types
   ```python
   # BAD
   raise Exception("Invalid input")
   
   # GOOD
   raise ValueError("Expected int, got str")
   ```

3. **Not handling file operations** - Always wrap I/O in try/except
   ```python
   # BAD
   content = path.read_text()  # May crash
   
   # GOOD
   try:
       content = path.read_text()
   except OSError:
       return None
   ```

4. **Swallowing exceptions silently** - Log or return sentinel, don't just pass
   ```python
   # BAD
   try:
       process()
   except:
       pass
   
   # GOOD
   try:
       process()
   except ValueError as e:
       log_warn(f"Processing failed: {e}")
       return None
   ```

---

## Examples from Codebase

### Git Operations

```python
def run_git(args: list[str], cwd: Path | None = None) -> tuple[int, str, str]:
    """Run a git command and return (returncode, stdout, stderr)."""
    try:
        git_args = ["git", "-c", "i18n.logOutputEncoding=UTF-8"] + args
        result = subprocess.run(
            git_args,
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        return result.returncode, result.stdout, result.stderr
    except Exception as e:
        return 1, "", str(e)
```

**Location**: `common/git.py:13-31`

### File Reading with Fallback

```python
def _load_config(repo_root: Path | None = None) -> dict:
    """Load and parse config.yaml. Returns empty dict on any error."""
    config_file = _get_config_path(repo_root)
    try:
        content = config_file.read_text(encoding="utf-8")
        return parse_simple_yaml(content)
    except (OSError, IOError):
        return {}
```

**Location**: `common/config.py:39-46`
