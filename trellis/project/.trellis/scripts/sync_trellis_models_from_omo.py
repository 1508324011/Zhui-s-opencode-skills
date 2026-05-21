#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

DEFAULT_PROJECT_ROOT = Path("/home/zhurui")
DEFAULT_OMO_CONFIG_PATH = DEFAULT_PROJECT_ROOT / ".config/opencode/oh-my-openagent.jsonc"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Deprecated compatibility shim for the old Trellis model sync command"
    )
    _ = parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_OMO_CONFIG_PATH,
        help="Unused compatibility flag from the old sync workflow",
    )
    _ = parser.add_argument(
        "--project-root",
        type=Path,
        default=DEFAULT_PROJECT_ROOT,
        help="Unused compatibility flag from the old sync workflow",
    )
    _ = parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Accepted for compatibility; no files are rewritten",
    )
    _ = parser.parse_args()

    print("Trellis model sync is deprecated in this repository.")
    print(
        "Vendored Trellis runtime files are no longer rewritten from ~/.config/opencode/oh-my-openagent.jsonc."
    )
    print("Review and edit runtime defaults directly when you intentionally want different local models:")
    print("- /home/zhurui/.opencode/agents/dispatch.md")
    print("- /home/zhurui/.opencode/agents/trellis-plan.md")
    print("- /home/zhurui/.opencode/commands/trellis/start.md")
    print("No changes were made.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
