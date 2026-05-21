#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import cast

DEFAULT_PROJECT_ROOT = Path("/home/zhurui")
DEFAULT_OMO_CONFIG_PATH = DEFAULT_PROJECT_ROOT / ".config/opencode/oh-my-openagent.jsonc"

CATEGORY_TO_AGENT = {
    "trellis-research": "research",
    "trellis-implement": "implement",
    "trellis-check": "check",
    "trellis-debug": "debug",
    "trellis-finish": "finish",
}

def strip_jsonc_comments(text: str) -> str:
    result: list[str] = []
    i = 0
    in_string = False
    escaped = False
    while i < len(text):
        ch = text[i]
        nxt = text[i + 1] if i + 1 < len(text) else ""
        if in_string:
            result.append(ch)
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                in_string = False
            i += 1
            continue
        if ch == '"':
            in_string = True
            result.append(ch)
            i += 1
            continue
        if ch == "/" and nxt == "/":
            while i < len(text) and text[i] != "\n":
                i += 1
            continue
        result.append(ch)
        i += 1
    return "".join(result)


def expect_object(value: object, label: str, config_path: Path) -> dict[str, object]:
    if not isinstance(value, dict):
        raise SystemExit(f"Expected {label} object in {config_path}")
    return cast(dict[str, object], value)


def load_trellis_models(config_path: Path) -> dict[str, str]:
    raw = config_path.read_text(encoding="utf-8")
    parsed_obj = cast(object, json.loads(strip_jsonc_comments(raw)))
    parsed = expect_object(parsed_obj, "top-level", config_path)
    categories_obj = expect_object(parsed.get("categories", {}), "categories", config_path)

    missing: list[str] = []
    resolved: dict[str, str] = {}
    for category_name, agent_name in CATEGORY_TO_AGENT.items():
        category_value = expect_object(categories_obj.get(category_name, {}), category_name, config_path)
        if not category_value:
            missing.append(category_name)
            continue

        model_obj = category_value.get("model")
        if not isinstance(model_obj, str) or not model_obj:
            missing.append(category_name)
            continue
        resolved[agent_name] = model_obj
    if missing:
        joined = ", ".join(missing)
        raise SystemExit(f"Missing Trellis model categories in {config_path}: {joined}")
    return resolved


TASK_BLOCK_RE = re.compile(
    r'(?P<prefix>Task\(\s*subagent_type:\s*"(?P<agent>[^"]+)".*?prompt:\s*"(?P<prompt>.*?)".*?model:\s*")(?P<model>[^"]+)(?P<suffix>"\s*[,)])',
    re.DOTALL,
)


def resolve_model(agent: str, prompt: str, models: dict[str, str]) -> str | None:
    if agent == "research":
        return models["research"]
    if agent == "implement":
        return models["implement"]
    if agent == "debug":
        return models["debug"]
    if agent == "check":
        if "[finish]" in prompt:
            return models["finish"]
        return models["check"]
    return None


def build_target_files(project_root: Path) -> dict[Path, dict[str, bool]]:
    return {
        project_root / ".opencode/agents/dispatch.md": {"replace_constraint": True},
        project_root / ".opencode/agents/trellis-plan.md": {},
        project_root / ".opencode/commands/trellis/start.md": {},
    }


def sync_file(path: Path, models: dict[str, str], replace_constraint: bool = False) -> tuple[bool, str]:
    text = path.read_text(encoding="utf-8")
    replacements = 0

    def replacer(match: re.Match[str]) -> str:
        nonlocal replacements
        agent = match.group("agent")
        prompt = match.group("prompt")
        model = resolve_model(agent, prompt, models)
        if model is None or match.group("model") == model:
            return match.group(0)
        replacements += 1
        return f'{match.group("prefix")}{model}{match.group("suffix")}'

    updated = TASK_BLOCK_RE.sub(replacer, text)

    if replace_constraint:
        updated = updated.replace(
            '3. **All subagents should use opus model for complex tasks**',
            '3. **All subagent models are synced from `~/.config/opencode/oh-my-openagent.jsonc` via `categories.trellis-*`**',
        )

    changed = updated != text
    if changed:
        _ = path.write_text(updated, encoding="utf-8")
    return changed, f"{path}: updated {replacements} task model reference(s)"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Sync active Trellis OpenCode model selections from oh-my-openagent.jsonc"
    )
    _ = parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_OMO_CONFIG_PATH,
        help="Path to oh-my-openagent.jsonc",
    )
    _ = parser.add_argument(
        "--project-root",
        type=Path,
        default=DEFAULT_PROJECT_ROOT,
        help="Project root containing .opencode/ and .trellis/",
    )
    _ = parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report intended changes without writing files",
    )
    args = parser.parse_args()

    project_root = cast(Path, args.project_root).resolve()
    config_path = cast(Path, args.config).resolve()
    target_files = build_target_files(project_root)

    models = load_trellis_models(config_path)
    print("Resolved Trellis subagent models:")
    for agent_name in ["research", "implement", "check", "debug", "finish"]:
        print(f"- {agent_name}: {models[agent_name]}")

    if cast(bool, args.dry_run):
        print("\nDry run mode - no files written.")
        return 0

    for path, options in target_files.items():
        changed, message = sync_file(path, models, **options)
        status = "changed" if changed else "unchanged"
        print(f"[{status}] {message}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
