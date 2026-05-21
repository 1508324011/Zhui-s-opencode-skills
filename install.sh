#!/bin/bash
# =============================================================================
# Zhui's OpenCode Skills Installer
# =============================================================================
# 双发布线安装器：
#   - superpowers（默认，兼容旧行为）
#   - trellis
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${BLUE}ℹ️  $1${NC}"; }
log_success() { echo -e "${GREEN}✅ $1${NC}"; }
log_warning() { echo -e "${YELLOW}⚠️  $1${NC}"; }
log_error() { echo -e "${RED}❌ $1${NC}"; }

usage() {
  cat <<'EOF'
用法：
  ./install.sh [superpowers|trellis] [--skip-npm-install] [--help]

说明：
  superpowers   安装 OMO + DCP + 官方 Superpowers 工作流（默认）
  trellis       安装 OpenAgent + DCP + Trellis skeleton

参数：
  --skip-npm-install   只验证/复制文件，不执行 npm install
  --help, -h           显示帮助

环境变量覆盖：
  OPENCODE_TARGET_DIR              默认 ~/.config/opencode
  OPENCODE_HOME_OPENCODE_DIR       默认 ~/.opencode
  OPENCODE_TRELLIS_TEMPLATE_DIR    默认 ~/.config/opencode/trellis-project-template/.trellis
  OPENCODE_SUPERPOWERS_REPO_DIR    默认 ~/.config/opencode/superpowers
  OPENCODE_SUPERPOWERS_GIT_URL     默认 https://github.com/obra/superpowers.git
  OPENCODE_SKIP_OPENCODE_CHECK=1   跳过 opencode 命令存在性检查

兼容性：
  不带模式参数时，默认按 superpowers 线安装，以兼容旧的 ./install.sh 用法。
EOF
}

MODE="superpowers"
SKIP_NPM_INSTALL=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    superpowers|trellis)
      MODE="$1"
      shift
      ;;
    --skip-npm-install)
      SKIP_NPM_INSTALL=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      log_error "未知参数：$1"
      echo
      usage
      exit 1
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="${OPENCODE_TARGET_DIR:-$HOME/.config/opencode}"
HOME_OPENCODE_DIR="${OPENCODE_HOME_OPENCODE_DIR:-$HOME/.opencode}"
TRELLIS_TEMPLATE_DIR="${OPENCODE_TRELLIS_TEMPLATE_DIR:-$TARGET_DIR/trellis-project-template/.trellis}"
SUPERPOWERS_REPO_DIR="${OPENCODE_SUPERPOWERS_REPO_DIR:-$TARGET_DIR/superpowers}"
SUPERPOWERS_GIT_URL="${OPENCODE_SUPERPOWERS_GIT_URL:-https://github.com/obra/superpowers.git}"
SKIP_OPENCODE_CHECK="${OPENCODE_SKIP_OPENCODE_CHECK:-0}"
TIMESTAMP="$(date +%Y%m%d-%H%M%S)"

echo "=========================================="
echo "🚀 Zhui's OpenCode 安装程序"
echo "=========================================="
echo
log_info "安装线路：$MODE"
log_info "源目录：$SCRIPT_DIR"
log_info "目标目录：$TARGET_DIR"
if [ "$MODE" = "trellis" ]; then
  log_info "~/.opencode 目标：$HOME_OPENCODE_DIR"
  log_info "Trellis 模板目录：$TRELLIS_TEMPLATE_DIR"
else
  log_info "Superpowers 仓库目录：$SUPERPOWERS_REPO_DIR"
fi
echo

ensure_opencode_installed() {
  if [ "$SKIP_OPENCODE_CHECK" = "1" ]; then
    log_warning "已跳过 opencode 命令检查（OPENCODE_SKIP_OPENCODE_CHECK=1）"
    return
  fi

  if ! command -v opencode >/dev/null 2>&1; then
    log_error "未检测到 OpenCode"
    echo
    echo "请先安装 OpenCode:"
    echo "  npm install -g opencode-ai"
    echo "  或访问：https://opencode.ai"
    echo
    exit 1
  fi

  log_success "检测到 OpenCode: $(command -v opencode)"
}

install_npm_package() {
  local package_name="$1"

  if [ "$SKIP_NPM_INSTALL" = "1" ]; then
    log_warning "跳过 npm 安装：$package_name"
    return
  fi

  log_info "安装 npm 包：$package_name"
  npm install -g "$package_name"
}

copy_if_missing() {
  local src="$1"
  local dest="$2"

  if [ -e "$dest" ] || [ -L "$dest" ]; then
    log_warning "已存在，保留原文件：$dest"
    return
  fi

  mkdir -p "$(dirname "$dest")"
  cp "$src" "$dest"
  log_success "已创建：$dest"
}

copy_with_backup() {
  local src="$1"
  local dest="$2"

  mkdir -p "$(dirname "$dest")"

  if [ -e "$dest" ] || [ -L "$dest" ]; then
    if cmp -s "$src" "$dest" 2>/dev/null; then
      log_info "无需更新：$dest"
      return
    fi

    local backup_path="${dest}.bak-${TIMESTAMP}"
    mv "$dest" "$backup_path"
    log_warning "已备份：$dest -> $backup_path"
  fi

  cp "$src" "$dest"
  log_success "已安装：$dest"
}

replace_dir_with_backup() {
  local src_dir="$1"
  local dest_dir="$2"

  if [ -e "$dest_dir" ] || [ -L "$dest_dir" ]; then
    local backup_path="${dest_dir}.bak-${TIMESTAMP}"
    mv "$dest_dir" "$backup_path"
    log_warning "已备份目录：$dest_dir -> $backup_path"
  fi

  mkdir -p "$dest_dir"
  cp -R "$src_dir/." "$dest_dir/"
  log_success "已同步目录：$dest_dir"
}

replace_with_symlink() {
  local source_path="$1"
  local target_path="$2"

  mkdir -p "$(dirname "$target_path")"

  if [ -L "$target_path" ]; then
    rm -f "$target_path"
  elif [ -e "$target_path" ]; then
    local backup_path="${target_path}.bak-${TIMESTAMP}"
    mv "$target_path" "$backup_path"
    log_warning "已备份：$target_path -> $backup_path"
  fi

  ln -s "$source_path" "$target_path"
  log_success "已创建符号链接：$target_path -> $source_path"
}

write_gitignore_if_missing() {
  local gitignore_path="$TARGET_DIR/.gitignore"

  if [ -e "$gitignore_path" ]; then
    log_warning "已存在，保留原 .gitignore：$gitignore_path"
    return
  fi

  mkdir -p "$(dirname "$gitignore_path")"

  cat > "$gitignore_path" <<'EOF'
# API Keys & Credentials
opencode.json
auth.json
*.env
*.env.*
.env.local
.env.*.local
*.local
*.local.*
*api*key*
*apikey*
*secret*
*credentials*
*password*

# Certificates & Keys
*.pem
*.key
*.p12
*.pfx
*.jks

# User-specific Config
*.local.jsonc
*.local.json
config.local.*

# IDE & Editor
.idea/
.vscode/
*.swp
*.swo
*~
.DS_Store

# Logs & Temp
*.log
logs/
tmp/
temp/
.cache/
.claude_resources.json
EOF

  log_success "已创建 .gitignore"
}

install_shared_skills() {
  mkdir -p "$TARGET_DIR/skills"

  if command -v rsync >/dev/null 2>&1; then
    rsync -a \
      --exclude 'superpowers' \
      --exclude '__pycache__/' \
      --exclude '*.pyc' \
      --exclude '*.pyo' \
      --exclude '.DS_Store' \
      "$SCRIPT_DIR/skills/" "$TARGET_DIR/skills/"
  else
    local entry
    for entry in "$SCRIPT_DIR"/skills/*; do
      local name
      name="$(basename "$entry")"
      [ "$name" = "superpowers" ] && continue
      if [ -d "$entry" ] && [ ! -L "$entry" ]; then
        mkdir -p "$TARGET_DIR/skills/$name"
        cp -R "$entry/." "$TARGET_DIR/skills/$name/"
      else
        cp -R "$entry" "$TARGET_DIR/skills/$name"
      fi
    done
  fi

  local skill_count
  skill_count="$(find "$TARGET_DIR/skills" -mindepth 1 -maxdepth 1 | wc -l | tr -d ' ')"
  log_success "已同步共享 skills：$skill_count 个顶层条目"
}

ensure_superpowers_repo() {
  if [ -d "$SUPERPOWERS_REPO_DIR/.git" ]; then
    log_info "使用已有 Superpowers 仓库：$SUPERPOWERS_REPO_DIR"
    return
  fi

  if [ -e "$SUPERPOWERS_REPO_DIR" ] || [ -L "$SUPERPOWERS_REPO_DIR" ]; then
    local backup_path="${SUPERPOWERS_REPO_DIR}.bak-${TIMESTAMP}"
    mv "$SUPERPOWERS_REPO_DIR" "$backup_path"
    log_warning "已备份现有路径：$SUPERPOWERS_REPO_DIR -> $backup_path"
  fi

  log_info "克隆官方 Superpowers 仓库..."
  git clone --depth 1 "$SUPERPOWERS_GIT_URL" "$SUPERPOWERS_REPO_DIR"
  log_success "已克隆 Superpowers 仓库"
}

install_superpowers_symlinks() {
  replace_with_symlink "$SUPERPOWERS_REPO_DIR/.opencode/plugins/superpowers.js" "$TARGET_DIR/plugins/superpowers.js"
  replace_with_symlink "$SUPERPOWERS_REPO_DIR/skills" "$TARGET_DIR/skills/superpowers"
}

install_home_opencode_dependencies() {
  if [ "$SKIP_NPM_INSTALL" = "1" ]; then
    log_warning "跳过 ~/.opencode 依赖安装"
    return
  fi

  if [ ! -f "$HOME_OPENCODE_DIR/package.json" ]; then
    log_warning "未找到 ~/.opencode/package.json，跳过本地依赖安装"
    return
  fi

  log_info "安装 ~/.opencode 本地依赖..."
  (
    cd "$HOME_OPENCODE_DIR"
    npm install
  )
  log_success "~/.opencode 本地依赖安装完成"
}

install_superpowers_line() {
  log_info "开始安装 Superpowers 线路..."

  mkdir -p "$TARGET_DIR" "$TARGET_DIR/skills" "$TARGET_DIR/plugins"
  install_npm_package "oh-my-opencode"
  install_npm_package "@tarquinen/opencode-dcp"

  write_gitignore_if_missing
  install_shared_skills
  copy_with_backup "$SCRIPT_DIR/config/oh-my-opencode.json" "$TARGET_DIR/oh-my-opencode.json"
  copy_with_backup "$SCRIPT_DIR/dcp.jsonc" "$TARGET_DIR/dcp.jsonc"
  copy_with_backup "$SCRIPT_DIR/config/opencode.json.example" "$TARGET_DIR/opencode.json"
  ensure_superpowers_repo
  install_superpowers_symlinks

  log_success "Superpowers 线路安装完成"
}

install_trellis_line() {
  log_info "开始安装 Trellis 线路..."

  mkdir -p "$TARGET_DIR"
  install_npm_package "oh-my-openagent"
  install_npm_package "@tarquinen/opencode-dcp"

  write_gitignore_if_missing
  install_shared_skills
  copy_with_backup "$SCRIPT_DIR/trellis/config/oh-my-openagent.jsonc" "$TARGET_DIR/oh-my-openagent.jsonc"
  copy_with_backup "$SCRIPT_DIR/trellis/config/dcp.jsonc" "$TARGET_DIR/dcp.jsonc"
  copy_with_backup "$SCRIPT_DIR/trellis/config/opencode.json.example" "$TARGET_DIR/opencode.json"

  replace_dir_with_backup "$SCRIPT_DIR/trellis/home/.opencode" "$HOME_OPENCODE_DIR"
  install_home_opencode_dependencies

  replace_dir_with_backup "$SCRIPT_DIR/trellis/project/.trellis" "$TRELLIS_TEMPLATE_DIR"

  log_success "Trellis 线路安装完成"
}

print_summary() {
  echo
  echo "=========================================="
  log_success "安装完成！"
  echo "=========================================="
  echo
  echo "安装线路：$MODE"
  echo "目标目录：$TARGET_DIR"
  if [ "$MODE" = "trellis" ]; then
    echo "~/.opencode：$HOME_OPENCODE_DIR"
    echo "Trellis 模板：$TRELLIS_TEMPLATE_DIR"
  else
    echo "Superpowers 仓库：$SUPERPOWERS_REPO_DIR"
  fi
  echo
  echo "后续步骤："
  echo "1. 填写 API Key：nano $TARGET_DIR/opencode.json"
  if [ "$MODE" = "superpowers" ]; then
    echo "2. 重启 OpenCode：opencode"
    echo "3. 在 OpenCode 中输入 /models 验证模型加载"
  else
    echo "2. 将 $TRELLIS_TEMPLATE_DIR 复制到项目根作为 .trellis 模板"
    echo "3. 重启 OpenCode：opencode"
    echo "4. 在 OpenCode 中输入 /models 验证模型加载"
  fi
  echo
}

ensure_opencode_installed
echo

case "$MODE" in
  superpowers)
    install_superpowers_line
    ;;
  trellis)
    install_trellis_line
    ;;
esac

print_summary
