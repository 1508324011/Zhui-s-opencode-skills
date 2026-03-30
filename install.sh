#!/bin/bash
# =============================================================================
# Zhui's OpenCode Skills Installer
# =============================================================================
# 一键安装所有 Skills、插件和配置文件到新设备
# =============================================================================

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${BLUE}ℹ️  $1${NC}"; }
log_success() { echo -e "${GREEN}✅ $1${NC}"; }
log_warning() { echo -e "${YELLOW}⚠️  $1${NC}"; }
log_error() { echo -e "${RED}❌ $1${NC}"; }

echo "=========================================="
echo "🚀 Zhui's OpenCode 安装程序"
echo "=========================================="
echo ""

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="$HOME/.config/opencode"

log_info "源目录：$SCRIPT_DIR"
log_info "目标目录：$TARGET_DIR"
echo ""

# 检查 OpenCode 是否已安装
if ! command -v opencode &> /dev/null; then
    log_error "未检测到 OpenCode"
    echo ""
    echo "请先安装 OpenCode:"
    echo "  npm install -g opencode-ai"
    echo "  或访问：https://opencode.ai"
    echo ""
    exit 1
fi

log_success "检测到 OpenCode: $(which opencode)"
echo ""

# 创建目标目录
log_info "创建配置目录..."
mkdir -p "$TARGET_DIR/skills"
mkdir -p "$TARGET_DIR/plugins"

for legacy_path in "$TARGET_DIR/plugins/superpowers.js" "$TARGET_DIR/skills/superpowers" "$TARGET_DIR/superpowers"; do
    if [ -e "$legacy_path" ] || [ -L "$legacy_path" ]; then
        rm -rf "$legacy_path"
        log_success "已清理旧的 Superpowers 残留：$legacy_path"
    fi
done

# =============================================================================
# 0. 安装核心 npm 插件
# =============================================================================
log_info "安装核心 npm 插件..."

# oh-my-opencode
log_info "安装 oh-my-opencode..."
npm install -g oh-my-opencode 2>/dev/null || npm install oh-my-opencode 2>/dev/null || log_warning "oh-my-opencode 安装失败，请手动安装"

# @tarquinen/opencode-dcp
log_info "安装 @tarquinen/opencode-dcp..."
npm install -g @tarquinen/opencode-dcp 2>/dev/null || npm install @tarquinen/opencode-dcp 2>/dev/null || log_warning "opencode-dcp 安装失败，请手动安装"

log_success "核心插件安装完成"
echo ""

# =============================================================================
# 1. 安装 Skills
# =============================================================================
log_info "安装 Skills (146 个科学和技术技能)..."

if [ -d "$SCRIPT_DIR/skills" ]; then
    # 复制所有 skills
    cp -r "$SCRIPT_DIR/skills/"* "$TARGET_DIR/skills/"
    
    SKILL_COUNT=$(ls -1 "$TARGET_DIR/skills/" | wc -l)
    log_success "已安装 $SKILL_COUNT 个 skills"
else
    log_error "未找到 skills 目录"
fi

echo ""

# =============================================================================
# 2. 安装配置文件
# =============================================================================
log_info "安装配置文件..."

# oh-my-opencode.json (可以安全覆盖)
if [ -f "$SCRIPT_DIR/config/oh-my-opencode.json" ]; then
    cp "$SCRIPT_DIR/config/oh-my-opencode.json" "$TARGET_DIR/"
    log_success "已安装 oh-my-opencode.json"
fi

if [ -f "$SCRIPT_DIR/config/oh-my-openagent.jsonc.example" ] && [ ! -f "$TARGET_DIR/oh-my-openagent.jsonc" ]; then
    cp "$SCRIPT_DIR/config/oh-my-openagent.jsonc.example" "$TARGET_DIR/oh-my-openagent.jsonc"
    log_success "已创建 oh-my-openagent.jsonc (Trellis 模型同步模板)"
fi

# opencode.json (仅当不存在时创建)
if [ -f "$SCRIPT_DIR/config/opencode.json.example" ] && [ ! -f "$TARGET_DIR/opencode.json" ]; then
    cp "$SCRIPT_DIR/config/opencode.json.example" "$TARGET_DIR/opencode.json"
    log_success "已创建 opencode.json (请填入 API Key)"
fi

# dcp.jsonc (如果存在)
if [ -f "$SCRIPT_DIR/dcp.jsonc" ] && [ ! -f "$TARGET_DIR/dcp.jsonc" ]; then
    cp "$SCRIPT_DIR/dcp.jsonc" "$TARGET_DIR/"
    log_success "已安装 dcp.jsonc"
fi

echo ""

# =============================================================================
# 3. 创建 .gitignore
# =============================================================================
log_info "创建 .gitignore..."

cat > "$TARGET_DIR/.gitignore" << 'EOF'
# API Keys & Credentials
opencode.json
opencode-mem.jsonc
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
EOF

log_success "已创建 .gitignore"

echo ""

# =============================================================================
# 4. 验证安装
# =============================================================================
log_info "验证安装..."

SKILL_COUNT=$(ls -1 "$TARGET_DIR/skills/" | wc -l)
log_success "Skills 数量：$SKILL_COUNT"

if [ -f "$TARGET_DIR/dcp.jsonc" ]; then
    log_success "DCP 配置：已安装"
else
    log_warning "DCP 配置：未找到"
fi

if [ -f "$TARGET_DIR/oh-my-openagent.jsonc" ]; then
    log_success "Trellis 模型模板：已安装"
else
    log_warning "Trellis 模型模板：未找到"
fi

if [ -f "$TARGET_DIR/oh-my-opencode.json" ]; then
    log_success "oh-my-opencode.json: 已安装"
fi

echo ""

# =============================================================================
# 5. 后续步骤
# =============================================================================
echo "=========================================="
log_success "安装完成！"
echo "=========================================="
echo ""
echo "后续步骤："
echo ""
echo "1. 配置 API Key:"
echo "   nano $TARGET_DIR/opencode.json"
echo "   (将 YOUR_API_KEY_HERE 替换为你的真实 API Key)"
echo ""
echo "2. 如需调整 DCP 行为："
echo "   nano $TARGET_DIR/dcp.jsonc"
echo ""
echo "3. 如果你使用 Trellis："
echo "   nano $TARGET_DIR/oh-my-openagent.jsonc"
echo "   # 编辑 categories.trellis-*，然后在 Trellis 环境中运行 /trellis:sync-models"
echo ""
echo "4. 如果旧设备曾经启用过 Superpowers，请确认 opencode.json 中已移除 ./superpowers"
echo ""
echo "5. 重启 OpenCode:"
echo "   opencode"
echo ""
echo "6. 验证安装:"
echo "   在 OpenCode 中输入 /models 和 /dcp stats 查看配置是否生效"
echo ""
echo "=========================================="
echo ""
