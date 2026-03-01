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
mkdir -p "$TARGET_DIR/superpowers"

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

# opencode-mem
log_info "安装 opencode-mem..."
npm install -g opencode-mem 2>/dev/null || npm install opencode-mem 2>/dev/null || log_warning "opencode-mem 安装失败，请手动安装"

log_success "核心插件安装完成"
echo ""

# =============================================================================
# 1. 安装 Skills
log_info "创建配置目录..."
mkdir -p "$TARGET_DIR/skills"
mkdir -p "$TARGET_DIR/plugins"
mkdir -p "$TARGET_DIR/superpowers"

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
# 2. 安装 Superpowers 插件
# =============================================================================
log_info "安装 Superpowers 插件..."

if [ -f "$SCRIPT_DIR/plugins/superpowers.js" ]; then
    cp "$SCRIPT_DIR/plugins/superpowers.js" "$TARGET_DIR/plugins/"
    log_success "已安装 superpowers.js"
elif [ -d "$SCRIPT_DIR/superpowers" ]; then
    # 如果没有 plugins/superpowers.js，尝试安装 superpowers 仓库
    log_info "克隆 Superpowers 仓库..."
    
    if [ ! -d "$TARGET_DIR/superpowers/.git" ]; then
        git clone https://github.com/obra/superpowers.git "$TARGET_DIR/superpowers" 2>/dev/null || {
            log_warning "无法克隆 superpowers，跳过..."
        }
    fi
    
    if [ -d "$TARGET_DIR/superpowers/.git" ]; then
        # 创建 symlinks
        ln -sf "$TARGET_DIR/superpowers/.opencode/plugins/superpowers.js" "$TARGET_DIR/plugins/superpowers.js" 2>/dev/null || true
        ln -sf "$TARGET_DIR/superpowers/skills" "$TARGET_DIR/skills/superpowers" 2>/dev/null || true
        log_success "已安装 Superpowers 插件和技能"
    fi
else
    log_warning "未找到 Superpowers，跳过..."
fi

echo ""

# =============================================================================
# 3. 安装配置文件
# =============================================================================
log_info "安装配置文件..."

# oh-my-opencode.json (可以安全覆盖)
if [ -f "$SCRIPT_DIR/config/oh-my-opencode.json" ]; then
    cp "$SCRIPT_DIR/config/oh-my-opencode.json" "$TARGET_DIR/"
    log_success "已安装 oh-my-opencode.json"
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
# 4. 创建 .gitignore
# =============================================================================
log_info "创建 .gitignore..."

cat > "$TARGET_DIR/.gitignore" << 'EOF'
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
EOF

log_success "已创建 .gitignore"

echo ""

# =============================================================================
# 5. 验证安装
# =============================================================================
log_info "验证安装..."

SKILL_COUNT=$(ls -1 "$TARGET_DIR/skills/" | wc -l)
log_success "Skills 数量：$SKILL_COUNT"

if [ -f "$TARGET_DIR/plugins/superpowers.js" ] || [ -L "$TARGET_DIR/plugins/superpowers.js" ]; then
    log_success "Superpowers 插件：已安装"
else
    log_warning "Superpowers 插件：未找到"
fi

if [ -f "$TARGET_DIR/oh-my-opencode.json" ]; then
    log_success "oh-my-opencode.json: 已安装"
fi

echo ""

# =============================================================================
# 6. 后续步骤
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
echo "2. 重启 OpenCode:"
echo "   opencode"
echo ""
echo "3. 验证安装:"
echo "   在 OpenCode 中输入 /models 查看可用模型"
echo ""
echo "=========================================="
echo ""
