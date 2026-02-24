#!/bin/bash
# OpenCode Skills Setup Script
# 此脚本设置新的 OpenCode 安装，包含所有科学技能和 Oh My OpenCode 配置

set -e

echo "=== OpenCode Skills Setup ==="
echo ""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # 无颜色

# 检查 opencode 是否已安装
if ! command -v opencode &> /dev/null; then
    echo -e "${RED}错误：OpenCode 未安装${NC}"
    echo "请先安装 OpenCode: https://opencode.ai"
    exit 1
fi

echo -e "${GREEN}✓ 找到 OpenCode${NC}"
OPENCODE_VERSION=$(opencode --version)
echo "  版本：$OPENCODE_VERSION"

# 创建配置目录
CONFIG_DIR="$HOME/.config/opencode"
SKILLS_DIR="$CONFIG_DIR/skills"
mkdir -p "$CONFIG_DIR" "$SKILLS_DIR"

echo ""
echo "配置目录：$CONFIG_DIR"
echo "技能目录：$SKILLS_DIR"
echo ""

# 克隆仓库（如果不存在）
REPO_DIR="$CONFIG_DIR/skills-repo"
if [ ! -d "$REPO_DIR/.git" ]; then
    echo "克隆技能仓库..."
    git clone https://github.com/1508324011/Zhui-s-opencode-skills.git "$REPO_DIR"
else
    echo "更新技能仓库..."
    cd "$REPO_DIR" && git pull
fi

echo -e "${GREEN}✓ 仓库就绪${NC}"
echo ""

# 安装 Oh My OpenCode
echo -e "${BLUE}=== 安装 Oh My OpenCode ===${NC}"
if npm list -g oh-my-opencode &> /dev/null; then
    echo -e "${GREEN}✓ oh-my-opencode 已安装${NC}"
else
    echo "安装 oh-my-opencode..."
    npm install -g oh-my-opencode
    echo -e "${GREEN}✓ oh-my-opencode 安装完成${NC}"
fi

# 配置 Oh My OpenCode
echo ""
echo "配置 Oh My OpenCode..."
if [ -f "$REPO_DIR/config/oh-my-opencode.json" ]; then
    cp "$REPO_DIR/config/oh-my-opencode.json" "$CONFIG_DIR/oh-my-opencode.json"
    echo -e "${GREEN}✓ Oh My OpenCode 配置完成${NC}"
fi

# 配置提供商（示例文件）
echo ""
echo -e "${YELLOW}=== 配置 AI 提供商 ===${NC}"
if [ -f "$REPO_DIR/config/opencode.json.example" ]; then
    if [ ! -f "$CONFIG_DIR/opencode.json" ]; then
        echo "复制提供商配置示例..."
        cp "$REPO_DIR/config/opencode.json.example" "$CONFIG_DIR/opencode.json"
        echo -e "${GREEN}✓ 配置示例已复制${NC}"
        echo -e "${YELLOW}⚠ 请编辑 ~/.config/opencode/opencode.json 填入你的 API Key${NC}"
    else
        echo -e "${GREEN}✓ opencode.json 已存在${NC}"
    fi
fi

# 复制所有技能
echo ""
echo "安装技能..."
if [ -d "$REPO_DIR/skills" ]; then
    cp -r "$REPO_DIR/skills"/* "$SKILLS_DIR/"
    SKILL_COUNT=$(find "$SKILLS_DIR" -name "SKILL.md" 2>/dev/null | wc -l)
    echo -e "${GREEN}✓ 已安装 $SKILL_COUNT 个技能${NC}"
else
    echo -e "${YELLOW}⚠ 技能目录未找到${NC}"
fi

# 安装 Superpowers 插件（可选）
echo ""
echo -e "${BLUE}=== 安装 Superpowers 插件（可选）${NC}"
SUPERPOWERS_DIR="$CONFIG_DIR/superpowers"

if [ ! -d "$SUPERPOWERS_DIR/.git" ]; then
    echo "克隆 superpowers 仓库..."
    git clone --depth 1 https://github.com/obra/superpowers.git "$SUPERPOWERS_DIR" 2>&1 || {
        echo -e "${YELLOW}⚠ 无法克隆 superpowers（可选组件）${NC}"
    }
else
    echo "更新 superpowers..."
    cd "$SUPERPOWERS_DIR" && git pull 2>&1 || true
fi

# 设置 superpowers 符号链接
if [ -d "$SUPERPOWERS_DIR/.opencode/plugins" ]; then
    mkdir -p "$CONFIG_DIR/plugins"
    
    # 删除旧的正确文件名
    if [ -f "$CONFIG_DIR/plugins/superpower.js" ]; then
        echo -e "${YELLOW}删除旧的插件文件 (superpower.js)${NC}"
        rm -f "$CONFIG_DIR/plugins/superpower.js"
    fi
    
    # 创建正确的符号链接
    ln -sf "$SUPERPOWERS_DIR/.opencode/plugins/superpowers.js" "$CONFIG_DIR/plugins/superpowers.js" 2>/dev/null || true
    ln -sf "$SUPERPOWERS_DIR/skills" "$SKILLS_DIR/superpowers" 2>/dev/null || true
    echo -e "${GREEN}✓ Superpowers 插件已链接${NC}"
    echo -e "${YELLOW}注意：插件文件名为 'superpowers.js'（带 's'）${NC}"
fi

# 创建 .gitignore（如果不存在）
if [ ! -f "$CONFIG_DIR/.gitignore" ] && [ -f "$REPO_DIR/config/.gitignore" ]; then
    cp "$REPO_DIR/config/.gitignore" "$CONFIG_DIR/.gitignore"
    echo -e "${GREEN}✓ .gitignore 已创建${NC}"
fi

# 完成
echo ""
echo -e "${GREEN}================================${NC}"
echo -e "${GREEN}     安装完成！${NC}"
echo -e "${GREEN}================================${NC}"
echo ""
echo "📁 配置位置：$CONFIG_DIR"
echo "📚 技能数量：$SKILL_COUNT"
echo ""
echo "⚠️  下一步操作："
echo "1. 编辑 ~/.config/opencode/opencode.json 填入你的 API Key"
echo "2. 运行 opencode auth login 配置 Kimi API（可选）"
echo "3. 重启 OpenCode 开始使用"
echo ""
echo -e "${BLUE}🚀 提示：在提示词中使用 'ultrawork' 或 'ulw' 激活全自动模式！${NC}"
echo ""
