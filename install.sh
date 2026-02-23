#!/bin/bash
# OpenCode Skills Setup Script
# This script sets up a new OpenCode installation with all 146 scientific skills

set -e

echo "=== OpenCode Skills Setup ==="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if opencode is installed
if ! command -v opencode &> /dev/null; then
    echo -e "${RED}Error: OpenCode is not installed${NC}"
    echo "Please install OpenCode first: https://opencode.ai"
    exit 1
fi

echo -e "${GREEN}✓ OpenCode found${NC}"

# Create skills directory
SKILLS_DIR="$HOME/.config/opencode/skills"
mkdir -p "$SKILLS_DIR"

echo "Setting up skills in: $SKILLS_DIR"
echo ""

# Clone this repository if not already present
REPO_DIR="$HOME/.config/opencode/skills-repo"
if [ ! -d "$REPO_DIR/.git" ]; then
    echo "Cloning skills repository..."
    git clone https://github.com/1508324011/Zhui-s-opencode-skills.git "$REPO_DIR"
else
    echo "Updating skills repository..."
    cd "$REPO_DIR" && git pull
fi

echo -e "${GREEN}✓ Repository ready${NC}"
echo ""

# Method 1: Direct copy of all skills (recommended for offline use)
echo "Method 1: Installing skills directly..."
if [ -d "$REPO_DIR/skills" ]; then
    cp -r "$REPO_DIR/skills"/* "$SKILLS_DIR/"
    echo -e "${GREEN}✓ Copied all skills to $SKILLS_DIR${NC}"
else
    echo -e "${YELLOW}⚠ Skills directory not found in repo${NC}"
    echo "Using Method 2 instead..."
fi

# Method 2: Setup superpowers plugin (optional, for extra features)
echo ""
echo "Method 2: Setting up Superpowers plugin (optional)..."
SUPERPOWERS_DIR="$HOME/.config/opencode/superpowers"

if [ ! -d "$SUPERPOWERS_DIR/.git" ]; then
    echo "Cloning superpowers repository..."
    git clone --depth 1 https://github.com/obra/superpowers.git "$SUPERPOWERS_DIR" 2>&1 || {
        echo -e "${YELLOW}⚠ Could not clone superpowers (optional)${NC}"
    }
else
    echo "Updating superpowers..."
    cd "$SUPERPOWERS_DIR" && git pull 2>&1 || true
fi

# Setup superpowers symlinks if superpowers is installed
if [ -d "$SUPERPOWERS_DIR/.opencode/plugins" ]; then
    mkdir -p "$HOME/.config/opencode/plugins"
    
    # Remove old incorrect plugin file if exists
    if [ -f "$HOME/.config/opencode/plugins/superpower.js" ]; then
        echo -e "${YELLOW}Removing old plugin file (superpower.js)${NC}"
        rm -f "$HOME/.config/opencode/plugins/superpower.js"
    fi
    
    # Create correct symlinks
    ln -sf "$SUPERPOWERS_DIR/.opencode/plugins/superpowers.js" "$HOME/.config/opencode/plugins/superpowers.js" 2>/dev/null || true
    ln -sf "$SUPERPOWERS_DIR/skills" "$SKILLS_DIR/superpowers" 2>/dev/null || true
    echo -e "${GREEN}✓ Superpowers plugin linked${NC}"
    echo -e "${YELLOW}Note: Plugin filename is 'superpowers.js' (with 's')${NC}"
fi

# Count installed skills
SKILL_COUNT=$(find "$SKILLS_DIR" -name "SKILL.md" 2>/dev/null | wc -l)
echo ""
echo -e "${GREEN}=== Setup Complete! ===${NC}"
echo "Skills installed: $SKILL_COUNT"
echo ""
echo "Restart OpenCode to use the new skills."
echo ""
echo "To verify, ask OpenCode: 'do you have superpowers?'"
