#!/bin/bash
# Installation script for Zhui's OpenCode Skills Collection
# This script sets up both scientific skills and superpowers plugin

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Zhui's OpenCode Skills Installer${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Check if OpenCode is installed
if ! command -v opencode &> /dev/null; then
    echo -e "${YELLOW}Warning: opencode command not found in PATH${NC}"
    echo "Please ensure OpenCode is installed: https://opencode.ai"
    echo ""
fi

# Configuration
REPO_URL="https://github.com/1508324011/Zhui-s-opencode-skills.git"
SUPERPOWERS_URL="https://github.com/obra/superpowers.git"
INSTALL_DIR="${HOME}/zhui-opencode-skills"
OPENCODE_CONFIG="${HOME}/.config/opencode"

echo -e "${YELLOW}Step 1: Cloning skills repository...${NC}"
if [ -d "$INSTALL_DIR" ]; then
    echo "Directory exists, updating..."
    cd "$INSTALL_DIR"
    git pull origin main
else
    git clone "$REPO_URL" "$INSTALL_DIR"
fi
echo -e "${GREEN}✓ Skills repository ready${NC}"
echo ""

echo -e "${YELLOW}Step 2: Installing scientific skills...${NC}"
mkdir -p "${OPENCODE_CONFIG}/skills"
cp -r "${INSTALL_DIR}/skills/"* "${OPENCODE_CONFIG}/skills/"
SKILL_COUNT=$(ls -1 "${OPENCODE_CONFIG}/skills" | wc -l)
echo -e "${GREEN}✓ Installed ${SKILL_COUNT} scientific skills${NC}"
echo ""

echo -e "${YELLOW}Step 3: Installing Superpowers plugin...${NC}"
if [ -d "${OPENCODE_CONFIG}/superpowers" ]; then
    echo "Superpowers exists, updating..."
    cd "${OPENCODE_CONFIG}/superpowers"
    git pull
else
    git clone "$SUPERPOWERS_URL" "${OPENCODE_CONFIG}/superpowers"
fi
echo -e "${GREEN}✓ Superpowers plugin ready${NC}"
echo ""

echo -e "${YELLOW}Step 4: Creating plugin symlink...${NC}"
mkdir -p "${OPENCODE_CONFIG}/plugins"
if [ -L "${OPENCODE_CONFIG}/plugins/superpowers.js" ]; then
    rm "${OPENCODE_CONFIG}/plugins/superpowers.js"
fi
ln -s "${OPENCODE_CONFIG}/superpowers/.opencode/plugins/superpowers.js" \
      "${OPENCODE_CONFIG}/plugins/superpowers.js"
echo -e "${GREEN}✓ Plugin symlink created${NC}"
echo ""

echo -e "${YELLOW}Step 5: Creating skills symlink...${NC}"
if [ -L "${OPENCODE_CONFIG}/skills/superpowers" ]; then
    rm "${OPENCODE_CONFIG}/skills/superpowers"
fi
ln -s "${OPENCODE_CONFIG}/superpowers/skills" \
      "${OPENCODE_CONFIG}/skills/superpowers"
echo -e "${GREEN}✓ Skills symlink created${NC}"
echo ""

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Installation Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "1. Restart OpenCode"
echo "2. Test with: 'do you have superpowers?'"
echo ""
echo -e "${YELLOW}Installed:${NC}"
echo "  • ${SKILL_COUNT} scientific skills"
echo "  • 12+ superpowers capabilities"
echo "  • Plugin integration"
echo ""
echo -e "${YELLOW}To update in the future:${NC}"
echo "  cd ${INSTALL_DIR} && git pull"
echo "  cd ${OPENCODE_CONFIG}/superpowers && git pull"
echo ""
echo -e "${GREEN}Happy researching!${NC}"
