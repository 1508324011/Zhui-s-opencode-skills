---
name: frontend-design
description: Create distinctive, production-grade frontend interfaces with high design quality. Use this skill when the user asks to build web components, pages, or applications. Generates creative, polished code that avoids generic AI aesthetics.
metadata:
    skill-author: Anthropic (adapted for local OpenCode installation)
    upstream-source: https://github.com/anthropics/claude-code/tree/main/plugins/frontend-design
---

# Frontend Design

This skill is adapted from Anthropic's `frontend-design` Claude Code plugin and installed as a native OpenCode skill. It guides creation of distinctive, production-grade frontend interfaces that avoid generic "AI slop" aesthetics.

Implement real working code with exceptional attention to aesthetic details and creative choices.

## When to Use This Skill

Use this skill when the user asks to build or refine:

- web components
- full pages
- frontend applications
- dashboards, landing pages, settings panels, marketing sites, or product UIs
- UI polish, visual direction, interaction quality, or frontend aesthetics

The user may include context about the purpose, audience, or technical constraints.

## Design Thinking

Before coding, understand the context and commit to a **bold** aesthetic direction:

- **Purpose**: What problem does this interface solve? Who uses it?
- **Tone**: Pick an extreme direction when appropriate: brutally minimal, maximalist chaos, retro-futuristic, organic/natural, luxury/refined, playful/toy-like, editorial/magazine, brutalist/raw, art deco/geometric, soft/pastel, industrial/utilitarian, and so on.
- **Constraints**: Respect technical requirements such as framework, performance, and accessibility.
- **Differentiation**: What makes this interface unforgettable? What is the one thing someone will remember?

**Critical**: Choose a clear conceptual direction and execute it with precision. Bold maximalism and refined minimalism can both work. The key is intentionality, not intensity.

Then implement working code (HTML/CSS/JS, React, Vue, and similar frontend stacks) that is:

- production-grade and functional
- visually striking and memorable
- cohesive with a clear aesthetic point of view
- meticulously refined in every detail

## Frontend Aesthetics Guidelines

Focus on the following areas.

### Typography

- Choose fonts that are beautiful, unique, and interesting.
- Avoid generic defaults like Arial, Inter, Roboto, and system fonts when the project allows stronger choices.
- Prefer distinctive display fonts paired with refined body fonts.

### Color and Theme

- Commit to a cohesive aesthetic.
- Use CSS variables or the equivalent theme system for consistency.
- Dominant colors with sharp accents usually outperform timid, evenly distributed palettes.

### Motion

- Use animation for impact and micro-interaction quality.
- Prefer CSS-only solutions for plain HTML where practical.
- Use the Motion library for React when it is already available or clearly appropriate.
- Focus on high-impact moments: orchestrated page loads, staggered reveals, scroll-triggered moments, and hover states that surprise.

### Spatial Composition

- Embrace unexpected layouts.
- Use asymmetry, overlap, diagonal flow, grid-breaking elements, generous negative space, or controlled density when the concept benefits from them.

### Backgrounds and Visual Details

- Create atmosphere and depth rather than defaulting to flat solid backgrounds.
- Use effects and textures that match the overall aesthetic: gradient meshes, noise textures, geometric patterns, layered transparency, dramatic shadows, decorative borders, custom cursors, and grain overlays.

## Avoid Generic AI Aesthetics

Never fall back to:

- overused font families and safe defaults without intent
- cliched color schemes, especially purple gradients on white backgrounds
- predictable layouts and component patterns
- cookie-cutter designs that ignore context
- repeatedly converging on the same aesthetic across different tasks

Interpret the request creatively and make unexpected choices that feel genuinely designed for the context. Vary light and dark themes, fonts, and visual systems. Do not converge on the same default choices across generations.

## Match Complexity to Vision

- Maximalist designs may need elaborate code, richer visuals, and more extensive animation.
- Minimalist or refined designs need restraint, precision, and careful control of spacing, typography, and subtle details.

Elegance comes from executing the chosen vision well.

## Final Principle

Do not hold back creatively, but keep the implementation grounded in real, working frontend code. Aim for interfaces that feel intentionally designed rather than statistically generated.
