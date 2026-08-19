# CODEX_PROJECT Target B/C Lessons Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reusable scientific-processing and figure-reproduction lessons from the completed GMSL/GIA/OBD work to targets B and C in `CODEX_PROJECT.md`.

**Architecture:** Insert one focused subsection directly below target B's Level-2 to Level-3 chain and one below target C's three comparison requirements. Keep the additions self-contained, operational, and independent of any out-of-bound code or output path.

**Tech Stack:** UTF-8 Markdown, PowerShell read-only checks, Git for the plan record.

## Global Constraints

- Modify only `D:\AAAA海平面变化\CODEX_PROJECT.md` target B and target C content.
- Do not change target A, target D, execution boundaries, roadmap, or existing acceptance thresholds.
- Do not reference or import any code, configuration, path, or output from excluded repositories.
- Do not invoke the L2toL3 skill.
- Preserve the existing Markdown hierarchy and UTF-8 encoding.
- Use direct project-manual Chinese; avoid promotional or generic summary language.

---

### Task 1: Add Target B Processing Constraints

**Files:**
- Modify: `CODEX_PROJECT.md:77-94`

**Interfaces:**
- Consumes: the existing target B Level-2 to Level-3 pipeline definition.
- Produces: a `#### 本次实践形成的处理约束` subsection immediately before target C.

- [ ] **Step 1: Record the pre-edit structure**

Run:

```powershell
rg -n "^### 目标 B|^### 目标 C|^#### 本次实践形成的处理约束" CODEX_PROJECT.md
```

Expected before editing: one target B heading, one target C heading, and no processing-constraints subsection.

- [ ] **Step 2: Insert the target B subsection with `apply_patch`**

Add nine numbered requirements covering:

```text
variable semantics and sign registry
GIA versus contemporary elastic loading versus OBD
signed upward-positive OBD on the budget side
shared Level-2 preprocessing for mass and OBD
separate observed and filled month masks
fill each reference center before averaging
separate main and sensitivity results
separate reference/adopted/computed trends
configuration and SHA-256 provenance
```

The OBD requirement must state this relation explicitly:

```text
SSH/GMSL = mass + steric + OBD
```

and must state that OBD is not added again to observed altimetry SSH/GMSL.

- [ ] **Step 3: Verify target B placement and prohibited references**

Run:

```powershell
rg -n -A 45 "^### 目标 B" CODEX_PROJECT.md
rg -n "reproduce_figure1|SaGEA_L2_to_L3_optimized|L2toL3 工作流" CODEX_PROJECT.md
```

Expected: the new subsection appears between targets B and C; the only excluded-repository references remain the pre-existing execution-boundary statements, not the new subsection.

---

### Task 2: Add Target C Reproduction and Acceptance Lessons

**Files:**
- Modify: `CODEX_PROJECT.md:95-102`

**Interfaces:**
- Consumes: the existing target C three-way comparison requirements.
- Produces: a `#### 本次实践形成的复现与验收经验` subsection immediately before target D.

- [ ] **Step 1: Insert the target C subsection with `apply_patch`**

Add eight numbered requirements covering:

```text
define each scientific curve before plotting
keep observed GMSL, mass, steric, and OBD closure diagnostics separate
separate scientific processing from display offsets
report source and smoothed valid-month counts
separate reference/adopted and recomputed trends with uncertainty method
write a complete figure bundle
verify equations, hashes, rerun identity, and rendered outputs
preserve old results when introducing a new scientific definition
```

The Figure 1 requirement must state that the mass-plus-steric curve excludes OBD and that OBD belongs in a separate closure diagnostic.

- [ ] **Step 2: Apply the Chinese humanization pass**

Read only the two new subsections and remove:

```text
empty transitions
promotional claims
repeated “此外/至关重要/确保” phrasing
formulaic three-part conclusions
```

Keep technical terms, equations, units, file-format names, and validation actions unchanged.

- [ ] **Step 3: Verify structure, content, and encoding**

Run:

```powershell
rg -n "^### 目标 [ABCD]|^#### 本次实践形成的" CODEX_PROJECT.md
rg -n "mass \+ steric \+ OBD|不得再次加到|质量海平面与比容海平面之和|SHA-256|PNG|PDF" CODEX_PROJECT.md
Get-Content -Raw -Encoding UTF8 CODEX_PROJECT.md | Out-Null
```

Expected: target order remains A, B, C, D; exactly two new fourth-level headings exist; all required phrases are present; UTF-8 reading succeeds.

- [ ] **Step 4: Inspect the final edited region**

Run:

```powershell
Get-Content CODEX_PROJECT.md | Select-Object -Skip 60 -First 120
```

Expected: the new requirements read as concise project instructions, with no change to target A or target D content.

- [ ] **Step 5: Report the edit without staging unrelated files**

Because `CODEX_PROJECT.md` is currently an untracked user file, do not stage or commit it automatically. Report the exact path and the two inserted subsection headings.
