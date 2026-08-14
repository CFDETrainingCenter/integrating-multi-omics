#!/usr/bin/env python3
"""Production text hygiene checker for the CFDE lung course package.

Two rule sets by file class (see build_instructions/WORKFLOW_POLISH_2026-08-10.md):

1. Code / config / machine-read data (.py, .yaml, .json, .tsv, .csv, .sh, and
   filenames): ASCII only (codepoint > 127 is a finding), with an exact-string
   allowlist for upstream values that must not be altered.
2. Learner-facing prose (.md, notebook markdown cells): flag typographic dashes,
   curly quotes, ellipsis, NBSP, arrows, and (optionally) Greek letters; also
   flag common AI-filler phrases for human review.

Exit code 1 if any hard findings (ASCII / typography) are reported.
Phrase findings are warnings and do not fail the process unless --strict-phrases.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Exact upstream strings that must remain unmodified (join keys / portal labels).
UPSTREAM_ALLOWLIST: frozenset[str] = frozenset(
    {
        "Interstitial Mφ perivascular",
        "Monocyte-derived Mφ",
        # Verbatim published title (Argelaguet et al., MOFA+, Genome Biol 2020).
        # A citation title is quoted, not authored, and must not be reworded.
        "comprehensive integration of multi-modal single-cell data",
        # Statistical sense of leverage (influence of an outlier on a fit),
        # not the corporate verb the phrase list targets.
        "outlier gene leverage",
    }
)

CODE_SUFFIXES = {".py", ".yaml", ".yml", ".json", ".tsv", ".csv", ".sh", ".txt", ".log", ".cfg", ".ini"}
PROSE_SUFFIXES = {".md"}

# Extensionless files that are still machine-read or shipped text (LICENSE, Makefile).
# Checked as code, so any codepoint above 127 is flagged.
EXTENSIONLESS_CHECKED = {"LICENSE-CODE", "Makefile", "Dockerfile"}

# LICENSE carries the complete unmodified CC BY 4.0 legal code, which contains two
# curly double quotes in the official text. Altering it would break the "complete and
# unmodified" requirement that makes the license valid, so it is exempt by path.
# Do not "fix" those characters.
PATH_EXEMPT = {"LICENSE"}

# Upstream Azimuth labels in the Module 3 cache use Greek small letter phi
# (codepoint U+03C6). Exempt that codepoint on the cache TSV only; do not
# transcribe the glyph in this file, and do not "fix" the labels.
PHI_CODEPOINT = 0x03C6
PHI_EXEMPT_RELPATHS = frozenset(
    {
        "outputs/tables/module3_label_cache_HBM828.GPVG.252.tsv",
    }
)

TYPOGRAPHY = {
    "\u2014": "em dash",
    "\u2013": "en dash",
    "\u2018": "curly single quote",
    "\u2019": "curly single quote",
    "\u201c": "curly double quote",
    "\u201d": "curly double quote",
    "\u2026": "ellipsis",
    "\u00a0": "non-breaking space",
    "\u2192": "arrow",
    "\u2194": "arrow",
    "\u00b7": "middle dot",
    "\u00b2": "superscript two",
    "\u00b3": "superscript three",
    "\u2260": "not-equal sign",
    "\u2212": "minus sign",
    "\u00b1": "plus-minus sign",
    "\u2264": "less-or-equal sign",
    "\u2265": "greater-or-equal sign",
}

PHRASE_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("delve into", re.compile(r"\bdelve into\b", re.I)),
    ("dive into", re.compile(r"\bdive into\b", re.I)),
    ("it is worth noting that", re.compile(r"\bit is worth noting that\b", re.I)),
    ("it should be noted that", re.compile(r"\bit should be noted that\b", re.I)),
    ("leverage (verb)", re.compile(r"\bleverage[sd]?\b", re.I)),
    ("seamless", re.compile(r"\bseamless\b", re.I)),
    ("showcase", re.compile(r"\bshowcase[sd]?\b", re.I)),
    ("underscore (verb)", re.compile(r"\bunderscore[sd]?\b", re.I)),
    ("in the realm of", re.compile(r"\bin the realm of\b", re.I)),
    ("landscape (metaphor)", re.compile(r"\blandscape\b", re.I)),
    ("pivotal", re.compile(r"\bpivotal\b", re.I)),
    ("vital", re.compile(r"\bvital\b", re.I)),
    ("comprehensive", re.compile(r"\bcomprehensive\b", re.I)),
    ("furthermore", re.compile(r"\bfurthermore\b", re.I)),
    ("moreover", re.compile(r"\bmoreover\b", re.I)),
    ("not only ... but also", re.compile(r"\bnot only\b.{0,80}\bbut also\b", re.I | re.S)),
]

SKIP_DIR_NAMES = {
    ".git",
    ".venv",
    "__pycache__",
    ".ipynb_checkpoints",
    "node_modules",
    # Large binary / source trees: do not scan contents (filenames still checked
    # when walked via parent listings for ship tree only).
    "source",
}


def _masked_for_allowlist(text: str) -> str:
    out = text
    for s in UPSTREAM_ALLOWLIST:
        if s in out:
            out = out.replace(s, " " * len(s))
    return out


def _iter_ship_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if any(part in SKIP_DIR_NAMES for part in p.parts):
            continue
        # Skip large binary / processed objects by suffix
        if p.suffix.lower() in {
            ".h5ad",
            ".h5mu",
            ".hdf5",
            ".h5",
            ".png",
            ".pdf",
            ".gz",
            ".zip",
            ".npy",
            ".npz",
            ".pkl",
            ".pickle",
            ".so",
            ".dylib",
        }:
            continue
        # Skip absolute-path-heavy run params and huge contrast tables by default?
        # Keep them: ASCII rule applies to machine-read TSVs.
        files.append(p)
    return files


def check_filename_ascii(path: Path, root: Path) -> list[str]:
    rel = str(path.relative_to(root))
    hits = []
    for i, ch in enumerate(rel):
        if ord(ch) > 127:
            hits.append(f"FILENAME non-ASCII U+{ord(ch):04X} in {rel!r} at pos {i}")
    return hits


def check_code_ascii(
    path: Path,
    text: str,
    *,
    allowed_codepoints: frozenset[int] | None = None,
) -> list[str]:
    masked = _masked_for_allowlist(text)
    allowed = allowed_codepoints or frozenset()
    hits = []
    for lineno, line in enumerate(masked.splitlines(), 1):
        for col, ch in enumerate(line, 1):
            if ord(ch) > 127 and ord(ch) not in allowed:
                hits.append(
                    f"{path}:L{lineno}:C{col} non-ASCII U+{ord(ch):04X} ({ch!r})"
                )
                if len(hits) >= 50:
                    hits.append(f"{path}: ... truncated after 50 ASCII findings")
                    return hits
    return hits


def check_prose_typography(path: Path, text: str, *, flag_greek: bool) -> list[str]:
    hits = []
    for lineno, line in enumerate(text.splitlines(), 1):
        for col, ch in enumerate(line, 1):
            if ch in TYPOGRAPHY:
                hits.append(
                    f"{path}:L{lineno}:C{col} {TYPOGRAPHY[ch]} U+{ord(ch):04X}"
                )
            elif flag_greek and 0x0370 <= ord(ch) <= 0x03FF:
                hits.append(
                    f"{path}:L{lineno}:C{col} Greek U+{ord(ch):04X} ({ch!r})"
                )
            if len(hits) >= 80:
                hits.append(f"{path}: ... truncated after 80 typography findings")
                return hits
    return hits


def check_phrases(path: Path, text: str) -> list[str]:
    # Mask upstream / quoted strings first, so the phrase list cannot flag a
    # verbatim citation title or a term used in its technical sense. Masking
    # substitutes equal-length spaces, so reported line numbers stay correct.
    text = _masked_for_allowlist(text)
    hits = []
    for label, pat in PHRASE_PATTERNS:
        for m in pat.finditer(text):
            # approximate line
            lineno = text.count("\n", 0, m.start()) + 1
            hits.append(f"{path}:L{lineno} phrase-review: {label}")
            if len(hits) >= 40:
                return hits
    return hits


def notebook_markdown_cells(path: Path) -> list[tuple[int, str]]:
    nb = json.loads(path.read_text(encoding="utf-8"))
    out = []
    for i, cell in enumerate(nb.get("cells") or []):
        if cell.get("cell_type") != "markdown":
            continue
        src = cell.get("source") or []
        if isinstance(src, list):
            text = "".join(src)
        else:
            text = str(src)
        out.append((i, text))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--root",
        type=Path,
        default=ROOT,
        help="Package root (default: the repository root)",
    )
    ap.add_argument(
        "--flag-greek",
        action="store_true",
        help="Flag Greek letters in learner prose (default: off; spell-out is preferred)",
    )
    ap.add_argument(
        "--strict-phrases",
        action="store_true",
        help="Treat phrase-review hits as failing findings",
    )
    ap.add_argument(
        "--skip-outputs",
        action="store_true",
        help="Skip outputs/ (generated tables may contain upstream non-ASCII labels)",
    )
    args = ap.parse_args()
    root: Path = args.root.resolve()

    hard: list[str] = []
    warn: list[str] = []

    for path in _iter_ship_files(root):
        rel = path.relative_to(root)
        if args.skip_outputs and rel.parts and rel.parts[0] == "outputs":
            continue
        hard.extend(check_filename_ascii(path, root))
        suf = path.suffix.lower()
        try:
            raw = path.read_bytes()
        except OSError as exc:
            warn.append(f"{path}: unreadable ({exc})")
            continue
        # Skip obvious binaries
        if b"\x00" in raw[:4096]:
            continue
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            hard.append(f"{path}: not valid UTF-8")
            continue

        if path.name in PATH_EXEMPT:
            # Verbatim upstream legal text: must not be altered. See PATH_EXEMPT.
            continue

        if suf in CODE_SUFFIXES or (not suf and path.name in EXTENSIONLESS_CHECKED):
            rel_posix = rel.as_posix()
            allowed = (
                frozenset({PHI_CODEPOINT}) if rel_posix in PHI_EXEMPT_RELPATHS else None
            )
            hard.extend(check_code_ascii(path, text, allowed_codepoints=allowed))
        elif suf in PROSE_SUFFIXES:
            hard.extend(check_prose_typography(path, text, flag_greek=args.flag_greek))
            phrases = check_phrases(path, text)
            (hard if args.strict_phrases else warn).extend(phrases)
        elif suf == ".ipynb":
            for cell_i, md in notebook_markdown_cells(path):
                label = f"{path}#cell{cell_i}"
                hard.extend(
                    check_prose_typography(
                        Path(label), md, flag_greek=args.flag_greek
                    )
                )
                phrases = check_phrases(Path(label), md)
                (hard if args.strict_phrases else warn).extend(phrases)

    print(f"check_text_hygiene: root={root}")
    print(f"hard findings: {len(hard)}")
    for h in hard[:200]:
        print("HARD", h)
    if len(hard) > 200:
        print(f"... {len(hard) - 200} more hard findings")
    print(f"phrase / soft findings: {len(warn)}")
    for w in warn[:80]:
        print("WARN", w)
    if len(warn) > 80:
        print(f"... {len(warn) - 80} more soft findings")
    return 1 if hard else 0


if __name__ == "__main__":
    raise SystemExit(main())
