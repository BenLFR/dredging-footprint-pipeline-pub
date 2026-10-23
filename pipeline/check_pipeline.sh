#!/bin/bash
# Pipeline guard: does PIPELINE.md still tell the truth?
#
# DESIGN: the guard checks the RESULT, never a hand-typed count.
#
# Everything is recomputed from disk and compared with what the manifest
# claims. No number is written in this script.
#
# Modifies NOTHING. Exits with 1 and names each gap.
set -uo pipefail

ROOT="${ROOT:-${DFP_REPO:?set DFP_REPO}}"
V6="$ROOT/pipeline"
MAN="${MAN:-$V6/PIPELINE.md}"
GAPS=0

flag_gap() { echo "  ✗ $*"; GAPS=$((GAPS + 1)); }

echo "=== pipeline guard — $(date '+%Y-%m-%d %H:%M') ==="

# --- 1. every Rscript <file> of a .sbatch points to an existing file ---------
echo
echo "[1] do the .sbatch files point to existing files?"
n_ok=0
while IFS= read -r f; do
  base=$(dirname "$f")
  cd_dir=$(grep -hoE '^[[:space:]]*cd[[:space:]]+"?[^"]+"?' "$f" 2>/dev/null | tail -1 \
           | sed -E 's/^[[:space:]]*cd[[:space:]]+//; s/"//g')
  cd_dir=$(eval echo "$cd_dir" 2>/dev/null)
  grep -hoE 'Rscript[[:space:]]+[^[:space:];|&]+' "$f" 2>/dev/null | while read -r line; do
    target=$(printf '%s' "${line#Rscript}" | tr -d " \"'")
    case "$target" in -*|''|*'$'*) continue ;; esac
    found=0
    for cand in "$target" "$cd_dir/$target" "$base/$target" "$V6/$target"; do
      if [ -n "$cand" ] && [ -f "$cand" ]; then found=1; break; fi
    done
    # an `exit` here would leave the pipe subshell and skip the following
    # Rscript calls of the same file: mark and continue.
    [ "$found" -eq 0 ] && echo "MISSING|${f#$HOME/}|$target"
  done
done < <(find "$ROOT" -name '*.sbatch' 2>/dev/null | sort) > /tmp/guard_$$.txt
while IFS='|' read -r _ launcher target; do
  [ -n "${launcher:-}" ] && flag_gap "$launcher runs \`$target\` — file not found"
done < /tmp/guard_$$.txt
rm -f /tmp/guard_$$.txt
[ "$GAPS" -eq 0 ] && echo "  ✓ every .sbatch points to yr existing file"

# --- 2. does the manifest exist, and are its fingerprints up to date? --------
echo
echo "[2] does the manifest match the disk?"
if [ ! -f "$MAN" ]; then
  flag_gap "$MAN missing — run make_manifest.sh"
else
  before=$GAPS
  # rows of the "By step" table: | stepN | `path` | `md5` | ...
  grep -E '^\| step[0-9] \| `' "$MAN" | while IFS='|' read -r _ _ path fingerprint _; do
    path=$(printf '%s' "$path" | tr -d ' `')
    # the column holds `md5` · size · date since the inner separator
    # is a middle dot: keep only the md5, otherwise the comparison always
    # fails and the guard reports "stale manifest" on intact files.
    fingerprint=$(printf '%s' "$fingerprint" | tr -d ' `' | sed 's/·.*$//')
    [ -n "$path" ] && [ "$path" != "—" ] || continue
    if [ ! -f "$V6/$path" ]; then
      echo "ABSENT|$path"
    else
      actual=$(md5sum "$V6/$path" | cut -c1-8)
      [ "$actual" = "$fingerprint" ] || echo "DIFFERS|$path|$fingerprint|$actual"
    fi
  done > /tmp/guard2_$$.txt
  while IFS='|' read -r what_kind path exp actual; do
    case "$what_kind" in
      ABSENT)  flag_gap "the manifest cites \`$path\`, missing from the disk" ;;
      DIFFERS) flag_gap "\`$path\`: manifest \`$exp\`, disk \`$actual\` — stale manifest" ;;
    esac
  done < /tmp/guard2_$$.txt
  rm -f /tmp/guard2_$$.txt
  [ "$GAPS" -eq "$before" ] && echo "  ✓ every fingerprint of the manifest matches the disk"
fi

# --- 3. is the live code in git? ------------------------------------------------
echo
echo "[3] is the live code committed?"
dirty=$(cd "$ROOT" && git status --porcelain pipeline 2>/dev/null)
if [ -n "$dirty" ]; then
  while IFS= read -r l; do flag_gap "not committed: $l"; done <<< "$dirty"
  echo "     → git cannot tell which version ran as long as this is the case"
else
  echo "  ✓ pipeline/ is clean with respect to git"
fi

echo
if [ "$GAPS" -eq 0 ]; then
  echo "=== PASS — the manifest tells the truth ==="
  exit 0
fi
echo "=== $GAPS gap(s) — see above ==="
exit 1
