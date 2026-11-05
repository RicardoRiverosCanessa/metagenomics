#!/usr/bin/env bash
# check_def_syntax.sh – deterministic static checks for .def files
# No container building; simply verifies required sections exist.

set -euo pipefail

def_dir="$(dirname "${0}")/../defs"

if [[ ! -d "$def_dir" ]]; then
  echo "[ERROR] defs directory not found: $def_dir" >&2
  exit 1
fi

missing=0
for def in "$def_dir"/*.def; do
  [[ -e "$def" ]] || continue
  # Required keywords (case‑insensitive)
  for kw in "Bootstrap" "From" "%post"; do
    if ! grep -i "^$kw" "$def" > /dev/null; then
      echo "[MISSING] $kw in $def"
      missing=$((missing+1))
    fi
  done
  # Disallow use of eval or backticks which imply dynamic execution
  if grep -E "eval|\`" "$def" > /dev/null; then
    echo "[DISALLOWED] dynamic eval/backticks in $def"
    missing=$((missing+1))
  fi
done

if (( missing > 0 )); then
  echo "[FAIL] $missing issues detected in .def files"
  exit 1
else
  echo "[OK] All .def files passed static checks"
fi
