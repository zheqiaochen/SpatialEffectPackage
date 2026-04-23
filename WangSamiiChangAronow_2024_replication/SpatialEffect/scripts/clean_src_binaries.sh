#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="${SCRIPT_DIR}/../src"

if [[ ! -d "${TARGET_DIR}" ]]; then
  echo "Error: target directory not found: ${TARGET_DIR}" >&2
  exit 1
fi

count=$(find "${TARGET_DIR}" -type f \( -name "*.so" -o -name "*.o" \) | wc -l | tr -d ' ')

if [[ "${count}" -eq 0 ]]; then
  echo "No .so or .o files found in ${TARGET_DIR}"
  exit 0
fi

echo "Deleting .so/.o files from ${TARGET_DIR}:"
find "${TARGET_DIR}" -type f \( -name "*.so" -o -name "*.o" \) -print -delete
echo "Done. Deleted ${count} file(s)."
