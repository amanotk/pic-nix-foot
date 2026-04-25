#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH="${SCRIPT_DIR}/patches/ion-bg-subtraction.patch"

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Clone pic-nix and apply the ion background B-field subtraction patch."
    echo ""
    echo "Options:"
    echo "  -r, --ref REF    pic-nix branch or tag to checkout (default: develop)"
    echo "  -d, --dir  DIR   pic-nix clone destination (default: ./pic-nix)"
    echo "  -h, --help       Show this help"
}

PICNIX_REF="develop"
PICNIX_DST=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--ref)
            PICNIX_REF="$2"
            shift 2
            ;;
        -d|--dir)
            PICNIX_DST="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# destination
PICNIX_DST="${PICNIX_DST:-${SCRIPT_DIR}/pic-nix}"

# clone if needed
if [ ! -d "${PICNIX_DST}" ]; then
    echo "Cloning pic-nix into ${PICNIX_DST} ..."
    git clone https://github.com/amanotk/pic-nix.git "${PICNIX_DST}"
fi

# checkout requested ref
echo "Checking out ${PICNIX_REF} ..."
(cd "${PICNIX_DST}" && git checkout "${PICNIX_REF}")

# apply patch
echo "Applying patch: ${PATCH}"
(cd "${PICNIX_DST}" && git apply --check "${PATCH}" 2>/dev/null && git apply "${PATCH}") \
    || echo "Warning: patch did not apply cleanly (may already be applied)"

echo ""
echo "Done. Patched pic-nix is at: ${PICNIX_DST}"
echo "Build the 'foot' target from there, e.g.:"
echo "  cmake -S ${PICNIX_DST} -B build && cmake --build build --target foot"
