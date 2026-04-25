#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH="${SCRIPT_DIR}/patches/ion-bg-subtraction.patch"
BUILD_DIR="${SCRIPT_DIR}/build"

# pic-nix source directory
if [ -n "${PICNIX_DIR:-}" ]; then
    PICNIX_SRC="$(cd "${PICNIX_DIR}" && pwd)"
else
    PICNIX_SRC="${SCRIPT_DIR}/pic-nix"
    if [ ! -d "${PICNIX_SRC}" ]; then
        echo "Cloning pic-nix into ${PICNIX_SRC} ..."
        git clone https://github.com/amanotk/pic-nix.git "${PICNIX_SRC}"
    fi
fi

# branch / tag to use (default: develop)
PICNIX_REF="${PICNIX_REF:-develop}"

# checkout requested ref
(cd "${PICNIX_SRC}" && git checkout "${PICNIX_REF}")

# apply patch
echo "Applying patch: ${PATCH}"
(cd "${PICNIX_SRC}" && git apply --check "${PATCH}" 2>/dev/null && git apply "${PATCH}" \
    || echo "Patch already applied or does not apply cleanly — continuing")

# build
CMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER:-mpicxx}"
CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS:--O2 -fopenmp}"

cmake -S "${PICNIX_SRC}" -B "${BUILD_DIR}" \
    -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}" \
    -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}"
cmake --build "${BUILD_DIR}" --target foot

echo ""
echo "Build complete. Executable at: ${BUILD_DIR}/pic/example/foot/alfven/main.out"
echo "Scenario directories: ${BUILD_DIR}/pic/example/foot/{alfven,buneman,weibel}"
