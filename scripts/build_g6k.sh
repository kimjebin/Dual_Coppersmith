#!/usr/bin/env bash
set -euo pipefail
# Usage: bash scripts/build_g6k.sh /absolute/path/to/experiment-work
task_work="${1:?An experiment work directory is required}"
task_python="${PYTHON:-python}"
task_prefix="$($task_python -c 'import sys; print(sys.prefix)')"
task_commit=c71e084dd003e1ffbfe8c91a7211b0ce9da22245
task_source="$task_work/third_party/g6k-$task_commit"
mkdir -p "$task_work/third_party"
if [ ! -d "$task_source/.git" ]; then
    git clone https://github.com/fplll/g6k.git "$task_source"
fi
git -C "$task_source" checkout --detach "$task_commit"
test "$(git -C "$task_source" rev-parse HEAD)" = "$task_commit"
export CC=/usr/bin/gcc CXX=/usr/bin/g++ LD=/usr/bin/ld
export PKG_CONFIG_PATH="$task_prefix/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="$task_prefix/lib:${LD_LIBRARY_PATH:-}"
export PATH="/usr/bin:/bin:$task_prefix/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$task_source"
if [ ! -e configure ]; then autoreconf -i; fi
./configure --disable-shared --enable-static --prefix="$task_source/local"
make -j1
