#!/usr/bin/env bash
# Build and run the unmodified RAPTOR (GRMHD image of the bundled dump040) twice:
# with its own makefile flags (-Ofast => -ffast-math) and with IEEE-conforming -O3.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
for variant in Ofast O3; do
  D=build/raptor_$variant
  rm -rf "$D"; mkdir -p "$D"
  cp third_party/raptor/*.c third_party/raptor/*.h third_party/raptor/model.in third_party/raptor/dump040 "$D"/
  ( cd "$D"
    # -fcommon: RAPTOR defines globals in headers (GCC >= 10 defaults to -fno-common)
    gcc -fopenmp -std=c99 -$variant -fcommon -w -o RAPTOR *.c -lm -lgsl -lgslcblas
    mkdir -p output
    start=$(date +%s.%N)
    ./RAPTOR model.in dump040 1e19 60 1 1 0 > log.txt 2> err.txt
    end=$(date +%s.%N)
    echo "RAPTOR -$variant: wall $(echo "$end - $start" | bc) s, $(grep -o 'Integrated flux density = .*' err.txt)" )
done
