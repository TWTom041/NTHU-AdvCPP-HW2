#!/usr/bin/env bash
# Run every experiment of the report; console output -> results/*.txt,
# raw data for the figures -> build/.   Total ~25 min on 4 cores.
#
#   make all cross raptor && ./scripts/run_all.sh [section ...]
# sections: env validate raptor q1 parallel cache q2 q3 q4 q5 q6 renders   (default: all)
set -uo pipefail
cd "$(dirname "$0")/.."
R=results
B=build
mkdir -p $R
QEMU="qemu-aarch64 -L /usr/aarch64-linux-gnu"
SECTIONS=${*:-"env validate raptor q1 parallel cache q2 q3 q4 q5 q6 renders"}
want() { [[ " $SECTIONS " == *" $1 "* ]]; }
# run a command, echo it like a terminal prompt, append output to a log
run() { local log=$1; shift; { echo "\$ $*"; "$@" 2>&1; echo; } | tee -a "$log"; }

if want env; then
  f=$R/q0_environment.txt; : > $f
  run $f lscpu
  run $f grep -m1 -o "flags.*" /proc/cpuinfo
  run $f uname -srmo
  run $f cat /etc/os-release
  run $f g++ --version
  run $f aarch64-linux-gnu-g++ --version
  run $f clang++ --version
  run $f ldd --version
  run $f qemu-aarch64 --version
  run $f valgrind --version
  run $f free -h
  run $f nproc
fi

if want validate; then
  f=$R/validate_shadow.txt; : > $f
  run $f $B/kerr_validate --angles 12
fi

if want raptor; then
  f=$R/raptor_compare.txt; : > $f
  run $f $B/raptor_compare --res 100 --method rk4 --step 0.03 --dump $B/raptor_cmp_rk4.bin
  run $f $B/raptor_compare --res 100 --method rk4 --step 0.01
  run $f $B/raptor_compare --res 100 --method rk2 --step 0.03
  run $f $B/raptor_compare_ofast --res 100 --method rk4 --step 0.03
  run $f $B/raptor_compare --res 100 --method rk4 --step 0.03 --margin 0.01
  f=$R/raptor_full.txt; : > $f
  run $f ./scripts/run_raptor_full.sh
  run $f python3 scripts/compare.py raptor $B/raptor_O3/output/img_data_0_1.000000e+11_60.00.dat \
                                         $B/raptor_Ofast/output/img_data_0_1.000000e+11_60.00.dat
fi

if want q1; then
  f=$R/q1_microbench.txt; : > $f
  run $f $B/q1_base
  run $f $B/q1_native
  run $f bash -c "objdump -d --no-show-raw-insn $B/q1_base | grep -m3 -E 'vaddsh|__addhf3|call.*hf' ; \
                   objdump -d --no-show-raw-insn $B/q1_native | grep -m3 -E 'vaddsh|vaddph'; \
                   objdump -d --no-show-raw-insn $B/q1_base | grep -m3 -E 'fadd|faddp'; \
                   objdump -d --no-show-raw-insn $B/q1_base | grep -m3 -E 'call.*__addtf3|call.*__multf3'"
  f=$R/q1_kerr_precision.txt; : > $f
  run $f $B/kerr_precision --res 48 --reps 3
  run $f $B/kerr_precision --res 48 --reps 3 --exact-ktheta
fi

if want parallel; then
  f=$R/parallel.txt; : > $f
  run $f $B/parallel_bench --res 256 --reps 3 --type double
  run $f $B/parallel_bench --res 256 --reps 3 --type float --false-sharing 0
fi

if want cache; then
  f=$R/cache_sweep.txt; : > $f
  run $f $B/cache_sweep 1024
  f=$R/cachegrind.txt; : > $f
  for m in row col; do for t in float double; do
    run $f $B/cache_patterns $m $t 4096
  done; done
  for m in row col; do for t in float double; do
    run $f bash -c "valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=/dev/null \
      $B/cache_patterns $m $t 2048 2>&1 | grep -E 'n=|D   refs|D1  miss|LLd miss|D1  miss rate|LLd miss rate'"
  done; done
  # the ray tracer itself: small working set -> almost no misses
  g++ -std=gnu++23 -O2 -Iinclude src/kerr_render.cpp -o $B/kerr_render_serial -lquadmath
  for t in float double; do
    run $f bash -c "valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=/dev/null \
      $B/kerr_render_serial --type $t --res 24 --disk 2>&1 | grep -E 'type=|I   refs|I1  miss rate|D   refs|D1  miss|LLd miss|D1  miss rate|LLd miss rate'"
  done
fi

if want q2; then
  f=$R/q2_sin.txt; : > $f
  run $f $B/q2_sin 200000
  run $f ./scripts/which_sin.sh $B/q2_ifunc
  run $f env GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA ./scripts/which_sin.sh $B/q2_ifunc
fi

if want q3; then
  f=$R/q3_fastmath.txt; : > $f
  run $f $B/q3_strict
  run $f $B/q3_fast
  for b in O3 fast; do for t in double float; do
    run $f env OMP_NUM_THREADS=1 $B/kerr_render_$b --type $t --res 192 --disk --out $B/fm_${b}_$t.bin
  done; done
  run $f python3 scripts/compare.py kimg $B/fm_O3_double.bin $B/fm_fast_double.bin "-O3" "-O3 -ffast-math"
  run $f python3 scripts/compare.py kimg $B/fm_O3_float.bin $B/fm_fast_float.bin "-O3" "-O3 -ffast-math"
fi

if want q4; then
  f=$R/q4_printf.txt; : > $f
  run $f $B/q4_printf
  run $f $B/q4_printf_musl
fi

if want q5; then
  f=$R/q5_env.txt; : > $f
  run $f $B/q5_x86_glibc $B/q5_x86_glibc.bin
  run $f env GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA $B/q5_x86_glibc $B/q5_x86_glibc_nofma.bin
  run $f $B/q5_x86_musl $B/q5_x86_musl.bin
  run $f $QEMU $B/q5_arm_glibc $B/q5_arm_glibc.bin
  run $f python3 scripts/compare.py libm $B/q5_x86_glibc.bin $B/q5_x86_glibc_nofma.bin "glibc(FMA ifunc)" "glibc(-AVX2,-FMA)"
  run $f python3 scripts/compare.py libm $B/q5_x86_glibc.bin $B/q5_x86_musl.bin "x86-64 glibc" "x86-64 musl"
  run $f python3 scripts/compare.py libm $B/q5_x86_glibc.bin $B/q5_arm_glibc.bin "x86-64 glibc" "aarch64 glibc"
  run $f $B/kerr_render --res 96 --disk --out $B/q5_kerr_x86.bin
  run $f $QEMU $B/kerr_render_arm --res 96 --disk --out $B/q5_kerr_arm.bin
  run $f python3 scripts/compare.py kimg $B/q5_kerr_x86.bin $B/q5_kerr_arm.bin x86-64 aarch64
  run $f bash -c "aarch64-linux-gnu-objdump -d $B/q5_arm_glibc | grep -A1 '<mul_sub>:'; objdump -d $B/q5_x86_glibc | grep -A3 '<mul_sub>:'"
  for t in 1 2 3 4; do run $f env OMP_NUM_THREADS=$t $B/q5_threads; done
fi

if want q6; then
  f=$R/q6_extra.txt; : > $f
  run $f $B/q6_extra
fi

if want renders; then
  f=$R/renders.txt; : > $f
  run $f $B/kerr_render --type double --res 512 --disk --out $B/render_double.bin
  run $f $B/kerr_render --type float --res 512 --disk --out $B/render_float.bin
  run $f $B/kerr_render --type double --res 512 --out $B/render_shadow.bin
  run $f $B/kerr_render --type double --res 512 --disk --spin 0 --out $B/render_schw.bin
fi
echo "done: $SECTIONS"
