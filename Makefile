# Kerr geodesic engine + floating-point experiments (Advanced C++ HW2)
#
#   make            build everything (x86-64, needs GSL, MPFR, TBB, libquadmath)
#   make cross      also build the aarch64 / musl variants used in Q5
#   make raptor     fetch RAPTOR (pinned commit) - required by raptor_compare
#   make run        run all experiments, logs -> results/
#
# Flags are spelled out per target on purpose: several experiments are about
# what the flags do.

CXX      ?= g++
CC       ?= gcc
STD      := -std=gnu++23
INC      := -Iinclude -Iraptor_bridge
B        := build
RAPTOR   := third_party/raptor

CXXFLAGS := $(STD) -O2 -Wall -Wextra -Wno-unused-parameter $(INC)

ENGINE_HDR := $(wildcard include/kerr/*.hpp)

TARGETS := $(B)/kerr_render $(B)/kerr_render_O3 $(B)/kerr_render_fast $(B)/kerr_validate \
           $(B)/kerr_precision $(B)/parallel_bench \
           $(B)/raptor_compare $(B)/raptor_compare_ofast \
           $(B)/q1_base $(B)/q1_native $(B)/cache_sweep $(B)/cache_patterns \
           $(B)/q2_sin $(B)/q2_ifunc $(B)/q3_strict $(B)/q3_fast $(B)/q4_printf $(B)/q5_x86_glibc $(B)/q5_threads \
           $(B)/q6_extra

all: $(TARGETS)

$(B):
	mkdir -p $(B)

# ---- engine tools ---------------------------------------------------------
$(B)/kerr_render: src/kerr_render.cpp $(ENGINE_HDR) | $(B)
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@ -lquadmath
$(B)/kerr_render_O3: src/kerr_render.cpp $(ENGINE_HDR) | $(B)
	$(CXX) $(STD) -O3 $(INC) -fopenmp $< -o $@ -lquadmath
$(B)/kerr_render_fast: src/kerr_render.cpp $(ENGINE_HDR) | $(B)
	$(CXX) $(STD) -O3 -ffast-math $(INC) -fopenmp $< -o $@ -lquadmath
$(B)/kerr_validate: src/kerr_validate.cpp $(ENGINE_HDR) | $(B)
	$(CXX) $(CXXFLAGS) $< -o $@ -lquadmath
$(B)/kerr_precision: src/kerr_precision.cpp $(ENGINE_HDR) | $(B)
	$(CXX) $(CXXFLAGS) $< -o $@ -lquadmath
$(B)/parallel_bench: src/parallel_bench.cpp $(ENGINE_HDR) | $(B)
	$(CXX) $(STD) -O2 -march=native $(INC) -fopenmp $< -o $@ -lquadmath -ltbb

# ---- RAPTOR comparison ------------------------------------------------------
$(RAPTOR)/metric.c:
	./scripts/fetch_raptor.sh
raptor: $(RAPTOR)/metric.c

$(B)/raptor_bridge_O2.o: raptor_bridge/raptor_bridge.c raptor_bridge/raptor_bridge.h $(RAPTOR)/metric.c | $(B)
	$(CC) -std=gnu99 -O2 -w -I$(RAPTOR) -Iraptor_bridge -c $< -o $@
# RAPTOR's own makefile compiles with -Ofast (= -O3 -ffast-math)
$(B)/raptor_bridge_Ofast.o: raptor_bridge/raptor_bridge.c raptor_bridge/raptor_bridge.h $(RAPTOR)/metric.c | $(B)
	$(CC) -std=gnu99 -Ofast -w -I$(RAPTOR) -Iraptor_bridge -c $< -o $@
$(B)/raptor_compare: src/raptor_compare.cpp $(B)/raptor_bridge_O2.o $(ENGINE_HDR)
	$(CXX) $(CXXFLAGS) -fopenmp $< $(B)/raptor_bridge_O2.o -o $@ -lquadmath -lgsl -lgslcblas -lm
$(B)/raptor_compare_ofast: src/raptor_compare.cpp $(B)/raptor_bridge_Ofast.o $(ENGINE_HDR)
	$(CXX) $(CXXFLAGS) -fopenmp $< $(B)/raptor_bridge_Ofast.o -o $@ -lquadmath -lgsl -lgslcblas -lm

# ---- experiments ------------------------------------------------------------
$(B)/q1_base: experiments/q1_microbench.cpp | $(B)
	$(CXX) $(STD) -O2 -DBUILD_FLAGS='"g++ -O2 (x86-64 baseline ISA)"' $< -o $@ -lquadmath
$(B)/q1_native: experiments/q1_microbench.cpp | $(B)
	$(CXX) $(STD) -O3 -march=native -DBUILD_FLAGS='"g++ -O3 -march=native"' $< -o $@ -lquadmath
$(B)/cache_sweep: experiments/cache_sweep.cpp | $(B)
	$(CXX) $(STD) -O3 -march=native $< -o $@
# valgrind cannot decode AVX-512, keep this one at the baseline ISA
$(B)/cache_patterns: experiments/cache_patterns.cpp | $(B)
	$(CXX) $(STD) -O2 $< -o $@
$(B)/q2_sin: experiments/q2_sin.cpp | $(B)
	$(CXX) $(STD) -O2 $< -o $@ -lmpfr -lgmp -lquadmath
$(B)/q2_ifunc: experiments/q2_ifunc.c | $(B)
	$(CC) -O2 $< -o $@ -lm -ldl
$(B)/q3_strict: experiments/q3_fastmath.cpp | $(B)
	$(CXX) $(STD) -O3 $< -o $@
$(B)/q3_fast: experiments/q3_fastmath.cpp | $(B)
	$(CXX) $(STD) -O3 -ffast-math $< -o $@
$(B)/q4_printf: experiments/q4_printf.cpp | $(B)
	$(CXX) $(STD) -O2 $< -o $@
$(B)/q5_x86_glibc: experiments/q5_env.c | $(B)
	$(CC) -O2 $< -o $@ -lm
$(B)/q5_threads: experiments/q5_threads.cpp | $(B)
	$(CXX) $(STD) -O2 -fopenmp $< -o $@
$(B)/q6_extra: experiments/q6_extra.cpp $(ENGINE_HDR) | $(B)
	$(CXX) $(CXXFLAGS) $< -o $@ -lquadmath

# ---- Q5 cross builds: same GCC 13.3, same -O2 --------------------------------
cross: $(B)/q5_x86_musl $(B)/q5_arm_glibc $(B)/kerr_render_arm $(B)/q4_printf_musl
$(B)/q5_x86_musl: experiments/q5_env.c | $(B)
	musl-gcc -O2 -static $< -o $@ -lm
$(B)/q5_arm_glibc: experiments/q5_env.c | $(B)
	aarch64-linux-gnu-gcc -O2 $< -o $@ -lm
$(B)/kerr_render_arm: src/kerr_render.cpp $(ENGINE_HDR) | $(B)
	aarch64-linux-gnu-g++ $(CXXFLAGS) -fopenmp $< -o $@
$(B)/q4_printf_musl: experiments/q4_printf.cpp | $(B)
	musl-gcc -O2 -static -x c $< -o $@ -lm

run: all cross
	./scripts/run_all.sh

clean:
	rm -rf $(B)

.PHONY: all cross raptor run clean
