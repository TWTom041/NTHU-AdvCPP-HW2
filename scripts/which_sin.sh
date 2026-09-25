#!/usr/bin/env bash
# Which glibc sin() variant did the IFUNC resolver pick on this CPU?
# Disassemble libm at the resolved address and count FMA instructions.
set -euo pipefail
BIN=${1:-build/q2_ifunc}
LIBM=$(ldd "$BIN" | awk '/libm.so/{print $3}')
line=$("$BIN")
off=$(echo "$line" | sed 's/.*+ 0x//')
n=$(objdump -d --no-show-raw-insn --start-address=0x$off --stop-address=$(printf '0x%x' $((16#$off + 1200))) "$LIBM" \
    | grep -cE 'vfmadd|vfmsub|vfnmadd|vfnmsub' || true)
echo "$line  -> $n FMA instructions in the first 1200 bytes (GLIBC_TUNABLES=${GLIBC_TUNABLES:-<unset>})"
