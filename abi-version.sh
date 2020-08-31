#!/bin/sh

CMDDIR=$(dirname $0)

CONVEX_H="${CMDDIR}/src/include/cvxc.h"
if [ ! -f "$CONVEX_H" ]; then
    echo "abi-version.sh: header cvxc.h not available"
    exit 1
fi

CURRENT=$(awk '/^#define +CVXC_ABI_CURRENT/ {print $3}' $CONVEX_H)
REVISION=$(awk '/^#define +CVXC_ABI_REVISION/ {print $3}' $CONVEX_H)
AGE=$(awk '/^#define +CVXC_ABI_AGE/ {print $3}' $CONVEX_H)

if [ -z "$CURRENT" -o -z "$REVISION" -o -z "$AGE" ]; then
    echo "abi-version.sh: abi version macros not found."
    exit 1
fi

case $1 in
    -libtool)
        printf '%s' "$CURRENT:$REVISION:$AGE"
        ;;
    *)
        printf '%s' "$CURRENT.$REVISION.$AGE"
        ;;
esac

