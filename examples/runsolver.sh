#!/bin/bash
DIRNAME=$(dirname $0)
SOLVER=$(find $DIRNAME/.. -name cvxsolver -executable | head -1)
if [ "$CVXSOLVER" != "" ]; then
    SOLVER=CVXSOLVER
fi
if [ "$SOLVER" = "" ]; then
    echo "No solver" >&2
    exit 1
fi

if [ ! -x $SOLVER ]; then
    echo "'$SOLVER' not excutable" >&2
    exit 1
fi

echo $SOLVER $* >&2
$SOLVER $*

