#!/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPT="../5+5/5_plus_5.py"
PARAMS="params.py"
REGEXP="\[Node 0\] Simulation lasted: \([0-9]\+\) seconds"

echo "Hybrid (matrix/Krylov) regime benchmarking script."
echo "See small_matrix_size.*.log for resulting execution times."

for small_matrix_size in 0 1 5 10 25 50 100
do
    cp "${PARAMS}.in" "${PARAMS}"
    echo "p[\"krylov_small_matrix_size\"] = ${small_matrix_size}" >> ${PARAMS}
    echo "results_file_name = \"small_matrix_size.${small_matrix_size}.h5\"" >> ${PARAMS}
    for run in `seq 0 5`; do
        ${SCRIPT} | grep "${REGEXP}" | sed -e "s/${REGEXP}/\1/" >> "small_matrix_size.${small_matrix_size}.log"
    done
    rm ${PARAMS} ${PARAMS}c
done
