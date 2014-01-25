#!/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPT="./5_plus_5.py"
PARAMS="./params.py"

cp "${PARAMS}" "${PARAMS}.bak"
for small_matrix_size in 0 1 5 10 25 50 100
do
    sed -i "s/^\(p\[\"krylov_small_matrix_size\"\] =\)[ 0-9]\+/\1 ${small_matrix_size}/" ${PARAMS}
    sed -i "s/^\(results_file_name =\).*/\1 \"small_matrix_size.${small_matrix_size}.h5\"/" ${PARAMS}
    for run in `seq 0 5`; do
        ${SCRIPT} | grep "Simulation lasted:" >> "small_matrix_size.${small_matrix_size}.log"
    done
done
mv -f "${PARAMS}.bak" "${PARAMS}"
