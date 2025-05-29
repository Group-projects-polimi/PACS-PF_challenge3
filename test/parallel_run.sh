#!/bin/bash

EXECUTABLE="./executable"
TEST_DATA="./test/testingdata.json"
OUTPUT_DIR="./test/output"
PLOT_DIR="./plots"
BASE_DATA_FILE="${OUTPUT_DIR}/paralleltime.txt"
PLOT_FILE="${PLOT_DIR}/time_parallel.png"

PROCESS_COUNTS=(2 4 8)
TOTAL_CORES=8 

for NPROCS in "${PROCESS_COUNTS[@]}"; do
    OMP_THREADS=$((TOTAL_CORES / NPROCS))
    export OMP_NUM_THREADS=$OMP_THREADS

    mpirun --oversubscribe -np "$NPROCS" "$EXECUTABLE" "$TEST_DATA"

    unset OMP_NUM_THREADS
done

PLOT_COMMAND="plot "
FIRST_PLOT=true

for NPROCS in "${PROCESS_COUNTS[@]}"; do
    OMP_THREADS=$((TOTAL_CORES / NPROCS))
    CURRENT_DATA_FILE="${OUTPUT_DIR}/paralleltime_${NPROCS}.txt"

    if [ -f "$CURRENT_DATA_FILE" ]; then
        if [ "$FIRST_PLOT" = true ]; then
            PLOT_COMMAND+="'${CURRENT_DATA_FILE}' using 1:2 with linespoints pt 7 lw 2 title '${NPROCS}p / ${OMP_THREADS}t'"
            FIRST_PLOT=false
        else
            PLOT_COMMAND+=", '${CURRENT_DATA_FILE}' using 1:2 with linespoints pt 7 lw 2 title '${NPROCS}p / ${OMP_THREADS}t'"
        fi
    else
        echo "Warning: Data file '$CURRENT_DATA_FILE' not found, skipping for plot." >&2
    fi
done

gnuplot << EOF
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output '${PLOT_FILE}'
set title 'Execution Time vs Grid Points (Parallel Cases)'
set xlabel 'Grid points'
set ylabel 'Execution Time (ms)'
set grid
set logscale y # Often useful for comparing times, remove if not needed
set key top left # Turn on the legend and place it
${PLOT_COMMAND}
EOF

if [ -f "$PLOT_FILE" ]; then
    echo "Plot successfully generated at '${PLOT_FILE}'"
else
    echo "Error: Gnuplot might have failed. Plot file not found."
fi

