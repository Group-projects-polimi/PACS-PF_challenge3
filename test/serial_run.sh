#!/bin/bash

./executable ./test/testingdata.json

DATA_FILE="./test/output/serialtime.txt"
OUTPUT_FILE="./plots/time_serial.png"

gnuplot << EOF
set terminal pngcairo size 800,600 enhanced font 'Verdana,10' 
set output '${OUTPUT_FILE}'
set title 'Execution Time, serial case'
set xlabel 'Grid points'
set ylabel 'Execution Time (ms)'
set grid
set key off
plot '${DATA_FILE}' using 1:2 with linespoints pt 7 lw 2
EOF
