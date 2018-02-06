#!/usr/bin/gnuplot
set term pdfcairo enhanced color fontscale 0.7 lw 1.5
set out "V.pdf"
set encoding iso_8859_1

set border 3
set xtics nomirror
set ytics nomirror

set xlabel "time, 10^3 steps (1 step = 2 fs)"
set ylabel "volume, 10^3 \305^3"
set xr [0:20.1]
set yr [26:32]
set ytics 26,1,32

p "V.dat" u ($1/1000):($2/1000.) not w l

