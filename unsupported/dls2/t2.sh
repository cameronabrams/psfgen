Vmax=`cat V.dat | awk 'BEGIN{m=0}{if (m<$2) m=$2}END{print m}'`
Vmin=`cat V.dat | awk 'BEGIN{m=99999999}{if (m>$2) m=$2}END{print m}'`
cat > tmp.gp << EOF
set term pdfcairo enhanced color fontscale 0.7 lw 1.5
set out "V.pdf"
set encoding iso_8859_1
set border 3
set xtics nomirror
set ytics nomirror
set xlabel "time, 10^3 steps (1 step = 2 fs)"
set ylabel "volume, 10^3 \305^3"
set xr [0:20.1]
ymin=int($Vmin/1000) - 1
ymax=int($Vmax/1000) + 1
set yr [ymin:ymax]
set ytics ymin,2,ymax
p "V.dat" u (\$1/1000):(\$2/1000.) not w l
EOF
gnuplot tmp.gp
rm tmp.gp

echo "Done."
