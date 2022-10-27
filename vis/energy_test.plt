set terminal pngcairo enhanced font'font,14'
set output 'vis/energy_test.png'


file = "'e.dat'"

set xlabel 'Time'
set ylabel 'Energy'
set key right center
plot @file using 1:2 title 'kinetic', @file using 1:3 title 'potential', @file using 1:4 title 'total' w lines
