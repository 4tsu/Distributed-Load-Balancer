set terminal pngcairo enhanced font'font,14'

file0 = "'time_net_0.dat'"
file1 = "'time_net_1.dat'"
file2 = "'time_net_2.dat'"
file3 = "'time_net_3.dat'"
file4 = "'time_net_4.dat'"
file5 = "'time_net_5.dat'"

set output 'calc.png'
set format y "%.1t{/Symbol \26410^{%T}"
set key left center
set xlabel 'step'
set ylabel 'calculation-time per step [ms]' offset -0.5, 0.0
set yrange [0:]
plot @file0 using 1:3 title 'w/o load-balancer' w linespoints, @file1 using 1:3 title 'global sort' w linespoints, @file2 using 1:3 title 'voronoi' w linespoints, @file3 using 1:3 title 'RCB' w linespoints, @file4 using 1:3 title '1D parallel' w linespoints, @file5 using 1:3 title 'skew boundary' w linespoints

# ----------------------------------------------------------------

file0 = "'time_gross_0.dat'"
file1 = "'time_gross_1.dat'"
file2 = "'time_gross_2.dat'"
file3 = "'time_gross_3.dat'"
file4 = "'time_gross_4.dat'"
file5 = "'time_gross_5.dat'"

set output 'calc_comm.png'
set xlabel 'step'
set ylabel 'execution-time per step [ms]' offset -0.5, 0.0
set key right center
set yrange [0:]
plot @file0 using 1:3 title 'w/o load-balancer' w linespoints, @file1 using 1:3 title 'global sort' w linespoints, @file2 using 1:3 title 'voronoi' w linespoints, @file3 using 1:3 title 'RCB' w linespoints, @file4 using 1:3 title '1D parallel' w linespoints, @file5 using 1:3 title 'skew boundary' w linespoints

# ----------------------------------------------------------------

file0 = "'time_comm_0.dat'"
file1 = "'time_comm_1.dat'"
file2 = "'time_comm_2.dat'"
file3 = "'time_comm_3.dat'"
file4 = "'time_comm_4.dat'"
file5 = "'time_comm_5.dat'"

set output 'comm.png'
set xlabel 'step'
set ylabel 'communication-time per step [ms]'
set key right center
set yrange [0:]
plot @file0 using 1:3 title 'w/o load-balancer' w linespoints, @file1 using 1:3 title 'global sort' w linespoints, @file2 using 1:3 title 'voronoi' w linespoints, @file3 using 1:3 title 'RCB' w linespoints, @file4 using 1:3 title '1D parallel' w linespoints, @file5 using 1:3 title 'skew boundary' w linespoints

# ----------------------------------------------------------------

file0 = "'load_balance_0.dat'"
file1 = "'load_balance_1.dat'"
file2 = "'load_balance_2.dat'"
file3 = "'load_balance_3.dat'"
file4 = "'load_balance_4.dat'"
file5 = "'load_balance_5.dat'"


set output 'load_balance.png'
set xlabel 'step'
set ylabel 'workload per step'
unset format y
set yrange [0:]
set key left center
plot @file0 using 1:3 title 'w/o load-balancer' w linespoints, @file1 using 1:3 title 'global sort' w linespoints, @file2 using 1:3 title 'voronoi' w linespoints, @file3 using 1:3 title 'RCB' w linespoints, @file4 using 1:3 title '1D parallel' w linespoints, @file5 using 1:3 title 'skew boundary' w linespoints

