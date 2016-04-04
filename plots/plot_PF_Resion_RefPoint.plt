#set terminal pdf
#set output "PF_Region_RefPoint.pdf"

set terminal eps enhanced
set output "PF_Region_RefPoint.eps"


#regions
#region  1
x1_l = 0.8
x1_r = 0.95
#region  2
x2_l = 0.4
x2_r = 0.5
#region 3
x3_l = 0.15
x3_r = 0.2

xl=0; xh=1; yl=0; yh=1;

set table 'shadowkey.dat'
splot [xl:xh] [yl:yh] x/(xh-xl)
unset table


probName = "ZDT1"

set xrange [0:1.1]
set yrange [0:1.3]

set title probName
set label "Region # 1" at 0.875, 0.65 rotate left front
set label "Region # 2" at 0.450, 0.65 rotate left front
set label "Region # 3" at 0.175, 0.65 rotate left front

set arrow from x1_l, 0 to x1_l, 1.0 nohead lc rgb 'red' front 
set arrow from x1_r, 0 to x1_r, 1.0 nohead lc rgb 'red' front

set arrow from x2_l, 0 to x2_l, 1.0 nohead lc rgb 'green' front 
set arrow from x2_r, 0 to x2_r, 1.0 nohead lc rgb 'green' front

set arrow from x3_l, 0 to x3_l, 1.0 nohead lc rgb 'blue' front 
set arrow from x3_r, 0 to x3_r, 1.0 nohead lc rgb 'blue' front

unset colorbox
set palette defined (0 "#8888ff", 1 "#ffffff")
#With JV reference points plot [0:1.0]  'shadowkey.dat' w ima notitle, "..\\paretoFronts\\ZDT1.pf" using 1:2 t "Pareto-front", "..\\paretoFronts\\RefPoints.txt" using 1:2 t "Reference Point" lc rgb "black", "..\\rNSGAII_FUN_seed_831779" using 1:2 t "rNSGAII result" lc rgb "green"  pt 7 ps 1  

plot [0:1.0]  'shadowkey.dat' w ima notitle, "..\\paretoFronts\\ZDT1.pf" using 1:2 t "Pareto-front", "..\\paretoFronts\\RefPoints.txt" using 1:2 t "Reference Point" lc rgb "black", "..\\rNSGAII_FUN_seed_831779" using 1:2 t "rNSGAII result" lc rgb "green"  
unset yrange
unset yrange
unset terminal
unset output
unset label