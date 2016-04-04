set terminal pdf
set output "pMOEA_results_test.pdf"

set multiplot

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

#xl=0; xh=1; yl=0; yh=1;

#set table 'shadowkey.dat'
#splot [xl:xh] [yl:yh] x/(xh-xl)
#unset table


probName = "ZDT1"

set xrange [0:1.0]
set yrange [0:1.2]

#set title probName
set label "Region # 1" at 0.875, 0.65 rotate left front
set label "Region # 2" at 0.450, 0.65 rotate left front
set label "Region # 3" at 0.175, 0.65 rotate left front

set arrow from x1_l, 0 to x1_l, 1.0 nohead lc rgb 'red' front 
set arrow from x1_r, 0 to x1_r, 1.0 nohead lc rgb 'red' front

set arrow from x2_l, 0 to x2_l, 1.0 nohead lc rgb 'green' front 
set arrow from x2_r, 0 to x2_r, 1.0 nohead lc rgb 'green' front

set arrow from x3_l, 0 to x3_l, 1.0 nohead lc rgb 'blue' front 
set arrow from x3_r, 0 to x3_r, 1.0 nohead lc rgb 'blue' front

#unset colorbox
#set palette defined (0 "#8888ff", 1 "#ffffff")
#With JV reference points plot [0:1.0]  'shadowkey.dat' w ima notitle, "..\\paretoFronts\\ZDT1.pf" using 1:2 t "Pareto-front", "..\\paretoFronts\\RefPoints.txt" using 1:2 t "Reference Point" lc rgb "black", "..\\rNSGAII_FUN_seed_831779" using 1:2 t "rNSGAII result" lc rgb "green"  pt 7 ps 1  

plot "..\\paretoFronts\\ZDT1.pf" using 1:2 t "Pareto-front" lc rgb "grey", "..\\experiments\\pAGEStudyZDT\\data\\pAGE.F\\ZDT1\\FUN.1" using 1:2 t "pAGEOffline" lc rgb "red", "..\\experiments\\pNSGAIIStudyZDT\\data\\pNSGAII\\ZDT1\\FUN.1" using 1:2 t "pNSGAII" lc rgb "green"  


#now set option for the smaller plot
unset arrow
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
set size 0.4,0.3
set origin 0.16, 0.5
#set xtics 0.02
#set ytics 0.05
set xrange [0.15:0.20]
set yrange[0.55:0.65]


plot "..\\experiments\\pAGEStudyZDT\\data\\pAGE.F\\ZDT1\\FUN.1" using 1:2 lc rgb "red" notitle, "..\\experiments\\pNSGAIIStudyZDT\\data\\pNSGAII\\ZDT1\\FUN.1" using 1:2  lc rgb "green"  notitle

unset multiplot
unset yrange
unset yrange
unset xlabel
unset ylabel

unset size
unset origin
unset xtics

unset output
unset terminal
