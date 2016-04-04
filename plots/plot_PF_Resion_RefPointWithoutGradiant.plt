set terminal pdf
set output "PF_Region_RefPoint.pdf"

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


probName = "ZDT1"

set xrange [0:1.1]
set yrange [0:1.2]

set title probName
set label "Region # 1" at 0.875, 0.65 rotate left
set label "Region # 2" at 0.450, 0.65 rotate left
set label "Region # 3" at 0.175, 0.65 rotate left

set arrow from x1_l, 0 to x1_l, 1.0 nohead lc rgb 'red'
set arrow from x1_r, 0 to x1_r, 1.0 nohead lc rgb 'red'

set arrow from x2_l, 0 to x2_l, 1.0 nohead lc rgb 'green'
set arrow from x2_r, 0 to x2_r, 1.0 nohead lc rgb 'green'

set arrow from x3_l, 0 to x3_l, 1.0 nohead lc rgb 'blue'
set arrow from x3_r, 0 to x3_r, 1.0 nohead lc rgb 'blue'

set style line 1 lt 1 lw 3 pt 3 linecolor rgb "red"

plot ".\\paretoFronts\\ZDT1.pf" using 1:2 t "Pareto-front", ".\\paretoFronts\\RefPoints.txt" using 1:2 t "Reference Point" lc rgb "blue"

unset yrange
unset yrange
unset terminal
unset output
unset label