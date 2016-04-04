#set terminal pdf enhanced
#set output "PF_Region_ES.pdf"

set terminal eps enhanced
set output "PF_Region_ES.eps"


#regions
#region  1
x1_l = 0.40
x1_r = 0.50
#region  2
x2_l = 0.0
x2_r = 0.15
#region 3
x3_l = -0.50
x3_r = -0.40


set xlabel "CO_2 emission (in million  tons)"
set ylabel "Annual cost (in  million DKK)"

set xrange [-0.6:0.65]
set yrange [2500:6500]

set label "Region # 1" at 0.455, 4500 rotate left
set label "Region # 2" at 0.075, 4500 rotate left
set label "Region # 3" at -0.45, 4500 rotate left

set arrow from x1_l, 2500 to x1_l, 6000 nohead lc rgb 'red'
set arrow from x1_r, 2500 to x1_r, 6000 nohead lc rgb 'red'

set arrow from x2_l, 2500 to x2_l, 6000 nohead lc rgb 'green'
set arrow from x2_r, 2500 to x2_r, 6000 nohead lc rgb 'green'

set arrow from x3_l, 2500 to x3_l, 6000 nohead lc rgb 'blue'
set arrow from x3_r, 2500 to x3_r, 6000 nohead lc rgb 'blue'

set style line 1 lt 1 lw 3 pt 3 linecolor rgb "red"

plot "C:\\Users\\mahbub\\Documents\\GitHub\\EnergyPLANDomainKnowledgeEAStep1\\IntegratedMOEAResults\\ParetoFront\\mergeFUN.pf" using 1:2 t "Pareto-front" lc rgb "gray", "..\\FUN_ES_seed_831779" using 1:2 t "Preferred regional Pareto-front" lc rgb "red"

unset yrange
unset yrange
unset terminal
unset output
unset label