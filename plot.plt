set terminal pdf
set output "pZDT6.pdf"

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


probName = "ZDT6"

set xrange [0:1.1]
set yrange [0:1.1]

set title probName

set arrow from x1_l, 0 to x1_l, 1.0 nohead lc rgb 'red'
set arrow from x1_r, 0 to x1_r, 1.0 nohead lc rgb 'red'

set arrow from x2_l, 0 to x2_l, 1.0 nohead lc rgb 'green'
set arrow from x2_r, 0 to x2_r, 1.0 nohead lc rgb 'green'

set arrow from x3_l, 0 to x3_l, 1.0 nohead lc rgb 'blue'
set arrow from x3_r, 0 to x3_r, 1.0 nohead lc rgb 'blue'

#plot "C:\\Users\\mahbub\\Documents\\GitHub\\jMetal\\jmetal-problem\\src\\test\\resources\\pareto_fronts\\ZDT3.pf" using 1:2 linetype rgb "#C0C0C0" title "True front", "FUNN.tsv" using 1:2 lt rgb "red" title "NSGAII", "FUN_PN.tsv" using 1:2 lt rgb "green" title "pNSGAII", "FUN_PS.tsv" using 1:2 lt rgb "blue" title "pSPEA2"

#plot "C:\\Users\\mahbub\\Documents\\GitHub\\jMetal\\jmetal-problem\\src\\test\\resources\\pareto_fronts\\ZDT6.pf" using 1:2 linetype rgb "#C0C0C0" title "True front", "FUN" using 1:2 lt rgb "red" title "NSGAII_4.5 version (600 evaluations)", "FUN_norN600.tsv" using 1:2 lt rgb "green" title "NSGAII_5.0 version (600 evaluations)","FUN.tsv" using 1:2 lt rgb "blue" title "NSGAII45 version (600 evaluations)" 

plot "C:\\Users\\mahbub\\Documents\\GitHub\\EnergyPLANDomainKnowledgeEAStep1\\paretoFronts\\ZDT6.pf" using 1:2 t "true front", "experiments\\pNSGAIIStudy\\data\\pNSGAII\\ZDT6\\FUN.0" using 1:2 t "pNSGAII", "experiments\\pNSGAIIStudy\\data\\NSGAII1\\ZDT6\\FUN.0" using 1:2 t "NSGAII1"

unset yrange
unset terminal
unset output