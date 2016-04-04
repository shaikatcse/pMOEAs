#set terminal pdf enhanced 
#set output "AssignedRegions.pdf"	

set terminal eps enhanced 
set output "AssignedRegions.eps"	

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

set bmargin at screen 0.08

set xrange [-0.05:1.0]
set yrange [0.0:2.0]

   
set label "Region # 1" at 0.875, 1.15 rotate left
set label "Region # 2" at 0.175, 1.15 rotate left

set label "{/Symbol a_1}=6, {/Symbol a_2}=6 " at 0.42, 1.40 

set arrow from x1_l, 0 to x1_l, 1.7 nohead lc rgb 'red'
set arrow from x1_r, 0 to x1_r, 1.7 nohead lc rgb 'red'
	
set arrow from x3_l, 0 to x3_l, 1.7 nohead lc rgb 'blue'
set arrow from x3_r, 0 to x3_r, 1.7 nohead lc rgb 'blue'
	
unset xtics
unset ytics

set label "R_l^1" at x1_l-0.01, -0.06 
set label "R_u^1" at x1_r-0.01, -0.06

set label "{R_{l}}^2" at x3_l-0.01, -0.06
set label "{R_{u}}^{2}" at x3_r-0.01, -0.06


plot "..\\assignRegion\\Region1_229" using 1:2 t "Solutions associated with Region # 1" lc rgb 'red', "..\\assignRegion\\Region2_229" using 1:2 t "Solutions associated with Region # 2" lc rgb 'blue'

unset bmargin
unset xrange
unset yrange
unset label
unset output
