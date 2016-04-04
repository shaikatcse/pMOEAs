set terminal png

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

set xrange [0:1.1]
set yrange [0:6]

probName = "ZDT1"
do for [t=2:400] {
    
	outfile = sprintf('Gen%d.png',t)
	title = sprintf('ZDT1: Gen: %d', t)
	
	set output outfile	
	set title title
	
	

	set label "Region # 1" at 0.875, 0.65 rotate left
	set label "Region # 2" at 0.175, 0.65 rotate left

	set arrow from x1_l, 0 to x1_l, 6.0 nohead lc rgb 'red'
	set arrow from x1_r, 0 to x1_r, 6.0 nohead lc rgb 'red'
	
	set arrow from x3_l, 0 to x3_l, 6.0 nohead lc rgb 'blue'
	set arrow from x3_r, 0 to x3_r, 6.0 nohead lc rgb 'blue'
	
	Region1 = sprintf('Region1_%d',t)
	Region2 = sprintf('Region2_%d',t)

	plot Region1 using 1:2 t "Region # 2" lc rgb 'red', Region2 using 1:2 t "Region # 2" lc rgb 'blue'

	
	unset label
	
	unset output
}