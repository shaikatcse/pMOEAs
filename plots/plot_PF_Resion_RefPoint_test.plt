set terminal png 
set output "PF_Region_RefPoint.png"

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
my_val=(x1_r-x1_l)/20
splot [x1_l+my_val:x1_r-my_val] [0:1] x/(xh-xl)
unset table


probName = "ZDT1"

set xrange [0:1.1]
set yrange [0:1.2]


set arrow from x1_l, 0 to x1_l, 1.0 nohead lc rgb 'red' front 
set arrow from x1_r, 0 to x1_r, 1.0 nohead lc rgb 'red' front

unset colorbox
set palette defined (0 "#8888ff", 1 "#ffffff")
plot [0.8:0.95] [0:1] 'shadowkey.dat' w ima, \
"\\Points.txt" using 1:2 t "Point" lc rgb "blue"

unset yrange
unset yrange
unset terminal
unset output
unset label