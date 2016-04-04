reset
f(v,T)=4/sqrt(pi*T**3)*exp(-v**2/T)*v**2
set sample 500
set title "Maxwell speed distribution"
set xlabel "Speed [(2kT_c/m)^{1/2}]"
set ylabel "Probability Intensity f(v)"
set style fill transparent solid 0.5
set key Left
set term png truecolor enhanced font "Times,15"
set output "maxwell_speed_distribution.png"
plot [0:5] f(x,0.1) w filledcurves x1 title "T=0.1T_c",\
           f(x,1)   w filledcurves x1 title "T=T_c",\
           f(x,5)   w filledcurves x1 title "T=5T_c"
set term pdf enhanced
set output "maxwell_speed_distribution.pdf"
replot