cd 'C:\Users\mao\Desktop\Numerical-Analysis'
set terminal png
set output 'FithHugh_Nagumo_param.png'
set xrange [0:64]
set yrange [0:64]
set view map
plot "julia.dat" matrix w image
unset output