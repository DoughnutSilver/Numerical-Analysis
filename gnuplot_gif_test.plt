set terminal gif animate delay 100 optimize size 640,480     #delay=100 is fps 1000/10
set output 'test.gif'
set xrange [-2*pi:2*pi] 
set yrange [-1:1] 
set samples 10000
do for [i=1:100]{
    plot cos(i*x)
}
unset output