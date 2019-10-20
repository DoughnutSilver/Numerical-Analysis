using Plots
using LinearAlgebra
gr()


β=10^10
spin=0.5
step=10^7

println("input size")
size=parse(Int16, readline())
println("size is :   ",size)


spin_range=-spin:1:spin

site=rand(spin_range,size,size)

function transition(energy,site,next)
    if(energy>0 && rand()>exp(-β*energy))
        return site 
    else
        return next
    end
end

for i in 1:step
    x=rand(1:size)
    y=rand(1:size)
    next=rand(spin_range)
    if(1<x<size)
        xl=x+1
        xr=x-1
    elseif(x==1)
        xl=2
        xr=size
    else
        xl=1
        xr=x-1
    end
    if(1<y<size)
        yl=y+1
        yr=y-1
    elseif(y==1)
        yl=2
        yr=size
    else
        yl=1
        yr=y-1
    end
    energy=-(next-site[x,y])*(site[xl,y]+site[xr,y]+site[x,yl]+site[x,yr])
    site[x,y]=transition(energy,site[x,y],next) 
end

x=1:1:size
y=1:1:size
heatmap(x,y,site)
