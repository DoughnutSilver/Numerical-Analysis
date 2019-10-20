using Plots
using LinearAlgebra
gr(
    legend=false,
)


β=10^10
spin=0.5
step=10^7

println("input size")
sys_size=parse(Int16, readline())
println("size is :   ",sys_size)


spin_range=-spin:1:spin
site=rand(spin_range,sys_size,sys_size)
x_range=1:1:sys_size
y_range=1:1:sys_size


function pick_up_site(sys_size)
    x=rand(1:sys_size)
    y=rand(1:sys_size)
    if(1<x<sys_size)
        xl=x+1
        xr=x-1
    elseif(x==1)
        xl=2
        xr=sys_size
    else
        xl=1
        xr=x-1
    end
    if(1<y<sys_size)
        yl=y+1
        yr=y-1
    elseif(y==1)
        yl=2
        yr=sys_size
    else
        yl=1
        yr=y-1
    end
    return x,xl,xr,y,yl,yr
end

function diff_energy(next,x,xl,xr,y,yl,yr)
    energy=-(next-site[x,y])*(site[xl,y]+site[xr,y]+site[x,yl]+site[x,yr])
    return  energy
end

function transition(energy,site,next)
    if(energy>0 && rand()>exp(-β*energy))
        return site 
    else
        return next
    end
end

function one_step_MC(sys_size,site)
    x,xl,xr,y,yl,yr=pick_up_site(sys_size)
    next=rand(spin_range)
    site[x,y]=transition(diff_energy(next,x,xl,xr,y,yl,yr),site[x,y],next) 
end

anim=Animation()

function make_anim(anim,sys_size,site)
    one_step_MC(sys_size,site)
    plt=heatmap(x_range,y_range,site,size=(6*sys_size,6*sys_size))   
    frame(anim,plt)
end
anim=Animation()
amin=@animate for i=1:step
    make_anim(anim,sys_size,site)
end

gif(anim, "C:\Users\mao\AppData\Local\Temp\anim_fps15.gif", fps = 15)
displayfile("C:\Users\mao\AppData\Local\Temp\","anim_fps15.gif")
