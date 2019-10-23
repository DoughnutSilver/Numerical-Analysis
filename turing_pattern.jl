using DifferentialEquations
using ParameterizedFunctions
using Plots
gr()

# FitzHugh-南雲##########################################

#########################################################


function sp_diff_2(f,dx,sys_size)
    nabla_square=-2*f
    inverse_dx=1/dx
    for site in 2:sys_size-1
        nabla_square[site] += f[site+1]+f[site-1]
    end
    nabla_square[1] += f[sys_size]+f[2]
    nabla_square[sys_size] += f[sys_size-1]+f[1]
    return inverse_dx*nabla_square
end

function diffusion(f,dx,sys_size)
    output=fill(0.0,size(f)) #output=fill(0,size(f))だとoutputが整数列になってしまってうまくいかなかった
    for site in 1:sys_size
        output[site,:]+=sp_diff_2(f[site,:],dx,sys_size)
    end
    for site in 1:sys_size
        output[:,site]+=sp_diff_2(f[:,site],dx,sys_size)
    end
    return output
end

function onestep(site,dx,sys_size)
    site+=0.01*diffusion(site,dx,sys_size)
    return site
end

function oneshot(site,interval,dx,sys_size)
    for t in 1:interval
        site=onestep(site,dx,sys_size)
    end
    return site
end

function make_anim(sys_size,site)
    for i in 1:shot
        oneshot(site,interval,dx,sys_size)
    end
    heatmap(x_range,y_range,site,size=(6*sys_size,6*sys_size))   
end

# 系の作成
println("input size")
sys_size=parse(Int16, readline())
println("size is :   ",sys_size)

x_range=1:1:sys_size
y_range=1:1:sys_size

u=fill(0,(sys_size,sys_size))
u[1,:]=fill(1000,(1,sys_size))


dt=2^(-10)
dx=2^(-10)
interval=1000
shot=1000



out=@animate for i=1:shot
    make_anim(sys_size,u)
end

gif(out, "test.gif", fps = 15)
