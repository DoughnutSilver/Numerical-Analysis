#using DifferentialEquations
#using ParameterizedFunctions
using Plots
gr()
ENV["PLOTS_TEST"]="true"

# FitzHugh-南雲##########################################

#########################################################
function diffusion(dx,dt,sys_size,variable)
    out=zeros(sys_size,sys_size)
    out=-4*copy(variable)
    #diff x
    out+=[variable[2:sys_size,:] ; variable[1,:]']
    out+=[variable[sys_size,:]' ; variable[1:sys_size-1,:]]
    #diff y
    out+=[variable[:,2:sys_size] variable[:,1]]
    out+=[variable[:,sys_size] variable[:,1:sys_size-1]]
    invdx=1.0/dx^2
    return out*invdx
end

function onetime(dx,dt,sys_size,variable)
    variable =variable+dt*diffusion(dx,dt,sys_size,variable)
    return variable
end

function oneshot(interval,dx,dt,sys_size,variable)
    for time in 1:interval
        variable+=dt*onetime(dx,dt,sys_size,variable)
    end
    return variable
end

function make_frame(x_range,y_range,interval,dx,dt,sys_size,variable)
    variable=oneshot(interval,dx,dt,sys_size,variable)
    heatmap(x_range,y_range,variable,size=(2*sys_size,2*sys_size))
end
####################################################################################
println("input size")
sys_size=parse(Int16, readline())
println("size is :   ",sys_size)

x_range=1:1:sys_size
y_range=1:1:sys_size

u=zeros(sys_size,sys_size)
u[1,:]=fill(10.0,(1,sys_size))
u[:,1]=fill(-10.0,(sys_size,1))

dx=1.0
dt=1.0
interval=100
shot=50
#######################################################################################

#=
anim = @animate for i=1:shot
    global u=oneshot(interval,dx,dt,sys_size,u)     #forの中では新しくlocal scopeを作るのでglobalにしないと値の更新がされない。（adapさんに教えてもらった）
    plt=heatmap(x_range,y_range,u,size=(5*sys_size,5*sys_size))
    frame(anim,plt)
end
=#

anim=@animate for i in 1:shot
    global u=oneshot(interval,dx,dt,sys_size,u)
end

#heatmap(x_range,y_range,u)
mp4(anim,"C:\\Users\\mao\\Desktop\\数値計算\\tmp\\diffusion.mp4",fps=15)
