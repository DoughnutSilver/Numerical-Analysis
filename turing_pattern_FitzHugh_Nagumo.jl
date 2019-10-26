#=----------------------------------------------------------------------------------
This program makes turing pattern by using FitzHugh-Nagumo-model
<Turing-Model>
    du/dt   =Du △^2u + f(u,v)
    dv/dt   =Dv △^2v + g(v,v)
<FitzHugh-Nagumo-model>
    f(u,v)=u-u^3-v
    g(u,v)=γ(u-αv-β)

<Stability condition>
    DΔt/Δx^2 <1/2    
----------------------------------------------------------------------------------=#
#using DifferentialEquations
#using ParameterizedFunctions
using Plots
using LinearAlgebra
gr()

#=----------------------------------------------------------------------------------
関数
----------------------------------------------------------------------------------=#
function f(u::Float64,v::Float64)
    return u-u^3-v
end

function g(u::Float64,v::Float64)
    α =0.5
    β =0.0
    γ =26.0
    return γ*(u-α*v-β)
end

function diffusion(system_parameter,distribution)
    out=zeros(system_parameter.size,system_parameter.size)
    out=-4*copy(distribution)
    #diff x
        @inbounds out+=[distribution[2:system_parameter.size,:]   ; distribution[1,:]']
        @inbounds out+=[distribution[system_parameter.size,:]'    ; distribution[1:system_parameter.size-1,:]]
    #diff y
        @inbounds out+=[distribution[:,2:system_parameter.size]   distribution[:,1]]
        @inbounds out+=[distribution[:,system_parameter.size]     distribution[:,1:system_parameter.size-1]]
    local   invdx=1.0/system_parameter.dx^2
    return  out*invdx
end

function onestep_u(model_parameter,system_parameter,variable)
    @inbounds variable.u+=system_parameter.dt*(model_parameter.Du*diffusion(system_parameter,variable.u)+f.(variable.u,variable.v))
    return variable.u
end

function onestep_v(model_parameter,system_parameter,variable)
    @inbounds variable.v+=system_parameter.dt*(model_parameter.Dv*diffusion(system_parameter,variable.v)+g.(variable.u,variable.v))
    return variable.v
end

function onestep(model_parameter,system_parameter,variable)
    variable.u=onestep_u(model_parameter,system_parameter,variable)
    variable.v=onestep_v(model_parameter,system_parameter,variable)
    return variable
end

function oneshot(x,y,model_parameter,system_parameter,variable)
    for i=1:system_parameter.interval
        onestep(model_parameter,system_parameter,variable)
    end
    heatmap(x,y,variable.u-variable.v,size=(400,400))
end

function make_mp4(x,y,model_parameter,system_parameter,variable)
    out=@animate for i=1:system_parameter.shot
        oneshot(x,y,model_parameter,system_parameter,variable)
    end
    mp4(out,"tmp//FitzHugh_Nagumo.mp4",fps=15)
end

function make_gif(x,y,model_parameter,system_parameter,variable)
    out=@animate for i=1:system_parameter.shot
        oneshot(x,y,model_parameter,system_parameter,variable)
    end
end

#=----------------------------------------------------------------------------------
構造体
----------------------------------------------------------------------------------=#

struct model_parameter
    Du :: Float64
    Dv :: Float64
end

struct system_parameter
    size :: Int16
    interval :: Int32
    shot :: Int32
    dx   :: Float64
    dt   :: Float64
end

mutable struct variable
    u   :: Array{Float64,2}
    v   :: Array{Float64,2}
end


model=model_parameter(1.0*10^(-4) ,1.0*10^(-2))
sys=system_parameter(128,10,75,10^(-2),10^(-4)) #dx^2>dt  #0.50,0.04がよかった
u=rand(sys.size,sys.size)
v=rand(sys.size,sys.size)
var=variable(u,v)

#=----------------------------------------------------------------------------------
計算
----------------------------------------------------------------------------------=#

x_range=1:1:sys.size
y_range=1:1:sys.size

#plt=heatmap(x_range,y_range,var.u)
#display(plt)

function calc(model_parameter,system_parameter,variable)
    for time in 1:1000000
        global variable.u=onestep_u(model_parameter,system_parameter,variable)
        global variable.v=onestep_v(model_parameter,system_parameter,variable)
    end
end

@timev calc(model,sys,var)

#make_gif(x_range,y_range,model,sys,var)

