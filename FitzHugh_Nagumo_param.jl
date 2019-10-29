#=----------------------------------------------------------------------------------
This program makes turing pattern by using FitzHugh-Nagumo-model
<Turing-Model>
    du/dt   =Du △^2u + f(u,v)
    dv/dt   =Dv △^2v + g(v,v)
<FitzHugh-Nagumo-model>
    f(u,v)=u-u^3-v
    g(u,v)=γ(u-αv-β)

<Stability condition>
    Let u_n=A^ne^{iknx} and substitute into Turing eq. Then we get
    Ae^{ikx}-1=DΔt(e^{ikx}-2+e^{-ikx})/Δx^2-A^{2n}e^{i2nkx}
    To make this program stable, we have to find |A|<1. Thus ignoring A^{2n} (<< A if |A|<1)
    Ae^{ikx}-1=DΔt(e^{ikx}-2+e^{-ikx})/Δx^2
    A=[DΔt(2cos(kx)-2)/Δx^2+1]e^{-ikx}
    A=[-D4sin^2(kx/2)Δt/Δx^2+1]e^{-ikx}
    |A|=[-D4sin^2(kx/2)Δt/Δx^2+1]
    if |A|<1 then 2>4Dsin^2(kx/2)Δt/Δx^2>0
    1/2>DΔt/Δx^2>0
    in tihs program 
    D=10^-2 so 
    Δt<Δx^2*50
    (Δt,Δx)=(1,0.25)0.125<0.5
----------------------------------------------------------------------------------=#
###using Plots
###using LinearAlgebra
###gr()

#=----------------------------------------------------------------------------------
関数
----------------------------------------------------------------------------------=#
function generate_g(model_parameter)
    func="function g(u ::Float64,v ::Float64) \n
            return $(model_parameter.γ)*(u-($model_parameter.α)*v-$(model_parameter.β))     \n
        end     \n"
    func=Meta.parse(func)
    eval(func)
    return g
end


function f(u::Float64,v::Float64)
    return u-u^3-v
end

function diffusion(system_parameter,distribution)
    out=similar(distribution)
    @. out=-4*distribution
    #diff x
        @. @inbounds out[1:system_parameter.size-1,:]+=@view distribution[2:system_parameter.size,:]
        @. @inbounds out[system_parameter.size,:]+=@view distribution[1,:]
        @. @inbounds out[2:system_parameter.size,:]+=@view distribution[1:system_parameter.size-1,:]
        @. @inbounds out[1,:]+=@view distribution[system_parameter.size,:]
    #diff y
        @. @inbounds out[:,1:system_parameter.size-1]+=@view distribution[:,2:system_parameter.size]
        @. @inbounds out[:,system_parameter.size]+=@view distribution[:,1]
        @. @inbounds out[:,2:system_parameter.size]+=@view distribution[:,1:system_parameter.size-1]
        @. @inbounds out[:,1]+=@view distribution[:,system_parameter.size]
    return  out*system_parameter.invdx
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
    out=@animate for i=1:system_parameter.shot i::Int32
        oneshot(x,y,model_parameter,system_parameter,variable)
    end
    mp4(out,"tmp//FitzHugh_Nagumo.mp4",fps=15)
end

function make_gif(x,y,model_parameter,system_parameter,variable)
    out=@animate for i=1:system_parameter.shot i::Int32
        oneshot(x,y,model_parameter,system_parameter,variable)
    end
end

#=----------------------------------------------------------------------------------
構造体
----------------------------------------------------------------------------------=#

struct model_parameter
    Du :: Float64
    Dv :: Float64
    α  :: Float64
    β  :: Float64
    γ  :: Float64
end

struct system_parameter
    size :: Int16
    interval :: Int32
    shot :: Int32
    invdx   :: Float64
    dt   :: Float64
end

mutable struct variable
    u   :: Array{Float64,2}
    v   :: Array{Float64,2}
end


model=model_parameter(1.0*10^(-4) ,1.0*10^(-2),0.5,0.05,26.0)
sys=system_parameter(64,10,75,10^4,10^(-4)) #dx^2>dt  #0.50,0.04がよかった
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

generate_g(model)

function calc(model_parameter,system_parameter,variable)
    for time in 1:500000 time::Int64
        variable.u=onestep_u(model_parameter,system_parameter,variable)
        variable.v=onestep_v(model_parameter,system_parameter,variable)
    end
end

function make_png(x_range,y_range,model_parameter,system_parameter,variable)
    calc(model_parameter,system_parameter,variable)
    contour(x_range,y_range,variable.u-variable.v,size=(400,400),fill=true)
end
function speed(model_parameter,system_parameter,variable)
    @time calc(model_parameter,system_parameter,variable)
end



#make_png(x_range,y_range,model,sys,var)

#speed(model,sys,var)

#make_gif(x_range,y_range,model,sys,var)

