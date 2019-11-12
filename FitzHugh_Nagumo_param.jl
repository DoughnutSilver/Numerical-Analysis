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
#function generate_g(model)  #meta                                      
#    func="function g(u ::Float64,v ::Float64) \n
#            return $(model.γ)*(u-($model.α)*v-$(model.β))     \n
#        end     \n"
#    func=Meta.parse(func)
#    eval(func)
#    eval(func)
#    eval(func)
#    eval(func)
#    eval(func)
#    eval(func)
#    return g
#end
#=----------------------------------------------------------------------------------
    Generating g in tihs way is not efficient, because jit compile doesn't work
----------------------------------------------------------------------------------=#
using  Base.Threads


function f(u::Float64,v::Float64)
    return u-u^3-v
end

function g(u ::Float64,v ::Float64, model)
    return model.γ*(u-model.α*v-model.β)     
end

function NL_f(u::Array{Float64,2},v::Array{Float64,2})
    out=similar(u)
    for i in  eachindex(u)  i:: Int64
        out[i]=f(u[i],v[i])
    end
    return out
end

function NL_g(u::Array{Float64,2},v::Array{Float64,2}, model)
    out=similar(u)
    for i in  eachindex(u)  i:: Int64
        out[i]=g(u[i],v[i],model)
    end
    return out
end

function diffusion(system,distribution)
    out=similar(distribution)
    @. @inbounds out=-4*distribution
    #diff x
        @. @inbounds out[1:system.size-1,:]+=@view distribution[2:system.size,:]
        @. @inbounds out[system.size,:]+=@view distribution[1,:]
        @. @inbounds out[2:system.size,:]+=@view distribution[1:system.size-1,:]
        @. @inbounds out[1,:]+=@view distribution[system.size,:]
    #diff y
        @. @inbounds out[:,1:system.size-1]+=@view distribution[:,2:system.size]
        @. @inbounds out[:,system.size]+=@view distribution[:,1]
        @. @inbounds out[:,2:system.size]+=@view distribution[:,1:system.size-1]
        @. @inbounds out[:,1]+=@view distribution[:,system.size]
    return  out*system.invdx
end

function onestep_u(model,system,variable)
    @inbounds variable.u+=system.dt*(model.Du*diffusion(system,variable.u)+NL_f(variable.u,variable.v))
    return variable.u
end

function onestep_v(model,system,variable)
    @inbounds variable.v+=system.dt*(model.Dv*diffusion(system,variable.v)+NL_g(variable.u,variable.v,model))
    return variable.v
end

function onestep(model,system,variable)
    variable.u=onestep_u(model,system,variable)
    variable.v=onestep_v(model,system,variable)
    return variable
end

#=----------------------------------------------------------------------------------
構造体
----------------------------------------------------------------------------------=#

struct model
    Du :: Float64
    Dv :: Float64
    α  :: Float64
    β  :: Float64
    γ  :: Float64
end

struct system
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

mdl=model(1.0*10^(-4) ,1.0*10^(-2),0.5,0.00,26.0)
sys=system(64,10,75,10^4,10^(-4)) #dx^2>dt  #0.50,0.04がよかった
u=rand(sys.size,sys.size)
v=rand(sys.size,sys.size)
var=variable(u,v)
#=----------------------------------------------------------------------------------
計算
----------------------------------------------------------------------------------=#

x_range=1:1:sys.size
y_range=1:1:sys.size

function calc(model,system,variable)
    for time in 1:500000 time::Int64
        variable.u=onestep_u(model,system,variable)
        variable.v=onestep_v(model,system,variable)
    end
end

function speed(model,system,variable)
    @time calc(model,system,variable)
end

function file_out(out,sys)
    file=open("julia.dat","w")
    for i in 1:sys.size
        for j in 1:sys.size
        print(file,out[i,j],"      ")
        end
        println(file,"")
    end
    close(file)
end

function gnu_plot()
    run(pipeline(`gnuplot`,stdin="FitzHugh_Nagumo_param.plt"))
end




#speed(mdl,sys,var)
#file_out(var.u-var.v,sys)
#gnu_plot()

#make_png(x_range,y_range,mdl,sys,var)

#speed(mdl,sys,var)

#make_gif(x_range,y_range,mdl,sys,var)

