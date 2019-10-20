using Plots
using LinearAlgebra
gr()

println("input dimentoin of system")
Sys_dim=parse(Int128, readline())
println("Dimention of the system is :   ",Sys_dim)
println("input size of BZ")
BZ_size=parse(Int128, readline())
println("size of BZ is :   ",BZ_size)

function BZ(n_st , size , dim)
    k_vector=zeros(Float64, dim, size )
    k_vector[1,:]=range(0,stop=2pi,length=size)
    for i in 2:dim
        k_vector[i,:]=k_vector[1,:]
    end
    k_vector=2pi*(n_st-1.5).+k_vector
    return k_vector
end

function Hamiltonian(wave_vector)
    wave_vector^2
end
first_BZ=BZ(1,BZ_size,Sys_dim)
second_BZ=BZ(2,BZ_size,Sys_dim)
third_BZ=BZ(0,BZ_size,Sys_dim)
plt=plot(first_BZ[1,:],Hamiltonian.(first_BZ[1,:])); plot!(first_BZ[1,:],Hamiltonian.(second_BZ[1,:]));plot!(first_BZ[1,:],Hamiltonian.(third_BZ[1,:]));

display(plt)