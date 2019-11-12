using Plots
gr()
using ImplicitEquations


function f(u,v)
    return u-u^3-v
end

function g(u,v)         #関数をパラメータ化できない
    return 3(u-4*v-4)
end



plot(f⩵0)
plot!(g⩵0)

