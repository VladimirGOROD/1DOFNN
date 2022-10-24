using Pkg;
if lowercase(Base.active_project()) != lowercase(@__DIR__()*"\\Project.toml")
    Pkg.activate(".");
end
using Jurvis;       
using PyPlot;
using LinearAlgebra
using Peaks;
dt = 1e-4;
t = (0.0:dt:5.0); # time frames

damp1(t) = 0.002 *(2- (t/30)^2); 
ζ1 = damp1.(t);

freq1(t)=(exp.((-0.7*t))+50);
f1 = freq1.(t);

ω1 = 2*π*f1;
ϕ1 = cumsum(ω1) * dt;

A1 = 100*exp.(-ζ1 .* ω1 .* t);

signal = A1.*sin.(ϕ1);


# plot(signal);
# plot(A1)


calc_env = envelope(signal)
# plot(calc_env)


# using BasicInterpolators;
max_pks, max_vals = findmaxima(signal)
min_pks, min_vals = findminima(signal)
# plot(t, signal); plot(t, A1)
# plot()
env = decomposeSSA(calc_env, 10, 1000)[:,1];
# plot(t, env)
phase = instphase(signal)
# phase = decomposeSSA(instphase(signal), 10, 1000)[:,1];
# plot(phase)


freq = decomposeSSA(findiff(phase, dt; order = 4), 10,1000)[:,1];
typeof(freq)

function calcfreq(env::AbstractVector{<:Number}, time::AbstractVector{<:Number}, step::Number)
    Ȧ = findiff(env, step; order = 2);
    g = - Ȧ ./ env;
    βₐₚ = Vector{Float64}(undef, length(time));
    for i in 2:length(time);
        βₐₚ[i] = (step * time[i] / ( step + time[i])) * ( βₐₚ[i-1] / step + g[i] / time[i]);
    end
    return βₐₚ
 end
 


B=calcfreq(env,t,dt) #произведение ω и ζ
typeof(B)
Da=B./freq;

plot(Da)
plot(ζ1)
# plot(damp)
# plot(ω1)