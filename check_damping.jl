using Pkg;
if lowercase(Base.active_project()) != lowercase(@__DIR__()*"\\Project.toml")
    Pkg.activate(".");
end
using Jurvis;       
using PyPlot;
using LinearAlgebra
dt = 1e-4;
t = (0.0:dt:5.0); # time frames

damp1(t) = 0.002 *(2- (t/30)^2); 
ζ1 = damp1.(t);

freq1(t)=(exp.((-0.7*t))+50);
f1 = freq1.(t);

ω1 = 2*π*g;
ϕ1 = cumsum(ω1) * dt;

A1 = 100*exp.(-ζ1 .* ω1 .* t);

signal = A1.*sin.(ϕ1);


plot(signal);
plot(A1)


calc_env = envelope(signal)
plot(calc_env)

# using Peaks;
# using BasicInterpolators;
max_pks, max_vals = findmaxima(signal)
min_pks, min_vals = findminima(signal)
plot(t, signal); plot(t, A1)

env = decomposeSSA(calc_env, 10, 1000)[:,1];
plot(t, env)
phase = instphase(signal)
# phase = decomposeSSA(instphase(signal), 10, 1000)[:,1];
plot(phase)


freq = decomposeSSA(findiff(phase, dt; order = 4), 10,1000)[:,1]

damp = instdamping(A1, f1)
plot(ζ1)
plot(damp ./ ζ1)
plot(ω1)