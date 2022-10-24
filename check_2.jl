using PyPlot;
dt = 0.001;
T = 10;
t = 0: dt: T;

ζ = 0.01;

damp1(t) = 0.002 *(2- (t/30)^2); 
ζ = damp1.(t);

ω = 80;

freq1(t)=(exp.((-0.7*t))/3+15);
f1 = freq1.(t);

ω = 2*π*f1;
plot(t, ω)

β = ζ .* ω;

A₀= 5;
A = A₀*exp.(-β.*t);

x = A.*sin.(ω.*t)

plot(t, x)
plot(t, A)

using Jurvis
Ȧ = findiff(A, dt; order = 2);
g = - Ȧ ./ A


βₐₚ = Vector{Float64}(undef, length(t));

βₐₚ[1] =  - Ȧ[1] / A[1];

for i in 2:length(t)
    βₐₚ[i] = (dt * t[i] / ( dt + t[i])) * ( βₐₚ[i-1] / dt + g[i] / t[i]);
end

plot(t, β); plot(t, βₐₚ)
plot(t, βₐₚ)