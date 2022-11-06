using Jurvis;
using PyPlot;
dt = 0.001; #step
T = 10;
t = 0: dt: T;#time

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
A = A₀*exp.(-β.*t);#огибающая

x = A.*sin.(ω.*t)# сигнал

plot(t, x)
plot(t, A)





function calcdecaycoeff(env::AbstractVector{<:Number}, time::AbstractVector{<:Number}, step::Number)
    Ȧ = findiff(env, step; order = 2);
    g = - Ȧ ./ env;
    βₐₚ = Vector{Float64}(undef, length(time));
    βₐₚ[1] =  - Ȧ[1] / A[1];
    for i in 2:length(time);
        βₐₚ[i] = (step * time[i] / ( step + time[i])) * ( βₐₚ[i-1] / step + g[i] / time[i]);
    end
    return βₐₚ
end

β_calc = calcdecaycoeff(A,t,dt);

plot(t, β_calc)
plot(t, β)