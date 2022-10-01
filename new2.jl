using Pkg;
if lowercase(Base.active_project()) != lowercase(@__DIR__()*"\\Project.toml")
    Pkg.activate(".");
end
using DifferentialEquations; 
using Jurvis;       
using PyPlot;
using LinearAlgebra

# dev https://github.com/kutsjuice/Jurvis.jl
function reduce_order(m_M, m_C, m_K)
    """
    function reduces order of system of 2nd-order ODE by creating matrix A wich define system of 1st-order ODE"""
    
    @assert size(m_M) == size(m_K) == size(m_C) "Matrices have different size"
    @assert size(m_M)[1] == size(m_M)[2] "Matrices are not square"
    N = size(m_M)[1]
    m_A = Matrix{Float64}(undef, 2N, 2N)
    m_A[1:N, 1:N] = fill!(m_A[1:N, 1:N], 0.);
    m_A[1:N, N+1:2N] = diagm(ones(N));
    m_A[N+1:2N, 1:N] = -m_M\m_K;
    m_A[N+1:2N, N+1:2N] = -m_M\m_C;        
    return m_A
end

# Load function definition
# in this function we take fisrt half of sin function and scale it.
function force1(t, ampl; start=0, duration=0.1)
    if start<t<(start+duration)
        return ampl*(sin((t-start)*π/duration));
    else
        return 0;
    end
end
# in this function we take full period of sin function, translate and scale it
function force2(t, ampl; start=0, duration=0.1)
    if start<t<(start+duration)
        return 0.5*ampl*(-cos((t-start)*2π/duration) + 1);
    else
        return 0;
    end
end
begin
    dur = 10e-3;
    f_lim = 100;
    tt = 0:dur/50:dur*10;#L"f(x) = \frac {1}{x+x^{2}}"

    plot_title = """Force plot with exitation amplitude 300N and duration 0.1s"""
    
    frc = force2.(tt, 2.0; duration = dur);

    fig, ax = subplots(2,1);
    ax[1].plot(tt, frc, linewidth=2);

    ax[1].set_xlabel("Time, [s]");
    ax[1].set_ylabel("Force, [N]");
    ax[1].set_xlim([0, 2*dur]);
    frc_spec = easyspectrum(frc; dt = tt[2]-tt[1])
    ax[2].plot(frc_spec[1], frc_spec[2])
    ax[2].set_xlabel("Frequency, [Hz]");
    ax[2].set_ylabel("Amplitude, [N]");
    ax[2].set_xlim([0,f_lim]); 
    fig.tight_layout()
    
end
# du - ветор производных, u - вектор состояний, p - вектор параметров.
# Делайм нелинейность в коэффициент p[12] и p[13] записаны коэффициенты а, для двух случаев нелинейности
# p[5] в функции изменяем на p[5]*(1+2*p[5]*(u[2]-u[1]))
#p[5] в функции изменяем на p[5]*(1+3*p[6]*(u[2]-u[1])^2)

# Define the parameters of the system and problem to solve
p = [0.1, 30, 0.05, 500, 0.3, 0.05,  0.001, 0.03, 0.2];
function problem_3DOF!(du, u, p, t)
    m_M = p[1];
    m_K = p[2]*(1+3*p[6]*(u[1])^2);
    m_C = p[3];
    # m_A = reduce_order(m_M, m_C, m_K); 
    du[1] = u[2];
    du[2]=1/m_M*(force2(t, p[4]; duration = dur)-m_K*u[1] - p[8]*u[2] - sign(u[2])*p[9]);
    # du[:] = m_A*u + vcat([0., 0., 0.],m_M\[0,force2(t, p[10]; duration = dur), 0.]);
end
# +p[8]*u[2]-sign(u[2])*p[9])
# function problem_3DOF!(du, u, p, t)
# m_M = diagm(p[1:3])
# m_K = [sum(vcat(p[4:5],p[7])) p[5]*(1+3*p[12]*(u[2]-u[1])^2) -p[7];
# p[5]*(1+3*p[12]*(u[2]-u[1])^2) sum(p[5:6]) -p[6];
        # -p[7] -p[6] sum(p[6:7])];
# m_C = p[8].*m_M .+ p[9].*m_K;
# m_A = reduce_order(m_M, m_C, m_K); 
# du[:] = m_A*u + vcat([0., 0., 0.],m_M\[0,force2(t, p[10]; duration = dur), 0.]);
# end

# Parameters for the solver 
u0 = zeros(Float64, 2); # initial state
ft = 5e2; # sampling rate
dt = 1/ft; # sampling period
# u0[1] = 0.02;  
tspan = (0.0, 20.0); # time frames

prob = ODEProblem(problem_3DOF!,u0,tspan, p);
# ODEProblem принимает функцию описания дифур, нач сост, промежуток, набор параметров
sol = solve(prob, AutoVern7(Rodas5()), reltol=1e-8, abstol=1e-8, dtmax=p[7]/10, saveat = dt, maxiters = 1e7);
# get the time vector and state vector of the solution
u = hcat(sol.u...)';  t = vcat(sol.t...);

ind=t.>0;

Data=MeasuredData("Data", t[ind], u[ind,:],"Time, sec", ["d1","v1"], sum(ind), 2, ["d1","v1"])

d = 1e3*Data.ydata[:, 1]; # get displacements and scale it to mm
# v = u[:, 4:6]; # get velocities
# calculate accelerations using finite difference scheme of 4th accuracy order

# a = Vector{Float64}(undef, size(d,1),1);
a=findiff(Data.ydata[:,2], t[2] - t[1]; order=4);
add_data!(Data, a);
# Теперь MD содержит d, v, a
 
Data.ydata_names=["d1","v1","a1"];
# Spec=spectrumall(Data);
Data.ydata_names=["d1","v1","a1"];


#add_data!- позволяет добавить канал в MD
#n - временные точки, m - кол-во измерительных каналов
# записать решение дифура в MD, найти перемещения, скорость, ускорения и записать в один MD
# в SPL сунуть матрицу со все, 


begin
    fig, ax = plt.subplots(2,2, figsize = [12,6])
    ax[1,1].sharex(ax[1,2])
    ch = 1;
    ax[1,1].plot(Data.xdata, Data.ydata[:,ch]);
    ax[1,1].set_xlabel("Time, "*L"[s]");
    ax[1,1].set_ylabel("Displacements, "*L"[mm]")
    d_spec = easyspectrum(Data.ydata[:,ch]; dt = dt);
    ax[2,1].plot(d_spec[1], d_spec[2]);
    ax[2,1].set_xlabel("Friquency, "*L"[Hz]");
    ax[2,1].set_ylabel("Displacements, "*L"[mm]")
    ax[2,1].set_yscale("log")

    ax[2,1].sharex(ax[2,2])
    ax[2,1].set_xlim([0,250])
    ax[1,2].plot(Data.xdata,Data.ydata[:,ch+2])
    ax[1,2].set_xlabel("Time, "*L"[s]");
    ax[1,2].set_ylabel("Acceleration, "*L"[m/s^2]");
    a_spec = easyspectrum(Data.ydata[:,ch+2]; dt = dt);
    ax[2,2].plot(a_spec[1], a_spec[2]);2
    ax[2,2].set_xlabel("Friquency, "*L"[Hz]");
    ax[2,2].set_ylabel("Amplitude, "*L"[m/s^2]");
    ax[2,2].set_yscale("log")
    fig.tight_layout()
end
##



##


env1=envelope(Data.ydata[:,1]); #строит огибающую
env=vec(env1);# сделал вектор что бы подставить в decomposeSSA
ogib = decomposeSSA(env, 10,1000);# в первом канале сглаженная огибающая 
# вычисление мгновенной фазы сигнала
phase1 = instphase(Data.ydata[:,1]);
# вычисление мгновенной частоты пример
inst_freq = findiff(phase1, t[2] - t[1]; order = 4) / (2*π);

plot(env1)
plot(ogib[:,1])
plot(Data.ydata[:,1])
# коэф демпфирования
freq_modes = decomposeSSA(inst_freq, 10,500);
 plot(ogib[:,1],freq_modes[:,1])
ζ = 100*instdamping(ogib[:,1], freq_modes[:,1]);

##
plot(ogib[:,1])
plot(ζ)

figure()
plot(freq_modes[:,1])
# построение графиков
begin
    fig, ax = plt.subplots(1,2, figsize = [12,6]);
    # ax[1].sharex(ax[2]);

    ax[1].plot(ogib[:,1],freq_modes[:,1]);
    ax[1].set_xlabel("Амплитуда, [м/с]");
    ax[1].set_ylabel("Частота, [Гц]");
    ax[1].set_ylim([2.65,3.1]);
    ax[1].set_xlim([0.2,1.25]);
    ax[1].set_title("Мгновенная собственная частота-амплитуда");

    ax[2].plot(ogib[:,1],ζ);
    ax[2].set_xlabel("Амплитуда, [м/с]");
    ax[2].set_ylabel("Демпфирование, [%]");
    ax[2].set_ylim([0,0.5]);
    ax[2].set_xlim([0.16,1.3]);
    ax[2].set_title("Демпфирование-амплитуда");
end







