using Pkg;
if lowercase(Base.active_project()) != lowercase(@__DIR__()*"\\Project.toml")
    Pkg.activate(".");
end
using Jurvis;       
using PyPlot;
using LinearAlgebra
dt = 1e-4;
t = (0.0:dt:5.0); # time frames
p=[-0.0030, -0.02375, 3.775, 0.5, 50, 100, 80 ];
damp1(t) = 0.002 *(2- (t/30)^2); 

ζ = damp1.(t);
damp2(t) = 0.008;
ζ2 = damp2.(t);
freq1(t)=(exp.((-0.7*t))+50);
g = freq1.(t);
freq2(t)=3*exp((-t*1.5)-p[4])+87;
g2 = freq2.(t);
# plot(t,g2);
ω1 = 2*π*g;
ϕ1 = cumsum(ω1) * dt;
ω2 = 2*π*g2;
ϕ2 = cumsum(ω2) * dt;
A1 = p[6]*exp.(-ζ.*ϕ1);
A2 =p[7].*exp.(-ζ2.*ϕ2);
θ1 = A1.*sin.(ϕ1);
θ2=A2.*sin.(ϕ2.+π/8);
θ=θ1+θ2;
plot(θ);


##

SSA_win = 2000;
signal = MeasuredData("test data", t, hcat(θ, θ1, θ2), "time", ["test signal", "mode1", "mode2"], size(θ,1), 3, ["", "", ""]);
decomposed_data = decomposeSSA(signal,1, 10, SSA_win);
plotall(spectrumall(decomposed_data)); yscale("log"); legend();xlim([0,200])
plotall(decomposed_data)
# figure();
# plot(decomposed_data[:,])
# Data=MeasuredData(50001, 0);
# begin
    # for i=1:10
    # add_data!(Data,decomposed_data[:,i]); 
    # end 
# end
# Data.xdata=t;
#sp_modes = spectrumall(Data);
#plotchannel(sp_modes;ch=1)
#plot(sp_modes.ydata[:,1])
# plotall(sp_modes)
# plot(sp_modes.ydata[:,1])
#  yscale("log")
# plot(Data.ydata[:,3])


groups = [[2,3,6,7,10, 11],
            [4,5,8,9]
            ];

grouped_modes = groupmodes(decomposed_data, groups)


gr_modes = copydata(grouped_modes, [1,2,3])
# gr_modes = groupmodes(Data, groups);#содержит исходный сигнал и две первые сгруп моды
mode1 = copydata(gr_modes, [2]);
#  plot(mode1.ydata[:,1]);
env1=envelope(mode1.ydata[:,1]);
env_modes1 = decomposeSSA(env1, 10, 1000);
add_data!(mode1, env_modes1[:,1]; ydata_name = "Огибающая первой моды"); #канал 2
phase1 = instphase(mode1.ydata[:,1]);
add_data!(mode1, phase1; ydata_name = "inst phase моды 1"); #канал 3
# plot(mode1.ydata[:,1])
inst_freq1 = findiff(mode1.ydata[:,3], t[2] - t[1]; order = 4) / (2*π);
freq_modes1 = decomposeSSA(inst_freq1, 10,1000);
plot(inst_freq1);plot(freq_modes1[:,1]); plot(g)
ζ1 = 100*instdamping(mode1.ydata[:,2], freq_modes1[:,1]);
add_data!(mode1, ζ1; ydata_name = "Коэффициент демпфирования первой моды"); #канал 4

begin
    fig, ax = plt.subplots(1,2, figsize = [12,6]);
    # ax[1].sharex(ax[2]);

    ax[1].plot(mode1.ydata[:,2],freq_modes1[:,1]);
    ax[1].set_xlabel("Амплитуда, [м/с]");
    ax[1].set_ylabel("Частота, [Гц]");
    # ax[1].set_ylim([1.335,1.38]);
    #   ax[1].set_xlim([0.2,5.25]);
    ax[1].set_title("Мгновенная собственная частота-амплитуда для моды 1");
  
    ax[2].plot(mode1.ydata[:,2],mode1.ydata[:,4]);
    ax[2].set_xlabel("Амплитуда, [м/с]");
    ax[2].set_ylabel("Демпфирование, [%]");
      # ax[2].set_ylim([0.0038,0.004]);
    # ax[2].set_xlim([0.5,5.2]);
    ax[2].set_title("Демпфирование-амплитуда");
end

mode2 = copydata(gr_modes, [3]);

env2=envelope(mode2.ydata[:,1]);
env_modes2 = decomposeSSA(env2, 10, 1000);
add_data!(mode2, env_modes2[:,1]; ydata_name = "Огибающая первой моды"); #канал 2
phase2 = instphase(mode2.ydata[:,1]);
add_data!(mode2, phase2; ydata_name = "inst phase моды 1"); #канал 3

inst_freq2 = findiff(mode2.ydata[:,3], t[2] - t[1]; order = 4) / (2*π);
freq_modes2 = decomposeSSA(inst_freq2, 10,1000);
ζ2_2 = 100*instdamping(mode2.ydata[:,2], freq_modes2[:,1]);
add_data!(mode2, ζ2_2; ydata_name = "Коэффициент демпфирования первой моды"); #канал 4
##
begin

plot(t,A1)
plotchannel(mode1;ch=2)
figure()
plot(t,A2)
plotchannel(mode2;ch=2)

figure()
plot(t,ζ)
plotchannel(mode1;ch=4)
ylim([0.0038,0.004]);
figure()
plot(t,ζ2)
plotchannel(mode2;ch=4)
ylim([0,0.01])
end


figure()
plot(t,ζ1)
plotchannel(mode1;ch=4)
ylim([0.0038,0.004]);


X=Vector{Float64}([])#вектор для точек времени перехода через 0
begin
  for i=1:(mode1.length-1)
    if mode1.ydata[i,1] < 0 && mode1.ydata[i+1,1] >0
      y=0;
      y_a=mode1.ydata[i,1];
      y_b=mode1.ydata[i+1,1];
      x_a=t[i];
      x_b=t[i+1];
      x_0=((y-y_a)/(y_b-y_a))*(x_b-x_a)+x_a;
      push!(X,x_0);

    end

  end
end

function calcfreq(sig::AbstractVector{<:Number}, time::AbstractVector{<:Number})
  X=Vector{Float64}([]);
  for i in 1: size(sig,1)-1
    if sig[i] < 0 && sig[i+1] >0
      y=0;
      y_a=sig[i];
      y_b=sig[i+1];
      x_a=time[i];
      x_b=time[i+1];
      x_0=((y-y_a)/(y_b-y_a))*(x_b-x_a)+x_a;
      push!(X,x_0);
    end
  end
  return X
end

t1 = calcfreq(θ1, t)
f1 = 1 ./findiff(t1,1; order = 4)
instfreq = 1 ./ findiff(X,1; order=4)

par = LinRange(0,1,size(t,1))
gg = 51 .* (1 .- par) + 40 .* (par)
ω = g*2*pi;
ϕ = cumsum(ω) * dt

testsig = sin.(ϕ);
t2 = calcfreq(testsig, t)
f2 = 1 ./ findiff(t2, 1; order = 1)
# plot(testsig)
# testsig2 = sin.(2*pi*51 .*t)
# t2 = calcfreq(testsig2, t)
# f2 = 1 ./ findiff(t2, 1; order = 4)

plot(t,g); plot(X,instfreq); plot(t1,f1); plot(t2  ,f2, "--"); ylim([49,52])

figure();plot(testsig)
plot(t,testsig)