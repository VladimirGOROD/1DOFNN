using Pkg;
if lowercase(Base.active_project()) != lowercase(@__DIR__()*"\\Project.toml")
    Pkg.activate(".");
end
using DifferentialEquations; 
using Jurvis;       
using PyPlot;
using LinearAlgebra

t = (0.0:1e-4:5.0); # time frames
p=[-0.0030, -0.02375, 3.775, 0.5, 50, 100, 80 ];
e=2.7182;
damp1(t) = 0.002 *(2- (t/30)^2); 
ζ = damp1.(t);
damp2(t) = 0.008;
ζ2 = damp2.(t);
damp3(t)=(exp.((-0.7*t))+50);
g = damp3.(t);
damp4(t)=3*exp((-t*1.5)-p[4])+87;
g2 = damp4.(t);
plot(t,g2);
ω1 = 2*π*g;
ω2 = 2*π*g2;
A1 = p[6]*exp.(-ζ.*ω1.*t)
θ1 = A1.*sin.(ω1.*t)
θ2=p[7]*sin.(ω2.*t.+π/8).*e.^(-ζ2.*ω2.*t)
θ=θ1+θ2
plot(t,A1)

L=size(t)
SSA_win = 2000;
decomposed_data = decomposeSSA(θ, 10, SSA_win);
# figure();
# plot(decomposed_data[:,])
Data=MeasuredData(50001, 0);
begin
    for i=1:10
    add_data!(Data,decomposed_data[:,i]); 
    end 
end
Data.xdata=t;
sp_modes = spectrumall(Data);
plotchannel(sp_modes;ch=1)
plot(sp_modes.ydata[:,1])
# plotall(sp_modes)
# plot(sp_modes.ydata[:,1])
  yscale("log")
# plot(Data.ydata[:,3])


groups = [[2,5,6,9],
            [3,4,7,8,10]
            ];

grsp_modes = groupmodes(sp_modes, groups);
plot(sp_modes.ydata[:,10])
begin
    figure()
    plotall(grsp_modes)
    # plot(sp_modes.ydata[:,1])
    yscale("log")
    title("");
end

gr_modes = groupmodes(Data, groups);#содержит исходный сигнал и две первые сгруп моды
mode1 = copydata(gr_modes, [2]);
#  plot(mode1.ydata[:,1]);
env1=envelope(mode1.ydata[:,1]);
env_modes1 = decomposeSSA(env1, 10, 1000);
add_data!(mode1, env_modes1[:,1]; ydata_name = "Огибающая первой моды"); #канал 2
phase1 = instphase(mode1.ydata[:,1]);
add_data!(mode1, phase1; ydata_name = "inst phase моды 1"); #канал 3

inst_freq1 = findiff(mode1.ydata[:,3], t[2] - t[1]; order = 4) / (2*π);
freq_modes1 = decomposeSSA(inst_freq1, 10,1000);
ζ1 = 100*instdamping(mode1.ydata[:,2], freq_modes1[:,1]);
add_data!(mode1, ζ1; ydata_name = "Коэффициент демпфирования первой моды"); #канал 4
plot(t,ζ1)
plot(freq_modes1[:,1])
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
      ax[2].set_ylim([0.0038,0.004]);
    # ax[2].set_xlim([0.5,5.2]);
    ax[2].set_title("Демпфирование-амплитуда");
end
plotchannel(mode1;ch=2)
plot(t,A1)
plot(ζ1)
figure()
plot(,g)