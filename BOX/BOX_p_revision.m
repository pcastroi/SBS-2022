%% BOX Exercise - Structure-Borne Sound
clear all; clc; close all;
addpath(genpath('data'))

% Read data
for d=1
ASFolder = 'data/Plate Average/';
IMFolder = 'data/InputMobility/';
VelocityFiles=dir([ASFolder,'*Velocity*.txt']);
ForceFiles=dir([ASFolder,'*Force*.txt']);
H1Files=dir([IMFolder,'*H1*.txt']);
IMobilFiles=dir([IMFolder,'*Auto*.txt']);
for i=1:length(VelocityFiles)
    [band_v(:,i),f_v(:,i),value_v(:,i)]=read_pulse_2021(VelocityFiles(i).name);
    [band_f(:,i),f_f(:,i),value_f(:,i)]=read_pulse_2021(ForceFiles(i).name);
    if i<=length(H1Files)
        [band_h(:,i),f_h(:,i),value_h(:,i)]=read_pulse_2021(H1Files(i).name);
        [band_i(:,i),f_i(:,i),value_i(:,i)]=read_pulse_2021(IMobilFiles(i).name);
    end
end
end

%% Calculate transfer mobility
for d=1
% Constants
rho = 7.8*10^3; % [kg/m^3], m/V, density (Table 5.1 Note 7016)
v = 0.28; % Poisson's Ratio (Table 5.1 Note 7016)
E = 2*10^11; % [N/m^2], Young's modulus (Table 5.1 Note 7016)
h = [3 3 3 3 1.5 3]*10^-3; % [m], thickness of plates
S = [0.109 0.109 0.176 0.176 0.246 0.246]; % [m], surface area of plates
M = rho.*h.*S;

% Average transfer mobility and system loss factor
idx_tm = 4;
idx_lf = 5;
NUM_tm = zeros(length(value_v),idx_tm);
DEN_tm = zeros(length(value_f),idx_tm);
DEN_lf = zeros(length(value_f),idx_lf);
for i=1:idx_lf
    DEN_lf(:,i)=value_v(:,i)*M(i);
    if i<=idx_tm
        NUM_tm(:,i)=value_v(:,i)*M(i);
        DEN_tm(:,i)=value_f(:,i)*M(i);
    end
end
Ps = value_h.*value_i;
Y_t = sum(NUM_tm,2)./sum(DEN_tm,2); % Eq. 12 in ProjectPlans System 1 transfer mobility
eta = Ps(:,2)./(2*pi*f_h(:,2).*sum(DEN_lf,2)); % Eq. 13 in ProjectPlans
Y_00 = sqrt(12*(1-v^2))/(8*h(1)^2*sqrt(rho*E)); % Eq. 15 in ProjectPlans
cl = sqrt(E/(rho*(1-v^2))); % cl = cl5
%Y_tc = Y_00./(2*pi.*f_h(:,2).*mean(M(1:4)).*eta.*(1+(0.003*S(5)*0.003)./(eta.*sum(S(1:4))*0.003))); % Eq. 16 in ProjectPlans, average transfer mobility

for d=1 %Transfer mobility revision
    eta_t=0.01; Cl_t= cl; M_t=sum(M(1:4)); S_t=sum(S(1:4)); h_t=mean(h(1:4)); 
    eta5_t= 0.003; Cl5_t = cl; M5_t=M(5); S5_t=S(5); h5_t=h(5);
    A5_t = (eta5_t*S5_t) / (h5_t*Cl5_t);
    A_t = (eta_t*S_t) / (h_t*Cl_t);
    w_t=2*pi.*f_h(:,2);
Y_tc = Y_00./(w_t.*M_t.*eta_t.*(1+A5_t/A_t)); % Eq. 16 in ProjectPlans, average transfer mobility
end

[third_freq,Y_00_third]=onethirdoctave_average(f_i(:,2),value_h(:,2));
[third_freq,Y_t_third]=onethirdoctave_average(f_v(:,idx_tm),Y_t);
[third_freq,eta_third]=onethirdoctave_average(f_h(:,2),eta);
[third_freq,Y_tc_third]=onethirdoctave_average(f_h(:,2),Y_tc);
end

%% System 1 transfer mobility plot
%{
figure
semilogx(f_v(:,idx_tm),mag2db(Y_t))
hold on
semilogx(third_freq,mag2db(Y_t_third),LineWidth=2)
% title('Average transfer mobility (panels 1 to 4), $\left\langle\left|Y_t\right|^2\right\rangle_{1-4}$',Interpreter='latex')
xlim([10^2, max(f_i(:,2))])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Broadband','1/3 octave band',Location='best')
title('System 1 transfer mobility')
%}

%% System damping factor plot
%{
figure
semilogx(f_h(:,2),abs(eta))
hold on
semilogx(third_freq,abs(eta_third),LineWidth=2)
%scatter([101, 171.25, 313.5, 339.375, 339.375, 739],[9.90E-03, 7.30E-03, 1.52E-02, 1.40E-02, 1.40E-02, 8.96E-03],Color='r',Marker='x') % 3dB bandwidth
scatter([101, 171.25, 313.5, 339.375, 339.375, 739],[9.90E-03, 7.30E-03, 1.52E-02, 1.40E-02, 1.40E-02, 8.96E-03])% 3dB bandwidth
% title('System loss factor (panels 1 to 5)')
xlim([10^2, max(f_i(:,2))])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Broadband','1/3 octave band',Location='best')
title('System damping factor')
%}

%% Input mobility plot
%{
figure
[minValue,closestIndex] = min(abs(f_i(:,2)-100)); % find idx of f = 100Hz
semilogx(f_i(:,2),mag2db(Y_00*ones(1,length(value_i(:,2)))),LineWidth=2,LineStyle="--")
hold on
semilogx(f_i(:,2),mag2db(value_h(:,2)))
semilogx(third_freq,mag2db(Y_00_third),LineWidth=2)
legend('Theory','Measurement','Measurement 1/3 octave band',Location='best')
% title('Input Mobility, Y_{00}')
xlim([10^2, max(f_i(:,2))]) 
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
tite('Input mobility')
%}

%% Average transfer mobility plot
%
figure
semilogx(f_v(:,idx_tm),mag2db(abs(Y_tc)),LineWidth=2)
hold on
%semilogx(third_freq,mag2db(abs(Y_tc_third)),LineWidth=2)
semilogx(f_v(:,idx_tm),mag2db(Y_t))
semilogx(third_freq,mag2db(Y_t_third),LineWidth=2)
xlim([10^2, max(f_i(:,2))])
xlabel('Frequency [Hz]','FontSize',14)
ylabel('Magnitude [dB]','FontSize',14)
legend('Theory','Measurement','Measurement 1/3 octave band',Location='best',FontSize=12)
title('Average transfer mobility',FontSize=14)
%}
