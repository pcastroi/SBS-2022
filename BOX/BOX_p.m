%% BOX Exercise - Structure-Borne Sound
clc; clear all;
addpath(genpath('data'))

% Read data
ASFolder = 'data\Plate Average\';
IMFolder = 'data\InputMobility\';
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
%%
close all;
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
Y_t = sum(NUM_tm,2)./sum(DEN_tm,2); % Eq. 12 in ProjectPlans
eta = Ps(:,2)./(2*pi*f_h(:,2).*sum(DEN_lf,2)); % Eq. 13 in ProjectPlans
Y_00 = sqrt(12*(1-v^2))/(8*h(1)^2*sqrt(rho*E)); % Eq. 15 in ProjectPlans
cl = sqrt(E/(rho*(1-v^2))); % cl = cl5
Y_tc = Y_00./(2*pi.*f_h(:,2).*mean(M(1:4)).*eta.*(1+(0.003*S(5)*0.003)./(eta.*sum(S(1:4))*0.003))); % Eq. 16 in ProjectPlans

[third_freq,Y_00_third]=onethirdoctave_average(f_i(:,2),value_h(:,2));
[third_freq,Y_t_third]=onethirdoctave_average(f_v(:,idx_tm),Y_t);
[third_freq,eta_third]=onethirdoctave_average(f_h(:,2),eta);
[third_freq,Y_tc_third]=onethirdoctave_average(f_h(:,2),Y_tc);

% Plots
figure
semilogx(f_v(:,idx_tm),mag2db(Y_t))
hold on
semilogx(third_freq,mag2db(Y_t_third),LineWidth=2)
% title('Average transfer mobility (panels 1 to 4), $\left\langle\left|Y_t\right|^2\right\rangle_{1-4}$',Interpreter='latex')
xlim([10^2, max(f_i(:,2))])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Broadband','1/3 octave band',Location='best')

figure
semilogx(f_h(:,2),abs(eta))
hold on
semilogx(third_freq,abs(eta_third),LineWidth=2)
scatter([101, 171.25, 313.5, 339.375, 339.375, 739],[9.90E-03, 7.30E-03, 1.52E-02, 1.40E-02, 1.40E-02, 8.96E-03],'filled','black','LineWidth',2) % 3dB bandwidth
% title('System loss factor (panels 1 to 5)')
xlim([10^2, max(f_i(:,2))])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Broadband','1/3 octave band','3-dB bandwidth',Location='best')

figure
[minValue,closestIndex] = min(abs(f_h(:,2)-100)); % find idx of f = 100Hz
semilogx(f_h(:,2),mag2db(Y_00*ones(1,length(value_h(:,2)))),LineWidth=2,LineStyle="--")
hold on
semilogx(f_h(:,2),mag2db(abs(value_h(:,2))))
semilogx(third_freq,mag2db(abs(Y_00_third)),LineWidth=2)
semilogx(f_h(:,2),mag2db(mean(abs(value_h(:,2)))*ones(1,length(value_h(:,2)))),LineWidth=2,LineStyle="--")
legend('Theory','Measurement','Measurement 1/3 octave band','Measurement average',Location='best')
% title('Input Mobility, Y_{00}')
xlim([10^2, max(f_h(:,2))]) 
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

figure
semilogx(f_v(:,idx_tm),mag2db(abs(Y_tc)))
hold on
semilogx(third_freq,mag2db(abs(Y_tc_third)),LineWidth=2)
semilogx(f_v(:,idx_tm),mag2db(abs(Y_t)))
semilogx(third_freq,mag2db(Y_t_third),LineWidth=2)
xlim([10^2, max(f_i(:,2))])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Theory','Theory 1/3 octave band','Measurement','Measurement 1/3 octave band',Location='best')

