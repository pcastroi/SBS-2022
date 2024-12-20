%% BEND Exercise - Structure-Borne Sound
close all; clc; clear all;
%% Part 1 - Flexural vibration of a beam
% Parameters
N = 10000;
f = linspace(0,12500,N); % [Hz]
t = 1./f; % [s]
w = 2*pi*f; % [rad/s]
m_u = 0.1877; % [kg], mass of undamped beam
m_d = 0.2083; % [kg], mass of damped beam
h_u = 0.003; % [m], thickness of undamped beam
h_d = 0.0048; % [m], thickness of damped beam
b = 0.025; % [m], width of beams
L = 0.3; % [m], lenght of beams
l = L/2; % [m], lenght of half-beams
mpu_u = m_u/L; % [kg/m], mass per unit length of undamped beam
mpu_d = m_d/L; % [kg/m], mass per unit length of damped beam
rho_u = m_u/(L*b*h_u); % [kg/m^3], m/V, density of undamped beam
rho_d = m_d/(L*b*h_d); % [kg/m^3], m/V, density of damped beam
eta = [0.1, 0.03, 0.01]; % Damping Loss Factor
I_u = b*h_u^3/12; % Area moment of inertia of undamped beam
I_d = b*h_u^3/12+2*(b*((h_d-h_u)/2)^3/12); % Area moment of inertia of damped beam
E_1 = 10.3*10^10; % Young's modulus of undamped beam
E_2 = E_1; % Young's modulus of damped beam
B_u = E_1*I_u; % [Pa], Bending stiffness of undamped beam
B_d = E_2*I_d; % [Pa], Bending stiffness of damped beam
C_n = [3.561, 9.815, 19.244, 31.808, 47.518]; % Constants (BC) for resonances
A_n = [2.240, 14.028, 39.280, 76.976, 127.231]; % Constants (BC) for anti-resonances
n = length(C_n);

% Calculating resonant and anti-resonant frequencies
for i=1:n
    f_n_u(i)=C_n(i)/(L^2)*sqrt(B_u/mpu_u);
    f_n_d(i)=C_n(i)/(L^2)*sqrt(B_d/mpu_d);
    f_a_n_u(i)=A_n(i)/(L^2)*sqrt(B_u/mpu_u);
    f_a_n_d(i)=A_n(i)/(L^2)*sqrt(B_d/mpu_d);
end

%% Calculate the mobilities
Y_ll_u=zeros(length(eta),N);
Y_ll_u_the=Y_ll_u;
Y_ll_d=Y_ll_u;
Y_ll_d_the=Y_ll_u;

for i=1:length(eta)
    % Undamped beam
    cb_u(i,:) = sqrt(w)*(B_u/mpu_u)^(1/4); % [m/s], Phase speed (bending) in undamped beam
    k_u(i,:) = w./cb_u(i,:)*(1-1i*eta(i)/4); % [1/m], complex wavenumber
    N_u_ll(i,:) = cos(k_u(i,:)*l).*cosh(k_u(i,:)*l)+1;
    N_u_0l(i,:) = cos(k_u(i,:)*l)+cosh(k_u(i,:)*l);
    D_u(i,:) = cos(k_u(i,:)*l).*sinh(k_u(i,:)*l)+sin(k_u(i,:)*l).*cosh(k_u(i,:)*l);
    Y_ll_u(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_u(i,:)*l*sqrt(B_u*mpu_u))).*N_u_ll(i,:)./D_u(i,:);
    Y_ll_u_the(i,:)=-((1i*l)./(2*k_u(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_1*I_u*(1+1i*eta(i))*mpu_u))).*N_u_ll(i,:)./D_u(i,:); % Theoretical form
    Y_0l_u(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_u(i,:)*l*sqrt(B_u*mpu_u))).*N_u_0l(i,:)./D_u(i,:);
    Y_0l_u_the(i,:)=-((1i*l)./(2*k_u(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_1*I_u*(1+1i*eta(i))*mpu_u))).*N_u_0l(i,:)./D_u(i,:); % Theoretical form
    % Damped beam
    cb_d(i,:) = sqrt(w)*(B_d/mpu_d)^(1/4); % [m/s], Phase speed (bending) in damped beam
    k_d(i,:) = w./cb_d(i,:)*(1-1i*eta(i)/4); % [1/m], complex wavenumber
    N_d_ll(i,:) = cos(k_d(i,:)*l).*cosh(k_d(i,:)*l)+1;
    N_d_0l(i,:) = cos(k_d(i,:)*l)+cosh(k_d(i,:)*l);
    D_d(i,:) = cos(k_d(i,:)*l).*sinh(k_d(i,:)*l)+sin(k_d(i,:)*l).*cosh(k_d(i,:)*l);
    Y_ll_d(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_d(i,:)*l*sqrt(B_d*mpu_d))).*N_d_ll(i,:)./D_d(i,:);
    Y_ll_d_the(i,:)=-((1i*l)./(2*k_d(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_2*I_d*(1+1i*eta(i))*mpu_d))).*N_d_ll(i,:)./D_d(i,:); % Theoretical form
    Y_0l_d(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_d(i,:)*l*sqrt(B_d*mpu_d))).*N_d_0l(i,:)./D_d(i,:);
    Y_0l_d_the(i,:)=-((1i*l)./(2*k_d(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_2*I_d*(1+1i*eta(i))*mpu_d))).*N_d_0l(i,:)./D_d(i,:); % Theoretical form
end

colors = ['#ff8809'; '#ff6cc7';'#00a6ff'];

%% Part 2 - Plot txt files
addpath BEND/data
Alltxt=dir('BEND/data/*.txt');
%Alltxt=dir('*.txt');
for i=1:length(Alltxt)
    [band(:,i),f_p2(:,i),value(:,i)]=read_pulse_2021(Alltxt(i).name);
%     plot_everything(Alltxt(i),f_p2(:,i),value(:,i)); % Uncomment if need
end

%% Plots


for i=1:length(eta)
    %|Y_{ll}|
    figure(1)
    plot(f,mag2db(abs(Y_ll_u(i,:))),Color=colors(i,:))
    hold on
    plot(f,mag2db(abs(Y_ll_d(i,:))),Color=colors(i,:),LineStyle="--")
    xlim([min(f) max(f)])
    title('|Y_{ll}|')
    ylabel('Magnitude [dB]')
    xlabel('Frequency [Hz]')
    if i==3;legend(['Undamped Beam: \eta = ' num2str(eta(1),'%.2e\n')],['Damped Beam: \eta = ' num2str(eta(1),'%.2e\n')],['Undamped Beam: \eta = ' num2str(eta(2),'%.2e\n')],['Damped Beam: \eta = ' num2str(eta(2),'%.2e\n')],['Undamped Beam: \eta = ' num2str(eta(3),'%.2e\n')],['Damped Beam: \eta = ' num2str(eta(3),'%.2e\n')],Location='best');end
    grid on

%     %'\angle Y_{ll}'
%     figure(2)
%     plot(f,angle(Y_ll_u(i,:)),Color=colors(i,:))
%     hold on
%     plot(f,angle(Y_ll_d(i,:)),Color=colors(i,:),LineStyle="--")
%     xlim([min(f) max(f)])
%     title('\angle Y_{ll}')
%     ylabel('Phase [rad]')
%     xlabel('Frequency [Hz]')
%     if i==3;legend(['Undamped Beam: \eta = ' num2str(eta(1),'%.2e\n')],['Damped Beam: \eta = ' num2str(eta(1),'%.2e\n')],['Undamped Beam: \eta = ' num2str(eta(2),'%.2e\n')],['Damped Beam: \eta = ' num2str(eta(2),'%.2e\n')],['Undamped Beam: \eta = ' num2str(eta(3),'%.2e\n')],['Damped Beam: \eta = ' num2str(eta(3),'%.2e\n')],Location='best');end
%     grid on

%     %|Y_{ll}| vs |Y_{0l}|
%     figure(3)
%     plot(abs(k_u*l),mag2db(abs(Y_ll_u(i,:))),Color=[0 0 0 1*i^-3])
%     hold on
%     plot(abs(k_u*l),mag2db(abs(Y_0l_u(i,:))),Color=[1 0 0 1*i^-3])
%     plot(abs(k_d*l),mag2db(abs(Y_ll_d(i,:))),Color=[0 0 0 1*i^-3],LineStyle="--")
%     plot(abs(k_d*l),mag2db(abs(Y_0l_d(i,:))),Color=[1 0 0 1*i^-3],LineStyle="--")
%     title('|Y_{ll}| vs |Y_{0l}|')
%     ylabel('Magnitude [dB]')
%     xlabel('Helmholtz Number kl')
%     if i==1;legend('|Y_{ll}| Undamped Beam','','','|Y_{0l}| Undamped Beam','','','|Y_{ll}| Mobility Damped Beam','','','|Y_{0l}| Damped Beam','','');end
%     grid on
    
    if i==3
        % 11 and 23 idx (undamped and damped measurements)
        figure(1)
        plot(f_p2(:,11),mag2db(abs(value(:,11))),Color='k')
        hold on
        plot(f_p2(:,23),mag2db(abs(value(:,23))),Color='k',LineStyle="--")
        xlim([min(f) max(f)])
        title('Measured |Y_{ll}|')
        ylabel('Magnitude [dB]')
        xlabel('Frequency [Hz]')
        legend('Undamped beam','Damped beam')
        grid on

%         %'\angle Y_{ll}'
%         figure(2)
%         plot(f_p2(:,12),angle(value(:,12)),Color='k')
%         hold on
%         plot(f_p2(:,25),angle(value(:,25)),Color='k',LineStyle="--")
%         xlim([min(f) max(f)])
%         title('\angle Y_{ll}')
%         ylabel('Phase [rad]')
%         xlabel('Frequency [Hz]')
%         legend('Undamped beam','Damped beam')
%         grid on
    end
end

xlim([0 4000])
set(gca,'FontSize',14)
%% Functions
function plot_everything(Alltxt,f_p2,value)
    if contains(Alltxt.name,'Coherence')
        figure
        plot(f_p2,value)
        title([ strrep(strrep(Alltxt.name,'_',' '),'.txt','') ': C_{xy}'])
        ylabel('Coherence')
        xlabel('Frequency [Hz]')
        xlim([min(f_p2) max(f_p2)])
        ylim([0 1])
        grid on
    else
        figure
        plot(f_p2,mag2db(abs(value)))
        ylabel('Magnitude [dB]')
        xlabel('Frequency [Hz]')
        xlim([min(f_p2) max(f_p2)])
        title([ strrep(strrep(Alltxt.name,'_',' '),'.txt','') ': |Y_{ll}|'])
        grid on
        figure
        plot(f_p2,angle(value))
        ylabel('Phase [rad]')
        xlabel('Frequency [Hz]')
        xlim([min(f_p2) max(f_p2)])
        title([ strrep(strrep(Alltxt.name,'_',' '),'.txt','') ': \angle Y_{ll}'])
        grid on
    end
end
