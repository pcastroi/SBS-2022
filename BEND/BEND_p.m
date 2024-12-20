%% BEND Exercise - Structure-Borne Sound
close all; clc; clear all;
%% Part 1 - Flexural vibration of a beam
% Parameters
N = 6401;
f = linspace(0,12500,N); % [Hz]
t = 1./f; % [s]
w = 2*pi*f; % [rad/s]
m_1 = 0.1877; % [kg], mass of undamped beam
m_2 = 0.2083; % [kg], mass of damped beam
h_1 = 0.003; % [m], thickness of undamped beam
h_2 = 0.0046; % [m], thickness of damped beam
h_3 = 0.00081; % [m], thickness of each damping layer
a = (h_1+h_3)/2; % [m], thickness of half damped beam
b = 0.025; % [m], width of beams
L = 0.3; % [m], lenght of beams
l = L/2; % [m], lenght of half-beams
mpu_1 = m_1/L; % [kg/m], mass per unit length of undamped beam
mpu_2 = m_2/L; % [kg/m], mass per unit length of damped beam
rho_1 = m_1/(L*b*h_1); % [kg/m^3], m/V, density of undamped beam
rho_2 = m_2/(L*b*h_2); % [kg/m^3], m/V, density of damped beam
rho_3 = 1670; % [kg/m^3], m/V, density of damping layer
eta = [0.1, 0.03, 0.01]; % Damping Loss Factor
I_1 = b*h_1^3/12; % Area moment of inertia of undamped beam
I_3 = b*h_3^3/12; % Area moment of inertia of damping layer
I_2 = b*h_1^3/12+h_3*a^2;
% I_2 = I_1+2*(I_3+h_3*b*a^2); % Area moment of inertia of damped beam
C_n = [3.561, 19.244, 47.518]; % Constants (BC) for resonances
A_n = [2.240, 39.280, 127.231]; % Constants (BC) for anti-resonances
n = length(C_n);
f_n_1 = [120.6, 649, 1608.06]; % Resonant f undamped - Given by the measurements
f_n_2 = [120, 639.25, 1564.12]; % Resonant f damped - Given by the measurements

% Calculating Young's modulus from measured resonant f
for i=1:n
    E_n_1(i)=(f_n_1(i)*L^2/C_n(i))^2*mpu_1/I_1*(1+1i*0.0007)^-1; % eta_1 = 0.0007 - 3dB method
    A(i)=(h_1^3*E_n_1(i))/(24*h_3*a^2); % ProjectPlans
    E_n_2(i)=A(i)*((f_n_2/f_n_1)^2*(1+(2*rho_3*h_3)/(rho_1*h_1))-1); % ProjectPlans

    B_n_1(i)=(f_n_1(i)*L^2/C_n(i))^2*mpu_1;
    B_n_2(i)=(f_n_2(i)*L^2/C_n(i))^2*mpu_2;
end
E_1 = mean(E_n_1); % Young's modulus of undamped beam
% E_2 = E_1; % Young's modulus of damped beam
E_2 = mean(E_n_2); % Young's modulus of damping layer
eta_2 = 0.021*(1+A/E_n_2) ; % eta_total = 0.021 - 3dB method
B_1 = E_1*I_1; % [Pa], Bending stiffness of undamped beam
B_2 = E_2*I_2; % [Pa], Bending stiffness of damped beam
B_2 = (E_1*h_1^3/12+E_2*h_2*a^2);
B_1 = mean(B_n_1);
B_2 = mean(B_n_2);

% Preallocate
Y_ll_1=zeros(length(eta),N);
Y_ll_1_the=Y_ll_1;
Y_ll_2=Y_ll_1;
Y_ll_2_the=Y_ll_1;

% eta=[0.0007 0.021];

for i=1:length(eta)
    % Undamped beam
    cb_1(i,:) = sqrt(w)*(B_1/mpu_1)^(1/4); % [m/s], Phase speed (bending) in undamped beam
    k_1(i,:) = w./cb_1(i,:)*(1-1i*eta(i)/4); % [1/m], complex wavenumber
    N_1_ll(i,:) = cos(k_1(i,:)*l).*cosh(k_1(i,:)*l)+1;
    N_1_0l(i,:) = cos(k_1(i,:)*l)+cosh(k_1(i,:)*l);
    D_1(i,:) = cos(k_1(i,:)*l).*sinh(k_1(i,:)*l)+sin(k_1(i,:)*l).*cosh(k_1(i,:)*l);
    Y_ll_1(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_1(i,:)*l*sqrt(B_1*mpu_1))).*N_1_ll(i,:)./D_1(i,:);
    Y_ll_1_the(i,:)=-((1i*l)./(2*k_1(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_1*I_1*(1+1i*eta(i))*mpu_1))).*N_1_ll(i,:)./D_1(i,:); % Theoretical form
    Y_0l_1(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_1(i,:)*l*sqrt(B_1*mpu_1))).*N_1_0l(i,:)./D_1(i,:);
    Y_0l_1_the(i,:)=-((1i*l)./(2*k_1(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_1*I_1*(1+1i*eta(i))*mpu_1))).*N_1_0l(i,:)./D_1(i,:); % Theoretical form
    % Damped beam
    cb_2(i,:) = sqrt(w)*(B_2/mpu_2)^(1/4); % [m/s], Phase speed (bending) in damped beam
    k_2(i,:) = w./cb_2(i,:)*(1-1i*eta(i)/4); % [1/m], complex wavenumber
    N_2_ll(i,:) = cos(k_2(i,:)*l).*cosh(k_2(i,:)*l)+1;
    N_2_0l(i,:) = cos(k_2(i,:)*l)+cosh(k_2(i,:)*l);
    D_2(i,:) = cos(k_2(i,:)*l).*sinh(k_2(i,:)*l)+sin(k_2(i,:)*l).*cosh(k_2(i,:)*l);
    Y_ll_2(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_2(i,:)*l*sqrt(B_2*mpu_2))).*N_2_ll(i,:)./D_2(i,:);
    Y_ll_2_the(i,:)=-((1i*l)./(2*k_2(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_2*I_2*(1+1i*eta(i))*mpu_2))).*N_2_ll(i,:)./D_2(i,:); % Theoretical form
    Y_0l_2(i,:)=-((0.25*eta(i)+1i)*l./(2*w./cb_2(i,:)*l*sqrt(B_2*mpu_2))).*N_2_0l(i,:)./D_2(i,:);
    Y_0l_2_the(i,:)=-((1i*l)./(2*k_2(i,:)*(1-1i*eta(i)/4)*l*sqrt(E_2*I_2*(1+1i*eta(i))*mpu_2))).*N_2_0l(i,:)./D_2(i,:); % Theoretical form
end


k_1=k_1(1,:);
k_2=k_2(1,:);

colors = ['#ff8809'; '#ff6cc7';'#00a6ff'];
%% Part 2 - Plot txt files
Alltxt=dir('*.txt');
for i=1:length(Alltxt)
    [band(:,i),f_p2(:,i),value(:,i)]=read_pulse_2021(Alltxt(i).name);
%     plot_everything(Alltxt(i),f_p2(:,i),value(:,i)); % Uncomment if needed
end

%% Plots
for i=1:length(eta)
    figure(1)
    plot(abs(k_1*l),mag2db(abs(Y_ll_1(i,:))),Color=colors(i,:))
    hold on
    title('|Y_{ll}|')
    ylabel('Magnitude [dB]')
    xlabel('Helmholtz Number kl')
    if i==length(eta);legend(['Undamped Beam: \eta = ' num2str(eta(1),'%.2f\n')],['Undamped Beam: \eta = ' num2str(eta(2),'%.2f\n')],['Undamped Beam: \eta = ' num2str(eta(3),'%.2f\n')],'AutoUpdate','off',Location='best');end
    grid on
    figure(2)
    plot(abs(k_2*l),mag2db(abs(Y_ll_2(i,:))),Color=colors(i,:),LineStyle="--")
    hold on
    title('|Y_{ll}|')
    ylabel('Magnitude [dB]')
    xlabel('Helmholtz Number kl')
    if i==length(eta);legend(['Damped Beam: \eta = ' num2str(eta(1),'%.2f\n')],['Damped Beam: \eta = ' num2str(eta(2),'%.2f\n')],['Damped Beam: \eta = ' num2str(eta(3),'%.2f\n')],'AutoUpdate','off',Location='best');end
    grid on

    figure(3)
    plot(abs(k_1*l),angle(Y_ll_1(i,:)),Color=colors(i,:))
    hold on
    title('\angle Y_{ll}')
    ylabel('Phase [rad]')
    xlabel('Helmholtz Number kl')
    if i==length(eta);legend(['Undamped Beam: \eta = ' num2str(eta(1),'%.2f\n')],['Undamped Beam: \eta = ' num2str(eta(2),'%.2f\n')],['Undamped Beam: \eta = ' num2str(eta(3),'%.2f\n')],'AutoUpdate','off',Location='best');end
    grid on
    figure(4)
    plot(abs(k_2*l),angle(Y_ll_2(i,:)),Color=colors(i,:),LineStyle="--")
    hold on
    title('\angle Y_{ll}')
    ylabel('Phase [rad]')
    xlabel('Helmholtz Number kl')
    if i==length(eta);legend(['Damped Beam: \eta = ' num2str(eta(1),'%.2f\n')],['Damped Beam: \eta = ' num2str(eta(2),'%.2f\n')],['Damped Beam: \eta = ' num2str(eta(3),'%.2f\n')],'AutoUpdate','off',Location='best');end
    grid on

    figure(5)
    plot(abs(k_1*l),mag2db(abs(Y_ll_1(i,:))),Color=[0 0 0 1*i^-3])
    hold on
    plot(abs(k_1*l),mag2db(abs(Y_0l_1(i,:))),Color=[1 0 0 1*i^-3])
    plot(abs(k_2*l),mag2db(abs(Y_ll_2(i,:))),Color=[0 0 0 1*i^-3],LineStyle="--")
    plot(abs(k_2*l),mag2db(abs(Y_0l_2(i,:))),Color=[1 0 0 1*i^-3],LineStyle="--")
    title('Input Mobility |Y_{ll}| vs Transfer Mobility |Y_{0l}|')
    ylabel('Magnitude [dB]')
    xlabel('Helmholtz Number kl')
    if i==1;legend('|Y_{ll}| Undamped Beam','','','|Y_{0l}| Undamped Beam','','','|Y_{ll}| Damped Beam','','','|Y_{0l}| Damped Beam','','','AutoUpdate','off');end
    grid on
    if i==length(eta)
        % 12 and 25 idx (undamped and damped measurements)
        figure(1)
        plot(abs(k_1*l),mag2db(abs(value(:,12))),Color='k')
        title('|Y_{ll}|')
        ylabel('Magnitude [dB]')
        xlabel('Helmholtz Number kl')
        if i==length(eta);legend(['Analytical data \eta = ' num2str(eta(1),'%.2f\n')],['Analytical data \eta = ' num2str(eta(2),'%.2f\n')],['Analytical data \eta = ' num2str(eta(3),'%.2f\n')],'Measured data \eta = 0.0007','AutoUpdate','off',Location='best');end
        grid on
        figure(2)
        plot(abs(k_2*l),mag2db(abs(value(:,25))),Color='k',LineStyle="--")
        title('|Y_{ll}|')
        ylabel('Magnitude [dB]')
        xlabel('Helmholtz Number kl')
        if i==length(eta);legend(['Analytical data \eta = ' num2str(eta(1),'%.2f\n')],['Analytical data \eta = ' num2str(eta(2),'%.2f\n')],['Analytical data \eta = ' num2str(eta(3),'%.2f\n')],'Measured data \eta = 0.0021','AutoUpdate','off',Location='best');end
        grid on
        
        figure(3)
        plot(abs(k_1*l),angle(value(:,12)),Color='k')
        title('\angle Y_{ll}')
        ylabel('Phase [rad]')
        xlabel('Helmholtz Number kl')
        if i==length(eta);legend(['Analytical data \eta = ' num2str(eta(1),'%.2f\n')],['Analytical data \eta = ' num2str(eta(2),'%.2f\n')],['Analytical data \eta = ' num2str(eta(3),'%.2f\n')],'Measured data \eta = 0.0007','AutoUpdate','off',Location='best');end
        grid on
        figure(4)
        plot(abs(k_2*l),angle(value(:,25)),Color='k',LineStyle="--")
        title('\angle Y_{ll}')
        ylabel('Phase [rad]')
        xlabel('Helmholtz Number kl')
        if i==length(eta);legend(['Analytical data \eta = ' num2str(eta(1),'%.2f\n')],['Analytical data \eta = ' num2str(eta(2),'%.2f\n')],['Analytical data \eta = ' num2str(eta(3),'%.2f\n')],'Measured data \eta = 0.0021','AutoUpdate','off',Location='best');end
        grid on
    end
end

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
