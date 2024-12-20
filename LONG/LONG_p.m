%% LONG Exercise - Structure-Borne Sound
close all; clc; clear all;
% Parameters
N = 10000;
f = linspace(0,12500,N); % [Hz]
t = 1./f; % [s]
w = 2*pi*f; % [rad/s]
d = 0.03; % [m], diameter of rod
L = 0.6; % [m], lenght of rod
S = pi*(d/2)^2; % [m], cross-section surface area
rho = 8490; % [kg/m^3], m/V, density
% eta = [0.1, 0.01, 0.001]; % Damping Loss Factor
E = 10.3*10^10; % Young's modulus

delta_f=[-2907.25+2908.75,-5808.38+5809.38,-8704.75+8706.5]; % 3-dB Bandwidth
f_n=[2908, 5808.88, 8705.75]; % Different center frequencies in exp
eta=sort(delta_f./f_n,'descend'); % Damping Loss Factor

%% Part 1 - Transfer and input mobilities (use estimated E and eta)
% Preallocate
Y_00=zeros(length(eta),N);
Y_0L=Y_00;
for i=1:length(eta)
    E_c(i) = E*(1+1i*eta(i)); % [N/m^2], Young's modulus (elastic)
    cl2(i) = sqrt(E_c(i)/rho); % [m/s], speed of sound in medium
    k(i,:) = w./cl2(i); % [1/m], wavenumber
    % Task 1 - Theoretical Input/Transfer Mobility
    Y_00(i,:)=-1i*cos(k(i,:)*L)./(S*sqrt(E_c(i)*rho)*sin(k(i,:)*L));
    Y_0L(i,:)=-1i./(S*sqrt(E_c(i)*rho)*sin(k(i,:)*L));
end
colors = ['#ff8809'; '#ff6cc7';'#00a6ff'];
figure
k=abs(k);
for i=1:length(eta)
    figure(1)
    plot(f,mag2db(abs(Y_00(i,:))),Color=colors(i,:))
    xlim([min(f) max(f)])
    title('|Y_{00}|')
    ylabel('Magnitude [dB]')
    xlabel('Frequency [Hz]')
    if i==3;legend(['\eta = ' num2str(eta(1),'%.2e\n')],['\eta = ' num2str(eta(2),'%.2e\n')],['\eta = ' num2str(eta(3),'%.2e\n')],Location='best');end
    grid on
    hold on
    figure(2)
    plot(f,angle(Y_00(i,:)),Color=colors(i,:))
    xlim([min(f) max(f)])
    title('\angle Y_{00}')
    ylabel('Phase [rad]')
    xlabel('Frequency [Hz]')
    if i==3;legend(['\eta = ' num2str(eta(1),'%.2e\n')],['\eta = ' num2str(eta(2),'%.2e\n')],['\eta = ' num2str(eta(3),'%.2e\n')],Location='best');end
    grid on
    hold on
    figure(3)
    plot(f,mag2db(abs(Y_0L(i,:))),Color=colors(i,:))
    xlim([min(f) max(f)])
    title('|Y_{0L}|')
    ylabel('Magnitude [dB]')
    xlabel('Frequency [Hz]')
    if i==3;legend(['\eta = ' num2str(eta(1),'%.2e\n')],['\eta = ' num2str(eta(2),'%.2e\n')],['\eta = ' num2str(eta(3),'%.2e\n')],Location='best');end
    grid on
    hold on
    figure(4)
    plot(f,angle(Y_0L(i,:)),Color=colors(i,:))
    xlim([min(f) max(f)])
    title('\angle Y_{0L}')
    ylabel('Phase [rad]')
    xlabel('Frequency [Hz]')
    if i==3;legend(['\eta = ' num2str(eta(1),'%.2e\n')],['\eta = ' num2str(eta(2),'%.2e\n')],['\eta = ' num2str(eta(3),'%.2e\n')],Location='best');end
    grid on
    hold on
end

%% Part 2 - Experiment
Alltxt=dir('*.txt');
Cohtxt=Alltxt(endsWith({Alltxt.name},{'Coherence.txt'}));
Freqtxt=Alltxt(~endsWith({Alltxt.name},{'Coherence.txt'}));
for i=1:length(Cohtxt)
    [c_band(:,i),c_f(:,i),c_value(:,i)]=read_pulse_2021(Cohtxt(i).name);
    [fq_band(:,i),fq_f(:,i),fq_value(:,i)]=read_pulse_2021(Freqtxt(i).name);
    figure
    plot(fq_f(:,i),mag2db(abs(fq_value(:,i)))) 
    if contains(Freqtxt(i).name,'Input') == 1
        title('|Y_{00}|')
    else
        title('|Y_{0L}|')
    end
    ylabel('Magnitude [dB]')
    xlabel('Frequency [Hz]')
    xlim([min(fq_f(:,i)) max(fq_f(:,i))])
    grid on
    figure
    plot(fq_f(:,i),angle(fq_value(:,i)))
    if contains(Freqtxt(i).name,'Input') == 1
        title('\angle Y_{00}')
    else
        title('\angle Y_{0L}')
    end
    ylabel('Phase [rad]')
    xlabel('Frequency [Hz]')
    xlim([min(fq_f(:,i)) max(fq_f(:,i))])
    grid on
    figure
    plot(c_f(:,i),c_value(:,i))
    title('C_{xy}')
    ylabel('Coherence')
    xlabel('Frequency [Hz]')
    xlim([min(c_f(:,i)) max(c_f(:,i))])
    grid on
end




