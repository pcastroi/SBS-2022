%% PRE Exercise - Structure-Borne Sound
close all; clc; clear all;
% Parameters
N = 10000;
f = linspace(0,50000,N); % [Hz]
t = 1./f; % [s]
w = 2*pi*f; % [rad/s]
m = 5.5; % [kg]
s = 450*10^3; % [N/m]
damp_ratio = [0.01 0.1 1];
w0 = sqrt(s/m);
r = damp_ratio*2*m*w0; % [Ns/m]


% Task 1
for i=1:length(r)

    % Displacement
%     xi(:,i)=exp(1i*w.*t)./(s-m*w.^2+1i*w*r(i));

    % Mobility
    Y(:,i)=(1i.*w)./(s-m*w.^2+1i*w*r(i));
    
    % Receptance
    H(:,i)=1./(s-m*w.^2+1i*w*r(i));

end

% Task 2
for i=1:length(r)

    % Displacement
%     xi_11(:,i)=-(exp(1i*w.*t).*(s-w.^2.*m+1i*r(i)))./(w.^2*m.*(s+1i*w.*r(i)));

    % Mobility
    Y_11(:,i)=-1i*(s-w.^2.*m+1i*r(i))./(w.*m.*(s+1i*w.*r(i)));
    
    % Receptance
    H_11(:,i)=-(s-w.^2.*m+1i*r(i))./(w.^2*m.*(s+1i*w.*r(i)));

end

%% Plots
figure;
subplot(2,2,1)
loglog(w,abs(Y))
xline(w0,':')
title('|Y|')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
subplot(2,2,3)
semilogx(w,angle(Y))
xline(w0,':')
title('\phi_Y')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
subplot(2,2,2)
loglog(w,abs(H))
xline(w0,':')
title('|H|')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
subplot(2,2,4)
semilogx(w,angle(H))
xline(w0,':')
title('\phi_H')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
sgtitle('Simple mass-spring-damper system')


figure;
subplot(2,2,1)
loglog(w,abs(Y_11))
xline(w0,':')
title('|Y_{11}|')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
subplot(2,2,3)
semilogx(w,angle(Y_11))
xline(w0,':')
title('\phi_{Y_{11}}')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
subplot(2,2,2)
loglog(w,abs(H_11))
xline(w0,':')
title('|H_{11}|')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
subplot(2,2,4)
semilogx(w,angle(H_11))
xline(w0,':')
title('\phi_{H_{11}}')
legend('\zeta = 0.01', '\zeta = 0.1', '\zeta = 1', 'Location','best')
set(gca, 'XTick', sort([10^5 w0]), 'XTickLabel', {'\omega_0', '10^5'})
grid on
sgtitle('Base driven free mass-spring-damper system')