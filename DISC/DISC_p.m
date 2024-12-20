%% DISC Exercise - Structure-Borne Sound
clc; clear all; close all;
addpath data

% Read data
DataFolder = 'data\';
H1Files=dir([DataFolder,'*H1*.txt']);
CohFiles=dir([DataFolder,'*Coherence*.txt']);
VelLasFiles=dir([DataFolder,'*m(Velocity (laser)*.txt']);
VelAccFiles=dir([DataFolder,'*m(Velocity (acc.)*.txt']);
for i=1:length(H1Files)
    [band_H1(:,i),f_H1(:,i),val_H1(:,i)]=read_pulse_2021(H1Files(i).name);
    [band_Coh(:,i),f_Coh(:,i),val_Coh(:,i)]=read_pulse_2021(CohFiles(i).name);
    [band_VelLas(:,i),f_VelLas(:,i),val_VelLas(:,i)]=read_pulse_2021(VelLasFiles(i).name);
    [band_VelAcc(:,i),f_VelAcc(:,i),val_VelAcc(:,i)]=read_pulse_2021(VelAccFiles(i).name);
end

%%

SNRLas = val_VelLas(:,2)./val_VelLas(:,3);
SNRAcc = val_VelAcc(:,2)./val_VelAcc(:,3);

%% Plots
close all;
CohTitles = ["Calibration";"Measurement with excitation";"Measurement w/o excitation"];
figure
i=0;
for j=1:2*size(val_Coh,2)
    subplot(3,2,j)
    if rem(j,2)==1
        i=i+1;
        plot(f_Coh(:,i),val_Coh(:,i))
        xlabel('Frequency [Hz]')
        ylabel('Coherence')
    else
        semilogx(f_VelLas(:,i),mag2db(val_VelLas(:,i)),f_VelAcc(:,i),mag2db(val_VelAcc(:,i)))
        xlabel('Frequency [Hz]')
        ylabel('Magnitude [dB]')
        legend('Laser', 'Accelerometer',Location='best')
    end
    title(CohTitles(i,:))
end

% figure
% for i=1:size(val_Coh,2)
%     subplot(3,1,i)
%     semilogx(f_VelLas(:,i),mag2db(val_VelLas(:,i)),f_VelAcc(:,i),mag2db(val_VelAcc(:,i)))
%     xlabel('Frequency [Hz]')
%     ylabel('Magnitude [dB]')
%     legend('Laser', 'Accelerometer',Location='best')
%     title(CohTitles(i,:))
% end