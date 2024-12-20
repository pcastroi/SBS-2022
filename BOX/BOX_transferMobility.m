close all; clc; clear all;
addpath(genpath('data'))
S='data/Plate average/';%*.txt";

addpath(S) %Get route
S=dir(S+"*.txt"); %Create folder

%load Frequency response
for d=1 %Foldable code
 SForce = S(contains({S.name},{'Force'})); %Alternative formulation
  SVelo = S(contains({S.name},'Velo'));

for i=1:numel(SForce)
    [idx(:,i), freq(:,i), SF(:,i)] = read_pulse_2021(SForce(i).name);
    SF(:,i) = SF(:,i)/max(abs(SF(:,i)));
    [fb3(:,i), SFb3(:,i)] = onethirdoctave_average(freq(:,i), SF(:,i));
    [idx(:,i), freq(:,i), SV(:,i)] = read_pulse_2021(SVelo(i).name);
    SV(:,i) = SV(:,i)/max(abs(SV(:,i)));
    [fb3(:,i), SVb3(:,i)] = onethirdoctave_average(freq(:,i), SV(:,i));
end
i = numel(SForce)+1;
[idx(:,i), freq(:,i), SF(:,i)] = deal(idx(:,i-1), freq(:,i-1), zeros(numel(SF(:,i-1),1)));
%SF(:,i) = SF(:,i)/max(abs(SF(:,i))); %Self evident ehh, jejeje
[fb3(:,i), SFb3(:,i)] = onethirdoctave_average(freq(:,i), SF(:,i));
[idx(:,i), freq(:,i), SV(:,i)] = deal(idx(:,i-1), freq(:,i-1), zeros(1,numel(SV(:,i-1))));
%SV(:,i) = SV(:,i)/max(abs(SV(:,i)));
[fb3(:,i), SVb3(:,i)] = onethirdoctave_average(freq(:,i), SV(:,i));
end

%% Calculate transfer mobility
for d=1
rho = 7.8E3;
h = [3 3 3 3 1.5 3];
S = [0.109 0.109 0.176 0.176 0.246 0.246];
M = rho.*h.*S;

%Squared magnitudes
DCskip = 10;
%V2 = sum(SV(DCskip:end,:),1); V2b3 = sum(SVb3,1); %Velocity 
%F2 = sum(SF(DCskip:end,:),1); F2b3 = sum(SFb3,1); %Force
F20 = mean(SF,2); %F20std = std(F2,0,2);
F20b3 = mean(SFb3,2); %F2b3std = std(F2b3,0,2);

%Transfer mobility
%Sum across 2nd dimension
Yt2 = M.*SV*[1 1 1 1 0 0]' ./ (F20*M*[1 1 1 1 0 0]');
Yt2b3 = M.*SVb3*[1 1 1 1 0 0]' ./ (F20b3*M*[1 1 1 1 0 0]');

%Loss factor
Ps = (SV+SF)*[1 1 1 1 0 0]'; Psb3 = (SVb3+SFb3)*[1 1 1 1 0 0]';
w = freq(:,1)/(2*pi);
wb3 = fb3(:,1)/(2*pi);
eta = Ps./(w*M.*SV*[1 1 1 1 1 0]');
etab3 = Psb3 ./ (wb3*M.*SVb3*[1 1 1 1 1 0]');
end
%% plot
%
for d=1
h = figure
plot(freq(:,1),mag2db(abs(Yt2)),'.','color','b','Linewidth',2)
%plot(freq(:,1),mag2db(abs(eta)),'.','color','b','Linewidth',2)
hold on
plot(fb3(:,1),mag2db(abs(Yt2b3)),'color','c','Linewidth',4)
%plot(fb3(:,1),mag2db(abs(etab3)),'color','c','Linewidth',4)
xlabel( 'Frequency (Hz)','FontSize',14);
ylabel('Magnitude','FontSize',14);
title("Transfer Mobility",'FontSize',1,'interpreter','none')
xlim([10 inf]); 
%ylim([0.8 inf])
end
%}
