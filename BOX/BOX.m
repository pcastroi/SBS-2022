
data='data/';
InpY='data/InputMobility/';%*.txt";
addpath(InpY) %Get route
InpY=dir(InpY+"*.txt"); %Create folder

%load Frequency response
autoSpec = InpY(startsWith({InpY.name},{'A'})); %Alternative formulation
Coherence = InpY(startsWith({InpY.name},'C'));
H1 = InpY(startsWith({InpY.name},'H'));

i = 2; %Frequency response index
disp(InpY(i).name(1:end-4));

[idx, freq, S] = read_pulse_2021(autoSpec(i).name);
S = S/max(abs(S));
[fb3, Sb3] = onethirdoctave_average(freq, S);
[idx, freq, C] = read_pulse_2021(Coherence(i).name);
C = C/max(abs(S));
[fb3, Cb3] = onethirdoctave_average(freq, C);
[idx, freq, H] = read_pulse_2021(H1(i).name);
H = H/max(abs(H));
[fb3, Hb3] = onethirdoctave_average(freq, H);
    
%% plot
%plot(freq,mag2db(abs(C)),'color','b','Linewidth',2)
plot(freq,abs(C),'color','b','Linewidth',2)
hold on
%plot(fb3,mag2db(abs(Cb3)),'color','c','Linewidth',4)
plot(freq,abs(C),'color','b','Linewidth',2)
xlabel( 'Frequency (Hz)','FontSize',14);
ylabel('Magnitude','FontSize',14);
title(InpY(i).name,'FontSize',1,'interpreter','none')
xlim([64 inf]); 
ylim([0.8 inf])