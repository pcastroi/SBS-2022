main='/'; %Data storage
data='/data/';
R='/data/Radiation/';%*.txt";

addpath(genpath('data')) %Get route
R=dir(R+"*.txt"); %Create folder

%load Frequency response
for d=1 %Foldable code
   front = R(contains({R.name},'front'));
 SFfront = front(contains({front.name},'Force'));
 SVfront = front(contains({front.name},'Velocity'));
  Ifront = front(contains({front.name},'Intensity'));
    back = R(contains({R.name},'back'));
  SFback = back(contains({back.name},'Force'));
  SVback = back(contains({back.name},'Velocity'));
   Iback = back(contains({back.name},'Intensity'));

for i=1:numel(SFfront) %Front force and velocity 
    [idx(:,i), freq(:,i), SFf(:,i)] = read_pulse_2021(SFfront(i).name);
    SFf(:,i) = SFf(:,i)/max(abs(SFf(:,i)));
    [fb3(:,i), SFfb3(:,i)] = onethirdoctave_average(freq(:,i), SFf(:,i));
    [idx(:,i), freq(:,i), SVf(:,i)] = read_pulse_2021(SVfront(i).name);
    SVf(:,i) = SVf(:,i)/max(abs(SVf(:,i)));
    [fb3(:,i), SVfb3(:,i)] = onethirdoctave_average(freq(:,i), SVf(:,i));
end
for i=1:numel(SFfront) %Back force and velocity 
    [idx(:,i), freq(:,i), SFb(:,i)] = read_pulse_2021(SFback(i).name);
    SFb(:,i) = SFb(:,i)/max(abs(SFb(:,i)));
    [fb3(:,i), SFbb3(:,i)] = onethirdoctave_average(freq(:,i), SFb(:,i));
    [idx(:,i), freq(:,i), SVb(:,i)] = read_pulse_2021(SVback(i).name);
    SVb(:,i) = SVb(:,i)/max(abs(SVb(:,i)));
    [fb3(:,i), SVbb3(:,i)] = onethirdoctave_average(freq(:,i), SVb(:,i));
end
for i=1:numel(Ifront) %Front intensity
    [idx(:,i), freq(:,i), If(:,i)] = read_pulse_2021(Ifront(i).name);
    If(:,i) = If(:,i)/max(abs(If(:,i)));
    [fb3(:,i), Ifb3(:,i)] = onethirdoctave_average(freq(:,i), If(:,i));
end
for i=1:numel(Iback) %Back intensity
    [idx(:,i), freq(:,i), Ib(:,i)] = read_pulse_2021(Iback(i).name);
    Ib(:,i) = Ib(:,i)/max(abs(Ib(:,i)));
    [fb3(:,i), Ibb3(:,i)] = onethirdoctave_average(freq(:,i), Ib(:,i));
end
end

%% Calculate transfer mobility
for d=1
S = 0.246; %m2
rhoc0 = 410; %kg/(m2*s)

%Pf = Ifront*S; Pb = Iback*S;
radf = If./(rhoc0*SVf);
radfb3 = Ifb3./(rhoc0*SVfb3);

end
%% plot
%
for d=1
h = figure
plot(freq(:,1),mag2db(abs(radf)),'.','color','b','Linewidth',2)
hold on
plot(fb3(:,1),mag2db(abs(radfb3)),'color','c','Linewidth',4)
xlabel( 'Frequency (Hz)','FontSize',14);
ylabel('Magnitude','FontSize',14);
title("Transfer Mobility",'FontSize',1,'interpreter','none')
xlim([10 inf]); 
%ylim([0.8 inf])
end
%}
