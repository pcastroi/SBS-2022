function [band,f,value]=read_pulse(fname)
%ex. [f,H]=trimBKfile('myout.txt')

fid=fopen(fname);
%read header
linecounter=0;
data = [];
myStop=[];

while isempty(myStop)
    tline = fgetl(fid);
    linecounter = linecounter+1;
    myStop =  str2num(tline(1)); % waits for a number
end
% Header is finished when first character in line is a number 

% %read data
while ~isempty(myStop)
    linedata=sscanf(tline,'%f %f %f %f',[4,1]);
    data= [data; linedata'];
    tline = fgetl(fid);
    linecounter = linecounter+1;
    if feof(fid)
        myStop=[];  % EOF
    else
        myStop =  feof(fid) & str2num(tline(1)); % waits for a number
    end
end
fclose(fid);
f=data(:,2);
realH=data(:,3);
if size(data,2) == 4
    imgH=data(:,4);
else
    imgH=zeros(length(data),1);
end
% added 0.25 gain to transfer functions (wrong measured sensitivity)
value = .25*(realH + j.*imgH);
band=data(:,1);
return;