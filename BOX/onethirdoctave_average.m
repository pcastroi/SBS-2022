% -------------------------------------------------
% AVERAGING OF DATA IN 1/3-OCTAVE BANDS
% CODE WRITTEN BY LARS FRIIS
% -------------------------------------------------
% 
% [mfreq, mdata] = onethirdoctave_average(freq, data)
% 
% function which averages the input data in 1/3-octave bands. The first
% 1/3-octave band has a centre frequency of 16 Hz. The last 1/3-octave band
% is automatically determined from the input data.
% 
% INPUT:
% freq:     a vector with the input frequencies
% data:     a vector with the input data corresponding to freq
% 
% OUTPUT:
% mfreq:    a vector with the centre frequencies of the 1/3-octave bands
% mdata:    a vector with the data averaged in 1/3-octave bands corresponding to mfreq


function [mfreq, mdata] = onethirdoctave_average(freq, data)

%resolution of data
resol=freq(2)-freq(1);

%band cross-over frequencies
midfreq=[16 20 25 32 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000];
crossfreq(1)=midfreq(1)/2^(1/6);
crossfreq(2:length(midfreq)+1)=midfreq*2^(1/6);

%cross-over indicies
crosselem=ceil(crossfreq/resol);

%averaging of data
nn=1;
while crossfreq(nn+1)<=freq(end)
    mdata(nn)=mean(data(crosselem(nn):crosselem(nn+1)-1));
    nn=nn+1;  
end

%center frequencies
mfreq=midfreq(1:length(mdata));

