function [stimon, stimoff] = findstims(stimdata,t,Fs)

stimdata(stimdata < 1) = 0;
stimdata(stimdata > 1) = 1;
convsig = ones(1,Fs);

stimcon = conv(stimdata,convsig,'same');
stimcon = [zeros(1,Fs/2) stimcon];

idx = find(stimcon == 0);
idr = diff(idx');