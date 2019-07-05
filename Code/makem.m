% This verion of makem supports LoadData v1.7.
function m = makem
m.outputs = 'odn'; % oxy deoxy neural
m.loadpct = [0 1]; % [0 1] is the whole run, while [0 .25] would be the first 1/4 of the run.
m.dsf = 1; 
m.dsfWebcam = 1;
m.dpf = [.7 .8];
m.Dr = .8; m.Dg = .4; % for jRGECO correction
m.baseline = 30:100;
m.nrot = 0; % n*90 degree rotation CW
m.noload = 0;
m.corr_flicker = [];
m.gettruestim = 0;
m.PCAcomps = 0;
m.bkgsub = 0;
m.loadrandstim = 0;
m.offsetfactor = 2;
m.smooth = [];
m.autocrop = 0;
m.greenfilter = 534;