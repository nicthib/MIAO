% m = makem creates the m struct used in LoadData. See below for a
% description of each field:
%
% m.outputs         : String with these characters: rgblodn.
%                     (Red, Green, Blue, Lime, hbO, hbR, Neural).
%                     Default is 'odn'.
%
% m.loadpct         : This loads a proportion of the frames, where 
%                     [0 1] loads the whole run. Default is [0 1].
%
% m.dsf             : Downsample factor (Integer). 
%                     Default is 1.
%
% m.dpf             : Array of two pathlengths for GCaMP correction 
%                     Default is [.7 .8]. 
%
% m.dr and m.dg     : Correction factors for jRGECO correction. 
%                     Defaults are Dr = .8 and Dg = .4.
%
% m.baseline        : Vector containing frames to use as a baseline for
%                     correction. Default is [30:100].
%
% m.nrot            : n*90 degree CW rotation of dataset.
%
%
% m.noload          : 1 if you only want metadata returned.
%                     Default is 0  
%
% m.corr_flicker    : Array of which LEDs to perform flicker correction.
%                     Default is empty []
%
% m.PCAcomps        : 1:n components to use for PCA denoising. 
%                     Default is 0
%
% m.smooth          : Array of which LEDsto perform a temporal smooth. 
%                     Default is empty []
%
% m.autocrop        : 1 if you want an automated cropmask.
%                     Default is 0.
%
% m.greenfilter     : Wavelength of green bandpass.
%                     Default is 534.
%
% This verion of makem supports LoadData v1.7+.
function m = makem
m.outputs = 'odn';
m.loadpct = [0 1]; 
m.dsf = 1; 
m.dpf = [.7 .8];
m.Dr = .8; m.Dg = .4; 
m.baseline = 30:100;
m.nrot = 0; 
m.noload = 0;
m.corr_flicker = [];
m.PCAcomps = 0;
m.bkgsub = 0;
m.offsetfactor = 2;
m.smooth = [];
m.autocrop = 0;
m.greenfilter = 534;
m.isgui = 0;