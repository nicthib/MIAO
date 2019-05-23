function m = makem
m.outputs = 'odn';
m.loadpct = 1;
m.dsf = 1; m.dsfWebcam = 1;
m.dpf = [.7 .8];
m.Dr = .8; m.Dg = .4;
m.baseline = [30:100];
m.nrot = 0;
m.noload = 0;
m.corr_flicker = [];
m.gettruestim = 0;
m.PCAcomps = 0;
m.bkgsub = 0;
m.loadrandstim = 0;
m.offsetfactor = 2;
m.smooth = [];
m.autocrop = 0;