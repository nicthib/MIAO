function data = convertNic(data,h)
if h.m.isgui
    h.status.String = sprintf('\n\nConverting data...'); drawnow;
end
load dpffsMW_400to700nm.mat
load Hb_spectra.mat
if h.m.greenfilter == 530
    load LED_spectra_new
elseif h.m.greenfilter == 534
    load LED_spectra_0711
end
x530 = spectra_green(:,1); y530 = spectra_green(:,2);
x630 = spectra_red(:,1);   y630 = spectra_red(:,2);

Lolim = 400; Hilim = 700;

splineyHb  = spline([Lolim:2:Hilim],Hb(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
splineyHbO = spline([Lolim:2:Hilim],Hb02(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
spliney530 = spline(x530,y530,[Lolim:0.5:Hilim]);
spliney630 = spline(x630,y630,[Lolim:0.5:Hilim]);

splineydpff_530 = spline([Lolim:2:Hilim],dpff_530(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
splineydpff_630 = spline([Lolim:2:Hilim],dpff_630(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);

EHb_green = sum((1/sum(spliney530))*spliney530.*splineyHb);
EHb_red = sum((1/sum(spliney630))*spliney630.*splineyHb);
EHbO_green = sum((1/sum(spliney530))*spliney530.*splineyHbO);
EHbO_red = sum((1/sum(spliney630))*spliney630.*splineyHbO);

DPF_g = sum((1/sum(spliney530))*spliney530.*splineydpff_530);
DPF_red = sum((1/sum(spliney630))*spliney630.*splineydpff_630);

clear mua chbo chb chbt

mua_g = -(1/DPF_g)*log(double(data.green)./repmat(mean(data.green(:,:,round(h.m.baseline)),3),[1,1,size(data.green,3)]));
mua_r = -(1/DPF_red)*log(double(data.red)./repmat(mean(data.red(:,:,round(h.m.baseline)),3),[1,1,size(data.red,3)]));
clear DPF1 DPF2
data.chbo = -(squeeze((EHb_green*mua_r-EHb_red*mua_g)/(EHbO_green*EHb_red-EHbO_red*EHb_green)));
data.chbr = (squeeze((EHbO_green*mua_r-EHbO_red*mua_g)/(EHb_red*EHbO_green-EHb_green*EHbO_red)));

h.status.String = sprintf('\n\nCorrecting gcamp...'); drawnow;
load('1003_spectra.mat'); % GFP spectrum
spectra_blue(:,2)=spectra_blue(:,2)./max(spectra_blue(:,2));
spectra_gfp = GFP(:,[1,3]);
spectra_gfp(:,2) = spectra_gfp(:,2)./max(spectra_gfp(:,2));
spectra_red(:,2) = spectra_red(:,2)./max(spectra_red(:,2));

x530g = spectra_gfp(:,1);
y530g = spectra_gfp(:,2)/2; %why?
x470 = spectra_blue(:,1);
y470 = spectra_blue(:,2)/2;

splineyHb = spline([Lolim:2:Hilim],Hb(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
splineyHbO = spline([Lolim:2:Hilim],Hb02(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
spliney470 = spline(x470,y470,[Lolim:0.5:Hilim]);
spliney530g = spline(x530g,y530g,[Lolim:0.5:Hilim]);

splineydpff_488 = spline([Lolim:2:Hilim],dpff_488(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
splineydpff_530 = spline([Lolim:2:Hilim],dpff_530(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);

EHb_blue = sum((1/sum(spliney470))*spliney470.*splineyHb);
EHb_gfp = sum((1/sum(spliney530g))*spliney530g.*splineyHb);
EHbO_blue = sum((1/sum(spliney470))*spliney470.*splineyHbO);
EHbO_gfp = sum((1/sum(spliney530g))*spliney530g.*splineyHbO);

DPF_blue = sum((1/sum(spliney470))*spliney470.*splineydpff_488);
DPF_gfp = sum((1/sum(spliney530g))*spliney530g.*splineydpff_530);
DPF_blue = h.m.dpf(1);
DPF_gfp = h.m.dpf(2);

mua_b(:,:,:)= EHb_blue*data.chbr+EHbO_blue*data.chbo; % estimate contamination from hb on gfp ex and em
mua_g(:,:,:) = EHb_gfp*data.chbr+EHbO_gfp*data.chbo;

blue_adjusted = exp(-DPF_blue*mua_b(:,:,:));% estimate contamination from hb on gfp ex and em
green_adjusted = exp(-DPF_g*mua_g(:,:,:));
drawnow

clear mua_b mua_g
data.blue = double(data.blue);
data.gcamp = -1+(data.blue./repmat(mean(data.blue(:,:, h.m.baseline),3),[1,1,size(data.blue,3)]))./(blue_adjusted.^.5)./((green_adjusted).^.5);
clear blue_adjusted green_adjusted
