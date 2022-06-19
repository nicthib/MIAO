function gcamp=GcampMcCorrection_MIAO(blue_avg,chbr_avg,chbo_avg,startrange,dpf1,dpf2)
addpath('/local_mount/space/juno/1/ConvertFiles/conversion_loadthesefiles');
load dpffsMW_400to700nm.mat
load Hb_spectra.mat
addpath('/local_mount/space/juno/1/ConvertFiles/LED_spectra')
%load 0930_spectra_mohammed
    %load matlab.mat
   % spectra_blue=mo_blue;
   load('1003_spectra.mat')
    load('LED_spectra_0711','spectra_blue','spectra_green','spectra_red');
    spectra_blue(:,2)=spectra_blue(:,2)./max(spectra_blue(:,2));
    spectra_green=GFP(:,[1,3]);
    spectra_green(:,2)=spectra_green(:,2)./max(spectra_green(:,2));
    spectra_red(:,2)=spectra_red(:,2)./max(spectra_red(:,2));
    %load('0930_spectra_mohammed');
    %load('1003_spectra.mat')
    %load('/local_mount/space/enterprise/4/Personal Folders/Teresa/m345_spectrum.mat')
    %spectra_red=RED;
    %load('spectra_mos','spectra_blue')
    %load('spectra_mos','spectra_blue','spectra_green','spectra_red');
    %load matlab.mat
    %spectra_blue=mo_blue_stroke;
    %spectra_green=mo_green_stroke;
    %spectra_red=mo_red_stroke;
% measured LED spectra: 
x530 = spectra_green(:,1);
y530 = spectra_green(:,2)/2;
x470 = spectra_blue(:,1);
y470 = spectra_blue(:,2)/2;
x630 = spectra_red(:,1);
y630 = spectra_red(:,2)/2;
Lolim = 400;
Hilim = 700;    

splineyHb = spline([Lolim:2:Hilim],Hb(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
splineyHbO = spline([Lolim:2:Hilim],Hb02(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
spliney470 = spline(x470,y470,[Lolim:0.5:Hilim]);
spliney530 = spline(x530,y530,[Lolim:0.5:Hilim]);
spliney630 = spline(x630,y630,[Lolim:0.5:Hilim]);

splineydpff_488 = spline([Lolim:2:Hilim],dpff_488(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
splineydpff_530 = spline([Lolim:2:Hilim],dpff_530(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
splineydpff_630 = spline([Lolim:2:Hilim],dpff_630(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);

% ~ area under curves for E ( do also for dpff)
% order is blue, green, red 
EHb(1) = sum((1/sum(spliney470))*spliney470.*splineyHb);
EHb(2) = sum((1/sum(spliney530))*spliney530.*splineyHb);
EHb(3) = sum((1/sum(spliney630))*spliney630.*splineyHb);
EHbO(1) = sum((1/sum(spliney470))*spliney470.*splineyHbO);
EHbO(2) = sum((1/sum(spliney530))*spliney530.*splineyHbO);
EHbO(3) = sum((1/sum(spliney630))*spliney630.*splineyHbO);

%still need to incorporate this into the DPF
%DPF(1) = sum((1/sum(spliney470))*spliney470.*splineydpff_488);
DPF(1) = sum((1/sum(spliney470))*spliney470.*splineydpff_488);
DPF(2) = sum((1/sum(spliney530))*spliney530.*splineydpff_530);
DPF(3) = sum((1/sum(spliney630))*spliney630.*splineydpff_630);
DPF(1) = dpf1;
DPF(2) = dpf2;
mua_b(:,:,:)= EHb(1)*chbr_avg+EHbO(1)*chbo_avg;
mua_g(:,:,:) = EHb(2)*chbr_avg+EHbO(2)*chbo_avg;
blue_adjusted = exp(-DPF(1)*mua_b(:,:,:));
green_adjusted = exp(-DPF(2)*mua_g(:,:,:));
clear mua_b mua_g
blue_avg = double(blue_avg);
%gcamp = (blue_avg./repmat(mean(blue_avg(:,:, startrange),3),[1,1,size(blue_avg,3)]))./(blue_adjusted.^.21)./((green_adjusted./repmat(mean(green_adjusted(:,:, startrange),3),[1,1,size(green_adjusted,3)])).^.21); 
gcamp = (blue_avg./repmat(mean(blue_avg(:,:, startrange),3),[1,1,size(blue_avg,3)]))./(blue_adjusted.^.5)./((green_adjusted).^.5); 
gcamp = (blue_avg./repmat(mean(blue_avg(:,:, startrange),3),[1,1,size(blue_avg,3)]))./((blue_adjusted.^.5).*((green_adjusted).^.5)); 
gcamp2 = (blue_avg./((blue_adjusted.^.5).*((green_adjusted).^.5)));
gcamp2 = gcamp2./repmat(mean(gcamp2(:,:, startrange),3),[1,1,size(blue_avg,3)]); 
%gcamp = ((blue_avg-repmat(min(blue_avg(:,:, startrange),[],3),[1,1,size(blue_avg,3)]))./(repmat(mean(blue_avg(:,:, startrange),3),[1,1,size(blue_avg,3)])-repmat(min(blue_avg(:,:, startrange),[],3),[1,1,size(blue_avg,3)])))./(blue_adjusted.^.5)./((green_adjusted).^.5); 
gcamp = gcamp-1;
%gcamp_avg=gcamp;
%edit('WorkFlow.m')
clear blue_adjusted green_adjusted