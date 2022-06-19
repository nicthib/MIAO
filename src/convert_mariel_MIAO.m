function [chbo,chbr,chbt] = convert_mariel_MIAO(data1,data2,c1,c2,baseinterval,greenfilter)
%--------------------------------------------------------------------------
% convert takes in [blue,green] or [green,red] data sets and converts to
% HbO HbR and HbT 
%
% chb     : [HbR] output
% chbo    : [HbO] output (calculate HbT = HbO + HbR)
% data1   : first input data set (x,y,t)
% data2   : second input data set (x,y,t)
% c1      : color of data1 'b', 'g', or 'r'
% c2      : color of data2 'b','g', or 'r'
% filter  : 530 or 534
%
% Files to load - 
% dpffsMW_400to700nm.mat  : waves, dpff_488, dpff_530, dpff_630
% Hb_spectra.mat          : Hb, Hb02, lambda
% LED_spectra_new         : spectra_green, spectra_blue, spectra_red
%--------------------------------------------------------------------------
addpath('/local_mount/space/juno/1/ConvertFiles/conversion_loadthesefiles');
load dpffsMW_400to700nm.mat
% figure; plot(waves,dpff_488)
load Hb_spectra.mat
    addpath('/local_mount/space/juno/1/ConvertFiles/LED_spectra')
if greenfilter == 530 
    %load LED_spectra_new
elseif greenfilter == 534 
    load('LED_spectra_0711','spectra_blue','spectra_green','spectra_red');
    %load('1003_spectra.mat')
    spectra_blue(:,2) = spectra_blue(:,2)./max(spectra_blue(:,2));
    spectra_green(:,2) = spectra_green(:,2)./max(spectra_green(:,2));
    spectra_red(:,2) = spectra_red(:,2)./max(spectra_red(:,2));
    %load('0930_spectra_mohammed');    
    %load('0930_spectra_mohammed');    
    %load('/local_mount/space/enterprise/4/Personal Folders/Teresa/m345_spectrum.mat')
    %spectra_red=RED;
    %load('spectra_mos','spectra_blue')
    %load matlab.mat
    %spectra_blue=mo_blue;
    %spectra_green=mo_green_stroke;
    %spectra_red=mo_red_stroke;
end
% measured LED spectra: 
x530 = spectra_green(:,1);
y530 = spectra_green(:,2);
 x470 = spectra_blue(:,1);
 y470 = spectra_blue(:,2);
x630 = spectra_red(:,1);
y630 = spectra_red(:,2);
%{
figure; 
subplot(1,2,1); plot(x530,y530,'b-x',x470,y470,'g-x',x630,y630,'r-x');
subplot(1,2,2); plot(Hb,'b'); hold on; plot(Hb02,'r');
%}
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
DPF(1) = sum((1/sum(spliney470))*spliney470.*splineydpff_488);
DPF(2) = sum((1/sum(spliney530))*spliney530.*splineydpff_530);
DPF(3) = sum((1/sum(spliney630))*spliney630.*splineydpff_630);

% 'b' = 98, 'g'= 103, 'r' = 114 
% 1: blue, 2: green, 3: red
if strcmp(c1,'b') && strcmp(c2,'g')
    c1d = 1; c2d = 2;
elseif strcmp(c1,'g') && strcmp(c2,'b')
    c1d = 2; c2d = 1;
elseif strcmp(c1,'g') && strcmp(c2,'r')
    c1d = 2; c2d = 3;
elseif strcmp(c1,'r') && strcmp(c2,'g')
    c1d = 3; c2d = 2;
elseif strcmp(c1,'r') && strcmp(c2,'b')
    c1d = 3; c2d = 1;
elseif strcmp(c1,'b') && strcmp(c2,'r')
    c1d = 1; c2d = 3;   
end

DPF1 = DPF(c1d);
DPF2 = DPF(c2d);
EHb1 = EHb(c1d);
EHb2 = EHb(c2d);
EHbO1 = EHbO(c1d);
EHbO2 = EHbO(c2d);

% make sure DPF and data are for same LED 
clear mua chbo chb chbt

mua(1,:,:,:)=-(1/DPF1)*log(squeeze(data1(:,:,:)./repmat(mean(data1(:,:,round(baseinterval)),3),[1,1,size(data1,3)])));
mua(2,:,:,:)=-(1/DPF2)*log(squeeze(data2(:,:,:)./repmat(mean(data2(:,:,round(baseinterval)),3),[1,1,size(data2,3)])));
clear data2 data1 DPF1 DPF2
chbo = squeeze((EHb2*mua(1,:,:,:)-EHb1*mua(2,:,:,:))/(EHbO1*EHb2-EHbO2*EHb1));
%save([rat run 'chbo.mat'],'chbo','-v7.3')
%clear chbo;
chbr = squeeze((EHbO2*mua(1,:,:,:)-EHbO1*mua(2,:,:,:))/(EHb1*EHbO2-EHb2*EHbO1));
%save([rat run 'chbr.mat'],'chb','-v7.3')
chbt = chbo+chbr;
%--------------------------------------------------------------------------

% data1 = blue;
% data2 = green;

% play movie: 

% figure; 
% for i = 2:size(chbo,3)
%     subplot(1,2,1); imagesc(data1(:,:,i)); axis image; 
%     if i==2; c = caxis; else; caxis(c); end
%     subplot(1,2,2); imagesc(data2(:,:,i)); axis image; 
%     if i==2; d = caxis; else caxis(d); end
%     pause(0.01); 
% end
% 
% chbt = chb + chbo;
% fr = 30;
% start = 90; ends = 180;
% 
% figure; plot([start:size(chbo,3)]./fr,(squeeze(mean(mean(chbo(5:end-5,5:end-5,start:end),1),2))-mean(mean(mean(chbo(5:end-5,5:end-5,start:ends))))),'r'); 
% hold on; plot([start:size(chbo,3)]./fr,(squeeze(mean(mean(chb(5:end-5,5:end-5,start:end),1),2))-mean(mean(mean(chb(5:end-5,5:end-5,start:ends))))),'b');
% plot([start:size(chbo,3)]./fr,(squeeze(mean(mean(chbt(5:end-5,5:end-5,start:end),1),2))-mean(mean(mean(chbt(5:end-5,5:end-5,start:ends))))),'g');
% 

% % for plotting purposes interpolate so Hb, Hb02, 530, and 470 are on same axes between 450 and
% 610./if strcmp(c1,'b') && strcmp(c2,'g')
%     c1d = 1; c2d = 2;
% elseif strcmp(c1,'g') && strcmp(c2,'b')
%     c1d = 2; c2d = 1;
% elseif strcmp(c1,'g') && strcmp(c2,'r')
%     c1d = 2; c2d = 3;
% elseif strcmp(c1,'r') && strcmp(c2,'g')
%     c1d = 3; c2d = 2;
% elseif strcmp(c1,'r') && strcmp(c2,'b')
%     c1d = 3; c2d = 1;
% elseif strcmp(c1,'b') && strcmp(c2,'r')
%     c1d = 1; c2d = 3;   
% end
% x530c = spectra_green(find(spectra_green(:,1)==459.96):find(spectra_green(:,1)==610.05),1);
% y530c = spectra_green(find(spectra_green(:,1)==459.96):find(spectra_green(:,1)==610.05),2);
% x470c = spectra_blue(find(spectra_blue(:,1)==459.96):find(spectra_blue(:,1)==610.05),1);
% y470c = spectra_blue(find(spectra_blue(:,1)==459.96):find(spectra_blue(:,1)==610.05),2);
% x630c = spectra_red(find(spectra_red(:,1)==459.96):find(spectra_red(:,1)==610.05),1);
% y630c = spectra_red(find(spectra_red(:,1)==459.96):find(spectra_red(:,1)==610.05),2);