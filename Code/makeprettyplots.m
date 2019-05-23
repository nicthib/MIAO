function [roi,roinames,tc] = makeprettyplots(h,data,fname)
if ~exist([h.outfolder.String,'/jpg_summaries']);
    try
        mkdir(fullfile(h.outfolder.String,'jpg_summaries'))
    catch
        h.status.String = sprintf('\n\nPermission denied. Cannot write jpg files to analysis folder.');
        drawnow
    end
end
filename = [h.outfolder.String,'/jpg_summaries/', ...
    h.m.stimlist{h.m.runselect},fname,'_summary.jpg'];
tname = strrep([h.m.mouse,' ',h.m.stimlist{h.m.runselect}, ' ', fname],'_',' ');
if ~h.chooseroi.Value
    roi = [30 54  67 90; 74 92 71 89; 13 36 67 90; 
           96 112 79 95; 37 53 35 51; 76 91 33 48];
else
    f = figure(1);
    imagesc(nanmean(data.blue,3));
    axis image; colormap gray; hold on
    for i = 1:6
       title(['Choose roi # ' mat2str(i)]);
       [tmpR,tmpC] = ginput(2);
       roi(i,:) = [tmpC' tmpR']; roi = round(roi);
       plot([roi(i,3) roi(i,3) roi(i,4),roi(i,4),roi(i,3)],...
         [roi(i,1) roi(i,2) roi(i,2),roi(i,1),roi(i,1)])
    end
    h.chooseroi.Value = 0;
    close(f)
end
roinames = {'top HP','bottom HP','top Whisker','bottom Whisker','top Motor', 'bottom Motor'};

numr = size(roi,1);
yy = fieldnames(data);
fr = h.m.framerate/h.m.nLEDs;
TT = [0:size(data.chbo,3)-1]/fr;
cols = 0.8*jet(numr);

f = figure(10); clf
set(gcf,'Position',[60 5 830 868])

if isfield(data,'blue');
    subplot(5,3,1);
    imagesc(data.blue(:,:,30));
    colormap gray;
    axis image;
else
    subplot(5,3,1);
    imagesc(data.chbo(:,:,30)+data.chbr(:,:,30))
    colormap gray;
    axis image;
end
title(tname);
set(gca,'FontSize',7)
hold on
data.chbo = data.chbo*1e6;
data.chbr = data.chbr*1e6;
data.gcamp = data.gcamp*200;
for i = 1:numr
    subplot(5,3,1)
    hold on
    plot([roi(i,3) roi(i,3) roi(i,4),roi(i,4),roi(i,3)],...
         [roi(i,1) roi(i,2) roi(i,2),roi(i,1),roi(i,1)],'color',cols(i,:))
    for k = 1:3;
        tc.(yy{k})(i,:) = squeeze(mean(mean(data.(yy{k})...
            (roi(i,1):roi(i,2),roi(i,3):roi(i,4),:),2),1));
        subplot(numr,3,(i-1)*3+2)
        plot(TT,tc.(yy{k})(i,:)-tc.(yy{k})(i,50),yy{k}(1)); hold on
    end
    axis tight
    a = axis;
    plot([h.m.tpre h.m.tpre],[a(3) a(4)],'k');
    plot(h.m.tstim+[h.m.tpre h.m.tpre],[a(3) a(4)],'k');
    set(gca,'FontSize',7)
    for k = 4:length(yy)
        tc.(yy{k})(i,:) = squeeze(mean(mean(data.(yy{k})...
            (roi(i,1):roi(i,2),roi(i,3):roi(i,4),:),2),1));
    end
    subplot(numr,3,(i-1)*3+3)
    plot(TT,tc.gcamp(i,:),'k'); hold on
    plot(TT,tc.chbo(i,:)+tc.chbr(i,:),'g');
    plot(TT,tc.chbr(i,:),'b');
    plot(TT,tc.chbo(i,:),'r');
    axis tight
    a = axis;
    plot([h.m.tpre h.m.tpre],[a(3) a(4)],'m');
    plot(h.m.tstim+[h.m.tpre h.m.tpre],[a(3) a(4)],'m');
    ylabel('Hb(uM), GCaMPx200')
    title(roinames{i},'color',cols(i,:));
    set(gca,'FontSize',7)
end

xlabel('time (s)')
for k = 4:length(yy)
    subplot(5,3,(k-4)*3+7);
    imagesc(real(mean(data.(yy{k})(:,:,round(h.m.tpre*fr+(0:h.m.tstim*fr))),3)...
        -mean(data.(yy{k})(:,:,round(1:20-21+(h.m.tpre*fr))),3)));
    axis image; colormap gray; set(gca,'FontSize',7)
    %colorbar
    ylabel(yy{k})
end
subplot(5,3,4); % Plotting chbt
tmphbo = (mean(data.chbo(:,:,round(h.m.tpre*fr+(0:h.m.tstim*fr))),3)...
    -mean(data.chbo(:,:,round(1:20-21+h.m.tpre*fr)),3));
tmphbr = (mean(data.chbr(:,:,round(h.m.tpre*fr+(0:h.m.tstim*fr))),3)...
    -mean(data.chbr(:,:,round(1:20-21+h.m.tpre*fr)),3));
imagesc(real(tmphbo+tmphbr));
axis image; colormap gray; %colorbar
ylabel('HbT (-20 to 20)'); set(gca,'FontSize',7)
try saveas(f,filename,'bmp'); catch; end
