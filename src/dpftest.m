function dpftest(data,dpfr1,dpfr2,len,baseline)
addpath /local_mount/space/juno/1/Software/MIAO/MIAO_v2
imagesc(mean(data.blue,3)); axis image;
[x,y] = ginput(2); x = round(x); y = round(y);
TCo = squeeze(mean(mean(data.chbo(y(1):y(2),x(1):x(2),1:len))));
TCr = squeeze(mean(mean(data.chbr(y(1):y(2),x(1):x(2),1:len))));
for i = 1:length(dpfr1)
    for j = 1:length(dpfr2)
        gcamp = GcampMcCorrection_MIAO(data.blue(:,:,1:len),data.chbr(:,:,1:len),data.chbo(:,:,1:len),baseline,dpfr1(i),dpfr2(j));
        stdmap{i,j} = std(gcamp,[],3);
        meanmap{i,j} = mean(gcamp,3);
        
        TCg{i,j} = squeeze(mean(mean(gcamp(y(1):y(2),x(1):x(2),1:len))));
        disp(['Converted [' mat2str(dpfr1(i)) ', ' mat2str(dpfr2(j)) ']'])
    end
end

a = 1;
for i = 1:length(dpfr1)
    for j = 1:length(dpfr2)
        subplot(numel(dpfr1),numel(dpfr2),a)
        imagesc(stdmap{i,j}); axis image; axis off
        a = a+1;
    end
end
a = 1;
for i = 1:length(dpfr1)
    for j = 1:length(dpfr2)
        %subplot(length(dpfr1),length(dpfr2),a)
        plotyy(1:364,[TCo';TCr'],1:364,TCg{i,j}')
        title(['dpf1: ' mat2str(dpfr1(i)) ' dpf2: ' mat2str(dpfr2(j))])
        axis off
        drawnow
        %hold on
        pause
        a = a+1;
    end
end