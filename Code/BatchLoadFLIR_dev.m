% BatchLoadFLIR(mouse,run,frames,framerate,savepath)
% mouse (string)
% run (string)
% frames to load (array)
% framerate (single value)
% savepath (string)
% flip (boolean)
%
function BatchLoadFLIR_dev(mouse,run,vf,fr,savepth,flip)

mdir = findmousefolder(mouse);
webdir = fullfile(mdir,'webcam',run);
vidObj = VideoWriter(fullfile(savepth,[mouse '_' run '_' mat2str(vf(1)) '_' mat2str(numel(vf)) '.avi']));
vidObj.FrameRate = fr;
open(vidObj)
outwidth = 720;
outheight = 540;
chk_sz = 1000;
chk_st = 1:chk_sz:numel(vf);
% determine order
w0 = LoadFLIR_dev(webdir,1,flip,outwidth);
w1 = LoadFLIR_dev(webdir,1,1-flip,outwidth);
tmp = cat(2,w0,w1);
imagesc(tmp); axis image; colormap gray
[x,y] = ginput(2);

% write chk_sz frame chunks, then append to video file
for c = 1:numel(chk_st)-1
    parfor k = 1:(chk_st(c+1)-chk_st(c))
        idx = vf(k+chk_st(c)-1);
        w0 = LoadFLIR_dev(webdir,idx,flip,outwidth);
        w1 = LoadFLIR_dev(webdir,idx,1-flip,outwidth);
        tmp = cat(2,w0,w1);
        tval = round(idx*10/fr)/10;
        txt = sprintf([mouse '_' run ' %.1f sec'],tval);
        a(k).cdata = uint8(mean(insertText(tmp,[10 10],txt,'FontSize',18),3));
        a(k).colormap = gray(256);
    end
    disp(['Writing chunk # ' mat2str(c) ' of ' mat2str(numel(chk_st))])
    writeVideo(vidObj,a)
end
clear a
% Last chunk
parfor k = 1:(numel(vf)-chk_st(end)+1)
    idx = vf(k+chk_st(end)-1);
    w0 = LoadFLIR(webdir,idx,0,outwidth);
    w1 = LoadFLIR(webdir,idx,1,outwidth);
    if isempty(w0)
        w0 = zeros(outheight, outwidth);
    end
    if isempty(w1)
        w1 = zeros(outheight, outwidth);
    end
    tmp = cat(2,w0,w1);
    tval = round(idx*10/fr)/10;
    txt = sprintf([mouse '_' run ' %.1f sec'],tval);
    a(k).cdata = uint8(mean(insertText(tmp,[10 10],txt,'FontSize',18),3));
    a(k).colormap = gray(256);
end
disp('Almost done...')
writeVideo(vidObj,a)
clear a
close(vidObj)
disp('FLIR video written!')
