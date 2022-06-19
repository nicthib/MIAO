function out = LoadFLIR_dev(path,i,cam,width,crop)
tmp = dir(fullfile(path,['*_' mat2str(i) '_cam' mat2str(cam) '.jpg']));
t = 200; h = 160;
l = 440; w = 200;
if ~isempty(tmp)
    im = uint8(imread(fullfile(path,tmp.name)));
    pupil = im(t:t+h-1,l:l+w-1);
    out = imresize(im,width/size(im,2));
    if crop
        out(1:h,end-w+1:end) = pupil;
    end
else
    disp(['FILE MISSING for frame # ' mat2str(i)])
    out = [];
end
