function out = LoadFLIR(path,i,cam,width)
tmp = dir(fullfile(path,['*_' mat2str(i) '_cam' mat2str(cam) '.jpg']));
im = uint8(imread(fullfile(path,tmp.name)));
out = imresize(im,width/size(im,2));