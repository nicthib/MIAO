function out = LoadFLIR(path,i,cam,dsf)
tmp = dir(fullfile(path,['*_' mat2str(i) '_cam' mat2str(cam) '.jpg']));
out = imresize(uint8(imread(fullfile(path,tmp.name))),1/dsf);
