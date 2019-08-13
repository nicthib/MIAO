function out = imresize_WFOM(in,sz)
out = squeeze(nanmean(nanmean(reshape(in,[sz,size(in,1)/sz,sz,size(in,2)/sz]),1),3));