% Takes an indexed map (IDX) and extracts timecourses from 3D datasets that
% correspond to those regions. Also erodes if desired
function H = getHfromKmeans(data,IDX,m,n_erode)
for i = 1:max(IDX(:))
    tmp = imfill(logical(reshape(IDX == i,[m.sz m.sz])),'holes');
    if n_erode>0
        for j = 1:n_erode
            tmp = imerode(tmp,ones(3));
        end
    end
    tmp = double(tmp);
    tmp(tmp == 0) = NaN;
    H(i,:) = squeeze(nanmean(nanmean(data.*repmat(tmp,[1 1 size(data,3)]))));
end