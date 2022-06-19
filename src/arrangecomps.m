% This code arranges spatial and temporal components W (x*y x n) and 
% H(n x t) into an order dictated by the vertical centroids of W. This
% assumes W is square when reshaped.
function I = arrangecomps(W)
sz = sqrt(size(W,1));
for i = 1:size(W,2)
    Wtmp = reshape(W(:,i),[sz sz]);
    a = regionprops(ones(sz,sz), Wtmp, {'Centroid','WeightedCentroid'});
    b(i,:) = a.WeightedCentroid;
end
[~,I] = sort(b(:,2),'descend');
a = 1;
end
