% Quick application of parallel-processed LSQ nonnegative analysis to a
% dataset.
function [H,W,resid] = LSQanalysis(data,H)
ss = size(data);
data = reshape(data,[prod(ss(1:2)),ss(3)]);
data(isnan(data)) = 0;
W = zeros([prod(ss(1:2)),size(H,1)]);
parfor i = 1:prod(ss(1:2))
    [W(i,:),~,resid(i,:)] = lsqnonneg(H',data(i,:)');
end