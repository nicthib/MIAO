function out = temporal_corr(H,win)
ncomps = size(H,1);
out = zeros(ncomps,ncomps,size(H,2)/10);
kk=0;
h = waitbar(0,'Calculating temporal correlation...');
for i = 1:10:size(H,2)-win
    kk=kk+1;
    for j = 1:ncomps
        for k = 1:ncomps
            uu = corrcoef(H(j,i:i+win-1),H(k,i:i+win-1));
            out(j,k,kk) = uu(1,2); 
        end
    end
    waitbar(i/size(H,2),h)
end
close(h)