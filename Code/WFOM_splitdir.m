function [mouse,run,stim,CCDdir] = WFOM_splitdir(fulldir)
fulldir_spl = strsplit(fulldir,'/');
CCD_idx = find(cellfun(@sum,strfind(fulldir_spl,'CCD')));
mouse = fulldir_spl{CCD_idx-1};
run = fulldir_spl{CCD_idx+1};
stim = str2double(char(regexp(fulldir_spl{CCD_idx+2},'\d*','Match')));
CCDdir = strjoin(fulldir_spl(1:CCD_idx),'/');
