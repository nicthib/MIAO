% h = GetMetaData(h) gets all relevant metadata for the mouse run
% specified.

function h = ReadInfoFile(h)
fid = fopen(fullfile(h.m.CCDdir, h.m.run, [h.m.run '_info.txt']),'r');
txt = fread(fid,'uint8=>char');
txt = lower(strcat(txt))';
nums = regexp(txt,'\d*\.?\d*','Match');
Fields = {'gain','bitdepth','binsize','height','width','framerate','','tpre','tstim','tpost','movielength','nstims'};
for i = [3:6 8:11]
    h.m.(Fields{i}) = str2double(nums{i});
end
h.m.height = h.m.height/h.m.binsize;
h.m.width = h.m.width/h.m.binsize;
h.m.LEDs = regexp(txt,{'blue[1-5]','lime[1-5]','green[1-5]','red[1-5]','speckle[1-5]','cyan[1-5]'},'Match');
h.m.LEDs = h.m.LEDs(~cellfun('isempty',h.m.LEDs)); h.m.LEDs =  [h.m.LEDs{:}];
h.m.nLEDs = numel(h.m.LEDs);
h.m.LEDs = regexprep(h.m.LEDs,'.$','');
if exist(fullfile(h.m.fulldir,'acquisitionmetadata.ini'))
    h.m.camera = 'zyla'; 
else 
    h.m.camera = 'ixon';
end

