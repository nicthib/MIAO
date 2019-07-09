% h = ReadInfoFile(h) gets all relevant metadata for the mouse run
% specified from info.txt

function m = ReadInfoFile(m)
fid = fopen(fullfile(m.CCDdir, m.run, [m.run '_info.txt']),'r');
txt = fread(fid,'uint8=>char');
txt = lower(strcat(txt))';
nums = regexp(txt,'\d*\.?\d*','Match');
Fields = {'gain','bitdepth','binsize','height','width','framerate','','tpre','tstim','tpost','movielength','nstims'};
for i = [3:6 8:11]
    m.(Fields{i}) = str2double(nums{i});
end
m.height = m.height/m.binsize;
m.width = m.width/m.binsize;
m.LEDs = regexp(txt,{'blue[1-5]','lime[1-5]','green[1-5]','red[1-5]','speckle[1-5]','cyan[1-5]'},'Match');
m.LEDs = m.LEDs(~cellfun('isempty',m.LEDs)); m.LEDs =  [m.LEDs{:}];
m.nLEDs = numel(m.LEDs);
m.LEDs = regexprep(m.LEDs,'.$','');
if exist(fullfile(m.fulldir,'acquisitionmetadata.ini'))
    m.camera = 'zyla'; 
else 
    m.camera = 'ixon';
end

