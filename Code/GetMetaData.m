% h = GetMetaData(h) gets all relevant metadata for the mouse run
% specified.

function h = GetMetaData(h)
if ~isfield(h.m,'run') %for GUI
    h.m.run = h.m.runNames{h.runs.Value};
end
ccd_dir = dir(fullfile(h.m.CCDdir, h.m.run));
h.m.date = ccd_dir(1).date;
ccd_dir(~[ccd_dir.isdir]) = []; ccd_dir(cellfun('prodofsize', {ccd_dir.name}) < 3) = [];
ccd_dir = {ccd_dir.name};
h.m.stimlist = ccd_dir(~cellfun(@isempty,regexp(ccd_dir,['(' h.m.run '_stim(_)?)[1-9][0-9]*'],'Match')));
stimnums = regexp(h.m.stimlist,'\d*$','Match'); stimnums = str2double([stimnums{:}]);
[~,idx] = sort(stimnums); h.m.stimlist = h.m.stimlist(idx);
fid = fopen(fullfile(h.m.CCDdir, h.m.run, [h.m.run '_info.txt']),'r');
txt = fread(fid,'uint8=>char');
txt = lower(strcat(txt))';
nums = regexp(txt,'\d*\.?\d*','Match');

Fields = {'gain','bitdepth','binsize','height','width','framerate','','tpre','tstim','tpost','movielength','nstims'};
for i = [1:6 8:12]
    h.m.(Fields{i}) = str2double(nums{i});
end
h.m.height = h.m.height/h.m.binsize;
h.m.width = h.m.width/h.m.binsize;

h.m.LEDs = regexp(txt,{'blue[1-5]','lime[1-5]','green[1-5]','red[1-5]','speckle[1-5]','cyan[1-5]'},'Match');
h.m.LEDs = h.m.LEDs(~cellfun('isempty',h.m.LEDs)); h.m.LEDs =  [h.m.LEDs{:}];
h.m.nLEDs = numel(h.m.LEDs);
h.m.LEDs = regexprep(h.m.LEDs,'.$','');

if h.m.isgui
    h.m.fulldir = fullfile(h.m.CCDdir, h.m.run, h.m.stimlist{h.m.runselect});
end

if exist(fullfile(h.m.fulldir,'acquisitionmetadata.ini'))
    h.m.camera = 'zyla'; 
    zylaInfoFilePath = fullfile(h.m.fulldir,'acquisitionmetadata.ini');
    FID = fopen(zylaInfoFilePath, 'r');
    zylaMetaData = fread(FID, '*char')'; fclose(FID);
    ImagesPerFile_start = strfind(zylaMetaData, 'ImagesPerFile = ');
    h.m.numFramesPerSpool = str2double(zylaMetaData(ImagesPerFile_start + length('ImagesPerFile = '):end)); 
    h.m.nFrames = (length(dir(fullfile(h.m.fulldir,'*.dat')))*h.m.numFramesPerSpool);
else 
    h.m.camera = 'ixon';
    h.m.nFrames = length(dir(fullfile(h.m.fulldir,'*.dat')));
end
h.m.t = linspace(0,h.m.movielength,h.m.nFrames);

% Get start frame and end frame dates
filenames = dir(h.m.fulldir);
filestart = filenames(3).name; % skip . and ..
fileend = filenames(end-1).name; % avoid reading acquisitioninfo.ini
[~,str1] = dos(['stat ' h.m.fulldir '/' filestart]);
[~,str2] = dos(['stat ' h.m.fulldir '/' fileend]);
startidx = regexp(str1,'(Modify:)');
endidx = regexp(str2,'(Modify:)');
try
    h.m.StartTime = str1(startidx+8:startidx+30);
    h.m.EndTime = str2(endidx+8:endidx+30);
catch
    disp('Could not get start and end times for this run.')
end

% Try getting unix start time
try
    load(fullfile(regexprep(h.m.CCDdir,'CCD','stimCCD'),[h.m.run '_info.mat']),'info')
    h.m.ZylaStartUnix = info.runstart(h.m.stim);
catch
    disp('No Zyla start time was found');
end

% Try getting stim info file
try
    load(fullfile(strrep(h.m.CCDdir,'CCD','stimCCD'),[h.m.run '_info.mat']));
    h.m.stiminfo = info;
catch 
    disp('no stim info file found')
end

% Info text for GUI
if h.m.isgui
    text1 = {'binsize','height','width','framerate','tpre','tstim','tpost','movielength','nstims','LEDs','camera'};
    for i = 1:length(text1)
        text2{i} = h.m.(text1{i});
        try
            text2{i} = strjoin(text2{i});
            text2{i} = num2str(text2{i});
        catch
        end
    end
    text1{end+1} = 'Date:';   text2{end+1} = h.m.date;
    h.text1.String = text1; h.text2.String = text2;
end
