
% [m,data] = LoadData(Path,...) executes for 1 stim run only.
%
% Version 1.4
%
% Inputs:
%
% Path - full path to stim run files (REQUIRED)
%
% Outputs - string with up to 8 characters: rgbodn12.
% Red, Green, Blue, hbO, hbR, Neural. 1 and 2 are for webcams in CCD and 
% stimCCD folders, respectively. Order doesn't matter. 
% Default is 'odn'.
%
% DSF - downsample factor (Integer). 
% Default is 1.
%
% Pathlengths - array of two pathlengths for GCaMP correction. 
% Default is [.5 .6].
%
% Rotation value equal to x*90 degrees ccw (so nrot of 2 will rotate 
% by 180 degrees). 
% Default is 0.
%
% Baseline: passed as a vector containing what frames
% you want to use as a baseline for corrrection (minimum 5 frames). 
% Default is [30:100].
%
% Percent to load: passed as a percentage string, e.g. '50%'. This loads
% only x% of the frames from the beginning of the run.
%
% m structure: pass it, and it will be preserved.
%
% Other optional string arguments are 'noload' to simply load the metadata,
% 'corr_flicker' to perform flicker correction, or 'gettruestim' to read the
% stimCCD bin files and get the real stim values.
%
% Outputs:
% m: metadata struct
%
% data: data struct that includes WFOM and webcam image matrices.

function [m,data] = LoadData_v14(varargin)
%           UPDATES
% 5-11-18 - added crop to zyla loading so that outputs are square.
%         - added flicker correction.
% 5-31-18 - removed large number of inputs and made most of them optional.
%         - added webcam load functionality.
%         - added stim file reading for accurate stimulus information.
% 6-1-18  - added functionality to load a percentage of the dataset.
% 6-8-18  - added real time of stim writing for first and last files in
%           directory

%           TO DO
%         - add blank frame deletion from zyla data
%         - add stim averaging functionality
%         - reintroduce rotation
%         - allow cropping from BW file

addpath(genpath('/local_mount/space/juno/1/Software/MEDM/MEDM_v2'))
if numel(varargin) == 1 && ~ischar(varargin{1})
    h = varargin{1};
    h = GetMetaData(h);
    h.m.nrot = 3;
    h.m.isgui = 1; noload = 0; corr_flicker = 0;
    % We use h.m instead of m in order to preserve compatibility with the GUI.
    h.m.dsf = str2double(h.downsample.String);
else
    h.m.isgui = 0;
    noload = 0; 
    corr_flicker = 0; 
    gettruestim = 0;
    % Input parsing
    for i = 1:numel(varargin)
        tmp = varargin{i};
        if ischar(tmp) && ~isempty(strfind(tmp,'/')) && isempty(strfind(tmp,'%'))
            h.m.fulldir = tmp;
            h.m.mouse = regexp(h.m.fulldir,'[cm]+[0-9]+(_)[0-9]+','Match');
            h.m.mouse = h.m.mouse{:};
            h.m.run = regexp(h.m.fulldir,'(run)[A-Z]*[0-9]?','Match');
            h.m.run = h.m.run{1};
            h.m.stim = str2double(regexp(h.m.fulldir,'[0-9]+$','Match'));
            h.m.CCDdir = regexprep(varargin{1},[h.m.run '.*'],'');
        elseif strcmp(tmp,'noload')
            noload = 1;
        elseif strcmp(tmp,'corr_flicker')
            corr_flicker = 1;
        elseif strcmp(tmp,'gettruestim')
            gettruestim = 1;
        elseif ischar(tmp) && length(tmp) <= 8 && isempty(strfind(tmp,'%'))
            h.m.outputs = tmp;
        elseif ischar(tmp) && ~isempty(strfind(tmp,'%'))
            h.m.loadpct = str2double(regexp(tmp,'[0-9]*','Match'))/100;
        elseif ~ischar(tmp) && numel(tmp) == 1 && ~isstruct(tmp)
            h.m.dsf = tmp;
        elseif ~ischar(tmp) && numel(tmp) == 2
            h.m.dpf = tmp;
        elseif ~ischar(tmp) && length(tmp) >= 5
            h.m.baseline = tmp;
        elseif isstruct(tmp)
            % This folds in the fieldnames of the imported m struct into
            % the already existing h.m struct. Imported m struct overrides
            % duplicate fields.
            f = fieldnames(tmp);
            for j = 1:length(f)
                h.m.(f{j}) = tmp.(f{j});
            end
        end
    end
    
    % Determine missing fields and substitute defaults
    if ~isfield(h.m,'fulldir')
        errordlg('You must provide the full path to the data as an argument')
        return
    end
    if ~isfield(h.m,'outputs')
       h.m.outputs = 'odn';
       disp(sprintf('No outputs given. Defaulting to HbO, HbR, and Neural...'))
    end
    if ~isfield(h.m,'dsf')
       h.m.dsf = 1;
       disp(sprintf('No downsample given. defaulting to 1...'))
    end
    if ~isfield(h.m,'dpf')
       h.m.dpf = [.5 .6];
       disp(sprintf(['No dpf''s given. Defaulting to ' mat2str(h.m.dpf) '...']))
    end
    if ~isfield(h.m,'loadpct')
       h.m.loadpct = 1;
    end
    if ~isfield(h.m,'baseline')
       h.m.baseline = 30:100;
       disp(sprintf('No baseline given. Defaulting to 30:100...'))
    end
end

h = GetMetaData(h); h.m.greenfilter = 534;

if noload
    data = 'No data loaded'; m = h.m;
    disp('Done')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% IXON LOAD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(h.m.camera,'ixon') && ~isempty(regexp(h.m.outputs,'[rgbodn]'))
    textprogressbar(sprintf(['Loading ' h.m.mouse ' ' h.m.run ' stim ' mat2str(h.m.stim) '\n']))
    
    newdim = [h.m.height/h.m.dsf h.m.width/h.m.dsf];
    for i = 1:h.m.nLEDs
        data.(h.m.LEDs{i}) = zeros(newdim(1), newdim(2), floor(h.m.nFrames*h.m.loadpct/h.m.nLEDs));
    end
    
    % Get stim folder names
    tmp = dir(fullfile(h.m.CCDdir, h.m.run));
    tmp(~[tmp.isdir]) = []; tmp(cellfun('prodofsize', {tmp.name}) < 3) = [];
    stimfolders = {tmp.name}; clear tmp; blankflag = 0;
    tmp = dir(fullfile(h.m.CCDdir, h.m.run,stimfolders{1}));
    tmp(cellfun('isempty', regexpi({tmp.name}, '\.dat$'))) = []; datfiles = {tmp.name};
    for i = 1:length(datfiles)
        fid = fopen(fullfile(h.m.CCDdir, h.m.run,stimfolders{1},datfiles{i}),'r','l');
        tmp = fread(fid,[h.m.height h.m.width],'uint16','l');
        data.blanks(:,:,i) = mean(mean(reshape(tmp,[h.m.dsf,newdim(1),h.m.dsf,newdim(2)]),3),1);
        fclose(fid);
    end
     databg = squeeze(nanmean(data.blanks,3));
%     if mean(mean(databg)) > 1000;
         databg = databg*0;
%     end
    n = 0;
    for i = 0:h.m.nLEDs:round((h.m.nFrames-1-h.m.nLEDs)*h.m.loadpct)
        n = n + 1;
        for j = 1:h.m.nLEDs
            fid = fopen(fullfile(h.m.fulldir,  [h.m.run repmat('0',[1,10-size(num2str(i+j-1),2)]) num2str(i+j-1) '.dat']),'r','l');
            tmp = fread(fid,[h.m.height h.m.width],'uint16','l');
            if h.m.dsf == 1
                data.(h.m.LEDs{j})(:,:,n) = tmp;
            else
                data.(h.m.LEDs{j})(:,:,n) = squeeze(mean(mean(reshape(tmp,[h.m.dsf,newdim(1),h.m.dsf,newdim(2)]),3),1))-databg;
            end
            fclose(fid);
        end
        if mod(i,10)
            if h.m.isgui
                h.status.String = sprintf('\n\n %i %% complete',round(i*100/h.m.nFrames)); drawnow
            else
                textprogressbar(round(i*100/(h.m.nFrames*h.m.loadpct)));
            end
        end
    end
    textprogressbar(sprintf('  Done'))
%%%%%%%%%%%%%%%%%%%%%%%%%%% ZYLA LOAD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(h.m.camera,'zyla') && ~isempty(regexp(h.m.outputs,'[rgbodn]'))
    textprogressbar(sprintf(['Loading ' h.m.mouse ' ' h.m.run ' stim ' mat2str(h.m.stim) '\n']))
    zylaInfoFilePath = fullfile(h.m.fulldir,'acquisitionmetadata.ini');
    FID = fopen(zylaInfoFilePath, 'r');
    zylaMetaData = fread(FID, '*char')';
    fclose(FID);
    AOIHeight_start = strfind(zylaMetaData, 'AOIHeight = ');
    AOIWidth_start = strfind(zylaMetaData, 'AOIWidth = ');
    AOIStride_start = strfind(zylaMetaData, 'AOIStride = ');
    PixelEncoding_start = strfind(zylaMetaData, 'PixelEncoding = ');
    ImageSizeBytes_start = strfind(zylaMetaData, 'ImageSizeBytes = ');
    ImageSizeBytes_end = strfind(zylaMetaData, '[multiimage]')-1;
    ImagesPerFile_start = strfind(zylaMetaData, 'ImagesPerFile = ');
    ImageSize = str2double(zylaMetaData(ImageSizeBytes_start+length('ImageSizeBytes = '):...
        ImageSizeBytes_end));
    numDepths = str2double(zylaMetaData(AOIHeight_start+length('AOIHeight = '):...
        AOIWidth_start-1));
    strideWidth = str2double(zylaMetaData(AOIStride_start+length('AOIStride = '):...
        PixelEncoding_start-1));
    numLatPix = strideWidth/2;
    numFramesPerSpool = str2double(zylaMetaData(ImagesPerFile_start+ length('ImagesPerFile = ')...
        :end));
    
    numColumns = numDepths + 2;
    numRows = numLatPix;    
    newdim = floor([numRows/h.m.dsf, numColumns/h.m.dsf]);
    
    for i = 1:length(dir(fullfile(h.m.fulldir,'*.dat')))-1;
        temp = i;
        for j = 1:10
            a(i,j) = mod(temp, 10^j)/(10^(j-1));
            temp = temp-mod(temp, 10^j);
        end
        tempName = mat2str(a(i, :));
        tempName = tempName(2:end-1);
        tempName = tempName(find(tempName ~= ' '));
        tempName = [tempName 'spool.dat'];
        namesOut{i} = tempName;
    end
    filesToLoad = [{'0000000000spool.dat'} namesOut];
    if h.m.loadpct ~= 1
       filesToLoad(round(numel(filesToLoad)*h.m.loadpct:end)) = []; 
    end
    
    j = 1;
    FID = fopen(fullfile(h.m.fulldir,filesToLoad{j}));
    rawData = fread(FID, 'uint16=>uint16');
    fclose(FID);
    numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
    numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
    rawData = (reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
    for i = 1:3
        data.([h.m.LEDs{i} '_fullgray']) = rawData(:,:,i);
    end
    
    for i = 1:h.m.nLEDs  % preallocate
        data.(h.m.LEDs{i}) = (zeros(h.m.height/h.m.dsf,h.m.width/h.m.dsf,ceil(numFramesPerSpool.*length(filesToLoad)/h.m.nLEDs)));
    end
    
    for k = 1:3;
        if strcmp(h.m.LEDs{k},'blue'); bluestart = k; end
        if strcmp(h.m.LEDs{k},'green'); greenstart = k; end
        if strcmp(h.m.LEDs{k},'red'); redstart = k; end
    end
    
    indexx.blue{1} = bluestart:h.m.nLEDs:numFramesPerSpool;
    dsindex.blue{1} = 1:size(indexx.blue{1},2);
    indnum.blue{1} = length(dsindex.blue{1});
    
    indexx.green{1} = greenstart:h.m.nLEDs:numFramesPerSpool;
    dsindex.green{1} = 1:size(indexx.green{1},2);
    indnum.green{1} = length(dsindex.green{1});
    
    indexx.red{1} = redstart:h.m.nLEDs:numFramesPerSpool;
    dsindex.red{1} = 1:size(indexx.red{1},2);
    indnum.red{1} = length(dsindex.red{1});
    
    for i = 2:length(filesToLoad)
        lastblue = max(indexx.blue{i-1});
        skipblue = numFramesPerSpool-lastblue;
        indexx.blue{i} = (h.m.nLEDs-skipblue-1)+[1:h.m.nLEDs:(numFramesPerSpool-(h.m.nLEDs-skipblue-1))]; % location in each batch of blue LEDs
        dsindex.blue{i} = max(dsindex.blue{i-1})+[1:size(indexx.blue{i},2)];
        indnum.blue{i} = length(dsindex.blue{i});
        
        lastgreen = max(indexx.green{i-1});
        skipgreen = numFramesPerSpool-lastgreen;
        indexx.green{i} = (h.m.nLEDs-skipgreen-1)+[1:h.m.nLEDs:(numFramesPerSpool-(h.m.nLEDs-skipgreen-1))]; % location in each batch of green LEDs
        dsindex.green{i} = max(dsindex.green{i-1})+[1:size(indexx.green{i},2)];
        indnum.green{i} = length(dsindex.green{i});
        
        lastred = max(indexx.red{i-1});
        skipred = numFramesPerSpool-lastred;
        indexx.red{i} = (h.m.nLEDs-skipred-1)+[1:h.m.nLEDs:(numFramesPerSpool-(h.m.nLEDs-skipred-1))]; % location in each batch of red LEDs
        dsindex.red{i} = max(dsindex.red{i-1})+[1:size(indexx.red{i},2)];
        indnum.red{i} = length(dsindex.red{i});
    end
    
    for j = 1:length(filesToLoad)
        FID = fopen(fullfile(h.m.fulldir,filesToLoad{j}));
        rawData = fread(FID, 'uint16=>uint16');
        fclose(FID);
        numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
        numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
        if j == 1
            try
                rawData=(reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
            catch
                numColumns = numDepths + 1;
                numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
                numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
                rawData=(reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
            end
        else
            rawData=(reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
        end
        % Crop rawdata to original resolution
        rawData = rawData(1:h.m.height,1:h.m.width,:);
        
        for k = 1:h.m.nLEDs
            if h.m.dsf > 1
                data.(h.m.LEDs{k})(:,:,dsindex.(h.m.LEDs{k}){j}) = squeeze(mean(mean(reshape(rawData(1:h.m.height,1:h.m.width,indexx.(h.m.LEDs{k}){j}),[h.m.dsf,h.m.height/h.m.dsf,h.m.dsf,h.m.height/h.m.dsf,indnum.(h.m.LEDs{k}){j}]),3),1));
            else
                data.(h.m.LEDs{k})(:,:,dsindex.(h.m.LEDs{k}){j}) = rawData(1:h.m.height,1:h.m.width,indexx.(h.m.LEDs{k}){j});
            end
        end
        
        if mod(j,10)
            if h.m.isgui
                h.status.String = sprintf('\n\n %i %% complete',round(j*100/length(filesToLoad))); drawnow
            else
                textprogressbar(round(j*100/length(filesToLoad)));
            end
        end
    end
    
    if size(data.red,3) < size(data.blue,3); data.red(:,:,end) = data.red(:,:,end-1); end
    if size(data.green,3) < size(data.blue,3); data.green(:,:,end) = data.green(:,:,end-1); end
    data = rmfield(data,{'blue_fullgray','green_fullgray','red_fullgray'}); % Need to implement databg here.
    textprogressbar(sprintf('  Done'))
elseif strcmp(h.m.camera,'none')
    data = 'No Data loaded';
    return
end

if corr_flicker
    data.red = remove_flickering(data.red,h.m.framerate);
    data.green = remove_flickering(data.green,h.m.framerate);
end

if ~isempty(regexp(h.m.outputs,'[1]'))
    [data.webcam1,data.webcam1_std] = LoadWebcam(fullfile(h.m.CCDdir,h.m.run,[h.m.run '_webcam' mat2str(h.m.stim)]),h.m.dsf,h.m,'CCD');
    disp('Done loading webcam 1')
end

if ~isempty(regexp(h.m.outputs,'[2]'))
    [data.webcam2,data.webcam2_std] = LoadWebcam(fullfile(strrep(h.m.CCDdir,'CCD','stimCCD'),[h.m.run 'stim' mat2str(h.m.stim) '_webcam']),h.m.dsf,h.m,'stimCCD');
    disp('Done loading webcam 2')
end

% try
%     data.blue = rot90(data.blue,h.m.nrot);
%     data.green = rot90(data.green,h.m.nrot);
%     data.red = rot90(data.red,h.m.nrot);
% catch
% end

if ~isempty(regexp(h.m.outputs,'[odn]'))
    disp(sprintf('Converting...'))
    ss = size(data.blue);
    
%    Flicker correct all
    
    disp('FLICKER BLUE')
    flkcor = squeeze(mean(mean((data.blue(:,:,:)),2),1));
    corfact(1,1,:) = mean(flkcor)+flkcor-smooth(flkcor,floor(h.m.framerate/12)*2+1);
    data.blue = mean(corfact)*data.blue./repmat(corfact,[ss(1),ss(2),1]);
    
    disp('FLICKER GREEN')
    flkcor = squeeze(mean(mean((data.green(:,:,:)),2),1));
    corfact(1,1,:) = mean(flkcor)+flkcor-smooth(flkcor,floor(h.m.framerate/6)*2+1);
    data.green = mean(corfact)*data.green./repmat(corfact,[ss(1),ss(2),1]);
    
    disp('FLICKER RED')
    flkcor = squeeze(mean(mean((data.red(:,:,:)),2),1));
    corfact(1,1,:) = mean(flkcor)+flkcor-smooth(flkcor,floor(h.m.framerate/6)*2+1);
    data.red = mean(corfact)*data.red./repmat(corfact,[ss(1),ss(2),1]);
    
    keep = 1:50;
    disp('PCA BLUE')
    [COEFF,SCORE,~,~,EXPLAINr] = pca(reshape(data.blue,[ss(1)*ss(2),ss(3)]),'Centered','off');
    data.blue = reshape(SCORE(:,keep)*COEFF(:,keep)',[ss(1),ss(2),ss(3)]);
    data.blue = data.blue-min(data.blue(:));
    
    disp('PCA GREEN')
    [COEFF,SCORE,~,~,EXPLAINg] = pca(reshape(data.green,[ss(1)*ss(2),ss(3)]),'Centered','off');
    data.green = reshape(SCORE(:,keep)*COEFF(:,keep)',[ss(1),ss(2),ss(3)]); 
    data.green= data.green-min(data.green(:));
    
    disp('PCA RED')
    [COEFF,SCORE,~,~,EXPLAINb] = pca(reshape(data.red,[ss(1)*ss(2),ss(3)]),'Centered','off');
    data.red = reshape(SCORE(:,keep)*COEFF(:,keep)',[ss(1),ss(2),ss(3)]); 
    data.red = data.red-min(data.red(:));
    
    disp('CONVERTING AND CORRECTING')
    [data.chbo,data.chbr,~] = convert_mariel_MIAO(data.green,data.red,'g','r',h.m.baseline,534);
    data.gcamp = GcampMcCorrection_MIAO(data.blue,data.chbr,data.chbo,h.m.baseline,h.m.dpf(1),h.m.dpf(2));
end

if gettruestim
    disp('Getting true stim info...')
    STIMpath = strrep(h.m.CCDdir,'CCD','stimCCD');
    fid = fopen([STIMpath '/' h.m.run 'stim' mat2str(h.m.stim) '.bin'],'r');
    if fid
        stimdata = fread(fid,[6 inf],'double');
        h.m.tpre = min(find(stimdata(2,:)>1))./max(find(stimdata(3,:)>1))*(h.m.nFrames/h.m.framerate);
        h.m.tstim = round((max(find(stimdata(2,:)>1))./max(find(stimdata(3,:)>1))*(h.m.nFrames/h.m.framerate))-h.m.tpre);
        h.m.tpost = floor(h.m.nFrames/h.m.framerate)-h.m.tstim-h.m.tpre;
        fclose(fid);
    else
        disp('No stim file found.')
    end
end

if ~h.m.isgui
    opts = 'rgbodn';
    Fields = {'red','green','blue','chbo','chbr','gcamp'};
    for i = 1:6
        if isempty(regexp(h.m.outputs,opts(i)))
            try
                data = rmfield(data,Fields{i});
            catch
            end
        end
    end
end

try data = rmfield(data,'blanks'); catch; end
m = h.m;
disp('Done')