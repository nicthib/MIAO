% [m,data] = LoadData(Path,m) loads a specified WFOM dataset.
%
% Version 1.6
%
% Inputs:
%
% Path - full path to stim run files (REQUIRED)
%
% NEW for v1.6: all options must be passed into the m structure. Use
% m = makem to create a default m structure, then alter your parameters as
% neccesary before running LoadData().
%
% m.outputs - string with these characters: rgblodn12R.
% Red, Green, Blue, Lime, hbO, hbR, Neural. 1 and 2 are for webcams 1 and 2.
% R is for rotary. Order doesn't matter. 
% Default is 'odn'.
%
% m.dsf - downsample factor (Integer). 
% Default is 1.
%
% m.dpf - array of two pathlengths for GCaMP correction (or for jRGECO
% correction)
% Default is [.5 .6]. 
%
% m.baseline - passed as a vector containing what frames
% you want to use as a baseline for correction. 
% Default is [30:100].
%
% m.loadpct - 0-1. This loads a ratio of the frames from the beginning of
% the run, where 1 is the whole run or .5 is half of it.
% Default is 1.
%
% m.PCAcomps - 0 if you don't want to perform PCA, and n if you want it
% reconstructed with components 1:n.
%
% Other optional arguments (placed into m struct): 
% m.noload = 1          : load the metadata only
% m.corr_flicker = 1    : perform flicker correction
% m.gettruestim = 1     : read the stimCCD bin files and get the real stim values
% m.bkgsub (zyla) = 1   : performs background subtraction using first 3
%                         frames from blank recording
% Outputs:
% m: metadata struct
% data: data struct that includes WFOM and webcam image matrices.

function [m,data] = LoadData(varargin)
%           UPDATES
% 5-11-18 - added crop to zyla loading so that outputs are square.
% v1.0    - added flicker correction.

% 5-31-18 - removed large number of inputs and made most of them optional.
% v1.1    - added webcam load functionality.
%         - added stim file reading for accurate stimulus information.

% 6-1-18  - added functionality to load a percentage of the dataset.
% v1.3

% 6-8-18  - added real time of stim writing for first and last files in
% v1.4      directory

% 8-9-18  - added flicker correction for all channels.
% v1.5    - added PCA via the 'do_PCA' flag

% 10-9-18 - added rotation
% v1.6    - updated input parsing to move all options to be placed into m
%           structure.
%         - added JRGECO1a analysis, as well as lime channel support
%         - added the ability to specify how many PCA components you want
%         - added new webcam .avi load functionality
%         - added rotary load code

%           TO DO
%         - add 'spool index' loading for random stim datasets
%         - add blank frame deletion from zyla data
%         - add stim averaging functionality
%         - allow cropping for BW file

disp('LoadData version 1.6')
addpath(genpath('/local_mount/space/juno/1/Software/MIAO/MIAO_v2'))
warning('off','all')
% if ~ischar(varargin{1}) || ~isstruct(varargin{2})
%     disp('Check your inputs. Make sure to send the data path, followed by m. Use m = makem to create a blank m struct.')
%     return
% end

h.m.fulldir = varargin{1};
h.m.mouse = regexp(h.m.fulldir,'[cm]+[0-9]+(_)[0-9]+','Match');
h.m.mouse = h.m.mouse{:};
h.m.run = regexp(h.m.fulldir,'(run)[A-Z]*[0-9]?','Match');
h.m.run = h.m.run{1};
h.m.stim = str2double(regexp(regexprep(h.m.fulldir,'/',''),'[0-9]+$','Match'));
h.m.CCDdir = regexprep(varargin{1},[h.m.run '.*'],'');

% This folds in the fieldnames of the imported m struct into
% the already existing h.m struct. Imported m struct overrides
% duplicate fields.
if numel(varargin) == 2
    f = fieldnames(varargin{2});
    for j = 1:length(f)
        h.m.(f{j}) = varargin{2}.(f{j});
    end
else
    h.m = makem;
    disp('No m found. using default m...')
end
% Determine missing fields and substitute defaults
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

h.m.greenfilter = 534; h.m.isgui = 0;
h = GetMetaData(h); 

if h.m.gettruestim
    disp('Getting true stim info...')
    STIMpath = strrep(h.m.CCDdir,'CCD','stimCCD');
    fid = fopen([STIMpath '/' h.m.run 'stim' mat2str(h.m.stim) '.bin'],'r');
    if fid
        tmp = fread(fid,[6 inf],'double');
        h.m.stimt = tmp(1,1:min(find(tmp(1,:)>h.m.movielength))-1);
        h.m.stimdata = tmp(2,1:numel(h.m.stimt));
        h.m.DAQfs = round(numel(h.m.stimdata)/h.m.movielength);
        [h.m.stimon, h.m.stimoff] = findstims(h.m.stimdata,h.m.stimt,h.m.DAQfs);
        fclose(fid);
    else
        disp('No stim file found.')
    end
end

if h.m.noload
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
    if mean(mean(databg)) > 300;
        databg = databg*0;
    end
    n = 0;
    for i = 0:h.m.nLEDs:round((h.m.nFrames-1-h.m.nLEDs)*h.m.loadpct)
        n = n + 1;
        for j = 1:h.m.nLEDs
            fid = fopen(fullfile(h.m.fulldir,  [h.m.run repmat('0',[1,10-size(num2str(i+j-1),2)]) num2str(i+j-1) '.dat']),'r','l');
            tmp = fread(fid,[h.m.height h.m.width],'uint16','l');
            if h.m.dsf == 1
                data.(h.m.LEDs{j})(:,:,n) = tmp-databg;
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
    FID = fopen(fullfile(regexprep(h.m.fulldir,'[_][0-9]$',''),filesToLoad{j}));
    rawData = fread(FID, 'uint16=>uint16');
    fclose(FID);
    numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
    numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
    rawData = (reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
    for i = 1:h.m.nLEDs
        tmp = rawData(:,:,i);
        tmp = tmp(1:h.m.height,1:h.m.width,:);
        if h.m.bkgsub
            data.bkg.([h.m.LEDs{i}]) = tmp;
        else
            data.bkg.([h.m.LEDs{i}]) = tmp*0;
        end
    end
    
    textprogressbar(sprintf(['Loading ' h.m.mouse ' ' h.m.run ' stim ' mat2str(h.m.stim) '\n']))
    for i = 1:h.m.nLEDs  % preallocate
        data.(h.m.LEDs{i}) = (zeros(h.m.height/h.m.dsf,h.m.width/h.m.dsf,ceil(numFramesPerSpool.*length(filesToLoad)/h.m.nLEDs)));
    end
    
    for k = 1:h.m.nLEDs
        start = k;
        indexx.(h.m.LEDs{k}){1} = start:h.m.nLEDs:numFramesPerSpool;
        dsindex.(h.m.LEDs{k}){1} = 1:size(indexx.(h.m.LEDs{k}){1},2);
        indnum.(h.m.LEDs{k}){1} = length(dsindex.(h.m.LEDs{k}));
    end
    
    for i = 2:length(filesToLoad)
        for k = 1:h.m.nLEDs
            last = max(indexx.(h.m.LEDs{k}){i-1});
            skip = numFramesPerSpool-last;
            indexx.(h.m.LEDs{k}){i} = (h.m.nLEDs-skip-1)+[1:h.m.nLEDs:(numFramesPerSpool-(h.m.nLEDs-skip-1))]; % location in each batch of LEDs
            dsindex.(h.m.LEDs{k}){i} = max(dsindex.(h.m.LEDs{k}){i-1})+[1:size(indexx.(h.m.LEDs{k}){i},2)];
            indnum.(h.m.LEDs{k}){i} = length(dsindex.(h.m.LEDs{k}){i});
        end
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
                data.(h.m.LEDs{k})(:,:,dsindex.(h.m.LEDs{k}){j}) = uint16(squeeze(mean(mean(reshape(rawData(1:h.m.height,1:h.m.width,indexx.(h.m.LEDs{k}){j}),[h.m.dsf,h.m.height/h.m.dsf,h.m.dsf,h.m.height/h.m.dsf,numel(indexx.(h.m.LEDs{k}){j})]),3),1))) - repmat(imresize(data.bkg.(h.m.LEDs{k}),1/h.m.dsf),[1 1 numel(indexx.(h.m.LEDs{k}){j})]);
            else
                data.(h.m.LEDs{k})(:,:,dsindex.(h.m.LEDs{k}){j}) = rawData(1:h.m.height,1:h.m.width,indexx.(h.m.LEDs{k}){j})-repmat(data.bkg.(h.m.LEDs{k}),[1 1 numel(indexx.(h.m.LEDs{k}){j})]);
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
    
    textprogressbar(sprintf('  Done'))
    
elseif strcmp(h.m.camera,'none')
    data = 'No Data loaded';
    return
end

% Rotate
if isfield(h.m,'nrot')
    for i = 1:h.m.nLEDs
        data.(h.m.LEDs{i}) = rot90(data.(h.m.LEDs{i}),h.m.nrot);
    end
    disp('Rotated data')
end

% apply mask
if isfield(h.m,'mask')
    for i = 1:h.m.nLEDs
        data.(h.m.LEDs{i}) = data.(h.m.LEDs{i}).*repmat(h.m.mask,[1 1 size(data.(h.m.LEDs{i}),3)]);
    end
    disp('Cropped data')
end

if h.m.corr_flicker
    ss = size(data.(h.m.LEDs{1}));
    for i = 1:h.m.nLEDs
        disp(['FLICKER ' h.m.LEDs{i}])
        flkcor = squeeze(mean(mean((data.(h.m.LEDs{i})),2),1));
        corfact(1,1,:) = mean(flkcor)+flkcor-smooth(flkcor,floor(h.m.framerate/12)*2+1);
        data.blue = mean(corfact)*data.(h.m.LEDs{i})./repmat(corfact,[ss(1),ss(2),1]);
    end
end

if h.m.PCAcomps > 0
    ss = size(data.(h.m.LEDs{1}));
    keep = 1:h.m.PCAcomps;
    for i = 1:h.m.nLEDs
        disp(['PCA ' h.m.LEDs{i}])
        [COEFF,SCORE,~,~,h.m.EXPLAIN.(h.m.LEDs{i})] = pca(reshape(data.(h.m.LEDs{i}),[ss(1)*ss(2),ss(3)]),'Centered','off');
        data.(h.m.LEDs{i}) = reshape(SCORE(:,keep)*COEFF(:,keep)',[ss(1),ss(2),ss(3)]);
        data.(h.m.LEDs{i}) = data.(h.m.LEDs{i})-min(min(min(data.(h.m.LEDs{i}))));
    end
end

if ~isempty(regexp(h.m.outputs,'[1]'))
    disp('Loading Webcam 1...')
    try
        data.webcam1 = LoadWebcam(fullfile(h.m.CCDdir,h.m.run,[h.m.run '_webcam' mat2str(h.m.stim)]),h.m.dsf,h.m,'CCD');
    catch
        [data.webcam1,h.m] = LoadWebcamAVI(fullfile(strrep(h.m.CCDdir,'CCD','webcam'),h.m.run,[h.m.run '_stim' mat2str(h.m.stim) '_cam1.avi']),h.m);
    end
    disp('Done loading webcam 1')
end

if ~isempty(regexp(h.m.outputs,'[2]'))
    disp('Loading Webcam 2...')
    try
        data.webcam2 = LoadWebcam(fullfile(strrep(h.m.CCDdir,'CCD','stimCCD'),[h.m.run 'stim' mat2str(h.m.stim) '_webcam']),h.m.dsf,h.m,'stimCCD');
    catch
        [data.webcam2,h.m] = LoadWebcamAVI(fullfile(strrep(h.m.CCDdir,'CCD','webcam'),h.m.run,[h.m.run '_stim' mat2str(h.m.stim) '_cam2.avi']),m);
    end
    disp('Done loading Webcam 2')
end

if ~isempty(regexp(h.m.outputs,'[R]'))
    [data.rotary,h.m] = LoadRotary(fullfile(h.m.CCDdir,h.m.run,[h.m.run '_stim' mat2str(h.m.stim) '_rotary.txt']),h.m);
    disp('Done loading rotary')
end

if ~isempty(regexp(h.m.outputs,'[od]'))
    disp('Converting Hemodynamics...')
    [data.chbo,data.chbr,~] = convert_mariel_MIAO(data.green,data.red,'g','r',h.m.baseline,534);
    if ~isempty(regexp(h.m.outputs,'[n]'))
        if isfield(data,'blue')
            disp('Converting GCaMP...')
            data.gcamp = GcampMcCorrection_MIAO(data.blue,data.chbr,data.chbo,h.m.baseline,h.m.dpf(1),h.m.dpf(2));
        elseif isfield(data,'lime')
            disp('Converting RCaMP...')
            h.m.bkg = 90; % PLACEHOLDER
            data.rcamp = (data.lime-h.m.bkg)./((abs(data.red-h.m.bkg).^h.m.dpf(1).*(abs(data.green.^h.m.dpf(2)))));
            bgGG=mean(data.rcamp(:,:,h.m.baseline),3);
            data.rcamp=data.rcamp./repmat(bgGG,[1 1 size(data.rcamp,3)]);
        end
    end
end

if ~h.m.isgui
    opts = 'rgblodnn';
    fields = {'red','green','blue','lime','chbo','chbr','gcamp','rcamp'};
    for i = 1:length(opts)
        if isempty(regexp(h.m.outputs,opts(i)))
            try
                data = rmfield(data,Fields{i});
            catch
            end
        end
    end
end
m = h.m;
disp('Done')