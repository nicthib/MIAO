% [m,data] = LoadData(Path,m) loads a specified WFOM dataset.
%
% Version 2 - fully compatible with all new runs taken on WFOM-odeaux,
% legacy code removed
%
% Inputs:
%
% Path - full path to stim run files (REQUIRED)
%
% NEW for v2: This version only works with newer datasets collected on the
% Zyla, using FLIR cameras and
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
% m.do_PCA = 1          : perform PCA noise reduction on the data before conversion
% m.bkgsub (zyla) = 1   : performs background subtraction using first 3
%                         frames from blank recording
% Outputs:
% m: metadata struct
% data: data struct that includes WFOM and webcam image matrices.

function [m,data] = LoadData(varargin)
% VERSION HISTORY
%
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

% 5-29-19 - added some error dlgs, ported v17 to the current version.
% v1.7    - added cropping for BW file
%         - renamed the output 'rcamp' to 'jrgeco'
%         - fixed the textprogressbar bug

% 7/1/19  - This version only supports datasets from June 2019 on, taken on
%           WFOM-odeaux and (hopefully soon) WFOM 1.
% v2.0

% TO DO
%         - add 'spool index' loading for random stim datasets
%         - add blank frame deletion from zyla data
%         - add stim averaging functionality
%         - FLIR camera support(?)

disp('LoadData version 2.0')
addpath(genpath('/local_mount/space/juno/1/Software/MIAO/'))
warning('off','all')

h.m.fulldir = varargin{1};
try
    h.m.mouse = regexp(h.m.fulldir,'[cm]+[0-9]+(_)[0-9]*','Match');
    h.m.mouse = h.m.mouse{:};
catch
    errordlg('ERROR: mouse number not parsed properly. Please send in mouse number via the m structure')
end
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
    fprintf('No outputs given. Defaulting to HbO, HbR, and Neural...')
end
if ~isfield(h.m,'dsf')
    h.m.dsf = 1;
    fprintf('No downsample given. defaulting to 1...')
end
if ~isfield(h.m,'dpf')
    h.m.dpf = [.5 .6];
    fprintf(['No dpf''s given. Defaulting to ' mat2str(h.m.dpf) '...'])
end

if ~isfield(h.m,'loadpct')
    h.m.loadpct = [0 1];
elseif numel(h.m.loadpct) == 1
    errordlg('Please input two elements for loadpct. [0 1] loads the whole run. Thanks!')
    data = []; m = h.m;
    return
end

if isfield(h.m,'PCAcomps')
    if numel(h.m.PCAcomps) > 1
        errordlg('Please input PCAcomps as a single value (aka the # of comps you want to represent your data with)')
        return
    end
end
if ~isfield(h.m,'baseline')
    h.m.baseline = 30:100;
    disp('No baseline given. Defaulting to 30:100...')
end

h.m.greenfilter = 534; h.m.isgui = 0;
h = GetMetaData(h);


if strcmp(h.m.camera,'zyla') && ~isempty(regexp(h.m.outputs,'[rgbodnl]'))
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
    
    numColumns = numDepths + h.m.offsetfactor;
    numRows = numLatPix;
    newdim = floor([numRows/h.m.dsf, numColumns/h.m.dsf]);
    
    for i = 1:length(dir(fullfile(h.m.fulldir,'*.dat')))-1
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
    if h.m.loadpct(1) == 0
        h.m.loadpct(1) = 1/numel(namesOut);
    end
    h.m.spoolsLoaded = round(numel(namesOut)*h.m.loadpct(1):round(numel(namesOut)*h.m.loadpct(2)));
    filesToLoad = namesOut(h.m.spoolsLoaded);
    FID = fopen(fullfile(regexprep(h.m.fulldir,'[_][0-9]$',''),'0000000000spool.dat'));
    rawData = fread(FID, 'uint16=>uint16');
    fclose(FID);
    numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
    numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
    try
        rawData = (reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
    catch
        disp([sprintf('reshaping failed. Please check offset factor and adjust to proper value.\n Height = %i, Width = %i, factors = ',numRows,numColumns) mat2str(divisor(numel(rawData)/numRows))])
        return
    end
    % This switches height and width is nrot is odd
    if mod(h.m.nrot,2) == 1
        height = h.m.height;
        h.m.height = h.m.width;
        h.m.width = height;
    end
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
                try
                    textprogressbar(round(j*100/length(filesToLoad)));
                catch
                    textprogressbar(sprintf(['Loading ' h.m.mouse ' ' h.m.run ' stim ' mat2str(h.m.stim) '\n']))
                end
            end
        end
    end
    ss = size(data.(h.m.LEDs{1}));
    textprogressbar(sprintf('  Done'))
    
elseif strcmp(h.m.camera,'none')
    data = 'No Data loaded';
    return
end
 
% Rotate
if isfield(h.m,'nrot') && ~isempty(regexp(h.m.outputs,'[rgblodn]'))
    for i = 1:h.m.nLEDs
        data.(h.m.LEDs{i}) = rot90(data.(h.m.LEDs{i}),h.m.nrot);
    end
    disp('Rotated data')
end

% apply mask
if isfield(h.m,'BW')  && ~isempty(regexp(h.m.outputs,'[rgblodn]'))
    for i = 1:h.m.nLEDs
        data.(h.m.LEDs{i}) = data.(h.m.LEDs{i}).*repmat(imresize(h.m.BW,ss(1:2)),[1 1 size(data.(h.m.LEDs{i}),3)]);
    end
    disp('Cropped data')
elseif ~isempty(regexp(h.m.outputs,'[rgblodn]'))
    qans = questdlg('No crop mask (BW) found. Would you like to create one?',...
        '???',...
        'Yes','No','Workspace');
    if strcmp(qans,'Yes')
        tmp = figure;
        if isfield(data,'blue')
            imagesc(std(data.blue,[],3))
        else
            imagesc(std(data.lime,[],3))
        end
        axis image
        h.m.BW = roipoly;
        close(tmp)
        for i = 1:h.m.nLEDs
            data.(h.m.LEDs{i}) = data.(h.m.LEDs{i}).*repmat(imresize(h.m.BW,ss(1:2)),[1 1 size(data.(h.m.LEDs{i}),3)]);
        end
    elseif strcmp(qans,'It''s in my workspace!')
        try        
            evalin('base','assignin(''caller'',''BW'',BW)');
        catch
            errordlg('No variable named BW! Skipiping crop...')
        end
    end
end

if ~isempty(h.m.corr_flicker)
    for i = h.m.corr_flicker
        disp(['Flicker ' h.m.LEDs{i}])
        flkcor = squeeze(nanmean(nanmean(data.(h.m.LEDs{i}),2),1));
        if strcmp(h.m.LEDs{i},'red') || strcmp(h.m.LEDs{i},'green')
            corfact(1,1,:) = nanmean(flkcor)+flkcor-smooth(flkcor,floor(h.m.framerate/6)*2+1);
        else
            corfact(1,1,:) = nanmean(flkcor)+flkcor-smooth(flkcor,floor(h.m.framerate/12)*2+1);
        end
        data.(h.m.LEDs{i}) = nanmean(corfact)*data.(h.m.LEDs{i})./repmat(corfact,[ss(1),ss(2),1]);
    end
end

if isfield(h.m,'corr_flicker_exp')
    for i = h.m.corr_flicker_exp  
        disp(['Flicker ' h.m.LEDs{i}])
        tc = squeeze(nanmean(nanmean(data.(h.m.LEDs{i}),2),1));
        [~,facty]= flkcor_beth_glob(tc',data.(h.m.LEDs{i}),0);
        facty = reshape(facty,[1 1 numel(facty)]);
        tc = reshape(tc,[1 1 numel(tc)]);
        data.(h.m.LEDs{i}) = data.(h.m.LEDs{i}).*repmat(tc,[ss(1) ss(1) 1]);
    end
end

if ~isempty(h.m.smooth)
    for i = h.m.smooth
        disp(['Smoothing ' h.m.LEDs{i}])
        data.(h.m.LEDs{i}) = smooth3(data.(h.m.LEDs{i}),'box',[1 1 3]);
    end
end

if h.m.PCAcomps > 0  && ~isempty(regexp(h.m.outputs,'[rgblodn]'))
    keep = 1:h.m.PCAcomps;
    for i = 1:h.m.nLEDs
        disp(['PCA ' h.m.LEDs{i}])
        [COEFF,SCORE,~,~,h.m.EXPLAIN.(h.m.LEDs{i})] = pca(reshape(data.(h.m.LEDs{i}),[ss(1)*ss(2),ss(3)]),'Centered','off');
        data.(h.m.LEDs{i}) = reshape(SCORE(:,keep)*COEFF(:,keep)',ss);
    end
end
clear COEFF SCORE

if ~isempty(regexp(h.m.outputs,'[R]'))
    [data.rotary,h.m] = LoadRotary(fullfile(h.m.CCDdir,h.m.run,[h.m.run '_stim' mat2str(h.m.stim) '_rotary.txt']),h.m);
    disp('Done loading rotary')
end

if ~isempty(regexp(h.m.outputs,'[odn]'))
    disp('Converting Hemodynamics...')
    [data.chbo,data.chbr,~] = convert_mariel_MIAO(data.green,data.red,'g','r',h.m.baseline,534);
end

if ~isempty(regexp(h.m.outputs,'[n]'))
    if isfield(data,'blue')
        disp('Converting GCaMP...')
        data.gcamp = GcampMcCorrection_MIAO(data.blue,data.chbr,data.chbo,h.m.baseline,h.m.dpf(1),h.m.dpf(2));
    elseif isfield(data,'lime')
        disp('Converting jRGECO...')
        h.m.bkg = 90; % PLACEHOLDER
        data.jrgeco = (data.lime-h.m.bkg)./((abs(data.red-h.m.bkg).^h.m.Dr).*(abs(data.green).^h.m.Dg));
        if ~isfield(h.m,'bgGG')
            h.m.bgGG = mean(data.jrgeco(:,:,h.m.baseline),3);
        end
        data.jrgeco = data.jrgeco./repmat(h.m.bgGG,[1 1 size(data.jrgeco,3)])-1;
    end
end

m = h.m;
disp('Done')