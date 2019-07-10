% [m,data] = LoadData(Path,m) loads WFOM datasets and obtains 
% metadata from the specified run.
%
% Version 2.0
%
% Inputs:
%
% Path - full path to stim run files (REQUIRED)
%
% m - struct containing load options (Optional). Use m = makem to create a
% default m struct, and edit your options (fields of m are explained in
% detail in the help file for makem())
%
% Outputs:
% m: metadata struct
% data: data struct that includes WFOM image matrices.

function [m,data] = LoadData_v2(varargin)
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
% v1.5    - added PCA via the 'do_PCA' flag (REMOVED)

% 10-9-18 - added rotation
% v1.6    - updated input parsing to move all options to be placed into m
%           structure.
%         - added JRGECO1a analysis, as well as lime channel support
%         - added the ability to specify how many PCA components you want
%         - added new webcam .avi load functionality
%         - added rotary load code
%         - all options must be passed into the m structure. Use
%           m = makem to create a default m structure, then alter your 
%           parameters as neccesary before running LoadData().

% 5-29-19 - added some error dlgs, ported v17 to the current version.
% v1.7    - added cropping for BW file
%         - renamed the output 'rcamp' to 'jrgeco'
%         - fixed the textprogressbar bug

% 7-9-19  - Removed ixon load code
% v2.0    - Removed default field assignment. makem replaces this
%         - Cleared out unused functionality (rotary, webcam, etc.). Rotary
%           will be reimplemented in a simpler way,and webcam will be a
%           seperate function (BatchConvertFLIR.m?).
%         - Improved mouse and run detection via strsplit()
%         - GetMetaData() replaced with ReadInfoFile()
%         - Help has been moved to makem.m
%         - Added auxillary and DAQ loading

% TO DO
%         - add 'spool index' loading for random stim datasets
%         - add blank frame deletion from zyla data
%         - add stim averaging functionality for random stim
%         - more comments

disp('LoadData version 2.0')
addpath(genpath('/local_mount/space/juno/1/Software/MIAO/'))
warning('off','all')

m.fulldir = varargin{1};
[m.mouse,m.run,m.stim,m.CCDdir] = WFOM_splitdir(m.fulldir);

% This folds in the fieldnames of the imported m struct into
% the already existing m struct. Imported m struct overrides
% duplicate fields.
if numel(varargin) == 2
    f = fieldnames(varargin{2});
    for j = 1:length(f)
        m.(f{j}) = varargin{2}.(f{j});
    end
else
    m = makem;
    disp('No m found. using default options...')
end

m = ReadInfoFile(m);
try
    m = ReadAuxillary(m);
    disp('aux and DAQ data loaded')
catch
    disp('Unable to load aux and DAQ data')
end

if m.noload
    disp('No data loaded'); 
    data = [];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% ZYLA LOAD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(m.camera,'zyla') && ~isempty(regexp(m.outputs,'[rgbodnl]','once'))
    zylaInfoFilePath = fullfile(m.fulldir,'acquisitionmetadata.ini');
    FID = fopen(zylaInfoFilePath, 'r');
    zylaMetaData = fread(FID, '*char')';
    fclose(FID);
    AOIHeight_start = strfind(zylaMetaData, 'AOIHeight = ');
    AOIWidth_start = strfind(zylaMetaData, 'AOIWidth = ');
    AOIStride_start = strfind(zylaMetaData, 'AOIStride = ');
    PixelEncoding_start = strfind(zylaMetaData, 'PixelEncoding = ');
    ImagesPerFile_start = strfind(zylaMetaData, 'ImagesPerFile = ');
    numDepths = str2double(zylaMetaData(AOIHeight_start+length('AOIHeight = '):...
        AOIWidth_start-1));
    strideWidth = str2double(zylaMetaData(AOIStride_start+length('AOIStride = '):...
        PixelEncoding_start-1));
    numLatPix = strideWidth/2;
    numFramesPerSpool = str2double(zylaMetaData(ImagesPerFile_start+ length('ImagesPerFile = ')...
        :end));
    numColumns = numDepths + m.offsetfactor;
    numRows = numLatPix;    
    for i = 1:length(dir(fullfile(m.fulldir,'*.dat')))-1
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
    if m.loadpct(1) == 0
        m.loadpct(1) = 1/numel(namesOut);
    end
    m.spoolsLoaded = round(numel(namesOut)*m.loadpct(1):round(numel(namesOut)*m.loadpct(2)));
    filesToLoad = namesOut(m.spoolsLoaded);
    FID = fopen(fullfile(regexprep(m.fulldir,'[_][0-9]$',''),'0000000000spool.dat'));
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
    % This switches height and width if nrot is odd
    if mod(m.nrot,2) == 1
        height = m.height;
        m.height = m.width;
        m.width = height;
    end
    for i = 1:m.nLEDs
        tmp = rawData(:,:,i);
        tmp = tmp(1:m.height,1:m.width,:);
        if m.bkgsub
            data.bkg.([m.LEDs{i}]) = tmp;
        else
            data.bkg.([m.LEDs{i}]) = tmp*0;
        end
    end
    textprogressbar(sprintf(['Loading ' m.mouse ' ' m.run ' stim ' mat2str(m.stim) '\n']))
    for i = 1:m.nLEDs  % preallocate
        data.(m.LEDs{i}) = (zeros(m.height/m.dsf,m.width/m.dsf,ceil(numFramesPerSpool.*length(filesToLoad)/m.nLEDs)));
    end
    
    for k = 1:m.nLEDs
        start = k;
        indexx.(m.LEDs{k}){1} = start:m.nLEDs:numFramesPerSpool;
        dsindex.(m.LEDs{k}){1} = 1:size(indexx.(m.LEDs{k}){1},2);
        indnum.(m.LEDs{k}){1} = length(dsindex.(m.LEDs{k}));
    end
    
    for i = 2:length(filesToLoad)
        for k = 1:m.nLEDs
            last = max(indexx.(m.LEDs{k}){i-1});
            skip = numFramesPerSpool-last;
            indexx.(m.LEDs{k}){i} = (m.nLEDs-skip-1)+[1:m.nLEDs:(numFramesPerSpool-(m.nLEDs-skip-1))]; % location in each batch of LEDs
            dsindex.(m.LEDs{k}){i} = max(dsindex.(m.LEDs{k}){i-1})+[1:size(indexx.(m.LEDs{k}){i},2)];
            indnum.(m.LEDs{k}){i} = length(dsindex.(m.LEDs{k}){i});
        end
    end
    
    for j = 1:length(filesToLoad)
        FID = fopen(fullfile(m.fulldir,filesToLoad{j}));
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
        rawData = rawData(1:m.height,1:m.width,:);
        for k = 1:m.nLEDs
            if m.dsf > 1
                data.(m.LEDs{k})(:,:,dsindex.(m.LEDs{k}){j}) = uint16(squeeze(mean(mean(reshape(rawData(1:m.height,1:m.width,indexx.(m.LEDs{k}){j}),[m.dsf,m.height/m.dsf,m.dsf,m.height/m.dsf,numel(indexx.(m.LEDs{k}){j})]),3),1))) - repmat(imresize(data.bkg.(m.LEDs{k}),1/m.dsf),[1 1 numel(indexx.(m.LEDs{k}){j})]);
            else
                data.(m.LEDs{k})(:,:,dsindex.(m.LEDs{k}){j}) = rawData(1:m.height,1:m.width,indexx.(m.LEDs{k}){j})-repmat(data.bkg.(m.LEDs{k}),[1 1 numel(indexx.(m.LEDs{k}){j})]);
            end
        end
        
        if mod(j,10)
            try
                textprogressbar(round(j*100/length(filesToLoad)));
            catch
                textprogressbar(sprintf(['Loading ' m.mouse ' ' m.run ' stim ' mat2str(m.stim) '\r']))
            end
        end
    end
    
    textprogressbar(sprintf('  Done'))
elseif strcmp(m.camera,'ixon')
    disp('LoadData_v2 is not compatible with ixon files. Please use LoadData. Thanks :)')
    data = [];
    m = m;
    return
else
    disp('Camera is unknown, no data loaded')
    data = [];
    m = m;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    ss = size(data.(m.LEDs{1}));
catch
end

% Rotate
if isfield(m,'nrot') && ~isempty(regexp(m.outputs,'[rgblodn]','once'))
    for i = 1:m.nLEDs
        data.(m.LEDs{i}) = rot90(data.(m.LEDs{i}),m.nrot);
    end
    disp('Rotated data')
end

% Autocrop
if m.autocrop
    disp('autocropping...')
    tmp = reshape(data.(m.LEDs{1}),[m.height*m.width/(m.dsf^2), size(data.(m.LEDs{1}),3)]);
    for i = 1:size(tmp,1)
        [~,croptst(i)] = runstest(tmp(i,:));
    end
    croptst(croptst < .00001) = 0;
    croptst(croptst > .00001) = 1;
    croptst = ~logical(reshape(croptst,[m.height/m.dsf m.width/m.dsf]));
    m.BW = imerode(croptst,ones(4));
end

% Apply mask
if isfield(m,'BW')  && ~isempty(regexp(m.outputs,'[rgblodn]','once'))
    for i = 1:m.nLEDs
        data.(m.LEDs{i}) = data.(m.LEDs{i}).*repmat(imresize(m.BW,ss(1:2)),[1 1 size(data.(m.LEDs{i}),3)]);
    end
    disp('Cropped data')
elseif ~isempty(regexp(m.outputs,'[rgblodn]','once'))
    qans = questdlg('No crop mask found. Would you like to create one?','Yes','No');
    if strcmp(qans,'Yes')
        tmp = figure;
        imagesc(std(data.(m.LEDs{1}),[],3))
        axis image
        m.BW = roipoly;
        close(tmp)
        for i = 1:m.nLEDs
            data.(m.LEDs{i}) = data.(m.LEDs{i}).*repmat(imresize(m.BW,ss(1:2)),[1 1 size(data.(m.LEDs{i}),3)]);
        end
    end
end

% Flicker correction
if ~isempty(m.corr_flicker)
    for i = m.corr_flicker
        disp(['Flicker ' m.LEDs{i}])
        flkcor = squeeze(nanmean(nanmean(data.(m.LEDs{i}),2),1));
        if strcmp(m.LEDs{i},'red') || strcmp(m.LEDs{i},'green')
            corfact(1,1,:) = nanmean(flkcor)+flkcor-smooth(flkcor,floor(m.framerate/6)*2+1);
        else
            corfact(1,1,:) = nanmean(flkcor)+flkcor-smooth(flkcor,floor(m.framerate/12)*2+1);
        end
        data.(m.LEDs{i}) = nanmean(corfact)*data.(m.LEDs{i})./repmat(corfact,[ss(1),ss(2),1]);
    end
end

% Smooth (in time)
if ~isempty(m.smooth)
    for i = m.smooth
        disp(['Smoothing ' m.LEDs{i}])
        data.(m.LEDs{i}) = smooth3(data.(m.LEDs{i}),'box',[1 1 3]);
    end
end

% PCA Denoising
if m.PCAcomps > 0  && ~isempty(regexp(m.outputs,'[rgblodn]'))
    keep = 1:m.PCAcomps;
    for i = 1:m.nLEDs
        disp(['PCA ' m.LEDs{i}])
        [COEFF,SCORE,~,~,m.EXPLAIN.(m.LEDs{i})] = pca(reshape(data.(m.LEDs{i}),[ss(1)*ss(2),ss(3)]),'Centered','off');
        data.(m.LEDs{i}) = reshape(SCORE(:,keep)*COEFF(:,keep)',ss);
    end
end
clear COEFF SCORE

% if ~isempty(regexp(m.outputs,'[R]'))
%     [data.rotary,m] = LoadRotary(fullfile(m.CCDdir,m.run,[m.run '_stim' mat2str(m.stim) '_rotary.txt']),m);
%     disp('Done loading rotary')
% end

if ~isempty(regexp(m.outputs,'[odn]','once'))
    disp('Converting Hemodynamics...')
    [data.chbo,data.chbr,~] = convert_mariel_MIAO(data.green,data.red,'g','r',m.baseline,m.greenfilter);
    m.conv_vars = {'chbo','chbr'};
end

if ~isempty(regexp(m.outputs,'[n]','once'))
    if isfield(data,'blue')
        disp('Correcting GCaMP...')
        data.gcamp = GcampMcCorrection_MIAO(data.blue,data.chbr,data.chbo,m.baseline,m.dpf(1),m.dpf(2));
        m.conv_vars = [m.conv_vars 'gcamp'];
    elseif isfield(data,'lime')
        disp('Correcting jRGECO...')
        m.bkg = 0; % PLACEHOLDER
        data.jrgeco = (data.lime-m.bkg)./((abs(data.red-m.bkg).^m.Dr).*(abs(data.green).^m.Dg));
        m.bgGG = mean(data.jrgeco(:,:,m.baseline),3);
        data.jrgeco = data.jrgeco./repmat(m.bgGG,[1 1 size(data.jrgeco,3)])-1;
        m.conv_vars = [m.conv_vars 'jrgeco'];
    end
end

opts = 'rgblodnn';
fields = {'red','green','blue','lime','chbo','chbr','gcamp','jrgeco'};
for i = 1:length(opts)
    if isempty(regexp(m.outputs,opts(i),'once'))
        try
            data = rmfield(data,fields{i});
        catch
        end
    end
end

disp('Done')