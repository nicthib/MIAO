function [m,data] = LoadData_v10(varargin)
% [m,data] = LoadData(Path,Outputs,DSF,Pathlengths,[optional],[optional]) executes for 1 stim run only.
%
% Version 1.0
%
%% Inputs:
% Path: full path to stim run files
% Outputs: string with up to 6 characters: rgbodn.
% Red, Green, Blue, hbO, hbR, Neural. Order doesn't matter
% DSF: downsample factor (Integer).
% Pathlengths: array of two pathlengths for GCaMP correction.
%
% An optional fifth argument can be 'noload' if you just want metadata for
% this run, or a value from 1-3 if you want to rotate the data by x*90 
% degrees ccw (so nrot of 2 will rotate by 180 degrees).
%
% If you want to pass an already created m struct, pass it as the 6th
% argument.

%% Outputs:
% m: the relevant metadata from this run
% data: a struct containing all data that was requested

addpath(genpath('/local_mount/space/juno/1/Software/MEDM/MEDM_v2'))
if numel(varargin) == 1
    h = varargin{1}; 
    % We use h.m instead of m in order to preserve compatibility with the GUI.
    h.m.dsfactor = str2double(h.downsample.String);
elseif numel(varargin) > 1
    if numel(varargin) == 6
        h.m = varargin{6};
    end
    h.m.isgui = 0;
    h.m.fulldir = varargin{1};
    h.m.run = regexp(h.m.fulldir,'(run)[A-Z]*[0-9]?','Match');
    h.m.run = h.m.run{1};
    h.m.CCDdir = regexprep(varargin{1},[h.m.run '.*'],'');
    h = GetMetaData(h);
    processing.downsample = varargin{3};
    h.m.dsfactor = varargin{3};
    processing.dpf = varargin{4};
    processing.greenfilter = 534;
end

try
    if strcmp(varargin{5},'noload')
        data = 0; m = h.m;
        return
    end
catch
end

textprogressbar(sprintf(['Loading ' h.m.run ' ']))

if strcmp(h.m.camera,'ixon')
    % IXON
    newdim = [h.m.height/h.m.dsfactor h.m.width/h.m.dsfactor];
    for i = 1:h.m.nLEDs
        data.(h.m.LEDs{i}) = zeros(newdim(1), newdim(2), floor(h.m.nFrames/h.m.nLEDs));
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
        data.blanks(:,:,i) = mean(mean(reshape(tmp,[h.m.dsfactor,newdim(1),h.m.dsfactor,newdim(2)]),3),1);
        fclose(fid);
    end
    
    databg = squeeze(nanmean(data.blanks,3));
    if mean(mean(databg)) > 1000;
        databg = databg*0;
    end
    n = 0;
    for i = 0:h.m.nLEDs:(h.m.nFrames-1-h.m.nLEDs)
        n = n + 1;
        for j = 1:h.m.nLEDs
            fid = fopen(fullfile(h.m.fulldir,  [h.m.run repmat('0',[1,10-size(num2str(i+j-1),2)]) num2str(i+j-1) '.dat']),'r','l');
            tmp = fread(fid,[h.m.height h.m.width],'uint16','l');
            data.(h.m.LEDs{j})(:,:,n) = squeeze(mean(mean(reshape(tmp,[h.m.dsfactor,newdim(1),h.m.dsfactor,newdim(2)]),3),1))-databg;
            fclose(fid);
        end
        if mod(i,10)
            if h.m.isgui
                h.status.String = sprintf('\n\n %i %% complete',round(i*100/h.m.nFrames)); drawnow
            else
                textprogressbar(round(i*100/h.m.nFrames));
            end
        end
    end
    
elseif strcmp(h.m.camera,'zyla')
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
    %     if (2 == h.m.binsize)
    %         numColumns = numDepths + 1;
    %         numRows = ImageSize/2/numColumns;
    %         if(mod(numRows, 1) ~= 0)
    %             numColumns = numDepths+2;
    %             numRows = ImageSize/2/numColumns;
    %         end
    %     else
    numColumns = numDepths + 2;
    numRows = numLatPix;
    %     end
    
    newdim = floor([numRows/h.m.dsfactor, numColumns/h.m.dsfactor]);
    numColcut = newdim(2)*h.m.dsfactor;
    numRowcut = newdim(1)*h.m.dsfactor;
    
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
    
    j=1;
    FID = fopen(fullfile(h.m.fulldir,filesToLoad{j}));
    rawData = fread(FID, 'uint16=>uint16');
    fclose(FID);
    numFramesInThisSpool = floor(length(rawData)/(numRows*numColumns));
    numPixelsToReshape = numRows * numColumns * numFramesInThisSpool;
    rawData = (reshape(rawData(1:numPixelsToReshape),[numRows,numColumns,numFramesPerSpool]));
    for i = 1:3
        data.([h.m.LEDs{i} '_fullgray']) = rawData(:,:,i);
        %         if get(h.verbose,'Value')
        %             figure(100)
        %             set(gcf,'Position',[17 32 1236 324]);
        %             subplot(1,3,i)
        %             imagesc(rawData(:,:,i)); axis image; colormap gray; colorbar; title(h.m.LEDs{i,ct});
        %         end
    end
    
    % preallocate
    for i = 1:h.m.nLEDs
        data.(h.m.LEDs{i}) = (zeros(newdim(1),newdim(2),ceil(numFramesPerSpool.*length(filesToLoad)/h.m.nLEDs)));
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
        for k = 1:h.m.nLEDs
            if h.m.dsfactor > 1
                data.(h.m.LEDs{k})(:,:,dsindex.(h.m.LEDs{k}){j}) = squeeze(mean(mean(reshape(rawData(1:numRowcut,1:numColcut,indexx.(h.m.LEDs{k}){j}),[h.m.dsfactor,newdim(1),h.m.dsfactor,newdim(2),indnum.(h.m.LEDs{k}){j}]),3),1));
            else
                data.(h.m.LEDs{k})(:,:,dsindex.(h.m.LEDs{k}){j}) = rawData(1:numRowcut,1:numColcut,indexx.(h.m.LEDs{k}){j});
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
elseif strcmp(h.m.camera,'none')
    data = [];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    data.blue = rot90(data.blue,varargin{5});
    data.green = rot90(data.green,varargin{5});
    data.red = rot90(data.red,varargin{5});
catch
end
if numel(varargin) > 1 && ~isempty(regexp(varargin{2},'[orn]'))
    [data] = convertNic(data,30:100,processing,h);
end

if numel(varargin) > 1
    opts = 'rgbodn';
    Fields = {'red','green','blue','chbo','chbr','gcamp'};
    for i = 1:6
        if isempty(regexp(varargin{2},opts(i)))
            data = rmfield(data,Fields{i});
        end
    end
    textprogressbar(sprintf(' Done\n'))
end
try data = rmfield(data,'blanks'); catch; end
m = h.m;