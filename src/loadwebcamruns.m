function loadwebcamruns(h,WebcamFileList,side)
for pp = 1:length(WebcamFileList)
    FileName = WebcamFileList{pp};
    stimname = regexprep(FileName,'webcam|_|stim','');
    webcam.info = h.m;
    webcam.info.dsf = str2double(h.downsample.String);
    if strcmp(side,'CCD')    
        WebcamPath = fullfile(h.m.CCDdir,h.m.run,FileName);
    elseif strcmp(side,'stimCCD')
        WebcamPath = fullfile(strrep(h.m.CCDdir,'CCD',side),FileName);
    end
    
    a = dir(fullfile(WebcamPath,'*.txt'));
    FID = fopen(fullfile(WebcamPath, a.name));
    tt = fgetl(FID);
    frw = str2double(tt((strfind(tt,':')+1):end));
    fclose(FID);
    % loading
    Files = dir(fullfile(WebcamPath,'*.jpg')); Files = {Files.name};
    h.status.String = sprintf('\n\nLoading Webcam Data...'); drawnow
    % preallocate
    if strcmp(side, 'CCD')
        tmp = rgb2gray(imresize(imread(fullfile(WebcamPath,[h.m.run,'0.jpg'])),1/webcam.info.dsf));
    elseif strcmp(side, 'stimCCD')
        tmp = rgb2gray(imresize(imread(fullfile(WebcamPath,'0.jpg')),1/webcam.info.dsf));
    end
    tmp = zeros([size(tmp) length(Files)]);
    numIm = length(Files);
    for i = 1:numIm-1;
        if strcmp(side, 'CCD')
            tmp(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[h.m.run num2str(i) '.jpg'])),1/webcam.info.dsf));
        elseif strcmp(side, 'stimCCD')
            tmp(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[num2str(i) '.jpg'])),1/webcam.info.dsf));
        end
        
        if mod(i,10) == 0 % progress
            h.status.String = sprintf(['\n\n' FileName ' ' side ', \n' num2str(round(100*i/numIm))  ' %% complete']); drawnow;
        end
    end
    ss = size(tmp);
    [X,Y,Z] = meshgrid(1:ss(2),1:ss(1),linspace(0,ss(3)/frw,ss(3)));
    [X1,Y1,Z1] = meshgrid(1:ss(2),1:ss(1),linspace(0,h.m.nFrames/h.m.framerate,h.m.nFrames));
    webcam.data = interp3(X,Y,Z,tmp,X1,Y1,Z1);
    clear tmp X Y Z X1 Y1 Z1
    ss = size(webcam.data);
    webcam.time_std = std(reshape(webcam.data(:,:,1:end-1)-webcam.data(:,:,2:end),[ss(1)*ss(2),ss(3)-1]),[],1);
    webcam.info.webcamfilename = FileName;

    % assigning to workspace
    h.status.String = sprintf('\n\nSaving Webcam...'); drawnow
    if h.workspace.Value;
        h.status.String = sprintf('\n\nAssigning to workspace...'); drawnow;
        assignin('base', ['webcam_' stimname '_' side], webcam);
    end
    
    % saving
    if h.save.Value;
        if ~exist([h.m.analysisdir '/webcam_proc/']);
            mkdir([h.m.analysisdir '/webcam_proc/']);
        end
        save([h.m.analysisdir '/webcam_proc/' stimname '_' side '.mat'],'webcam')
    end
    h.status.String = (sprintf('\n\nDone loading webcam.')); drawnow;
end
