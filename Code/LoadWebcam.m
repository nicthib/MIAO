function [webcam,time_std] = LoadWebcamJPG(WebcamPath,m)
% [webcam,time_std] = LoadWebcam(WebcamPath,m.dsfWebcam,m,side)
%
% Inputs:
% WebcamPath: full path to webcam images.
% m: Metadata from this run
% side: location of the webcam files, either 'CCD' or 'stimCCD'
%
% Outputs:
% webcam: the webcam, arranged as a 3D matrix interpolated to the framerate
% of the WFOM camera.
% time_std: a measure of the std of the webcam images over time, indicating
% movement.

a = dir(fullfile(WebcamPath,'*.txt'));
FID = fopen(fullfile(WebcamPath, a.name));
tt = fgetl(FID);
frw = str2double(tt((strfind(tt,':')+1):end));
fclose(FID);
% loading
Files = dir(fullfile(WebcamPath,'*.jpg')); Files = {Files.name};
% preallocate
try
    tmp = rgb2gray(imresize(imread(fullfile(WebcamPath,[m.run,'0.jpg'])),1/m.dsfWebcam));
catch
    tmp = rgb2gray(imresize(imread(fullfile(WebcamPath,'0.jpg')),1/m.dsfWebcam));
end
tmp = zeros([size(tmp) length(Files)]); tmp = uint8(tmp);
numIm = length(Files);
for i = 1:numIm-1
    try
        tmp(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[m.run num2str(i) '.jpg'])),1/m.dsfWebcam));
    catch
        tmp(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[num2str(i) '.jpg'])),1/m.dsfWebcam));
    end
end
ss = size(tmp);
[X,Y,Z] = meshgrid(1:ss(2),1:ss(1),linspace(0,ss(3)/frw,ss(3)));
X = single(X); Y = single(Y); Z = single(Z);
[X1,Y1,Z1] = meshgrid(1:ss(2),1:ss(1),linspace(0,m.nFrames/m.framerate,m.nFrames));
X1 = single(X1); Y1 = single(Y1); Z1 = single(Z1);
webcam = interp3(X,Y,Z,single(tmp),X1,Y1,Z1);
time_std = squeeze(mean(std(double(webcam(:,:,1:end-1))-double(webcam(:,:,2:end)),[],1)));
time_std(end) = time_std(end-1);
clear tmp X* Y* Z*
webcam(:,:,1) = [];
