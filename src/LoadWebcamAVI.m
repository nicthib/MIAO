function [webcam_i,m] = LoadWebcamAVI(WebcamPath,m)
% [webcam,time_std] = LoadWebcam(WebcamPath,m)
%
% Inputs:
% WebcamPath: full path to webcam AVI.
% m: Metadata from this run
%
% Outputs:
% webcam: the webcam, arranged as a 3D matrix interpolated to the framerate
% of the WFOM camera.

t1 = csvread(strrep(WebcamPath,'.avi','.txt')); t1 = t1(1:end-1);
numFrames = length(t1); 
if isfield(m,'ZylaStartUnix')
    loadFrames = find(m.loadpct(2) * m.movielength + m.ZylaStartUnix > t1 & t1 > m.loadpct(1) * m.movielength + m.ZylaStartUnix);
else
    t1 = t1-t1(1);
    loadFrames = find(m.loadpct(2) * m.movielength > t1 & t1 > m.loadpct(1) * m.movielength);
end
vidObj = VideoReader(WebcamPath);
webcam = zeros(vidObj.Height,vidObj.Width,numel(loadFrames),'uint8');
webcam = read(vidObj,[loadFrames(1) loadFrames(end)]);
webcam = squeeze(webcam(:,:,1,:));
clear vidObj

m.WebcamStartUnix = t1(loadFrames(1));
try
    t2 = linspace(m.ZylaStartUnix,m.ZylaStartUnix+m.movielength,m.nFrames);
catch
    disp('No zyla start time found')
    t2 = linspace(m.WebcamStartUnix,m.WebcamStartUnix+m.movielength,m.nFrames);
end
ss = size(webcam);
webcam = reshape(webcam,[ss(1)*ss(2),numel(loadFrames)]);
webcam_i = zeros([ss(1)*ss(2),numel(t2)]);
disp('Interpolating to WFOM dataset...')
for i = 1:(ss(1)*ss(2))
    webcam_i(i,:) = uint8(interp1(t1(loadFrames),double(webcam(i,:)),t2));
end
webcam_i = reshape(webcam_i,[ss(1:2),numel(t2)]);
webcam_i = imresize(webcam_i,m.dsfWebcam);