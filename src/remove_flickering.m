% remove_flickering.m
% Mary Kate Montgomery
% April 2018
%
% Function to remove effect of flickering in red LED in zyla data. 

function out = remove_flickering(in,fps)
% Function to remove effect of flickering in raw zyla data. 

% Take timecourse of entire image
red_tc = squeeze(mean(mean(in,1),2));

% Smooth timecourse (1 sec smoothing factor)
sm = smooth(red_tc,1*fps);

% Multiply red by smoothed timecourse and divide by unsmoothed
out = zeros(size(in));
for i = 1:size(in,3)
    out(:,:,i) = in(:,:,i)*sm(i)/red_tc(i);
end