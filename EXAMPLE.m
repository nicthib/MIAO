%% How to use LoadData to load and preprocess WFOM datasets
% First, let's add the required path that has all the functions we need to use LoadData().
addpath(genpath('/local_mount/space/juno/1/Software/MIAO'))

% Let's setup a default m structure by calling makem().
m = makem;

% Now, we need to determine the path to the data we want to load. We can
% use the function findmousefolder() to quickly find where the mouse folder
% is located on the server. For this example, let's use cm62 day 2, run C stim 1.
fulldir = findmousefolder('cm62_2','runF','2');
%% Example 1: Default loading
% Using the path to the data as the only input loads HbO, HbR, and GCaMP
% using standard options. These options will be indicated in the command
% line.
[m,data] = LoadData(fulldir,m);

%% Example 2: Customized loading
% Here's an example of a more customized use of the command. It uses custom
% dpf's, a downsample factor of 2, and corrects for flicker.
% It also loads blue in addition to the converted variables.
m.flk_corr = 1; m.dpf= [.5 .6]; m.sf = 2; m.outputs = 'bodn';
[m,data] = LoadData(fulldir,m);

% For more usage information, use 
% > help LoadData 
% to view the help file, or read the README.md file.