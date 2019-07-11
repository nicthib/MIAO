function m = ReadAuxillary(m)
rot_dir = dir(fullfile(strrep(m.CCDdir,'CCD','auxillary'),m.run,'*_b.mat'));
DAQ_dir = dir(fullfile(strrep(m.CCDdir,'CCD','auxillary'),m.run,'*_DAQ.mat'));
load(fullfile(rot_dir.folder,rot_dir.name))
load(fullfile(DAQ_dir.folder,DAQ_dir.name))
m.aux = aux';
m.DAQ = DAQdata;

% Determine stim times
w_stim = m.DAQ(1,:);
if sum(find(w_stim > 2.5)) == 0
    m.stimon = 0;
else
    m.stimon = 1;
    [~,LOCS] = findpeaks(diff(w_stim),'minpeakdistance',2e4,'threshold',2.5);
    m.stimtimes = LOCS/1e4;
end