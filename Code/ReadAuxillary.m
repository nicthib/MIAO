function m = ReadAuxillary(m)
rot_dir = dir(fullfile(strrep(m.CCDdir,'CCD','auxillary'),m.run,'*_b.mat'));
DAQ_dir = dir(fullfile(strrep(m.CCDdir,'CCD','auxillary'),m.run,'*_DAQ.mat'));
load(fullfile(rot_dir.folder,rot_dir.name))
load(fullfile(DAQ_dir.folder,DAQ_dir.name))
m.aux = aux';
m.DAQ = DAQdata;