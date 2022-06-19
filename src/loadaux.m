function aux = loadaux(m)
f = dir(fullfile(strrep(m.CCDdir,'CCD','auxillary'),m.run,'*r.txt'));
raw = importdata(fullfile(f.folder,f.name));
for i = 1:numel(raw)
    raw_row = strsplit(raw{i},' ');
    if length(raw_row) == 3
        aux.audio(i) = str2num(raw_row{1});
        aux.rotary(i) = str2num(raw_row{2});
        aux.bb(i) = str2num(raw_row{3});
    else
        aux.audio(i) = NaN;
        aux.rotary(i) = NaN;
        aux.bb(i) = NaN;
    end
end 