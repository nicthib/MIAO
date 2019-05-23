function [ang_i,m] = LoadRotary(file,m)
FID = fopen(file);
i = 1;
while 1
    tmp = fgetl(FID);
    if ~ischar(tmp), break, end
    comma = strfind(tmp,',');
    if numel(comma) > 1
        a = str2num(tmp(1:comma-1)); if isempty(a); a = NaN; end
        b = str2num(tmp(comma(1)+1:comma(2)-1)); if isempty(b); b = NaN; end
        c = str2num(tmp(comma(2)+1:end)); if isempty(c); c = NaN; end
        ang(i) = a; sound(i) = b; t1(i) = c;
    else
        t1(i) = str2num(tmp(comma+1:end));
        ang(i) = str2num(tmp(1:comma-1));
    end
    i = i+1;
end
[~,idx] = unique(t1);
t1 = t1(idx); ang = ang(idx);
try
    t2 = linspace(m.ZylaStartUnix,m.ZylaStartUnix+m.movielength,m.nFrames);
    ang_i = interp1(t1,ang,t2);
catch
    ang_i = ang;
    m.t_rot = t1;
end
m.RotaryStartUnix = t1(1);

