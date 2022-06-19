function mouseid = getrecentdataset(initDir)
dirs = dir(initDir);
mousedirs_tmp = {dirs(:).name};
mousedates_tmp = {dirs(:).date};
mousedirs = {};
mousedates = {};
for i = 1:numel(mousedirs_tmp)
    mousedirs{end+1} = mousedirs_tmp{i};
    mousedates{end+1} = mousedates_tmp{i};
end
mousedates = cellfun(@datenum,mousedates);
mouseid = mousedirs{mousedates==max(mousedates)};