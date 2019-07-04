function listDir(path,targetCallback,is_dir)
files = dir(path);
tmp = 1; files_exist = 0;
filenames = {};
for i = 1:length(files)
    if is_dir == 1
        if files(i).isdir && ~strcmp(files(i).name,'..') && ~strcmp(files(i).name,'.')
            filenames{tmp} = files(i).name;
            tmp=tmp+1;
            files_exist = 1;
        end
    elseif is_dir == 0
        filenames{tmp} = files(i).name;
        tmp=tmp+1;
    end
end
if files_exist
    targetCallback.String = filenames;
    targetCallback.Value = 1;
end
