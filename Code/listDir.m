function listDir(path,targetCallback,dirbool)
files = dir(path);
tmp = 1;
filenames = {};
for i = 1:length(files)
    if dirbool == 1
        if files(i).isdir && ~strcmp(files(i).name,'..') && ~strcmp(files(i).name,'.')
            filenames{tmp} = files(i).name;
            tmp=tmp+1;
        end
    elseif dirbool == 0
        filenames{tmp} = files(i).name;
        tmp=tmp+1;
    end
end
targetCallback.String = filenames;
targetCallback.Value = 1;