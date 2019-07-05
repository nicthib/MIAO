function outpath = findmousefolder(varargin)

outpath = [];
roots   =   {'S:/'; '';              '/local_mount/space'
            };         
servers =   {'enterprise/1/';        'enterprise/2/';
             'lfoivault/lfoivault1'; 'lfoivault/lfoivault2';
             'revault/revault1';     'revault/revault2';
             'voyager/1';             'dingus/1'; ''
            };
exts =      {'_stroke';              '_glioma'; ''
            };
mtypes =    {'cmdata';               'mdata'  ; ''
            };
for h = 1:numel(roots)
    for i = 1:numel(servers)
        for j = 1:numel(exts)
            for k = 1:numel(mtypes)
                tmppath = fullfile(roots{h},servers{i},mtypes{k},[varargin{1} exts{j}]);
                if(exist(tmppath))
                    outpath = tmppath;
                    if numel(varargin) == 3
                        tmppath1 = fullfile(tmppath,'CCD',varargin{2},[varargin{2} '_stim' varargin{3}]);
                        tmppath2 = fullfile(tmppath,'CCD',varargin{2},[varargin{2} '_stim_' varargin{3}]);
                        if(exist(tmppath1))
                            outpath = tmppath1;
                        elseif(exist(tmppath2))
                            outpath = tmppath2;
                        else
                            outpath = [];
                        end
                    end
                end
            end
        end
    end
end

if isempty(outpath)
    disp('Sorry, couldn''t find that directory!')
    return
else
    disp(['Found Directory: ' outpath])
end