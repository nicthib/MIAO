function outpath = findmousefolder(varargin)
outpath =   [];
roots   =   {'/local_mount/space'; 'S:/'; 'D:/'; 'F:/'; '' };         
servers =   {'enterprise/1';         'enterprise/2';
             'lfoivault/lfoivault1'; 'lfoivault/lfoivault2';
             'revault/revault1';     'revault/revault2';
             'voyager/1';            'dingus/1'; '' };
mtypes  =   {'cmdata'; 'mdata'; ''};
exts    =   {'_stroke'; '_glioma'; ''};
for a = 1:numel(roots)
    for b = 1:numel(servers)
        for c = 1:numel(mtypes)
            for d = 1:numel(exts)
                tmppath = fullfile(roots{a},servers{b},mtypes{c},[varargin{1} exts{d}]);
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
else
    disp(['Found Directory: ' outpath])
end