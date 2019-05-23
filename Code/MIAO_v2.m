% Mouse Experiment Data Manager v2
% Last updated 3/23/2018 by Nic Thibodeaux
% 5/3/18 -updated all uses of '/' to fullfile() to prevent slash direction issues.
%        -Added the ability to not make summary figures.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = MIAO_v2(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MIAO_v2_OpeningFcn, ...
    'gui_OutputFcn',  @MIAO_v2_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function MIAO_v2_OpeningFcn(hObject, ~, h, varargin)
% MUST CHANGE TO WORKING MEDM FOLDER
addpath(genpath('/local_mount/space/juno/1/Software/MIAO/MIAO_v2')) 
h.initDir = '/local_mount/space/'; h.m.isgui = 1;
h.output = hObject;
set(0,'DefaultFigureColor',[.7 .7 .7])
advancedmode = 1;
if advancedmode
   
end
guidata(hObject, h);

function varargout = MIAO_v2_OutputFcn(hObject, eventdata, h)
varargout{1} = h.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function browse_Callback(hObject, eventdata, h)
h.m.mainDir = uigetdir(h.initDir,'Select cmdata directory');
if length(h.m.mainDir) ~= 2
    listDir(h.m.mainDir,h.exps,1);
    set(h.foldername,'String',h.m.mainDir);
end
guidata(hObject,h);

function popupmenu1_Callback(hObject, eventdata, h)
set(h.stimaverage,'Enable','off')
set(h.webcam,'Enable','off')
servroot = {'lfoivault','revault' 'revault'};
h.m.mainDir = fullfile(h.initDir, servroot{ceil(h.popupmenu1.Value/4)}, h.popupmenu1.String{h.popupmenu1.Value});
set(h.foldername,'String',h.m.mainDir);
listDir(h.m.mainDir,h.exps,1);
guidata(hObject,h);

function exps_Callback(hObject, eventdata, h)
h.figure1.Pointer = 'circle'; h.status.String = ''; drawnow;
try
    h.stimaverage.Enable = 'off';
    h.webcam.Enable = 'off';
    expNames = get(h.exps,'String');
    h.m.mouse = expNames{get(h.exps,'Value')};
    h.m.analysisDir = strrep([h.m.mainDir h.m.mouse],'mdata','mdata_CCD_analysis');
    %listDir(h.m.analysisDir,h.m.analysisDir,0);
    h.m.analysisDir = strrep(h.m.analysisDir,'/glioma',''); % to prevent unneccesary folder path
    h.m.CCDdir = fullfile(h.m.mainDir, h.m.mouse, 'CCD');
    listDir(h.m.CCDdir,h.runs,1);
    if ~isdir(h.m.analysisDir)
        mkdir(h.m.analysisDir)
        h.save.Enable = 'on';
        h.status.String = sprintf('\n\nAnalysis directory created'); drawnow
    end
    set(h.outfolder,'String',h.m.analysisDir)
    listDir(h.m.analysisDir,h.m.analysisDir,0);
    set(h.runs,'Enable','on')
catch guierror
    if strcmp(guierror.message,'Permission denied')
        h.status.String = sprintf('\n\n Permission denied to write to analysis folder. Writing to disk is disabled.'); drawnow
        h.save.Enable = 'off'; h.save.Value = 0;
    elseif strcmp(guierror.identifier,'MATLAB:UndefinedFunction')
        h.status.String = sprintf('\n\n There''s something wrong with this mouse folder. Please fix it! Thanks.'); drawnow
    end
end
h.figure1.Pointer = 'arrow'; drawnow;
guidata(hObject,h);

function runs_Callback(hObject, eventdata, h)
h.figure1.Pointer = 'circle'; drawnow;
h.m.runNames = get(h.runs,'String');
try h.m = rmfield(h.m,'run'); catch; end
h.m.runselect = 1; % Placeholder to get a metadata file.
h = GetMetaData(h);

listDir(fullfile(h.m.CCDdir,h.m.run),h.stims,1);
h.m.extrafilesdir1 = fullfile(strrep(h.m.CCDdir,'CCD','stimCCD'), [h.m.run '*webcam*']);
h.m.extrafilesdir2 = fullfile(h.m.CCDdir, h.m.run,  [h.m.run '*webcam*']);

if size(dir(h.m.extrafilesdir1),1)>0
    listDir(h.m.extrafilesdir1,h.extrafiles1,1);
else
    set(h.extrafiles1,'String','');
    h.extrafiles1.Value = [];
end
if size(dir(h.m.extrafilesdir2),1)>0
    listDir(h.m.extrafilesdir2,h.extrafiles2,1);
else
    set(h.extrafiles2,'String','');
    h.extrafiles2.Value = [];
end

set(h.stimaverage,'Enable','on')
set(h.webcam,'Enable','on')
h.figure1.Pointer = 'arrow'; drawnow;
guidata(hObject,h);

function stimaverage_Callback(hObject, eventdata, h)
h.figure1.Pointer = 'circle'; drawnow;
h.m.dpf = [.6 .7];
h.m.greenfilter = 534;
h.m.dsfactor = str2double(h.downsample.String);
h.m.baseline = 30:100;

datatypes = {'red','green','blue','chbr','chbo','gcamp'};

if strcmp(h.stimrunsave.String,'all');
    stimstoload = 1:h.m.nstims;
else
    stimstoload = str2num(h.stimrunsave.String);
end

info.m = h.m;
info.stimnumsaveraged = stimstoload;
info.processing = processing;
yesno = questdlg(sprintf('%s: average stims: %s, \ntstim = %d s, tpre = %d s downsample = %d',h.m.mainDir,num2str(stimstoload),h.m.tstim, h.m.tpre,h.m.dsfactor),'Yes','No');
if strcmp(yesno,'Yes')
    for i = stimstoload
        h.m.runselect = i;
        [~,data] = LoadData(h);
        if (h.chbo.Value + h.chbr.Value + h.gcamp.Value > 0);
            [data] = convertNic(data,h);
        end; 
        try data = rmfield(data,'blanks'); catch; end
        if h.sumfigbox.Value
            h.status.String = sprintf('\n\n Making summary figure..'); drawnow;
            [roi,roinames,tcs] = makeprettyplots(h,data,'');
            tc_all(i) = tcs;
        end
        
        if h.AverageData.Value
            yy = fieldnames(data);
            if ~exist('data_ave')
                for j = 1:length(yy)
                    data_ave.(yy{j}) = data.(yy{j})/length(stimstoload);
                end
            else
                for j = 1:length(yy)
                    data_ave.(yy{j}) = data.(yy{j})/length(stimstoload) + data_ave.(yy{j});
                end
            end
        else
            data_out = data;
            for j = 1:length(fieldnames(data))
                if ~h.(datatypes{j}).Value;
                    data_out = rmfield(data_out,datatypes(j));
                end
            end
            if h.workspace.Value
                h.status.String = sprintf('\n\nassigning to workspace...'); drawnow;
                assignin('base', ['data_' h.m.stimlist{i}], data_out);
            end
            
            if h.save.Value
                h.status.String = sprintf('\n\nsaving to disk...'); drawnow; 
                filename = fullfile(h.m.analysisDir, [h.m.mouse,'_', h.m.stimlist{i}, '.mat']);
                if h.sumfigbox.Value
                    save(filename,'data_out', 'tcs', 'roi', 'roinames','info');%,'-v7.3'
                else
                    save(filename,'data_out','info');
                end
            end
        end
    end
    if h.AverageData.Value
        h.status.String = (sprintf('\n\nMaking summary figure..')); drawnow;
        if h.sumfigbox.Value
            [roi,roinames,tcs] = makeprettyplots(h,data_ave,'ave');
        end
        data_out_ave = data_ave;
        for j = 1:length(fieldnames(data_out_ave))
            if ~h.(datatypes{j}).Value;
                data_out_ave = rmfield(data_out_ave,datatypes(j));
            end
        end
    end
    
    if h.workspace.Value
        h.status.String = sprintf('\n\nassigning to workspace...'); drawnow;
        if h.AverageData.Value
            assignin('base', 'data_out_ave', data_ave);
        end
        if h.sumfigbox.Value
            assignin('base','tc_all', tc_all);
            assignin('base','roinames', roinames);
        end
        assignin('base','info', info);
    end
    
    if h.save.Value && h.AverageData.Value;
        h.status.String = sprintf('\n\nsaving to disk...'); drawnow;
        if h.sumfigbox.Value
            save(fullfile(h.m.analysisDir, [h.m.mouse,'_', h.m.run,'_stim' strrep(h.stimrunsave.String,' ','') '_ave.mat']),'data_out_ave', 'tcs', 'roi', 'roinames','info');%,'-v7.3'
        else
            save(fullfile(h.m.analysisDir, [h.m.mouse,'_', h.m.run,'_stim' strrep(h.stimrunsave.String,' ','') '_ave.mat']),'data_out_ave','info');%,'-v7.3'
        end
    end
    clear data
end
h.figure1.Pointer = 'arrow'; drawnow;
h.status.String = sprintf('\n\nDone.'); drawnow;
guidata(hObject,h);

function webcam_Callback(hObject, eventdata, h)
h.figure1.Pointer = 'circle'; drawnow;
WebcamFileList1 = h.extrafiles1.String;
WebcamFileList2 = h.extrafiles2.String;

if h.WebcamSingleRun.Value
    WebcamFileList1 = WebcamFileList1(h.extrafiles1.Value);
    WebcamFileList2 = WebcamFileList2(h.extrafiles2.Value);
end
if ~isempty(WebcamFileList1)
    loadwebcamruns(h,WebcamFileList1,'stimCCD')
end
if ~isempty(WebcamFileList2)
    loadwebcamruns(h,WebcamFileList2,'CCD')
end
h.figure1.Pointer = 'arrow'; drawnow;

function workspace_Callback(hObject, eventdata, h)
if ~h.workspace.Value
    h.save.Value = 1;
end

function save_Callback(hObject, eventdata, h)
if ~h.save.Value
    h.workspace.Value = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%% EXTERNAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HAVE BEEN MOVED TO MEDM/mfiles

%%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function secs_CreateFcn(hObject, eventdata, h)
function stims_CreateFcn(hObject, eventdata, h)
function runs_CreateFcn(hObject, eventdata, h)
function exps_CreateFcn(hObject, eventdata, h)
function downsample_CreateFcn(hObject, eventdata, h)
function foldername_CreateFcn(hObject, eventdata, h)
function edit4_CreateFcn(hObject, eventdata, h)
function runselect_CreateFcn(hObject, eventdata, h)
function analysisdir_CreateFcn(hObject, eventdata, h)
function status_CreateFcn(hObject, eventdata, h)
function stimrunsave_CreateFcn(hObject, eventdata, h)
function popupmenu1_CreateFcn(hObject, eventdata, h)
function text1_CreateFcn(hObject, eventdata, h)
function text2_CreateFcn(hObject, eventdata, h)
function extrafiles1_CreateFcn(hObject, eventdata, h)
function extrafiles2_CreateFcn(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%% UNUSED CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimrunsave_Callback(hObject, eventdata, h)
function red_Callback(hObject, eventdata, h)
function downsample_Callback(hObject, eventdata, h)
function verbose_Callback(hObject, eventdata, h)
function AverageData_Callback(hObject, eventdata, h)
function analysisdir_Callback(hObject, eventdata, h)
function status_Callback(hObject, eventdata, h)
function extrafiles1_Callback(hObject, eventdata, h)
function stims_Callback(hObject, eventdata, h)
function extrafiles2_Callback(hObject, eventdata, h)
function chooseroi_Callback(hObject, eventdata, handles)
function sumfigbox_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%% DOUBLECLICK CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%

function exps_ButtonDownFcn(hObject, eventdata, h)
if strcmp(get(gcf,'selectiontype'),'alt')
    try cd(h.m.mainDir); catch; end
    return
end

function runs_ButtonDownFcn(hObject, eventdata, h)
if strcmp(get(gcf,'selectiontype'),'alt')
    try cd(h.m.CCDdir); catch; end
    return
end

function stims_ButtonDownFcn(hObject, eventdata, h)
if strcmp(get(gcf,'selectiontype'),'alt')
    try cd(fullfile(h.m.CCDdir, h.m.run)); catch; end
    return
end

function analysisdir_ButtonDownFcn(hObject, eventdata, h)
if strcmp(get(gcf,'selectiontype'),'alt')
    try cd(h.m.analysisDir); catch; end
    return
end
