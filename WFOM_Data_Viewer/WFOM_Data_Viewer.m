% New WFOM Data Viewer
% Updated 08/06/2018 by Teresa Zhao
% Last updated 07/04/2019 by Nic Thibodeaux

%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = WFOM_Data_Viewer(varargin)
warning off
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WFOM_Data_Viewer_OpeningFcn, ...
    'gui_OutputFcn',  @WFOM_Data_Viewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function WFOM_Data_Viewer_OpeningFcn(hObject, eventdata, h, varargin)
% MUST CHANGE TO WORKING MIAO FOLDER
addpath(genpath('C:\MIAO'));
h.initDir = 'S:\'; h.m.isgui = 1;
h.output = hObject;
set(0,'DefaultFigureColormap',jet)

guidata(hObject, h);

function varargout = WFOM_Data_Viewer_OutputFcn(hObject, eventdata, h)
varargout{1} = h.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%% CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init_Callback(hObject, eventdata, h)
h.chooseROI.Enable = 'off';
h.plotTCs.Enable = 'off';
h.loadmovie.Enable = 'off';
h.m = makem;
h.m.run = h.run.String;
h.m.stimnum = 1;
h.m.isgui = 0;
if isfield(h,'data')
    h = rmfield(h,'data');
    clear h.data;
end
h.m.CCDdir = fullfile(findmousefolder(h.mouseid.String),'CCD');
h.m.fulldir = fullfile(h.m.CCDdir,h.m.run,[h.m.run '_stim_1']);

h = GetMetaData(h);
h.m.dsf = h.m.height/128;
h.m.nrot = 2; h.m.autocrop = 1;

set(h.tpre,'String',h.m.tpre);
set(h.tstim,'String',h.m.tstim);
set(h.tpost,'String',h.m.tpost);
set(h.movieblue,'Value',0,'Enable','on'); 
set(h.moviered,'Value',0,'Enable','on'); 
set(h.moviegreen,'Value',0,'Enable','on'); 
set(h.moviegcamp,'Value',0,'Enable','on'); set(h.moviehbt,'Value',0,'Enable','on');
set(h.moviehbo,'Value',0,'Enable','on'); set(h.moviehbr,'Value',0,'Enable','on');
set(h.xrange,'String','');
set(h.yrange,'String','');   
guidata(hObject, h);

function loadgo_Callback(hObject, eventdata, h)
set(h.movieblue,'Value',0,'Enable','on'); 
set(h.moviered,'Value',0,'Enable','on'); 
set(h.moviegreen,'Value',0,'Enable','on'); 
set(h.moviegcamp,'Value',0,'Enable','on'); set(h.moviehbt,'Value',0,'Enable','on');
set(h.moviehbo,'Value',0,'Enable','on'); set(h.moviehbr,'Value',0,'Enable','on');
h.figure1.Pointer = 'circle'; drawnow;

% Read in which runs to load/average
h.m.stimstoload = str2num(h.stimstoavg.String);
if isempty(h.m.stimstoload), h.m.stimstoload = 1:h.m.nstims; end

h.m.outputs = '';
if h.loadraw.Value
   h.m.outputs = [h.m.outputs 'rgbl'];
end
if h.loadconverted.Value
   h.m.outputs = [h.m.outputs 'odn'];
end

% Read in user input for flicker correction
if h.corr_flk.Value==1, h.m.corr_flicker = 1:h.m.nLEDs; else h.m.corr_flicker = []; end

% Read in user input for % of run to load
if ~isempty(h.percenttoload.String)
    h.m.loadpct = [0 str2num(h.percenttoload.String)/100];
end

info.m = h.m;
info.stimnumsaveraged = h.m.stimstoload;
yesno = questdlg(sprintf('Filepath: %s \nStims to average: %s \nPre-stim = %ds | Stim = %ds | Post-stim = %ds | downsample = %d',...
    h.m.fulldir,num2str(h.m.stimstoload),h.m.tpre,h.m.tstim,h.m.tpost,h.m.dsf),...
    'Are these the settings you want?',...
    'Yes','No','No');
if strcmp(yesno,'Yes')
    [h.m,h.data] = LoadData(h.m.fulldir,h.m);
    disp('~~~~~~~~~~~~~~~ DONE LOADING ~~~~~~~~~~~~~~~~~~~')
    h.figure1.Pointer = 'arrow'; drawnow;
    h.chooseROI.Enable = 'on';
    set(h.xrange,'String','');
    set(h.yrange,'String','');
else
    h.figure1.Pointer = 'arrow';
end
guidata(hObject, h);

function chooseROI_Callback(hObject, eventdata, h)
h.pkfr = round((h.m.tpre*h.m.framerate/h.m.nLEDs)+(h.m.framerate/h.m.nLEDs));
bl_im = mean(h.data.(h.m.LEDs{1})(:,:,round(((h.m.tpre-3)*h.m.framerate/h.m.nLEDs)):round(((h.m.tpre-1)*h.m.framerate/h.m.nLEDs))),3);
roifig = figure; imagesc(h.data.(h.m.LEDs{1})(:,:,h.pkfr)-bl_im); colormap gray; axis image;
[xroi,yroi,~] = genchoose(1);
h.xroi = sprintf('%d:%d',min([round(xroi(1,1)),round(xroi(1,2))]), max([round(xroi(1,1)),round(xroi(1,2))]));
h.yroi = sprintf('%d:%d',min([round(yroi(1,1)),round(yroi(1,2))]), max([round(yroi(1,1)),round(yroi(1,2))]));
set(h.xrange,'String',h.xroi);
set(h.yrange,'String',h.yroi);
close(roifig)
h.plotTCs.Enable = 'on';
guidata(hObject, h);

function plotTCs_Callback(hObject, eventdata, h)
h.m.stimstoload = str2num(h.stimstoavg.String);
xroi = str2num(h.xroi); xroi = [xroi(1), xroi(end)];
yroi = str2num(h.yroi); yroi = [yroi(1), yroi(end)];
rpr = 3; rpc = 3; cpr = 3; cpc = 3;
tm = linspace(0,h.m.movielength,size(h.data.(h.m.LEDs{1}),3));
rawplots = figure; pos=get(gcf,'Position'); set(gcf,'unit','normalized','Position',[0.1 0.1 .4 .65]);

for i = 1:h.m.nLEDs
    subplot(rpr,rpc,i*3-2); imagesc(h.data.(h.m.LEDs{i})(:,:,100)); colormap gray; axis image off; title(h.m.LEDs{i})
    genchoose_PLOT(xroi,yroi,1,'cyan',1)
    
    subplot(rpr,rpc,i*3-1:i*3); hold on
    tc = squeeze(mean(mean(h.data.(h.m.LEDs{i})(yroi(1):yroi(2),xroi(1):xroi(2),:))));
    plot(tm,tc); ylim([mean(tc)-std(tc)*2 mean(tc)+std(tc)*2])
    hold on; plot([str2double(h.tpre.String), str2double(h.tpre.String)],get(gca,'YLim'),'k','LineWidth',0.5)
    plot([str2double(h.tpre.String)+str2double(h.tstim.String), ...
        str2double(h.tpre.String)+str2double(h.tstim.String)],get(gca,'YLim'),'k','LineWidth',0.5)
    title([h.m.LEDs{i}]); ylabel('Intensity'); xlabel('Time (sec)')
end

convplots = figure;
pos=get(gcf,'Position'); set(gcf,'unit','normalized','Position',[0.5 0.1 .4 .65]);
for i = 1:numel(h.m.conv_vars)
    subplot(cpr,cpc,(i-1)*3+1);  imagesc(h.data.(h.m.conv_vars{i})(:,:,h.pkfr) - mean(h.data(1).(h.m.conv_vars{i})(:,:,30:100),3));...
        axis image; title(h.m.conv_vars{i});
    genchoose_PLOT(xroi,yroi,1,'cyan',1)
    
    subplot(cpr,cpc,i*3-1:i*3); hold on
    tc = squeeze(mean(mean(h.data.(h.m.conv_vars{i})(yroi(1):yroi(2),xroi(1):xroi(2),:))));
    tc(isinf(tc)) = 0;
    plot(tm,tc); ylim([nanmean(tc)-nanstd(tc)*2 nanmean(tc)+nanstd(tc)*2])
    hold on; plot([str2double(h.tpre.String), str2double(h.tpre.String)],get(gca,'YLim'),'k','LineWidth',0.5)
    plot([str2double(h.tpre.String)+str2double(h.tstim.String), ...
        str2double(h.tpre.String)+str2double(h.tstim.String)],get(gca,'YLim'),'k','LineWidth',0.5)
    title([h.m.conv_vars{i}]);
end
guidata(hObject, h);

function loadmovie_Callback(hObject, eventdata, h)
clear moviedata
% get proper number of subplots based on user inputs
if h.movieblue.Value==1 || h.moviegreen.Value==1 || h.moviered.Value==1
    rawr = 1; rawc = h.movieblue.Value + h.moviegreen.Value + h.moviered.Value;
else 
    rawr = 0; rawc = 0;
end
if h.moviegcamp.Value==1 || h.moviehbt.Value==1 || h.moviehbo.Value==1 || h.moviehbr.Value==1
    convr = 1; convc = h.moviegcamp.Value + h.moviehbt.Value + h.moviehbo.Value + h.moviehbr.Value;
else 
    convr = 0; convc = 0;
end
rnum = rawr + convr;    cnum = max(rawc,convc);

titles = {'Blue','Green','Red','GCaMP','HbT','HbO','HbR'};
if     isfield(h.Mdata,'blue')==1, ss = size(h.Mdata.blue);
elseif isfield(h.Mdata,'gcamp')==1, ss = size(h.Mdata.gcamp);
end

ct=1; % get list of variables to make movie of
for i = 1:length(h.datatype)
    if i<=4, eval(['ynmovie = h.movie' h.datatype{i} '.Value;'])
    else eval(['ynmovie = h.movie' h.datatype{i}(2:end) '.Value;'])
    end
    if isfield(h.Mdata,h.datatype{i}) && ynmovie==1
        moviedata{ct} = h.datatype{i};
        movietitles{ct} = titles{i};
        ct=ct+1;
    end
end

% movie of selected variables
moviefig = figure(201); set(moviefig,'unit','normalized','Position',[0.2 0.2 0.6 0.6]);
spct = 1; rct = 1; cct = 1;
for k=1:ss(3)
    hcounter = 1;
    figure(201)
    for i=1:length(moviedata)
        if strcmp(moviedata{i},'blue')==1 || strcmp(moviedata{i},'green')==1 || strcmp(moviedata{i},'red')==1
            if k==1
                eval(['h' num2str(hcounter) '=subplot(rnum,cnum,rct); rct=rct+1;']); 
            else
                eval(['subplot(h' num2str(hcounter) ');']); 
            end
            eval(['imagesc(h.Mdata.' moviedata{i} '(:,:,' num2str(k) ')); colormap(h' num2str(hcounter) ',''gray''); axis image off; title(''' movietitles{i} ''')']);
            eval(['caxis([0 max(max(max(h.Mdata.' moviedata{i} ')))/1.75]);'])
            hcounter=hcounter+1;
        elseif strcmp(moviedata{i},'gcamp')==1
            if k==1
                if rnum>1,      eval(['h' num2str(hcounter) '=subplot(rnum,cnum,cnum+cct); cct=cct+1;']);
                elseif rnum==1, eval(['h' num2str(hcounter) '=subplot(rnum,cnum,cct); cct=cct+1;']);
                end
            else
                eval(['subplot(h' num2str(hcounter) ');']);
            end
            eval(['imagesc(h.Mdata.' moviedata{i} '(:,:,' num2str(k) ')); colormap(h' num2str(hcounter) ',''gray''); axis image off; title(''' movietitles{i} ''')']);
            eval(['caxis([min(min(min(h.Mdata.' moviedata{i} '(:,:,30:end-30))))/5, max(max(max(h.Mdata.' moviedata{i} '(:,:,30:end-30))))/3]);'])
            hcounter=hcounter+1;
        elseif strcmp(moviedata{i},'chbt')==1 || strcmp(moviedata{i},'chbo')==1 || strcmp(moviedata{i},'chbr')==1
            if k==1
                if rnum>1,      eval(['h' num2str(hcounter) '=subplot(rnum,cnum,cnum+cct); cct=cct+1;']);
                elseif rnum==1, eval(['h' num2str(hcounter) '=subplot(rnum,cnum,cct); cct=cct+1;']);
                end
            else
                eval(['subplot(h' num2str(hcounter) ');']);
            end
            eval(['imagesc(h.Mdata.' moviedata{i} '(:,:,' num2str(k) ')); colormap(h' num2str(hcounter) ',''jet''); axis image off; title(''' movietitles{i} ''')']);
            eval(['c = max([abs(min(min(min(h.Mdata.' moviedata{i} '(:,:,30:end-30))))), abs(max(max(max(h.Mdata.' moviedata{i} '(:,:,30:end-30)))))]);'])
            caxis([-c/4 c/4])
            hcounter=hcounter+1;
        end
    end
    if k==1
        tmptitle = get(h1,'title'); tmptitle = get(tmptitle,'String');
        subplot(rnum,cnum,1); title([tmptitle, ', t = ', num2str(round(k/(h.m.framerate/h.m.nLEDs),2)) ,' s']);
    else
        subplot(rnum,cnum,1); title([tmptitle, ', t = ', num2str(round(k/(h.m.framerate/h.m.nLEDs),2)) ,' s']);
    end
    drawnow
end
guidata(hObject, h);

function pausebtn_Callback(hObject, eventdata, h)
button_state = get(h.pausebtn,'Value');
if button_state == get(h.pausebtn,'Max')
    disp('yes')
elseif button_state == get(h.pausebtn,'Min')
    disp('no')
end
guidata(hObject, h)

function closeTC_Callback(hObject, eventdata, h)
h.m.figIDs = [99,101,102,201];

for i = [2,3]
    if ishandle(h.m.figIDs(i))
        close(num2str(h.m.figIDs(i)))
    end
end
clc
guidata(hObject, h);

function closeall_Callback(hObject, eventdata, h)
set(h.output, 'HandleVisibility', 'off');
close all
clc
set(h.output, 'HandleVisibility', 'on');
guidata(hObject, h);


function mouseid_Callback(hObject, eventdata, h)
clc
try
    init_Callback(hObject, eventdata, h)
    h.loadgo.Enable = 'on';
catch
    h.loadgo.Enable = 'off';
end

function run_Callback(hObject, eventdata, h)
clc
try
    init_Callback(hObject, eventdata, h)
    h.loadgo.Enable = 'on';
catch
    h.loadgo.Enable = 'off';
end

%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tpre_CreateFcn(hObject, eventdata, h)
function tstim_CreateFcn(hObject, eventdata, h)
function tpost_CreateFcn(hObject, eventdata, h)
function text2_CreateFcn(hObject, eventdata, h)
function loadraw_CreateFcn(hObject, eventdata, h)
function FigTitle_CreateFcn(hObject, eventdata, h)
function stimstoavg_CreateFcn(hObject, eventdata, h)
function percenttoload_CreateFcn(hObject, eventdata, h)
function xrange_CreateFcn(hObject, eventdata, h)
function yrange_CreateFcn(hObject, eventdata, h)

%%%%%%%%%%%%%%%%%%%%%%%%%% UNUSED CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadraw_Callback(hObject, eventdata, h)
function loadconverted_Callback(hObject, eventdata, h)
function corr_flk_Callback(hObject, eventdata, h)
function stimstoavg_Callback(hObject, eventdata, h)
function percenttoload_Callback(hObject, eventdata, h)
function xrange_Callback(hObject, eventdata, h)
function yrange_Callback(hObject, eventdata, h)
function movieblue_Callback(hObject, eventdata, h)
function moviegreen_Callback(hObject, eventdata, h)
function moviered_Callback(hObject, eventdata, h)
function moviegcamp_Callback(hObject, eventdata, h)
function moviehbt_Callback(hObject, eventdata, h)
function moviehbo_Callback(hObject, eventdata, h)
function moviehbr_Callback(hObject, eventdata, h)
function tpre_Callback(hObject, eventdata, handles)
function tstim_Callback(hObject, eventdata, handles)
function tpost_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%% OTHER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadraw_DeleteFcn(hObject, eventdata, h)
function loadraw_ButtonDownFcn(hObject, eventdata, h)
function loadraw_KeyPressFcn(hObject, eventdata, h)
function run_CreateFcn(hObject, eventdata, handles)
function mouseid_CreateFcn(hObject, eventdata, handles)
