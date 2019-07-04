% New WFOM Data Viewer
% Last updated 08/06/2018 by Teresa Zhao

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
% End initialization code - DO NOT EDIT

function WFOM_Data_Viewer_OpeningFcn(hObject, eventdata, h, varargin)
% MUST CHANGE TO WORKING WFOM_DATA_VIEWER FOLDER
addpath(genpath('C:\WFOM_Data_Viewer'));
h.initDir = 'G:\'; h.m.isgui = 1;
h.output = hObject;
set(0,'DefaultFigureColormap',jet)
guidata(hObject, h);


% --- Outputs from this function are returned to the command line.
function varargout = WFOM_Data_Viewer_OutputFcn(hObject, eventdata, h)
varargout{1} = h.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%% CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function browse_Callback(hObject, eventdata, h)
clc;
h.m.mainDir = uigetdir(h.initDir,'Select cmdata directory');
if length(h.m.mainDir)~=2
    set(h.pathname,'String',h.m.mainDir);
end

if isfield(h,'data')
    h = rmfield(h,'data');
    clear h.data;
end

% Read path to info file to get stim info
[a, b, c] = fileparts(h.m.mainDir);
infofile = fullfile(a,b,[b,'_info.txt']);
FID = fopen(infofile, 'r');
infofiledata = fread(FID,'*char')';
fclose(FID);

% Fill in values for pre-stim, stim, and post-stim
stimtmp = strfind(infofiledata,'Stim (s)'); stimtmp = stimtmp(2);
h.m.prestim = str2double(infofiledata(strfind(infofiledata, 'Pre-Stim')+length('Pre-Stim (s)  '):stimtmp-1));
h.m.stim = str2double(infofiledata(stimtmp+length('Stim (s)  '):strfind(infofiledata, 'Post-Stim')-1));
h.m.poststim = str2double(infofiledata(strfind(infofiledata, 'Post-Stim')+length('Post-Stim (s)  '):strfind(infofiledata, 'Movielength')-1));
set(h.prestim,'String',h.m.prestim);
set(h.stim,'String',h.m.stim);
set(h.poststim,'String',h.m.poststim);

% Get framerate, movielength, #stims, and nLEDs
h.m.framerate = str2double(infofiledata(strfind(infofiledata, 'Framerate')+length('Framerate (hz)  '):strfind(infofiledata, 'Exp Time')-1));
h.m.movielength = str2double(infofiledata(strfind(infofiledata, 'Movielength')+length('Movielength (s)  '):strfind(infofiledata, 'NStims')-1));
h.m.nstims = str2double(infofiledata(strfind(infofiledata, 'NStims')+length('NStims   '):strfind(infofiledata, 'LED 1')-1));
set(h.stimstoavg,'String',['1:',num2str(h.m.nstims)]);

tmpb = str2double(infofiledata(strfind(infofiledata, 'LED 1')+length('LED 1 is  Blue ')));
tmpg = str2double(infofiledata(strfind(infofiledata, 'LED 2')+length('LED 2 is  Green ')));
tmpr = str2double(infofiledata(strfind(infofiledata, 'LED 3')+length('LED 3 is  Red ')));
tmps = str2double(infofiledata(strfind(infofiledata, 'LED 4')+length('LED 4 is  Speckle ')));
h.m.nLEDs = (tmpb >0)+ (tmpg>0) + (tmpr>0) + (tmps>0);

% Get binning and downsample size
h.m.binsz = str2double(infofiledata(strfind(infofiledata, 'Bin Size')+length('Bin Size  '):strfind(infofiledata, 'Height')-1));
h.m.height = str2double(infofiledata(strfind(infofiledata, 'Height')+length('Height  '):strfind(infofiledata, 'Width')-1));
h.m.width= str2double(infofiledata(strfind(infofiledata, 'Width')+length('Width  '):strfind(infofiledata, 'Framerate')-1));
h.m.pixres = h.m.height/h.m.binsz;
h.m.dsfactor = h.m.pixres/128;

% reset buttons
set(h.movieblue,'Value',0,'Enable','on'); 
set(h.moviered,'Value',0,'Enable','on'); 
set(h.moviegreen,'Value',0,'Enable','on'); 
set(h.moviegcamp,'Value',0,'Enable','on'); set(h.moviehbt,'Value',0,'Enable','on');
set(h.moviehbo,'Value',0,'Enable','on'); set(h.moviehbr,'Value',0,'Enable','on');

set(h.xrange,'String','');
set(h.yrange,'String','');   

guidata(hObject, h);

function loadgo_Callback(hObject, eventdata, h)
if isfield(h,'data') ==1, h=rmfield(h,'data');  end
if isfield(h,'Mdata')==1, h=rmfield(h,'Mdata'); end

% reset buttons
set(h.movieblue,'Value',0,'Enable','on'); 
set(h.moviered,'Value',0,'Enable','on'); 
set(h.moviegreen,'Value',0,'Enable','on'); 
set(h.moviegcamp,'Value',0,'Enable','on'); set(h.moviehbt,'Value',0,'Enable','on');
set(h.moviehbo,'Value',0,'Enable','on'); set(h.moviehbr,'Value',0,'Enable','on');

h.figure1.Pointer = 'circle'; drawnow;
h.m.dpf = [0.6 0.7];
h.m.greenfilter = 534;
h.m.baseline = round(3*h.m.framerate/h.m.nLEDs):round(9.5*h.m.framerate/h.m.nLEDs);

h.datatype = {'blue','green','red','gcamp','chbt','chbo','chbr'};

% Read in which runs to load/average
h.m.stimstoload = str2num(h.stimstoavg.String);
if isempty(h.m.stimstoload), h.m.stimstoload = 1:h.m.nstims; end

% Read in user input for what to load
if h.loadraw.Value==1 % Raw data
    tmprgb = 'rgb';
	set(h.movieblue,'Value',1); set(h.moviered,'Value',1); set(h.moviegreen,'Value',1); 
else
    tmprgb = ''; 
	set(h.movieblue,'Value',0,'Enable','off'); set(h.moviered,'Value',0,'Enable','off'); set(h.moviegreen,'Value',0,'Enable','off'); 
end 
if h.loadconverted.Value==1 % Converted data
    tmpodn = 'odn'; 
    set(h.moviegcamp,'Value',1); set(h.moviehbt,'Value',1);
    set(h.moviehbo,'Value',1); set(h.moviehbr,'Value',1);
else
    tmpodn = ''; 
    set(h.moviegcamp,'Value',0,'Enable','off'); set(h.moviehbt,'Value',0,'Enable','off');
    set(h.moviehbo,'Value',0,'Enable','off'); set(h.moviehbr,'Value',0,'Enable','off');
end

h.m.varstoload = cat(2,tmprgb,tmpodn);

% Read in user input for flicker correction
if h.corr_flk.Value==1, h.m.flk = 'corr_flicker'; else h.m.flk = []; end

% Read in user input for % of run to load
if ~isempty(h.percenttoload.String)
    h.m.loadpct = str2double(regexp(h.percenttoload.String,'[0-9]*','Match'))/100;
end

info.m = h.m;
info.stimnumsaveraged = h.m.stimstoload;
yesno = questdlg(sprintf('Filepath: %s \nStims to average: %s \nPre-stim = %ds | Stim = %ds | Post-stim = %ds | downsample = %d',...
        h.m.mainDir,num2str(h.m.stimstoload),h.m.stim,h.m.prestim,round(h.m.poststim),h.m.dsfactor),...
        'Are these the settings you want?',...
        'Yes','No','No');

if strcmp(yesno,'Yes')
	[a, b, c] = fileparts(h.m.mainDir);
    for i = h.m.stimstoload
        info.m.runselect = i;
        info.m.filepath = fullfile(a,b,[b,'_stim_',num2str(info.m.runselect)]);
        [~,tmpdata] = LoadData(info.m.filepath,info.m.varstoload,info.m.dsfactor,info.m.dpf,info.m.flk,info.m.baseline);
        
        if exist('nstim')==0, nstim = 1; else nstim = nstim + 1; end
            
        for j = 1:length(h.datatype)
            if isfield(tmpdata,h.datatype{j})
                eval(['data(' num2str(nstim) ').' h.datatype{j} '= tmpdata.' h.datatype{j} '(:,:,1:end-3);']);
                eval(['data(' num2str(nstim) ').' h.datatype{j} '= rot90(data(' num2str(nstim) ').' h.datatype{j} ',2);']);
            end
        end
        clear tmpdata
        
        eval(['data(' num2str(nstim) ').chbt  = data(' num2str(nstim) ').chbo + data(' num2str(nstim) ').chbr;']);
        
        try data = rmfield(data,'blanks'); catch; end
    end

    for j = 1:length(h.datatype)
        if isfield(data,h.datatype{j})     
            eval(['Mdata.' h.datatype{j} ' = mean(cat(4, data.' h.datatype{j} '),4);']);
        end
    end

disp('~~~~~~~~~~~~~~~ DONE LOADING ALL RUNS ~~~~~~~~~~~~~~~~~~~')
h.figure1.Pointer = 'arrow'; drawnow;

h.data  = data;
h.Mdata = Mdata;
set(h.xrange,'String','');
set(h.yrange,'String','');   

clear data Mdata
else
    h.figure1.Pointer = 'arrow'; 
end
guidata(hObject, h);


function chooseROI_Callback(hObject, eventdata, h)
% data = h.data;
% Mdata = h.Mdata;

% select a rectangular region
h.pkfr = round((h.m.prestim*h.m.framerate/h.m.nLEDs)+(h.m.framerate/h.m.nLEDs));

if isfield(h.Mdata,'gcamp'), tmpvar = h.Mdata.gcamp;
else tmpvar = h.Mdata.blue; end

bl = mean(tmpvar(:,:,round(((h.m.prestim-3)*h.m.framerate/h.m.nLEDs)):round(((h.m.prestim-1)*h.m.framerate/h.m.nLEDs))),3);
figure(99); imagesc(tmpvar(:,:,h.pkfr)-bl); colormap gray; axis image;
[xroi,yroi,~] = genchoose(1);
h.xroi = sprintf('%d:%d',min([round(xroi(1,1)),round(xroi(1,2))]), max([round(xroi(1,1)),round(xroi(1,2))]));
h.yroi = sprintf('%d:%d',min([round(yroi(1,1)),round(yroi(1,2))]), max([round(yroi(1,1)),round(yroi(1,2))]));

% fill in xrange and yrange with selected ROI coordinates
set(h.xrange,'String',h.xroi);
set(h.yrange,'String',h.yroi);   
guidata(hObject, h);


function plotTCs_Callback(hObject, eventdata, h)
h.m.stimstoload = str2num(h.stimstoavg.String);
cmap = flipud(jet(length(h.m.stimstoload)));

xroi = str2num(h.xroi); xroi = [xroi(1), xroi(end)];
yroi = str2num(h.yroi); yroi = [yroi(1), yroi(end)];

rpr  = 3; rpc  = 3;
cpr = 4; cpc = 3;

if ishandle(99)==1
    close(num2str(99))
end

if isfield(h.data,'blue') % raw channels
    tm = linspace(0,h.m.movielength,size(h.data(1).blue,3));
    figure(101); pos=get(gcf,'Position'); set(gcf,'unit','normalized','Position',[0.1 0.1 .4 .65]);
    % brain stills
    subplot(rpr,rpc,1); imagesc(h.data(1).blue (:,:,100)); colormap gray; axis image off; title('Blue')
                        genchoose_PLOT(xroi,yroi,1,'cyan',1)
    subplot(rpr,rpc,4); imagesc(h.data(1).green(:,:,100)); colormap gray; axis image off; title('Green')
                        genchoose_PLOT(xroi,yroi,1,'cyan',1)
    subplot(rpr,rpc,7); imagesc(h.data(1).red  (:,:,100)); colormap gray; axis image off; title('Red')
                        genchoose_PLOT(xroi,yroi,1,'cyan',1)
    
    % timecourses
    for i=1:length(h.data)
        subplot(rpr,rpc,[2:3]); hold on
            plot(tm,squeeze(mean(mean(h.data(i).blue(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'color',cmap(i,:))
%             plot(
        subplot(rpr,rpc,[5:6]); hold on
            plot(tm,squeeze(mean(mean(h.data(i).green(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'color',cmap(i,:))
        subplot(rpr,rpc,[8:9]); hold on
            plot(tm,squeeze(mean(mean(h.data(i).red(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'color',cmap(i,:))
        disp(i)
    end
    
    % average timecourse
    subplot(rpr,rpc,[2:3]); hold on; plot(tm,squeeze(mean(mean(h.Mdata.blue(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'k','LineWidth',2)
    subplot(rpr,rpc,[5:6]); hold on; plot(tm,squeeze(mean(mean(h.Mdata.green(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'k','LineWidth',2)
    subplot(rpr,rpc,[8:9]); hold on; plot(tm,squeeze(mean(mean(h.Mdata.red(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'k','LineWidth',2)
    
    % stim on/off lines
    hax=subplot(rpr,rpc,[2:3]); hold on; plot([str2double(h.prestim.String), str2double(h.prestim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
                                         plot([str2double(h.prestim.String)+str2double(h.stim.String), ...
                                               str2double(h.prestim.String)+str2double(h.stim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
    hax=subplot(rpr,rpc,[5:6]); hold on; plot([str2double(h.prestim.String), str2double(h.prestim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
                                         plot([str2double(h.prestim.String)+str2double(h.stim.String), ...
                                               str2double(h.prestim.String)+str2double(h.stim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
    hax=subplot(rpr,rpc,[8:9]); hold on; plot([str2double(h.prestim.String), str2double(h.prestim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
                                         plot([str2double(h.prestim.String)+str2double(h.stim.String), ...
                                               str2double(h.prestim.String)+str2double(h.stim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
    
    % titles and labels
    subplot(rpr,rpc,[2:3]); title('Blue timecourses'); ylabel('Intensity')
    subplot(rpr,rpc,[5:6]); title('Green timecourses'); ylabel('Intensity')
    subplot(rpr,rpc,[8:9]); title('Red timecourses'); ylabel('Intensity'); xlabel('Time (sec)')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(h.data,'gcamp') % converted channels
	tm = linspace(0,h.m.movielength,size(h.data(1).gcamp,3));
    figure(102); pos=get(gcf,'Position'); set(gcf,'unit','normalized','Position',[0.5 0.1 .4 .65]);
    % brain stills
    h1=subplot(cpr,cpc,1);  imagesc(h.data(1).gcamp(:,:,h.pkfr) - mean(h.data(1).gcamp(:,:,30:100),3));...
                                    colormap(h1,'gray'); axis image off; title('GCaMP'); caxis([-0.02 0.04]);
                            genchoose_PLOT(xroi,yroi,1,'cyan',1)
    h2=subplot(cpr,cpc,4);  imagesc(h.data(1).chbt (:,:,h.pkfr) - mean(h.data(1).chbt(:,:,30:100),3));...
                                    colormap(h2,'jet'); axis image off; title('HbT'); caxis([-1 1]*1e-5);
                            genchoose_PLOT(xroi,yroi,1,'k',1)
    h3=subplot(cpr,cpc,7);  imagesc(h.data(1).chbo (:,:,h.pkfr) - mean(h.data(1).chbo(:,:,30:100),3));...
                                    colormap(h3,'jet'); axis image off; title('HbO'); caxis([-1 1]*1e-5);
                            genchoose_PLOT(xroi,yroi,1,'k',1)
    h4=subplot(cpr,cpc,10); imagesc(h.data(1).chbr (:,:,h.pkfr) - mean(h.data(1).chbr(:,:,30:100),3));...
                                    colormap(h4,'jet'); axis image off; title('HbR'); caxis([-1 1]*1e-5);
                            genchoose_PLOT(xroi,yroi,1,'k',1)
    
    % timecourses
    for i=1:length(h.data)
        subplot(cpr,cpc,[2:3]); hold on
            plot(tm,squeeze(mean(mean(h.data(i).gcamp(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'color',cmap(i,:))
        subplot(cpr,cpc,[5:6]); hold on
            plot(tm,squeeze(mean(mean(h.data(i).chbt(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'color',cmap(i,:))
        subplot(cpr,cpc,[8:9]); hold on
            plot(tm,squeeze(mean(mean(h.data(i).chbo(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'color',cmap(i,:))
        subplot(cpr,cpc,[11:12]); hold on
            plot(tm,squeeze(mean(mean(h.data(i).chbr(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'color',cmap(i,:))
    end
    
    % average timecourse
    subplot(cpr,cpc,[2:3]);   hold on; plot(tm,squeeze(mean(mean(h.Mdata.gcamp(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'k','LineWidth',2)
    subplot(cpr,cpc,[5:6]);   hold on; plot(tm,squeeze(mean(mean(h.Mdata.chbt(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'k','LineWidth',2)
    subplot(cpr,cpc,[8:9]);   hold on; plot(tm,squeeze(mean(mean(h.Mdata.chbo(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'k','LineWidth',2)
    subplot(cpr,cpc,[11:12]); hold on; plot(tm,squeeze(mean(mean(h.Mdata.chbr(yroi(1):yroi(2),xroi(1):xroi(2),:)))),'k','LineWidth',2)
    
    % stim on/off lines
    hax=subplot(cpr,cpc,[2:3]); hold on; plot([str2double(h.prestim.String), str2double(h.prestim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
                                         plot([str2double(h.prestim.String)+str2double(h.stim.String), ...
                                               str2double(h.prestim.String)+str2double(h.stim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
    hax=subplot(cpr,cpc,[5:6]); hold on; plot([str2double(h.prestim.String), str2double(h.prestim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
                                         plot([str2double(h.prestim.String)+str2double(h.stim.String), ...
                                               str2double(h.prestim.String)+str2double(h.stim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
    hax=subplot(cpr,cpc,[8:9]); hold on; plot([str2double(h.prestim.String), str2double(h.prestim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
                                         plot([str2double(h.prestim.String)+str2double(h.stim.String), ...
                                               str2double(h.prestim.String)+str2double(h.stim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
    hax=subplot(cpr,cpc,[11:12]); hold on; plot([str2double(h.prestim.String), str2double(h.prestim.String)],get(hax,'YLim'),'k','LineWidth',0.5)
                                         plot([str2double(h.prestim.String)+str2double(h.stim.String), ...
                                               str2double(h.prestim.String)+str2double(h.stim.String)],get(hax,'YLim'),'k','LineWidth',0.5)

	% titles and labels
    subplot(cpr,cpc,[2:3]);   title('GCaMP timecourses'); ylabel('\DeltaF/F')
    subplot(cpr,cpc,[5:6]);   title('HbT timecourses'); ylabel('M')
    subplot(cpr,cpc,[8:9]);   title('HbO timecourses'); ylabel('M')
    subplot(cpr,cpc,[11:12]); title('HbR timecourses'); ylabel('M'); xlabel('Time (sec)')
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

% h.pausebtn = uicontrol(moviefig,'Style','togglebutton','Callback',@pausebtn_Callback,...
%                        'units','normalized','Position',[0.92 0.65 0.05 0.05],...
%                        'Interruptible','off','String','Pause');
% % h.stopbtn  = uicontrol(moviefig,'Style','pushbutton','Callback',@stopbtn_Callback, 'units','normalized','Position',[0.92 0.55 0.05 0.05],'String','Stop');
% h.rewbtn   = uicontrol(moviefig,'Style','pushbutton','Callback',@rewbtn_Callback,  'units','normalized','Position',[0.92 0.45 0.05 0.05],'String','Rew');
% h.fwdbtn   = uicontrol(moviefig,'Style','pushbutton','Callback',@fwdbtn_Callback,  'units','normalized','Position',[0.92 0.35 0.05 0.05],'String','Fwd');

spct = 1; rct = 1; cct = 1;

for k=1:ss(3)
    hcounter = 1;
    figure(201)
    for i=1:length(moviedata)
                        %----------------------------%
        if strcmp(moviedata{i},'blue')==1 || strcmp(moviedata{i},'green')==1 || strcmp(moviedata{i},'red')==1
            if k==1
                eval(['h' num2str(hcounter) '=subplot(rnum,cnum,rct); rct=rct+1;']); 
            else
                eval(['subplot(h' num2str(hcounter) ');']); 
            end
            eval(['imagesc(h.Mdata.' moviedata{i} '(:,:,' num2str(k) ')); colormap(h' num2str(hcounter) ',''gray''); axis image off; title(''' movietitles{i} ''')']);
            eval(['caxis([0 max(max(max(h.Mdata.' moviedata{i} ')))/1.75]);'])
            hcounter=hcounter+1;
                        %----------------------------%
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
                        %----------------------------%
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
            
                        %----------------------------%
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
button_state = get(h.pausebtn,'Value')
    if button_state == get(h.pausebtn,'Max')
        disp('yes')
%         set(h.pausebtn,'String','Play');
    elseif button_state == get(h.pausebtn,'Min')
        disp('no')
%         set(h.pausebtn,'String','Pause');
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


%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pathname_CreateFcn(hObject, eventdata, h)
function prestim_CreateFcn(hObject, eventdata, h)
function stim_CreateFcn(hObject, eventdata, h)
function poststim_CreateFcn(hObject, eventdata, h)
function text2_CreateFcn(hObject, eventdata, h)
function loadraw_CreateFcn(hObject, eventdata, h)
function FigTitle_CreateFcn(hObject, eventdata, h)
function stimstoavg_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function percenttoload_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xrange_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function yrange_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


%%%%%%%%%%%%%%%%%%%%%%%%%% OTHER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadraw_DeleteFcn(hObject, eventdata, h)
function loadraw_ButtonDownFcn(hObject, eventdata, h)
function loadraw_KeyPressFcn(hObject, eventdata, h)
