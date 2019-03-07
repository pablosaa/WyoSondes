function [varargout] = RASOBS_GRAPH_DATA_RAW(varargin)
% Function to create basic visualization from the radiosonde data
% downloaded from the Wyoming University and fetched by the Matlab/Octave
% function  RASOBS_DOWNLOAD_DATA_RAW.
% This function is mostly useful for large dataset rather than for simgle
% profiles.
% (c) 2018 P. Saavedra Garfias, UNIVERSITY OF BERGEN
% Email: pablo.saa@uib.no
% See: LICENSE.TXT
% ---------------------------------------------------------------
% TODO:
% * QI control flag 
% * 
% ---------------------------------------------------------------    
if nargin==0,
    % open file browser to load .mat file
    [data0, metvar, metadata] = update_matfile();
elseif nargin==2,
    data = varargin{1};
    metvar = varargin{2};
    metadata = WyoRS_metadata();    
elseif nargin==3,
    data = varargin{1};
    metvar = varargin{2};
    metadata = varargin{3};
else
    error('two arguments needed: data_profile and radiosonde_indexes');
end

% ----------------------------------------------
% Variables default Initialization
nsonde = length(data0);  % PSG
isonde = min(5,nsonde);  % initial radiosonde to display
span = 2;  % +/- two time steps, it can be increased
rspan = min(nsonde,max(1,[isonde-span:isonde+span]));
nspan = length(rspan);
if nspan>nsonde,
    disp('Too few Profiles found!, try with at least 5.');
end
% Here are defined fine-tune parameters Octave/Matlab dependent:
if exist('OCTAVE_VERSION','builtin'),
    % Fosi: FontSize. 
    FoSi = 16;
    colorlin = summer(nspan);
else
    FoSi = 11;
    colorlin = summer(nspan);
end
i = min(find(rspan==isonde));
colorlin(i,:) = [0 0 0];

%% Including the Profile derivatives for T, Theta, Theta_v
data = Calculate_Derivatives(data0);  % PSG

maxz = 0;
xvar = 'D1TH'; %'TEMP';
yvar = 'HGHT';
names = fieldnames(data);
nvar = length(names);
tmp = strfind(names,xvar);
idxvar = find(arrayfun(@(x) ~isempty(tmp{x}),...
                       [1:nvar]));

% Converting the OBST variable to Matlab datenum format:
Nobs = length(metvar.OBST);
tmp = cell2mat(textscan(sprintf('20%8.2f',metvar.OBST),...
                        '%4f%2f%5f',Nobs,'CollectOutput',1));


%% Convert the series of profiles into a homogeneus matrix:
[PROFILER, TimePROF, Tnumobs] = InitializeProfiler(data,metvar.OBST,12000);


%% Variables to plot in the time series of the RS indices:
tsvar = 13; % 'CAPE' default RS indices to show;
xvar = getfield(metvar,'OBST');
yvar = getfield(metvar, metadata.RSINDICES{tsvar,1});

%% Plotting figure, axes and GUIs 
% creating the default figure:
figure;
set(gcf,'Position',[20 5 700 850],'PaperPositionMode','auto');
clf;
%% Sub-plot 1 for the individual profiles
ha.ax(1) = subplot(6,1,1:3);
hold on;
ha.hh = arrayfun(@(i){plot(data(i).D1TH,1e-3*data(i).HGHT,...
                              'DisplayName',sprintf('T0 %+2d',i-isonde),...
                              'LineWidth',2,'Tag','TEMP',...
                              'Color',colorlin(i-rspan(1)+1,:));...
                  },...
                 rspan,...
                 'UniformOutput',true);

maxz=1e-3*max(arrayfun(@(i) max([data(i).HGHT]),rspan));
set(ha.hh{span+1},'LineWidth',3);  % {isonde-span}
xlabel(sprintf('%s [%s]',metadata.PROFILE{idxvar,2:3}),'FontSize',FoSi);
ylabel(sprintf('%s [K%s]',metadata.PROFILE{2,2:3}),'FontSize',FoSi);
tmp = get(ha.ax(1),'Position').*[1 1 0.75 1];
i = [nspan:-1:1];  % to get the middle line to the top
i(span+1) = '';
i = [span+1, i];   % putting the middle line to the top
if exist('OCTAVE_VERSION','builtin'),
  objarray = cell2mat(ha.hh(i));
else
  objarray = [ha.hh{i}]';
end

set(ha.ax(1),'YLim',[-.05 maxz],'Box','on','YGrid','on',...
             'FontSize',FoSi,'Position',tmp,'Color',0.7*[1 1 1],...
             'Children',objarray);
legend('show');
set(findobj('Tag','legend'),'FontSize',FoSi,'Color','none');

%% Sub-plot 2 is for the 2D time series o RS profiles, default TEMP
ha.ax(2) = subplot(6,1,4:5);
if Nobs==1,
    disp('With only one profile no time-series needed');
    ha.p2d = pcolor([]);
else
    ha.p2d = pcolor(TimePROF,1e-3*PROFILER(:,:,2),PROFILER(:,:,3));
end
shading flat;
axis tight;
tmp = get(ha.ax(2),'Position').*[1 0.9 1 1.04];
ha.cob = colorbar;
set(ha.ax(2),'Position',tmp,'XLim',Tnumobs([1 end]),'XTickLabel',[]);
set(ha.cob,'Position',[tmp(1)+tmp(3)+0.005 tmp(2) 0.025 tmp(4)],...
           'FontSize',FoSi);
caxis([-30 25]);
ylabel('Height a.g.l. [km]','FontSize',FoSi);

%% Sub-plot 3 is for time series of Radiosonde Indices, default CAPE
ha.ax(3) = subplot(6,1,6);
hold on;
ha.htl = plot(Tnumobs(isonde)*[1 1],[min(yvar) max(yvar)],...
              '--k','LineWidth',2);
ha.hts = plot(Tnumobs,yvar,'-','LineWidth',1.5,...
              'Tag', metadata.RSINDICES{tsvar,1});
datetick('x','mm/dd');
tmp = get(ha.ax(3),'Position').*[1 1 1 1.1];
set(ha.ax(3),'Box','on','YGrid','on','FontSize',FoSi,'Position',tmp);
ha.tsl = legend(ha.hts,metadata.RSINDICES{tsvar,2});
set(ha.tsl,'FontSize',FoSi,'Color','none');
ylabel(sprintf('%s [%s]',metadata.RSINDICES{tsvar,[1,3]}),...
       'FontSize',FoSi);
xlabel('Observation date Month/Day','FontSize',FoSi);
tmp = get(ha.ax(3),'XTick');
set(ha.ax(2),'XTick',tmp,'XLim',Tnumobs([1 end]));
set(ha.ax,'TickDir','out','Box','on','FontSize',FoSi,'LineWidth',2);
clear tmp;

%xlabel('Observation date Month/Day','FontSize',15);
% ------------------------------------------------------
% creating the basic GUI
% ------------------------------------------------------
% text with the number of profiles:
uicontrol('Style','text','Units','normalized',...
          'Position',[0.1 0.95 0.2 0.0234],...
          'String',sprintf('N = %d',nsonde),'FontSize',FoSi);
% Text indicating the actual showing time:
ha.tdate = uicontrol('style','text','Units','normalized',...
                    'Position',[0.3 0.95 0.37 0.0234],...
                    'String',['Date T0: ' datestr(Tnumobs(isonde),...
                                                  'dd.mm.yyyy @ HHZ')],...
                    'FontSize',FoSi);

% Slider for changing the showing time:
ha.SliderH = uicontrol('style','slider','Units','normalized',...
                       'Position',[0.68 0.95 0.3 0.0234],...
                       'min',1, 'max', length(xvar),'Value',isonde,...
                       'SliderStep',[0.005 0.05],...
                       'Callback',{@update_time0,ha,Tnumobs,data}); 


% Button to select time point from time-series:
ha.SpotH = uicontrol('style','pushbutton','Units','normalized',...
                     'Position',[0.75 0.65 0.2 0.0234],...
                     'String','pick-up T0','Tag','pushme',...
                     'Callback',{@pick_time0,ha,Tnumobs,data});

% Button to open a new file:
ha.OpenfH = uicontrol('style','pushbutton','Units','normalized',...
                      'Position',[0.75 0.60 0.2 0.0234],...
                      'String','Open MAT',...
                      'Callback',{@update_matfile,ha});

% Button to close the GUI:
uicontrol('style','pushbutton','Units','normalized',...
          'Position',[0.75 0.55 0.2 0.0234],...
          'String','Close','Callback','close(gcf)');
% addlistener(SliderH, 'Value', 'PostSet', @callbackfn);
% Menu with the Profile variables
uicontrol('style','text','String','Profile Variable:','FontSize',FoSi,...
          'Units','normalized','Position',[0.73 0.88 0.26 0.0234]);
ha.ScrollH = uicontrol('style','popupmenu','String',names,...
                       'Units','normalized','Value',3,'Tag','ScrollH',...
                       'Position',[0.75 0.85 0.15 0.0234],...
                       'Callback',{@update_profile,ha,data,PROFILER,metadata.PROFILE});

% Menu with the Radiosonde Indices variables
uicontrol('style','text','String','Radiosonde Indices:','FontSize',FoSi,...
          'Units','normalized','Position',[0.73 0.8 0.26 0.0234]);
ha.SIndexH = uicontrol('Style','popupmenu',...
                       'String',fieldnames(metvar),...
                       'Units','normalized','Value',13,...
                       'Position',[0.75 0.77 0.15 0.0234],...
                       'Callback',{@update_rsindices,ha,metvar,metadata.RSINDICES}); ...


if nargout==1,
    % Optional output variable: the structure with handlers:
    varargout{1}=ha;
end
return;
end % end of main function


%% Here come the callback sub-functions for the GUIs:

% --------------------------------------------------------------
% Sub-function to change the time for the displayed Profile
function update_time0(ho,eventdata,ha,xvar,data)

    if ~ishandle(ho),
        [ho,fo] = gcbo;
        if isnumeric(fo), tt=fo; else tt=fo.Number; end;
    end
    k = round(get(ho,'Value'));
    % changing the title:
    set(ha.tdate,'String',['Date T0: ' datestr(xvar(k),'dd.mm.yyyy @ HHZ')]);

    % changing the profile data
    j = unique(min(max([k-2:k+2],1),length(data)));
    ispan = [1:length(j)];
    fieldvar = get(ha.hh{1},'Tag');
    arrayfun(@(i) {set(ha.hh{i},...
                       'XData',getfield(data(j(i)),fieldvar),...
                       'YData',1e-3*data(j(i)).HGHT);},...
             ispan,'UniformOutput',false);

    % changing vertical line indicating time
    set(ha.htl,'XData',[1 1]*xvar(k));
    return;
end
   
% -------------------------------------------------------------
% Sub-function to change the variable for the displayed Profile
function update_profile(h,eventdata,ha,data,PROFILER,meta)
    [ho,fo] = gcbo;
    k=get(ho,'Value');
    if isnumeric(fo), tt=fo; else tt=fo.Number; end;
    fnames = get(ho,'String');
    Nobs = length(data);
    nspan = fix(length(ha.hh)/2);
    j=min(Nobs,max(1,[-nspan:nspan]+round(get(ha.SliderH,'Value'))));
    % changing profile data:
    arrayfun(@(i) {set(ha.hh{i},'XData',...
                                getfield(data(j(i)),fnames{k}),...
                                'Tag',fnames{k});},...
             [1:length(j)],'UniformOutput',false);

    % updating the x-axis legend:
    %idx = arrayfun(@(i) strcmp(meta{i,1},fnames{k}),...
    %               [1:length(meta)],'UniformOutput',true);
    idx = strcmp(fnames{k},meta);
    set(get(ha.ax(1),'xlabel'),...
        'String',sprintf('%s [%s]',meta{idx,2:3}));

    axis tight;
    % updating the 2D profiles:
    tmp = PROFILER(:,:,k);
    set(ha.p2d,'cdata',tmp);
    tmp = quantile(tmp(:),[0.05 0.93]);
    set(ha.ax(2),'CLim',tmp);

    return;
end

% ----------------------------------------------------------------
% Sub-function to change the variable for the displayed RS indices
function update_rsindices(h,eventdata,ha,metvar,indimeta)
    [ho,fo] = gcbo;
    k = get(ho,'Value');
    fname = get(ho,'String');
    if strcmp(fname{k},'PROFILE_META'),
        return;
    end
    yvar = getfield(metvar,fname{k});
    if all(isnan(yvar)), yvar = []; end
    set(ha.hts,'YData',yvar);
    idx = strcmp(fname{k},indimeta);
    set(get(ha.ax(3),'ylabel'),'String',...
                      sprintf('%s [%s]',indimeta{idx,[1 3]}));
    set(ha.tsl,'String',indimeta{idx,2});
    if isempty(yvar),
        [vmin,vmax] = deal(-1,1);
    else
        [vmin,vmax] = deal(0.98*min(yvar),1.02*max(yvar));
    end
    set(ha.ax(3),'YLim',[vmin vmax]);
    
    return;
end

% ---------------------------------------------------------
% Sub-function to load MAT file with structure variables
function [data, metvar, metadata] = update_matfile(varargin);
    
    [matfile,matpath] = uigetfile({'*.mat','*.MAT'},...
                                  'Select MAT file to open',...
                                  'MultiSelect','off');
    load(fullfile(matpath,matfile),'data','metvar','metadata');
    if ~exist('data','var') | ~exist('metvar','var'),
        error('Structures "data" and/or "metvar" not found in MAT file');
    elseif ~exist('metadata','var'),
        % Checking whether MEATDATA has been assigned, 
        % otherwise read from function
        warning('METADATA isn''t assigned, default metadata loaded.');
        metadata = WyoRS_metadata();
        %[ProfileMeta,NameIndices,DescriptionIndices,UnitsIndices] = 
    
    else
        disp([matfile ' loaded ;D']);
    end

    if nargin>0,
        RASOBS_GRAPH_DATA_RAW(data,metvar,metadata);
    end
    
    return;
end

% ----------------------------------------------------------
% Sub-function to pick-up a Time0 directly from axes:
function []=pick_time0(h,eventdata,ha,Tnum,data)
    [x,y] = ginput(1);
    if x<min(Tnum) | x>max(Tnum),
        disp('Selection out of scope!');
        return;
    end
    [tmp,idx] = min(abs(x-Tnum));
    disp([datestr(Tnum(idx),'dd.mm.yyyy HH Z') '..selected!']);
    set(ha.SliderH,'Value',idx);
    update_time0(ha.SliderH,eventdata,ha,Tnum,data);
    return;
end


%% HERE STARTS AUXILIARY FUNCTIONS 

% ---------------------------------------------------------------
% Function to INITIALIZE the time-series Profiler:
function [PROFILER, TimePROF, Tnumobs] = InitializeProfiler(data,OBST,H_top);

names = fieldnames(data);
nvar = length(names);
Nobs = length(OBST);
% Converting the OBST variable to Matlab datenum format:
tmp = cell2mat(textscan(sprintf('20%8.2f',OBST),...
                        '%4f%2f%5f',Nobs,'CollectOutput',1));

Tnumobs = datenum(tmp);

% index of the highest layer below 15km
hmax = arrayfun(@(i) max(find(data(i).HGHT<H_top)),[1:Nobs]);
PROFILER = zeros(Nobs,max(hmax),nvar)*NaN;
for i=1:Nobs,
        for j=1:nvar,
            eval(['PROFILER(i,1:hmax(i),j) = data(i).'...
                  names{j} '(1:hmax(i));']);
        end
end
TimePROF = repmat(Tnumobs,1,max(hmax));

return;
end


% --------------------------------------------------------------
% Function to make the profiles derivatives respecto to Z
function data = Calculate_Derivatives(data0)
    
% new fields to include:
    d1dz_names = {'D1TE','D1TH','D1TV'};
    d2dz_names = {'D2TE','D2TH','D2TV'};
    in_names = {'TEMP','THTA','THTV'};
    
    Nobs = length(data0);
    
    tmp1 = arrayfun(@(i) cellfun(@(x) diff(getfield(data0(i),x))./diff(data0(i).HGHT),...
                                 in_names,'UniformOutput',false),...
                    [1:Nobs],'UniformOutput',false);
    
    tmp2 = arrayfun(@(i) cellfun(@(x) diff(getfield(data0(i),x),2)./diff(data0(i).HGHT,2),...
                                 in_names,'UniformOutput',false),...
                    [1:Nobs],'UniformOutput',false);


    %%qq=cell2struct(tmp1{1},d1dz_names,2);
    
    for i=1:Nobs,
        for j=1:3,
            data0(i).(d1dz_names{j}) = [tmp1{i}{j}; NaN];
        end
    end
    data = data0;
    return;
end
