% *****************************************************************
% Function to plot Wind Hodograph 
% USAGE:
% To generate a new Hodograph:
% > h = Wind_hodograph(WS, WD, H);
% > h = Wind_hodograph(WS, WD, H, 'axes',gca,...
%                      'vectors',0,'colored',1,'heights',0);
%
% To update existing Hodograph with new data:
% > Wind_hodograph(WS, WD, H, 'update',h);
%
% WHERE:
% INPUTS
% * WS Array with Wind Speed [m/s],
% * WD Array with Wind Direction [deg],
% * H  Array with Altitudes [km],
% OUTPUTS
% * h: structure with handler to plots and axis
%
% PAIR-OPTIONS:
% * 'axes': (gca default) axes handle where to plot,
% * 'update': (output handle) update handle h with new data,
% * 'vectors': (true default) plot WS vector field,
% * 'colored': (true default) wind shear line color-coded,
% * 'heights': (false default) show height km in windshear,
%
% ---
% (c) 2019, Pablo Saavedra Garfias
% pablo.saa@uib.no
% Geophysical Institute, University of Bergen
% SEE LICENSE.TXT
% *****************************************************************
function varargout = Wind_hodograph(WS, WD, H, varargin)
    
    if nargin<3,
        error('Three input arguments needed!');
    end
% define default parameters:
    HFLAG = true;
    CFLAG = true;
    VFLAG = true;
    % Checkign for input arguments:
    for i=1:2:length(varargin),
        VAL = varargin{i+1};
        switch varargin{i},
          case 'axes',
            if ishandle(VAL),
                ax = VAL;
            else
                error('"axes" needs to be a axes handle!');
            end
          case 'update',
            if isstruct(VAL),
                hodog = VAL;
                CFLAG = hodog.CFLAG;
                VFLAG = hodog.VFLAG;
                HFLAG = hodog.HFLAG;
            else
                error('Handler argument seems corrupted!');
            end        
          case 'vectors',
            VFLAG = logical(VAL);
          case 'colored',
            CFLAG = logical(VAL);
          case 'heights',
            HFLAG = logical(VAL);
          otherwise,
            error([varargin{i} ' Option not valid!']);
        end
    end
    

% definition of default units
ws_unit = 'm/s';
wh_unit = 'km';

% definition of labels
rosetext = {'N','NE','E','SE','S','SW','W','NW'};

%load('../../data/RASOBS/enzv/2008/RS_Y2008-2008_M01-12_D01-31_H00-18.mat');

Htop = 15.1;
H  = 1e-3*(H); %data(i).HGHT);
% considering by default 10km
imax = find(H<Htop);

WS = 0.51444*WS(imax); %data(i).SKNT(imax);
WD = deg2rad(WD(imax)); %data(i).DRCT(imax));
H = H(imax);

% converting to vector components:
Vy = WS.*cos(WD);
Ux = WS.*sin(WD);

U0 = zeros(size(Ux));
V0 = zeros(size(Vy));

Rmax = max(WS);

R = [0:5:10*floor(Rmax/10)];

if length(R)>6,
    R = R(1:2:end);
end
A = linspace(0,2*pi,360);
yy = cos(A')*R;
xx = sin(A')*R;

% Drawing in canvas
if exist('hodog','var'),
    % extracting handles to update data:
    hodog.VFLAG = VFLAG;
    hodog.CFLAG = CFLAG;
    hodog.HFLAG = HFLAG;
    resampling_hodograph(Ux, Vy, H, hodog);
    return;
else
    ax = gca;
end

% Creating canvas for Polar coordinates:
GridColor = [.5 .5 .5];
GridLine = ':';
TextColor = [.3 .3 .3];
TextSize = 12;

% creating concentric grid circles:
ch = plot(xx,yy,'LineStyle',GridLine,'Color',GridColor);
axis equal;
hold on;
[NN, iamin] = min(hist(WD, [0, pi/4, pi/2, 3/4*pi, pi,...
                    5/4*pi, 3/2*pi, 3/2*pi, 2*pi]));

% indexes for the wind-rose octants 
iradial = [1, 45, 90, 135, 180, 225, 270, 315];

% creating radial grid lines:
rh = plot(xx(iradial,:)', yy(iradial,:)','LineStyle',GridLine,'Color',GridColor);
txtlabel = ceil(rad2deg(A(iradial)));
% Angular labels:
Alabel = text(1.1*xx(iradial,end),1.1*yy(iradial,end),rosetext,...
              'FontSize',TextSize,'Color',TextColor);
% Radial labels:
Rlabel = text(1.01*xx(iradial(iamin),[1:end-1]),...
              1.1*yy(iradial(iamin),[1:end-1]),...
              num2cell(R(1:end-1)),...
              'FontSize',TextSize,'Color',TextColor);

set(ax, 'Color','none','Box','off','FontSize',13,...
        'XColor','none','YColor','none',...
        'XLimMode','manual','YLimMode','manual');

%% data products:
% Drawing vectors:
vh = quiver(U0, V0, Ux, Vy, 0, '-','Color',[.7 .7 1]);
set(vh, 'Visible',VFLAG);

% showing alitude color-code line
sh = surface([Ux Ux],[Vy Vy],[U0 V0],[H H],...
             'LineWidth',2,'EdgeColor','interp');

caxis([0 min(Htop,max(H))]);
if CFLAG,
    cmap = winter;
else
    cmap = [0 0 0];
end
colormap(ax, cmap);
hbar = colorbar('westoutside'); 
set(hbar,'Position',get(hbar,'Position').*[1 1.1 .5 .8],...
         'FontSize',TextSize,'Visible',CFLAG);
title(hbar,'km','FontSize',TextSize,'Visible',CFLAG);

% Writing in canvas selected altitudes:
[tmp idxh] = arrayfun(@(x) min(abs(x-H)), linspace(.2,max(H),10));
Hlabel= cellfun(@(x) sprintf('%2.1f',x),num2cell(H(idxh)),'UniformOutput',0);
htxt  = text(Ux(idxh), Vy(idxh), Hlabel);
set(htxt,'Visible',HFLAG,'Color','r');

% Returning structure handles:
if nargout == 1,
    hodog.ax = ax;
    hodog.ch = ch;
    hodog.rh = rh;
    hodog.vh = vh;
    hodog.sh = sh;
    hodog.hbar = hbar;
    hodog.Hdata = htxt;
    hodog.CFLAG = CFLAG;
    hodog.HFLAG = HFLAG;
    hodog.VFLAG = VFLAG;
    varargout{1} = hodog;
end
return;
end


%% HOW to update data:
% 0) update canvas (change radial scope and labels?)
% 1) calculate Ux Vy U0 V0
% 2) update new vectors:
%    > set(h,'XData',U0,'YData',V0,'UData',Ux,'VData',Vy)
% 3) update altitude line:
%    > set(sh,'XData',[Ux Ux],'YData',[Vy Vy],'Zdata',[U0 V0],'CData',[H H])
% 4) update altitude lavels:
%    > arrayfun(@(i) set(Hdata(i),'String',Hlabel{i},'Position',[Ux(idxh(i)) Vy(idxh(i))]),[1:10])
%
function resampling_hodograph(Ux, Vy, H, h),
% 
    U0 = zeros(size(Ux));
    V0 = zeros(size(Vy));

    % Updating the Vector Fields:
    set(h.vh,'XData',U0,'YData',V0,'UData',Ux,'VData',Vy,'Visible',h.VFLAG);
    
    % Updating the Shear line:
    set(h.sh,'XData',[Ux Ux],'YData',[Vy Vy],...
             'Zdata',[U0 V0],'CData',[H  H]);

    if h.CFLAG,
        cmap = 'winter';
    else
        cmap = [0 0 0];
    end
    set(h.hbar,'Visible',h.CFLAG);
    title(h.hbar,'km','Visible',h.CFLAG);
    colormap(h.ax,cmap);
    
    % Updating the Altitude labels along the line
    [tmp idxh] = arrayfun(@(x) min(abs(x-H)), linspace(.2,max(H),10));
    Hlabel= cellfun(@(x) sprintf('%2.1f',x),num2cell(H(idxh)),'UniformOutput',0);
        
    arrayfun(@(i) set(h.Hdata(i),'String',Hlabel{i},...
                                 'Visible', h.HFLAG,...
                                 'Position',[Ux(idxh(i)) Vy(idxh(i))]),[1:10])

end

% end of function