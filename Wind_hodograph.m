% *****************************************************************
% Function to plot Wind Hodograph 
% USAGE:
% > h = Wind_hodograph(WS,WD);
% 
% WHERE:
% * WS Wind Speed array [m/s]
% * WD Wind Direction [deg]
% * unit: (optional) alternative unit for WS
% * h: the plot handler
%
% (c) 2019, Pablo Saavedra Garfias
% pablo.saa@uib.no
% Geophysical Institute, University of Bergen
% SEE LICENSE.TXT
% *****************************************************************
%function h = Wind_hodograph(WS, WD, varargin)

% definition of default units
ws_unit = 'm/s';
wh_unit = 'km';

% definition of labels
rosetext = {'N','NE','E','SE','S','SW','W','NW'};

load('../../data/RASOBS/enzv/2008/RS_Y2008-2008_M01-12_D01-31_H00-18.mat');

i = 500;

Htop = 15.1;
H  = 1e-3*(data(i).HGHT);
% considering by default 10km
imax = find(H<Htop);

WS = 0.51444*data(i).SKNT(imax);
WD = deg2rad(data(i).DRCT(imax));
H = H(imax);

% converting to vector components:
Vy = WS.*cos(WD);
Ux = WS.*sin(WD);

U0 = zeros(size(Ux));
V0 = zeros(size(Vy));

% creating canvas:
ax = gca;
GridColor = [.5 .5 .5];
GridLine = ':';
TextSize = 12;

Rmax = max(WS);
%if Rmax<20,
%    R = unique(10*ceil([1:Rmax]/10));
%else
    R = [1:Rmax];%unique(10*floor([1:Rmax]/10));
    %end
if length(R)>6,
    R = R(1:2:end);
end
A = linspace(0,2*pi,360);
yy = cos(A')*R;
xx = sin(A')*R;


% creating lines:
ch = plot(xx,yy,'LineStyle',GridLine,'Color',GridColor);
axis equal;
hold on;
[NN, iamin] = min(hist(WD,[0 pi/2  pi 3/2*pi 2*pi]));

% indexes for the wind-rose octants 
iradial = [1, 45, 90, 135, 180, 225, 270, 315];

% plotting polar grid:
rh = plot(xx(iradial,:)', yy(iradial,:)','LineStyle',GridLine,'Color',GridColor);
txtlabel = ceil(rad2deg(A(iradial)));
Alabel = text(1.1*xx(iradial,end),1.1*yy(iradial,end),rosetext,'FontSize',TextSize);
Rlabel = text(1.01*xx(iradial(2*iamin),[1:end-1]),...
              1.1*yy(iradial(2*iamin),[1:end-1]),...
              num2cell(R(1:end-1)),'FontSize',TextSize);

set(ax, 'Color','none','Box','off','FontSize',13,...
        'XColor','none','YColor','none');

%% data products:
% showing vectors:
h = quiver(U0, V0, Ux, Vy, 0, '-','Color',[1 .7 .7]);

% showing alitude color-code line
sh = surface([Ux Ux],[Vy Vy],[U0 V0],[H H],'LineWidth',2,'EdgeColor','interp');

% writing selected altitudes:
[tmp idxh] = arrayfun(@(x) min(abs(x-H)), linspace(.2,max(H),10));
Hlabel= cellfun(@(x) sprintf('%2.1f',x),num2cell(H(idxh)),'UniformOutput',0);
Hdata = text(Ux(idxh), Vy(idxh), Hlabel);
%view([90 -90]);

hbar = colorbar('westoutside'); 
set(hbar,'Position',get(hbar,'Position').*[1.05 1.1 .6 .8],'FontSize',TextSize);
title(hbar,'km','FontSize',TextSize);
caxis([0 min(15,max(H))]);
colormap('winter');
%end
% end of function