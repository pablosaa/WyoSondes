% Function to plot the Skew-T_logP diagram

clear all;
close all;

if exist('OCTAVE_VERSION','builtin'),
  graphics_toolkit('qt');
  LiT = '-';
  LiW = '--';
  FS = 17;
  LW = 0.5;
else
  LiT = '-';
  LiW = '--';
  FS = 13;
  LW = 1.5;
endif
a = 6.112;
b = 17.67;
c = 243.5;
k = 0.286;
Cpv = 1.87;  % Head capacity of vapour [kJ/kg/K]
Cp  = 1.005; % Head capacity of air [kJ/kg/K]
hv  = 2502;  % specific entalphy of vapour [kJ/kg]
epsi = 0.622;  % R_dry/R_v
Rp = 0.287;  % Rv*epsilon [kJ/kg/K]
A = 2.53e9; % hPa
B = 5.42e3; % K
P0 = 1000;  % Reference pressure [hPa]
T=[-50:10:40];
P=[1000:-100:100];
THETA = [260:10:350];  % adiabatic
T_K = T + 273.15;
ws = 1e-3*[0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10];  % saturated mixing ration g/kg

% DRY ADIABAT [C]:
T_K = -273.15 + ((P'/P0).^k)*THETA;
%P_THTA = P0*(T_K'*(1./THETA)).^(1/k);

% ISOPLETS:
Td = B./log(A*epsi./(P(1:6)'*ws)) - 273.15;

% MOIST ADIABATIC
ep = @(T) a*exp(b*T./(T+c)); %6.112*exp(17.67*T./(T+243.5));

%ep = epsi*P'*(ws./(ws+epsi));
dws = @(T,P) epsi*b*c*P'.*ep(T)./((P'-ep(T)).*(T+c)).^2;
L = fliplr(ws)'.*(hv + Cpv*T);  % Latent Head [kJ/kg]
T_wet = T_K./(1+L.*dws(T,P)/Cp);
wsat = @(T,P) ep(T)./(P-ep(T));
%T_wet = T_K.*(1+L.*wsat./(T+273.15)/Rp)./(1+(L.^2).*epsi*wsat/Rp/Cp./(T+273.15).^2);

% Plotting parameters:
skew_angle = 45;
slope = @(P) -skew_angle*(P'-P0)/P0;
%m = tand(skew_angle);

% Starting Figure drawing
figure;

ha.isobar = plot(T,repmat(P,10,1), 'k', 'LineWidth', LW, 'LineStyle', LiT);
hold on;
ha.isoter = plot(repmat(T,10,1) + slope(P), P,'r','LineWidth', LW, 'LineStyle', LiT);

ha.isoadi = plot(T_K + slope(P), P, 'c', 'LineWidth', LW, 'LineStyle', LiT);
ha.isowet = plot(T_wet + slope(P), P, 'y', 'LineWidth', LW, 'LineStyle', LiW, 'Color',[0.4 .7 .3]);

ha.isoplet= plot(Td + slope(P(1:6)), P(1:6), 'g', 'LineWidth', LW, 'LineStyle', LiW);
tmptxt = cellfun (@num2str,num2cell(1e3*ws), 'UniformOutput',false);
ha.txtplet= text(1+Td(6,:)+ slope(P(6)), repmat(P(6)-20,1,length(ws)),...
         tmptxt,'Color','g','FontSize',11);

tmptxt = cellfun(@num2str, num2cell(THETA), 'UniformOutput',false);
ha.txtadi = text(T_K(5,:)+slope(P(5))-1.5, repmat(P(5)+30,1,length(THETA)),...
tmptxt,'Color','c','FontSize',11,'Rotation', 30);
% Default Axis configuration
set(gca, 'YDir', 'reverse', 'YLim',[100 1020], 'XLim', [-41 T(end)], ...
'Box','on','XGrid','on','XColor','r','FontSize',FS,'XTick',[-40:10:40],'YTick',[100:100:1000]);
xlabel('Temperature [C]', 'FontSize', FS);
ylabel('Pressure [hPa]', 'FontSize', FS);



%% --- end of main function
