% -------------------------------------------------------------------------------------
% Script to homogenize the Radiosonde Profiles downloaded from the Wyoming University.
% The objective is to have profiles at the same levels and with the same quality to
% apply radiative transfer simulations to them.
%
% (c) 2018 Pablo Saavedra Garifas, Geophysical Institute, University of Bergen
% See LICENSE.TXT
%
% --------------------------------------------------------------------------------------

clear all;
close all;

% Parameter definitions:
topH = 9000;   % max height to consider [m]
MinLev = 40;   % minimum number of continues layers the RS must have

% Default layer altitudes (if not specified as input)
H0 = [[10:20:1500],...
      [1550:50:3000],...
      [3200:200:topH]];    % layer to homoginize [m]

% P0 = 1024; % Reference Pressure [kPa]
% M = 29; % Molar mass [gr/mol]
% g = 9.807; % gravity acceleration [m/s^2]
% R = 8.31446; % Molar gas constant [J/mol/K]
% T0 = 280;   % Initial Ambient Temperature  [K]
% Lr = 6.5;   % Temperature lapse-rate [K/km]
%
% H0 = R*T0/M/g*log(P0./Pa)./(1+R*Lr/M/g*log(P0./P));

nlay = length(H0); % Number of layers for the homoginized grid
% Variables needs to be included in the ASCII file
surface_names = {'PRES','TEMP','RELH'};
layers_var = {'HGHT','PRES','TEMP','RELH','CLOUD','RAIN'};

% metadata for the NetCDF file (in the storing order):               
nc_metadata = {
    'HGHT',    'Geopotential Height','km';
    'PRES',    'Atmospheric Pressure','hPa';
    'TEMP',    'Temperature','K';
    'RELH',    'Relative Humidity','%';
    'CLOUD',   'Cloud Water Content','g/m^3';
    'RAIN','Rain Water Content','g/m^3';
    'ICE',  'Ice Content','g/m^3';
    'SNOW', 'Snow Content','g/m^3';
    'GRAUPEL','Graupel Content','g/m^3'
};
QC = [];   % quality factor to storage in mat-file
y = 0;

% Setting up the files to run:
%% INPATH = '~/GFI/data/RASOBS/norderney/2015/';
INPATH = '/home/pga082/GFI/data/RASOBS/polargmo/';
years = [2013:2013];

for x = 1:length(years),
    A = dir(fullfile(INPATH, num2str(years(x)), '*.mat'));
    n_y = length(A);
    for f = 1:n_y,
        filen = fullfile(INPATH,num2str(years(x)), A(f).name);
        %filen = 'RS_Y2015-2015_M09-11_D01-31_H00-12.mat';

        % loading the data:
        load(filen);
        names = regexprep(fieldnames(data),'HGHT','CLOUD');
        nobs  = length(data);

        QCflag = int8(zeros(nobs,1));   % Quality Control flag 8 bits:
                                        % bit 1: height reaches at-least 10km
                                        % bit 2: height increases monotonilcally
                                        % bit 3: no gaps in leyers
                                        % bit 4: pressure always decrease with height increase

        % Bit_1: Checking minimum height (e.g. 10km)
        tmp = arrayfun(@(i) max(data(i).HGHT)>topH & length(data(i).HGHT<topH)>=MinLev,[1:nobs]);
        QCflag(tmp) = bitset(QCflag(tmp),1);

        % Bit_2: Checking that height increses monotonically
        tmp = arrayfun(@(i) ~any(diff(data(i).HGHT)<=0),[1:nobs]);
        QCflag(tmp) = bitset(QCflag(tmp),2);

        % Bit_3: Checking no gaps on layers
        tmp = arrayfun(@(i) ~any(isnan(data(i).HGHT) & data(i).HGHT<topH),[1:nobs]);
        QCflag(tmp) = bitset(QCflag(tmp),3);

        % Bit_4: Checking the consistent decrease in Pressure
        tmp = arrayfun(@(i) ~any(diff(data(i).PRES)>0),[1:nobs]);
        QCflag(tmp) = bitset(QCflag(tmp),4);

        QC = [QC;QCflag];
        % ---------------------------------------------------------------
        % Open ASCII file to store the data for the RT code:

        if ~exist('fp','var'),
            % Saving the Header for ASCII file (only once at the
            % beggining)
            [tmppath,tmpfile,ext] = fileparts(filen);
            fp = fopen(fullfile(INPATH,[tmpfile '.dat']),'w');
            n_x = 0;
            month = 99;
            day = 99;
            hour = 99;
            %n_y = 0;
            % following is a fake line which will be updated at the end
            fprintf(fp,'%2d %2d %2d %3d %3d %3d 2.5 2.5\n',month,day,hour,n_x,n_y,nlay);
        end

        % Homogenazing the layes for all profiles:
        idxobs=0;
        for i=1:nobs,
            if QCflag(i) ~= 15,  % 15 corresponds to all Quality test passed
                continue;
            end
            idxobs = idxobs + 1;
            % Extracting the date for the specific Observation:
            yy = fix(metvar.OBST(i)/1e4);
            month = fix(metvar.OBST(i)/1e2) - yy*1e2;
            day   = fix(metvar.OBST(i)) - yy*1e4 - month*1e2;
            hour = 24*(metvar.OBST(i)-fix(metvar.OBST(i)));

            n_x = n_x + 1;
            % Converting Temperature units from C to Kelvin:
            data(i).TEMP = data(i).TEMP+273.15;

            % getting the Surface Data (extrapolating to Height 0):
            for k=1:length(surface_names),
                eval(['tmp = find(data(i).HGHT<=topH & ~isnan(data(i).'...
                      surface_names{k} '));']);
                eval(['SURFVars(idxobs,k) = interp1(data(i).HGHT(tmp),data(i).'...
                      surface_names{k} '(tmp),0,''linear'',''extrap'');']);
            end

            % Interpolating the profile variables to fixed layers below topH:
            tmp = find(data(i).HGHT<=topH);
            TEMPVars = structfun(@(x) (interp1(data(i).HGHT(tmp & ~isnan(x(tmp))),...
                                               x(tmp & ~isnan(x(tmp))),H0,'linear','extrap')),...
                                 data(i), 'UniformOutput',false);

            % Converting the Height units from m to km and from
            % station level to a.s.l.:
            TEMPVars.HGHT = (TEMPVars.HGHT+metvar.SELV(i))/1e3;

            % --------------------------------------------------------
            [TEMPVars.CLOUD, TEMPVars.RAIN, TEMPVars.IWC,...
             TEMPVars.SNOW, TEMPVars.GRAUPEL] = cloud_modell(TEMPVars.TEMP,...
                                                             TEMPVars.RELH,TEMPVars.PRES);

            names = fieldnames(TEMPVars);
            for j=1:length(names),
                eval(['ALLVARS(:,j,idxobs) = TEMPVars.' names{j} ';']);
                %eval([names{j} '(i,:) = TEMPVars.' names{j} ';']);
            end

            % --------------------------------------------------------------------
            % Saving data into a ASCII file (TODO: mat file?)
            fprintf(fp,'%2d %2d\n',day,hour); %n_x,1);
            fprintf(fp,'%6.3f %10.3f %10.3f %10.3f\n',...
                    metvar.SELV(i)/1e3,SURFVars(idxobs,:));
            
            idxcol = [2,1,3,5,12,13,14,15,16];
            fprintf(fp,['%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f ' ...
            '%10.3f %10.3f\n'],transpose(ALLVARS(:,idxcol,idxobs)));
            %cellfun(@(k) fprintf(fp,'%10.3f',transpose(k(i,:))),layers_var);

        end  % end over i: number of observations

            % Writing the first line in ASCII file
            frewind(fp);
            fprintf(fp,'%02d %02d %02d %3d %3d %3d 2.5 2.5\n',month,day,hour,n_x,n_y,nlay);
            % Closing ASCII file with "grid" data
            fclose(fp);
            clear fp;
            % Writing the whole data in a NetCDF file
            ncfile =  fullfile(INPATH,[tmpfile '.nc']);
            if exist(ncfile,'file'), delete(ncfile); end
            for j=1:length(idxcol),
                nccreate(ncfile,names{idxcol(j)},'Datatype','single',...
                         'Dimensions',{'lev',nlay, 'xn',nobs}, ...
                         'Format','classic');
                if j<=length(surface_names),
                    nccreate(ncfile,[surface_names{j} '_surf'],'Datatype','single',...
                             'Dimensions',{'xn',nobs},'Format', ...
                             'classic');
                    ncwrite(ncfile,[surface_names{j} '_surf'],SURFVars(:,j));    
                end
                ncwrite(ncfile,names{idxcol(j)},squeeze(ALLVARS(:,idxcol(j),:)));
                ncwriteatt(ncfile,names{idxcol(j)},'Short Name',nc_metadata{j,1});
                ncwriteatt(ncfile,names{idxcol(j)},'Long Name',nc_metadata{j,2});
                ncwriteatt(ncfile,names{idxcol(j)},'Units',nc_metadata{j,3});
            end
            ncwriteatt(ncfile,'/','Origin',filen);
            ncwriteatt(ncfile,'/','Station a.s.l [km]',metvar.SELV(1)*1e-3);
            ncwriteatt(ncfile,'/','Contact','pablo.saa@uib.no');

    end  % end over number of files (index f)
end  % end over years

    % ============ End of main function RS_homogenize_profiles =================

    % **************************************************************************
    % Cloud model function for estimation of liquid water cloud
    % profile
    % Including Mätzler cloud model (to be reviewed)
    %
    function [LWC,LWR,IWC,SWC,GWC] = cloud_modell(T,RH,P)
    % INPUT:
    % T-> Temperature [K]
    % RH -> Relative Humidity [%]
    % P -> Pressure [hPa]
    % OUTPUT:
    % LWC -> Cloud Liquid Water Content [g/m3]
    % LWR -> Rain Liquid Water Content [g/m3]
    % LWG -> Graupel Liquid Water Content [g/m3]
    % NOTE: This is only valid for Clouds! for Rain, Graupel, Ice
    % and Snow an assignation function
    % still needs to be implemented.

        b0 = 90;   % parameter b0 range from 85 to 90
        tmp = find(RH>=b0 & T>=240);
        % Cloud water content [g/m^3]:
        LWC(1,:) = zeros(1,length(T));
        LWC(1,tmp) = 2*((RH(tmp)-b0)/30).^2;  
        % Rain Water content [g/m^3]:
        LWR(1,:) = zeros(1,length(T));
        % Cloud Ice content [g/m^3]:
        IWC(1,:) = zeros(1,length(T));
        % Show content [g/m^3]:
        SWC(1,:) = zeros(1,length(T));
        % Graupel content [g/m^3]:
        GWC(1,:) = zeros(1,length(T));
    end
    % End of script
