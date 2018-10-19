% -------------------------------------------------------------------------------------
% Script to homogenize the Radiosonde Profiles downloaded from the Wyoming University.
% The objective is to have profiles at the same levels and with the same quality to
% apply radiative transfer simulations to them.
% The output is a plain ASCII file as the usual input for the
% RT3/RT4 model. For improved input version on RT3/RT4, a NetCDF is
% created too. The AIM is to mimic a NetCDF file alike the output of
% the WRF model therefore integration can be done seamlessly. 
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

% Definition of Variable names and attributes for NetCDF file:
% Variables needs to be included in the ASCII file
surface_names = {'PRES','TEMP','RELH'};
layers_var = {'HGHT','PRES','TEMP','RELH','CLOUD','RAIN'};

% metadata for the NetCDF file (in the storing order):
wrf_surfname = {'PSFC','Surface Pressure','hPa';
                'T2','2-m Temperature','K';
                'Q2','Surface Relative Humidity','%'};
wrf_varname = {'PHB','P','T','QVAPOR','QCLOUD','QRAIN','QICE', ...
               'QSNOW','QGRAUP'};

wrf_locname = {'HGT','Station Altitude','km';
              'LAT','Sation Latitude','deg';
              'LON','Station Longitude','deg'};
date_name = {'year','month','day','hour'};
% Radiosonde variables to load:
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


% Indexing for NetCDF files:
n_x = 0; % this is index for different stations
n_y = 0; % this index is for different years

y = 0;

% Setting up the files to run:
%% INPATH = '~/GFI/data/RASOBS/norderney/2015/';
OUTPATH = '/home/pga082/GFI/data/RT/';
STATIONS = {'polargmo','norderney'};
mxN_x = length(STATIONS);
years = [2015,2015];
mxN_y = length(unique(years));
outfilen = sprintf('RS_Y%04d-%04d_4RT',years(1),years(end));
if exist([OUTPATH outfilen '.nc']),
    disp('Deleting... NetCDF file');
    delete([OUTPATH outfilen '.nc']);
end

for n_x=1:mxN_x,
INPATH = ['/home/pga082/GFI/data/RASOBS/' STATIONS{n_x} '/'];
for n_y = 1:mxN_y,
    A = dir(fullfile(INPATH, num2str(years(n_y)), '*.mat'));
    N_f = length(A);
    for f = 1:N_f,
        filen = fullfile(INPATH,num2str(years(n_y)), A(f).name);

        % loading the data:
        disp(['Loading... ' filen]);
        load(filen);
        % names = regexprep(fieldnames(data),'HGHT','CLOUD');
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

            fp = fopen(fullfile(OUTPATH,[outfilen '.dat']),'w');
            month = 99;
            day = 99;
            hour = 99;

            % following is a fake line which will be updated at the end
            fprintf(fp,'%2d %2d %2d %3d %3d %3d 2.5 2.5\n',month,day,hour,n_x,n_y,nlay);
        end

        % Temporal variable for the Station coordinates:
        LOCVAR(1) = metvar.SELV(1)*1e-3;
        LOCVAR(2) = metvar.SLAT(1);
        LOCVAR(3) = metvar.SLON(1);
        
        % Homogenazing the layes for all profiles:
        idxobs=0;      % index for only quality passed profiles
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
            ndate(idxobs,:) = [2000+yy month day hour];

            %n_x = n_x + 1;
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
            TEMPVars.HGHT = (TEMPVars.HGHT)/1e3;  %+metvar.SELV(i)

            % --------------------------------------------------------
            [TEMPVars.CLOUD, TEMPVars.RAIN, TEMPVars.IWC,...
             TEMPVars.SNOW, TEMPVars.GRAUPEL] = cloud_modell(TEMPVars.TEMP,...
                                                             TEMPVars.RELH,TEMPVars.PRES);

            names = fieldnames(TEMPVars);
            for j=1:length(names),
                eval(['ALLVARS(:,j,idxobs) = TEMPVars.' names{j} ';']);
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
            
            %% Writing the whole data in a NetCDF file
            ncfile =  fullfile(OUTPATH,[outfilen '.nc']);
            
            % Writing 4-D variables in NetCDF file:
            for j=1:length(idxcol),
                if n_x==1 & n_y==1,
                    nccreate(ncfile,wrf_varname{j},'Datatype','single',...
                             'Dimensions',{'xn',mxN_x, 'yn', mxN_y,'lev',nlay,'time',Inf}, ...
                             'Format','netcdf4','Deflate',true,'DeflateLevel',9,'FillValue',NaN);
                end
                ncwrite(ncfile,wrf_varname{j},...
                        shiftdim(squeeze(ALLVARS(:,idxcol(j),:)),-2),...
                        [n_x n_y 1 1]);
                ncwriteatt(ncfile,wrf_varname{j},'short_name',nc_metadata{j,1});
                ncwriteatt(ncfile,wrf_varname{j},'long_name',nc_metadata{j,2});
                ncwriteatt(ncfile,wrf_varname{j},'units',nc_metadata{j,3});
            end
            
            % Writting 3-D variables in NetCDF file:
            for j=1:length(surface_names),
                if n_x==1 & n_y==1,
                    nccreate(ncfile,wrf_surfname{j,1},'Datatype','single',...
                             'Dimensions',{'xn',mxN_x,'yn',mxN_y,'time',Inf},'Format', ...
                             'netcdf4','Deflate',true,'DeflateLevel',9);
                end
                ncwrite(ncfile,wrf_surfname{j,1},shiftdim(SURFVars(:,j),-2),...
                        [n_x n_y 1]);
                ncwriteatt(ncfile,wrf_surfname{j,1},'short_name',wrf_surfname{j,1});
                ncwriteatt(ncfile,wrf_surfname{j,1},'long_name',wrf_surfname{j,2});
                ncwriteatt(ncfile,wrf_surfname{j,1},'units',wrf_surfname{j,3});
            end
            
            % Writting 2-D variables in NetCDF file:
            for j=1:length(wrf_locname),
                if n_x==1 & n_y==1,
                    nccreate(ncfile,wrf_locname{j,1},'Datatype','single',...
                             'Dimensions',{'xn',mxN_x,'yn',mxN_y});
                end
                ncwrite(ncfile,wrf_locname{j,1},LOCVAR(j),[n_x,n_y]);
                ncwriteatt(ncfile,wrf_locname{j,1},'short_name',wrf_locname{j,1});
                ncwriteatt(ncfile,wrf_locname{j,1},'long_name',wrf_locname{j,2});
                ncwriteatt(ncfile,wrf_locname{j,1},'units',wrf_locname{j,3});
            end
            
            % Writting 1-D variables in NetCDF file:
            for j=1:4,
                if n_x==1 & n_y==1,
                    nccreate(ncfile,date_name{j},'Datatype','single', ...
                             'Dimensions',{'time',nobs});
                end
                ncwrite(ncfile,date_name{j},ndate(:,j));
            end
            ncwriteatt(ncfile,'/','grid_x',NaN);
            ncwriteatt(ncfile,'/','grid_y',NaN);
            ncwriteatt(ncfile,'/','Origin',filen);
            ncwriteatt(ncfile,'/','Creation',datestr(today));
            ncwriteatt(ncfile,'/','Contact','pablo.saa@uib.no');

    end  % end over number of files (index f)
end  % end over years (index n_y)
end % end over stations (index n_x)

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
