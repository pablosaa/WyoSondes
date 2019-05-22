% -------------------------------------------------------------------------------------
% Script to homogenize the Radiosonde Profiles downloaded from the Wyoming University.
% The objective is to have profiles at the same levels and with the same quality to
% apply radiative transfer simulations to them.
% The output is a plain ASCII file as the usual input for the
% RT3/RT4 model. For improved input version on RT3/RT4, a NetCDF is
% created too. The AIM is to mimic a NetCDF file alike the output of
% the WRF model therefore integration can be done seamlessly. 
%
% (c) 2018 Pablo Saavedra Garfias, Geophysical Institute, University of Bergen
% See LICENSE.TXT
%
% --------------------------------------------------------------------------------------

clear all;
close all;

% boolean parameters
TXTF = false;    % set TRUE if a ASCII file as output is wanted
                 % (old version)


%% Parameter definitions:
topH = 10100;   % top profile height to consider [m]
MinLev = 40;   % minimum number of continues layers the RadioSonde must have

% Default layer altitudes (if not specified as input)
H0 = [[10:20:1500],...
      [1550:50:3000],...
      [3200:200:topH]];    % layer to homoginize [m]


NH = 145;
H0 = 10*unique(floor(logspace(1,4,NH)/10));
tmp = smooth(diff(H0));
H0 = cumsum(tmp);   % standard Atmosphere altitude [m]

%% Definition of Pressure layers where to homogenize RS profiles

P0 = 1013; % Surface Reference Pressure [hPa]
M = 28.8; % Molar mass [gr/mol]
g = 9.807; % gravity acceleration [m/s^2]
R = 8.31446; % Molar gas constant [J/mol/K]
T0 = 290;   % Initial Ambient Temperature  [K]
% Lr = 6.5;   % Temperature lapse-rate [K/km]
%Pa = [[P0:-1.6:710],[680:-30:560],[510:-50:300]];
Pa = P0*exp(-M*g/R/T0*H0*1e-3);
%
%H0 = 1e3*R*T0/M/g*log(P0./Pa);  

nlay = length(H0); % Number of layers for the homogenized layers

%% Definition of Variable names and attributes for NetCDF file:
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

%************************************************************
% Setting up the files to run:
% INPATH = '~/GFI/data/RASOBS/norderney/2015/';
OUTPATH = '/home/pga082/GFI/data/RT/';


%% The cell string 'STATIONS' must be organized according to
% the (Longitude,Latitude) grid cell desired into the NetCDF file.
% e.g. When only one file is specified, then the grid as only one point (1,1).
STATIONS = {'enzv'}; %{'polargmo';'enas';'enbj';'enan'}; %{'norderney'}; %;
[mxN_x, mxN_y] = size(STATIONS);

origin_str = '';   % string containing information about the origin RS
                   
SimTag = '_01COTUR';  % any simulation tag (short string max 8 char) 'fino1';
                   
%% The Database is orginized following the WRF-grid as close as possible:
% * (x_n, y_n): are the coordinates of the specified stations, when only
% one station is introduced then the dimension xn and yn is 1.
% * time: is the time converted to year,month,day,hour and it is
% sorted according to the readed from the directory there the files
% are located: A = dir([fullfile(INPATH, num2str(years)), '*.mat'])

years  = [2008:2013];
N_year = length(unique(years));

% Defining NetCDF and ASCII output file:
outfilen = sprintf('WyoRS_Y%04d-%04d_x%03dy%03dST%s.nc',...
years(1),years(end),...
mxN_x,mxN_y,SimTag(1:min(end,10)));

if exist(fullfile(OUTPATH, outfilen), 'file'),
    disp(['Deleting... NetCDF file: ' OUTPATH outfilen]);
    delete(fullfile(OUTPATH, outfilen) );
end

% Starting to read over number of Stations and years:
for n_x = 1:mxN_x,    % number of stations
    for n_y = 1:mxN_y,    % number of years
        Last_obs = 0;
        for i_year=1:N_year,
            INPATH = ['/home/pga082/GFI/data/RASOBS/'...
                      STATIONS{n_x,n_y} '/' num2str(years(i_year))];
            A = dir(fullfile(INPATH, 'RS_*.mat'));
            N_f = length(A);
            for f = 1:N_f,       % number of files in directory
                filen = fullfile(INPATH, A(f).name);

                % loading the data:
                disp(['Loading... ' filen]);
                clear data station stationname;
                load(filen);
                % checking if variables 'station' and 'stationname' are
                % present to create 'origin_str' for NetCDF output file 
                if ~exist('station','var'),
                    station = '???';
                end
                if ~exist('stationname','var'),
                    stationname = STATIONS{n_x,n_y};
                end
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
                
                % Bit_3: Checking no gaps on layers of Height,
                % Temperature and Relative Humidity
                tmp = arrayfun(@(i) ~any(isnan(data(i).HGHT) & data(i).HGHT<topH),[1:nobs]);
                tmp = and(tmp, arrayfun(@(i) ~any(isnan(data(i).TEMP) & data(i).HGHT<topH),[1:nobs]));
                tmp = and(tmp, arrayfun(@(i) ~any(isnan(data(i).RELH) & data(i).HGHT<topH),[1:nobs]));
                QCflag(tmp) = bitset(QCflag(tmp),3);

                % Bit_4: Checking the consistent decrease in Pressure
                tmp = arrayfun(@(i) ~any(diff(data(i).PRES)>0),[1:nobs]);
                QCflag(tmp) = bitset(QCflag(tmp),4);

                %% QC = [QC;QCflag];

                % ---------------------------------------------------------------
                % Open ASCII file to store the data for the RT code:
                if ~exist('fp','var') & TXTF,
                    % Saving the Header for ASCII file (only once at the
                    % beggining)
                    [tmppath,tmpfile,ext] = fileparts(filen);

                    fp = fopen(fullfile(OUTPATH,[outfilen(1:end-3) '.dat']),'w');
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
        
                idxobs=0;      % index for only quality passed profiles
                               
                % Homogenazing the layes for all profiles:
                for i=1:nobs,
                    idxobs = idxobs + 1;
                    if QCflag(i) ~= 15,  % 15 corresponds to all Quality test passed
                                         %continue; Including even
                                         %bad profiles!
                        QC(n_x,n_y,idxobs) = QCflag(i);
                    else
                        QC(n_x,n_y,idxobs) = QCflag(i);
                    end

                    % Extracting the date for the specific Observation:
                    yy = fix(metvar.OBST(i)/1e4);
                    month = fix(metvar.OBST(i)/1e2) - yy*1e2;
                    day   = fix(metvar.OBST(i)) - yy*1e4 - month*1e2;
                    hour = 24*(metvar.OBST(i)-fix(metvar.OBST(i)));
                    ndate(idxobs,:) = [2000+yy, month, day, hour];

                    % Converting Temperature units from C to Kelvin:
                    data(i).TEMP = data(i).TEMP+273.15;

                    % getting the Surface Data (extrapolating to Height 0):
                    [Zunq, Iunq] = unique(data(i).HGHT);
                    for k=1:length(surface_names),
                        eval(['tmp = find(Zunq<=topH & ~isnan(data(i).'...
                              surface_names{k} '(Iunq)));']);
                        
                        % In case no profile is available in the
                        % lower Atmospheric layers:
                        if isempty(tmp),
                            SURFVars(idxobs,k) = NaN;
                            continue;
                        end
                        % For relative humidity, the extrapolation
                        % is done via nearest value, since it has
                        % been found that linear extrapolation can
                        % produce negative of higher than 100%
                        % values:
                        if strcmp(surface_names{k},'RELH'),
                            interpmethod = 'nearest';
                        else
                            interpmethod = 'linear';
                        end
                                
                        eval(['SURFVars(idxobs,k) =interp1(Zunq(tmp),data(i).'...
                        surface_names{k} '(Iunq(tmp)),0,interpmethod,''extrap'');']);
                        
                    end

                    % Interpolating the profile variables to fixed layers below topH:
                    tmp = Iunq(find(Zunq <= topH)); %find(data(i).HGHT<=topH);
                    
                    h_tmp = data(i).HGHT(tmp);
                    Nh_tmp = length(h_tmp);
                    hi = structfun(@(x) h_tmp(tmp(~isnan(x(tmp)))),data(i),'UniformOutput',0);
                    vi = structfun(@(x) x(tmp(~isnan(x(tmp)))),data(i),'UniformOutput',0);
                    tt = structfun(@(l) length(l), vi);
                    names = fieldnames(vi);
                    for k=1:length(names),
                        Xin = getfield(hi,names{k});
                        Yin = getfield(vi,names{k});
                        % in case of Relative Humidity, the surface
                        % values is included for the interpolation
                        % to avoid high values or negarive values:
                        if strcmp(names{k},'RELH'),
                            Xin = [0; Xin];
                            Yin = [SURFVars(idxobs,3); Yin];
                        end
                        if tt(k)>fix(Nh_tmp/2),
                            Vinter = interp1(Xin,Yin,H0,'linear','extrap');
                        else
                            Vinter = NaN*ones(nlay,1);
                            for badin=1:tt(k),
                                [dummy, Idxin] = min(abs(Xin(badin)-H0));
                                Vinter(Idxin) = Yin(badin);
                            end
                       end
                       eval(['TEMPVars.' names{k} '=Vinter;']);
                       
                    end

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
                    clear TEMPVars;
                    % --------------------------------------------------------------------
                    idxcol = [2,1,3,5,12,13,14,15,16];   % index of variables to store

                    % Saving data into a ASCII file (TODO: mat
                    % file?)
                    if TXTF,
                        fprintf(fp,'%2d %2d\n',day,hour); %n_x,1);
                        fprintf(fp,'%6.3f %10.3f %10.3f %10.3f\n',...
                                metvar.SELV(i)/1e3,SURFVars(idxobs,:));
            
                        fprintf(fp,['%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f ' ...
                                    '%10.3f %10.3f\n'],transpose(ALLVARS(:,idxcol,idxobs)));
                        %cellfun(@(k)
                        %fprintf(fp,'%10.3f',transpose(k(i,:))),layers_var);
                    end
                    
                    % --
                    % end over i: number of observations per File 
                end
                
                if TXTF,
                    % Writing the first line in ASCII file
                    frewind(fp);
                    fprintf(fp,'%02d %02d %02d %3d %3d %3d 2.5 2.5\n',month,day,hour,n_x,n_y,nlay);
                    % Closing ASCII file with "grid" data
                    fclose(fp);
                    clear fp;
                end
                
                %% Writing the whole data in a NetCDF file
                ncfile =  fullfile(OUTPATH, outfilen);
            
                % Writing 4-D variables in NetCDF file:
                for j=1:length(idxcol),
                    if n_x==1 & n_y==1 & i_year==1,
                        nccreate(ncfile,wrf_varname{j},'Datatype','single',...
                                 'Dimensions',{'xn',mxN_x, 'yn', mxN_y,'lev',nlay,'time',Inf}, ...
                                 'Format','netcdf4','Deflate',true,'DeflateLevel',9,'FillValue',NaN);
                    end
                    ncwrite(ncfile,wrf_varname{j},...
                            shiftdim(squeeze(ALLVARS(:,idxcol(j),:)),-2),...
                            [n_x n_y 1 1+Last_obs]);
                    ncwriteatt(ncfile,wrf_varname{j},'short_name',nc_metadata{j,1});
                    ncwriteatt(ncfile,wrf_varname{j},'long_name',nc_metadata{j,2});
                    ncwriteatt(ncfile,wrf_varname{j},'units',nc_metadata{j,3});
                end
            
                % Writting 3-D variables in NetCDF file:
                for j=1:length(surface_names),
                    if n_x==1 & n_y==1 & i_year==1,
                        nccreate(ncfile,wrf_surfname{j,1},'Datatype','single',...
                                 'Dimensions',{'xn',mxN_x,'yn',mxN_y,'time',Inf},'Format', ...
                                 'netcdf4','Deflate',true,'DeflateLevel',9);
                    end
                    ncwrite(ncfile,wrf_surfname{j,1},shiftdim(SURFVars(:,j),-2),...
                            [n_x, n_y, 1+Last_obs]);
                    ncwriteatt(ncfile,wrf_surfname{j,1},'short_name',wrf_surfname{j,1});
                    ncwriteatt(ncfile,wrf_surfname{j,1},'long_name',wrf_surfname{j,2});
                    ncwriteatt(ncfile,wrf_surfname{j,1},'units',wrf_surfname{j,3});
                end

                if n_x==1 & n_y==1 & i_year==1,
                    nccreate(ncfile,'QIDX','Datatype','int8', ...
                             'Dimensions',{'xn',mxN_x,'yn',mxN_y,'time',Inf},...
                             'FillValue',NaN);
                end
                ncwrite(ncfile,'QIDX',shiftdim(QCflag,-2),[n_x, n_y, 1+Last_obs]);
                ncwriteatt(ncfile,'QIDX','long_name','Profile Quality index');
                ncwriteatt(ncfile,'QIDX','note:','15=good quality');
                
                % Writting 2-D variables in NetCDF file:
                for j=1:length(wrf_locname),
                    if n_x==1 & n_y==1 & i_year==1,
                        nccreate(ncfile,wrf_locname{j,1},'Datatype','single',...
                                 'Dimensions',{'xn',mxN_x,'yn',mxN_y});
                    end
                    ncwrite(ncfile,wrf_locname{j,1},LOCVAR(j),[n_x, n_y]);
                    ncwriteatt(ncfile,wrf_locname{j,1},'short_name',wrf_locname{j,1});
                    ncwriteatt(ncfile,wrf_locname{j,1},'long_name',wrf_locname{j,2});
                    ncwriteatt(ncfile,wrf_locname{j,1},'units',wrf_locname{j,3});
                end
            
                % Writting 1-D variables in NetCDF file:
                for j=1:4,
                    if n_x==1 & n_y==1 & i_year==1,
                        nccreate(ncfile,date_name{j},'Datatype','single',...                        
                        'Dimensions',{'xn',mxN_x,'yn',mxN_y,'time',Inf});
                    end
                    ncwrite(ncfile,date_name{j},shiftdim(ndate(:,j),-2),...
                    [n_x, n_y, 1+Last_obs]);
                end
            
                Last_obs = Last_obs + idxobs;
                clear ALLVARS SURFVars ndate;
                
                %--
                % end over number of files (index f)
            end
            
            % --
            % end over years (index i_year)
        end 
            
        % Writting global variables:
        origin_str = [origin_str, sprintf('[%03d,%03d]->%s(%s);',n_x,n_y,stationname,station)];
        ncwriteatt(ncfile,'/','grid_x',NaN);
        ncwriteatt(ncfile,'/','grid_y',NaN);
        ncwriteatt(ncfile,'/','origin',origin_str);
        ncwriteatt(ncfile,'/','Creation',datestr(today));
        ncwriteatt(ncfile,'/','Contact','Pablo.Saavedra@uib.no');
        ncwriteatt(ncfile,'/','Institution',['Geophysical Institute, ' ...
                            'Uni-Bergen']);
        % --
        % end over y station (index n_y)
    end 
    
    % --
    % end over x station (index n_x)
end 
return;

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

% *************************************************************************
% End of script
