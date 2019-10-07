% -------------------------------------------------------------------------------------
% Function to homogenize altitude levels from the Radiosonde
% Profiles downloaded from the Wyoming University by means of
% the WyoSonde Toolbox.
% 
% The objective is to have profiles at the same levels and with the same quality to
% apply radiative transfer simulations to them.
% The output is a plain ASCII file as the usual input for the
% RT3/RT4 model. For improved input version on RT3/RT4, a NetCDF is
% created too. The AIM is to mimic a NetCDF file alike the output of
% the WRF model therefore integration can be done seamlessly. 
%
% USAGE:
% > status = RS_homogenize_profiles(mat-file);
% > status = RS_homogenize_profiles('station',{'home','work'},'years',[1990:1999],...);
%
% WHERE the input arguments are given in pair as...
% ... mandatory input arguments:
% * 'station' = {'enzv'};  or {'polargmo';'enas';'enbj';'enan'}; or %{'norderney'};
% * 'years'   = [2014:2018]; % range of years for data to use;
%
% ... optional input arguments:
% * H: Array with height levels [Km] to use for homoginization,
% * OUTPATH: string with the path where to store the file,
% * SimTag = '_02COTUR';  % simulation tag (short char) 'fino1';
% * BASEPATH = '/home/username/data/RASOBS';  % the base folder there the RS data should be stores by stations;
% 
% The default values for the optional input pair arguments:
% * H = [10, 20, 30... 10000];
% * OUTPATH = '~/GFI/data/RT';
% * SimTag = '';
% * BASEPATH = '~/GFI/data/RASOBS'
%
% Given the above input arguments, the radiosonde MAT files will be searched at the following path:
% BASEPATH/[station]/[years]/*.mat
% 
% (c) 2018 Pablo Saavedra Garfias, Geophysical Institute, University of Bergen
% See LICENSE
% --------------------------------------------------------------------------------------
function [] = RS_homogenize_profiles(varargin)

	% HOME base path:
	HOME_BASE = getenv('HOME');

	% checking input arguments:
  switch(nargin)
    case 0,
      warning('Input arguments needed!');
      return;
    case 1,
      if ischar(varargin{1}), % & exist(varargin{1},'file'),
        matfile = varargin{1};            
      else
        error('When one input argument, it must be a string .MAT file name!');
      end
    otherwise,
      for i=1:2:nargin,
        MYARG = varargin{i};
        MYVAL = varargin{i+1};
        if ~ischar(MYARG),
          error('Pair arguments need to be string!');
        end
        switch(MYARG),
          case 'station',
            if ~iscell(MYVAL),
              error([MYARG ': needs cell array!']);
            end
            STATIONS = MYVAL;
          case 'years',
            if ~isnumeric(MYVAL),
              error([MYARG ': needs to be vector!']);
            end
            years = MYVAL;
          case 'tag',
            if ~ischar(MYVAL),
              error([MYARG ': needs a string!']);
            end
            SimTag = MYVAL;
          case 'BASEPATH',
            if ~ischar(MYVAL),
              error([MYARG ': needs characters!']);
            end
            BASEPATH = MYVAL;
          case 'OUTPUTPATH',
            if ~ischar(MYVAL),
              error([MYARG ': needs characters!']);
            end
            OUTPATH = MYVAL;
          case 'TXT',
            if ~islogical(MYVAL)
              error([MYARG ': needs false or true!']);
            end
            TXTF = MYVAL;
          otherwise,
            error(['Input arg ' MYARG ' not recognized!']);
        end
      end
  end
    
  % ASSIGNUNG PARAMETERS with default values:
  % boolean parameters:
  if ~exist('TXTF','var'),
    TXTF = false;    % TRUE for ASCII output file (old version)
  end
  % Output path:
  if ~exist('OUTPATH','var'),
    OUTPATH = [HOME_BASE '/GFI/data/RT/'];        
  end
  % Input path:
  if ~exist('BASEPATH','var'),
    BASEPATH = [HOME_BASE '/GFI/data/RASOBS/'];
  end
  % Stations names:
  if ~exist('STATIONS','var'),
    warning('No Station cell has been specified!');
  end
  % years:
  if ~exist('years','var'),
    warning('No year range has been specified!');
  end
  % Simulations tag:
  if ~exist('SimTag','var'),
    SimTag = '';
  end

  %% Parameter definitions:
  topH = 10100;   % top profile height [m]
  MinLev = 40;   % minimum continues layers the RadioSonde must have
	
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

  nlay = length(H0); % Number of layers for the homogenized layers

  %% Definition of Variable names and attributes for NetCDF file:
  % Variables needs to be included in the ASCII file
  surface_names = {'PRES','TEMP','RELH'};
  layers_var = {'HGHT','PRES','TEMP','RELH','CLOUD','RAIN'};

  % metadata for the NetCDF file (in the storing order):
  wrf_surfname = {'PSFC','Surface Pressure','hPa';...
									'T2','2-m Temperature','K';...
									'Q2','Surface Relative Humidity','%'};

  wrf_varname = {
								 'PHB','P','T','RH','QVAPOR','QCLOUD','QRAIN','QICE', ...
								 'QSNOW','QGRAUP','WINDDIR','WINDVEL'};

  wrf_locname = {'HGT','Station Altitude','km';...
								 'LAT','Sation Latitude','deg';...
								 'LON','Station Longitude','deg'};
    
  date_name = {'year','month','day','hour'};
    
  % Radiosonde variables to load:
  nc_metadata = {
								 'HGHT',    'Geopotential Height','km';...
								 'PRES',    'Atmospheric Pressure','hPa';...
								 'TEMP',    'Temperature','K';...
								 'RELH',    'Relative Humidity','%';...
								 'MIXR',    'Mixing Ratio', 'g/kg';...
								 'CLOUD',   'Cloud Water Content','g/m^3';...
								 'RAIN',    'Rain Water Content','g/m^3';...
								 'IWC',     'Ice Content','g/m^3';...
								 'SNOW',    'Snow Content','g/m^3';...
								 'GRAUPEL', 'Graupel Content','g/m^3';...
								 'DRCT',    'Wind Direction','deg';...
								 'SKNT',    'Wind Speed','knot'...
  };

  % quality factor of the profiles to storage in mat-file
  QC = [];
    
  % weather flag of the profiles: 1=clear sky, 2=cloudy, 3=rainy
  WC = [];
    
  % Indexing for NetCDF files:
  % * this is index for different stations
  n_x = 0;
  % * this index is for different years
  n_y = 0;
  
  y = 0;

    %************************************************************
    %% The cell string 'STATIONS' must be organized according to
    % the (Longitude,Latitude) grid cell desired into the NetCDF file.
    % e.g. When only one file is specified, then the grid as only one point (1,1).
    
    [mxN_x, mxN_y] = size(STATIONS);

    % string containing information about the origin RS
    origin_str = '';

    %% The Database is orginized following the WRF-grid as close as possible:
    % * (x_n, y_n): are the coordinates of the specified stations, when only
    % one station is introduced then the dimension xn and yn is 1.
    % * time: is the time converted to year,month,day,hour and it is
    % sorted according to the readed from the directory there the files
    % are located: A = dir([fullfile(INPATH, num2str(years)), '*.mat'])
    
    N_year = length(unique(years));

    % Defining NetCDF and ASCII output file:
    outfilen = sprintf('WyoRS_Y%04d-%04d_x%03dy%03dST%s.nc',...
    years(1),years(end),mxN_x,mxN_y,SimTag(1:min(end,10)));

    if exist(fullfile(OUTPATH, outfilen), 'file'),
        disp(['Deleting... NetCDF file: ' OUTPATH outfilen]);
        delete(fullfile(OUTPATH, outfilen) );
    end

    % Starting to read over number of Stations and years:
    % number of stations
    for n_x = 1:mxN_x,
        % number of years
        for n_y = 1:mxN_y,
            Last_obs = 0;
            for i_year=1:N_year,
                INPATH = [BASEPATH '/'...
                STATIONS{n_x,n_y} '/' num2str(years(i_year))];
                A = dir(fullfile(INPATH, 'RS_*.mat'));
                N_f = length(A);
		
				% number of files in directory:
                for f = 1:N_f,
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
                
                    % ----------------------------------------
                    % ** Quality Control flag 8 bits:
                    % bit 1: height reaches at-least 10km
                    % bit 2: height increases monotonilcally
                    % bit 3: no gaps in leyers
                    % bit 4: pressure always decrease with height increase
                    QCflag = int8(zeros(nobs,1));
                    
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
                    for i = 1:nobs,
                      idxobs = idxobs + 1;
                      % 15 corresponds to all Quality test passed
                      %continue; Including even bad profiles!
                      if QCflag(i) ~= 15,
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

											% ----
											% CONVERTION OF VARIABLES FOR BETTER TREATMENT
                      % Converting Temperature units from C to Kelvin:
                      data(i).TEMP = data(i).TEMP+273.15;
											% Converting wind speed and direction to vector:
											data(i).WVECX = data(i).SKNT.*cosd(data(i).DRCT);
											data(i).WVECY = data(i).SKNT.*sind(data(i).DRCT);
											
											% getting the Surface Data (extrapolating to Height 0):
                      [Zunq, Iunq] = unique(data(i).HGHT);
                      for k=1:length(surface_names),
                        eval(['tmp = find(Zunq<=topH & ~isnan(data(i).'...
																surface_names{k} '(Iunq)));']);
                        
                        % In case no profile is available or less
                        % than 3 points are available in the
                        % lower Atmospheric layers:
                        if isempty(tmp) | length(tmp)<3,
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

												% end over k, number of surface_names
                      end

											% Interpolating the profile variables to fixed layers below topH:
                      INidx = Iunq(find(Zunq <= topH)); % indices for unique H below topH
                      
                      h_tmp = data(i).HGHT(INidx);
                      Nh_tmp = length(h_tmp);

											% getting the Inter-quantile-region for profile variables:
											IQRvar = structfun(@(x) quantile(x(INidx), [.25 .75]), data(i), 'UniformOutput', 0);
											% getting the limits inside non-outliers values:
											OUTlie = structfun(@(x) [x(1) x(2)]+[-1 1]*1.3*diff(x), IQRvar, 'UniformOutput', 0);

											% getting the altitues and profiles only without NaNs:
                      hi = structfun(@(x) h_tmp( ~isnan(x(INidx))), data(i), 'UniformOutput', 0);
                      vi = structfun(@(x) x(INidx( ~isnan(x(INidx)) )), data(i), 'UniformOutput', 0);
                      tt = structfun(@(x) length(x), vi);
                      names = fieldnames(vi);
                      for k=1:length(names),
												% ---
												% Before interpolation, take out possible profile outlayers 
                        Xin = getfield(hi, names{k});
                        Yin = getfield(vi, names{k});
												YLIM = getfield(OUTlie, names{k});
												noliers = find(Yin>YLIM(1) & Yin<YLIM(2));
												Xin = Xin(noliers);
												Yin = Yin(noliers);
												NYin = length(Yin);
                        % For Relative Humidity variable, the surface
                        % values is included for the interpolation
                        % to avoid high values or negarive values:
                        if strcmp(names{k},'RELH'),
                          Xin = [0; Xin];
                          Yin = [SURFVars(idxobs,3); Yin];
                        end

												
                        if NYin>fix(Nh_tmp/2) && Nh_tmp>10,
													Vinter = interp1(Xin, Yin, H0, 'linear', 'extrap');
                        else
                          Vinter = NaN*ones(nlay,1);
                          for badin=1:NYin,
                            [dummy, Idxin] = min(abs(Xin(badin)-H0));
                            Vinter(Idxin) = Yin(badin);
                          end
                        end

												% checking for outliers after inter-/ extra-polation again:
												noliers = find(Vinter>YLIM(1) & Vinter<YLIM(2));
												liers = find(Vinter<=YLIM(1) | Vinter>=YLIM(2));
												if ~isempty(liers),
													% correcting the outliers:
													Xin = H0(noliers);
													Yin = Vinter(noliers);
													Vinter(liers) = interp1(Xin, Yin, H0(liers),'nearest');
												end
                        eval(['TEMPVars.' names{k} ' = Vinter;']);
                       
                      end

											% Converting the Height units from m to km and from
											% station level to a.s.l.:
                      TEMPVars.HGHT = (TEMPVars.HGHT)/1e3;

											% Converting and replacing the Wind Vector back to Wind Speed and Wind Direction
											TEMPVars.SKNT = sqrt( TEMPVars.WVECX.^2 + TEMPVars.WVECY.^2 );
											TEMPVars.DRCT = mod( 360 + atan2d(TEMPVars.WVECY, TEMPVars.WVECX), 360);
											
											% **********************************************
                      [TEMPVars.CLOUD, TEMPVars.RAIN, TEMPVars.IWC,...
                       TEMPVars.SNOW, TEMPVars.GRAUPEL] = cloud_modell(TEMPVars.TEMP,...
																																			 TEMPVars.RELH,TEMPVars.PRES);

                      names = fieldnames(TEMPVars);
                      for j=1:length(names),
                        eval(['ALLVARS(:,j,idxobs) = TEMPVars.' names{j} ';']);
                      end

                      % ---------------------------------------------
                      % ** Weather flag for the profiles:
                      % bit 1: clear sky
                      % bit 2: cloudy sky
                      % bit 3: rainy
                      % bit 4: no specification
                      WCflag(i,1) = int8(0);
                      if any( find(TEMPVars.CLOUD>0)),
												% cloudy
												WCflag(i) = bitset(WCflag(i), 2);
											else
												% clear
												WCflag(i) = bitset(WCflag(i), 1);
											end

											% ad-doc estimation of rainy profile:
											if WCflag(i) == 2 && all( TEMPVars.RELH(1:15) > 98),
												WCflag(i) = bitset(WCflag(i), 3);
											end
											
											clear TEMPVars;
                        
											% -----------------------------------------------
											% Finding the indexes for variables to store
											% from nc_metadata cell array (defined up)
											%idxcol = [2,1,3,5,12,13,14,15,16];   % index of variables to store
                      idxcol = cellfun(@(x) find(strcmp(x,names)), nc_metadata(:,1));
                    
                      % Saving data into a ASCII file:
                      if TXTF,
                        fprintf(fp,'%2d %2d\n',day,hour);
                        fprintf(fp,'%6.3f %10.3f %10.3f %10.3f\n',...
																metvar.SELV(i)/1e3,SURFVars(idxobs,:));
            
                        fprintf(fp,'%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',...
																transpose(ALLVARS(:,idxcol,idxobs)));

											end
                    
											% --
											% end over i: number of observations per File 
										end
										% ------------------------------------------------------------
		    
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

										% -- Quality Control index
                    if n_x==1 & n_y==1 & i_year==1,
                      nccreate(ncfile,'QIDX','Datatype','int8', ...
                               'Dimensions',{'xn',mxN_x,'yn',mxN_y,'time',Inf},...
                               'FillValue',NaN);
                    end
                    ncwrite(ncfile,'QIDX',shiftdim(QCflag,-2),[n_x, n_y, 1+Last_obs]);
                    ncwriteatt(ncfile,'QIDX','long_name','Profile Quality index');
                    ncwriteatt(ncfile,'QIDX','note:','15=good quality');

										% -- Weather Control index
										if n_x==1 & n_y==1 & i_year==1,
                      nccreate(ncfile,'WEIDX','Datatype','int8', ...
                               'Dimensions',{'xn',mxN_x,'yn',mxN_y,'time',Inf},...
                               'FillValue',NaN);
                    end
                    ncwrite(ncfile,'WEIDX',shiftdim(WCflag,-2),[n_x, n_y, 1+Last_obs]);
                    ncwriteatt(ncfile,'WEIDX','long_name','Weather Type index');
                    ncwriteatt(ncfile,'WEIDX','note:','1=clear, 2=cloudy, 6=rain');
		    
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
end
% ============ End of main function RS_homogenize_profiles =================

    
% **************************************************************************
% Cloud model function for estimation of liquid water cloud
% profile
% Including Mätzler cloud model (to be reviewed)
%

function [LWC, LWR, IWC, SWC, GWC] = cloud_modell(T, RH, P)
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
  LWC(1,:) = zeros(1, length(T));
  LWC(1,tmp) = 2*((RH(tmp)-b0)/30).^2;
  
	% Rain Water content [g/m^3]:
  LWR(1,:) = zeros(1, length(T));
  
	% Cloud Ice content [g/m^3]:
  IWC(1,:) = zeros(1, length(T));
  
	% Show content [g/m^3]:
  SWC(1,:) = zeros(1, length(T));
  
	% Graupel content [g/m^3]:
  GWC(1,:) = zeros(1, length(T));
end

% *************************************************************************
% End of script
