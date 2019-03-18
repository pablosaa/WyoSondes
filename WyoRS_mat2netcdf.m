% [] = WyoRS_mat2netcdf(data, metvar,station,station_number);
% [] = WyoRS_mat2netcdf(mat-file);
%
% Function to convert a Matlab mat-file containing Wyoming
% Radiosonde data fetched and stored by RASOBS_DOWNLOAD_DATA_RAW.m
% 
% (c) 2018 P. Saavedra Garfias, UNIVERSITY OF BERGEN
% Email: pablo.saa@uib.no
% See: LICENSE.TXT
% ---------------------------------------------------------------

function [] = WyoRS_mat2netcdf(varargin)

% defining default values:
    locationname = 'Unknown';
    stationname  = 'Unknown';
    station      = 'Unknown';
    
    % checking for input parameters:
    if nargin==1 && exist(varargin{1},'file'),
        matfile = varargin{1};
        load(matfile);
        [ncpath, ncfile, tmp] = fileparts(matfile);
        ncfile = fullfile(ncpath,[ncfile '.nc']);
    elseif nargin>1 && ~logical(mod(nargin,2)),
        for i=1:2:nargin,
            switch(lower(varargin{i})),
              case {'data'},
                data = varargin{i+1};
              case {'metvar'},
                metvar = varargin{i+1};
              case {'stationname'},
                stationname = varargin{i+1};
              case {'locationname'},
                locationname = varargin{i+1}
              case {'station'},
                station = num2str(varargin{i+1});
              case {'ncdffile'},
                ncfile = varargin{i+1};
              otherwise,
                error(['Argument "' varargin{i} '" not available!']);
            end
        end  % end loop over input arguments
    else
        error('Wrong number of input argumments. See help!');
    end
    
    if exist(ncfile,'file'),
        warning(['netCDF file already exist: ' ncfile]);
    end
    if isempty(ncfile),
       error('Introduce an Output netCDF file where to store!');
    end 
    if ~exist('data','var') || ~exist('metvar','var'),
        error('the variables "data" and "metvar" must be specified!');
    end
    % For GNU/Octave compatibility:
    if exist('OCTAVE_VERSION','builtin'),
        pkg load netcdf;
        disp('NETCDF package loaded... ');
    end
    
    nsonde = length(data);
    names{1} = fieldnames(data);
    
    metadata = WyoRS_metadata();
    ProfileMeta = metadata.PROFILE;
    NameIndices = metadata.RSINDICES(:,1);
    DescriptionIndices = metadata.RSINDICES(:,2);
    UnitsIndices = metadata.RSINDICES(:,3);
    
    for i=1:length(names{1}),
        nccreate(ncfile,names{1}{i},...
                 'Dimensions',{'nsonde',nsonde,'levels',Inf},...
                 'FillValue',NaN,'Format','netcdf4');
        tmp = strfind(ProfileMeta,names{1}{i});
        idxtmp = find(arrayfun(@(k) ~isempty(tmp{k,1}),...
                               [1:length(ProfileMeta)]));
        for j=1:nsonde,
            eval(['ncwrite(ncfile,names{1}{i},transpose(data(j).'...
                  names{1}{i} '),[j,1]);']);
        end
        ncwriteatt(ncfile,names{1}{i},'ShortName', ProfileMeta{idxtmp,1});
        ncwriteatt(ncfile,names{1}{i},'LongName', ProfileMeta{idxtmp,2});
        ncwriteatt(ncfile,names{1}{i},'units', ProfileMeta{idxtmp,3});

    end  % end over variable names
        
        for i=1:length(NameIndices),
            if any([strcmp(NameIndices{i},'SLAT'),...
                    strcmp(NameIndices{i},'SLON'),...
                    strcmp(NameIndices{i},'SELV')])
                continue;
            end
            nccreate(ncfile,NameIndices{i},...
                     'Dimensions',{'nsonde',nsonde},...
                     'FillValue',NaN,'Format','netcdf4');
            eval(['ncwrite(ncfile,NameIndices{i},'...
                  'transpose(metvar.' NameIndices{i} '));']);
            ncwriteatt(ncfile,NameIndices{i},'units',...
                       UnitsIndices{i});
            ncwriteatt(ncfile,NameIndices{i},'LongName',...
                       DescriptionIndices{i});
            
            
        end
        ncwriteatt(ncfile,'/','Observatory Name',[stationname '_' locationname]);
        ncwriteatt(ncfile,'/','Station Code Number',station);
        ncwriteatt(ncfile,'/','Observatory Latitude [deg]',metvar.SLAT(1));
        ncwriteatt(ncfile,'/','Observatory Longitude [deg]',metvar.SLON(1));
        ncwriteatt(ncfile,'/','Observatory Altitude [m]',metvar.SELV(1));
        ncwriteatt(ncfile,'/','Source','http://weather.uwyo.edu/cgi-bin/sounding');
        ncwriteatt(ncfile,'/','Indices Description',...
                   'http://weather.uwyo.edu/upperair/indices.html');
        ncwriteatt(ncfile,'/','Contact','Pablo.Saavedra@uib.no');
        ncwriteatt(ncfile,'/','Institution','Geophysical Institute, Uni-Bergen');


        return;
end

% end of function
