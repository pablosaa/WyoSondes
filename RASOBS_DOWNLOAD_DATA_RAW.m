% [data, metinfo, metadata] = RASOBS_DOWNLOAD_DATA_RAW(station,year,month,day,hour, [OPTIONALS]);
% This function gets the 'DATA_raw' of the soundings (raob) from the
% Wyoming University internet site (http://weather.uwyo.edu);
% INPUT
    % station-> (string) Code for the station to download;
    % year  --> (numeric) Range for Year of the date we want to download;
    % month --> (numeric) Range for Month of the date we want to dowload;
    % day   --> (numeric) Range for Day of the date we want to download;
    % hour  --> (numeric) Range for Hour of the date we want to download [00,12];
    % -- FOLLOWING OPTION INPUTS (keyword value pairs)
    % 'netcdf', true/false --> whether to storage as NetCDF file;
    % 'csvfile', true/false --> whether to storage as CSV files;
    % 'matfile', true/false --> whether to storage as v7 MATLAB files;
    % 'outputpath','/path_to/storage/' --> Dir where files are storaged;
    % -- Default options:
    % 'matfile'->true, 'outputpath'->'../data/RASOBS/station/yyyy/', 'csvfile'->false, 'netcdf'->false.
% OUTPUT
    % data --> Structure data containing the RS data and information.
    % In case of unsuccessful downloading it retures empty variable.
    % OPTIONAL OUTPUT
    % metinfo --> Sounding Station Parameters and Indices for every
    % profile (see HTML http://weather.uwyo.edu/upperair/indices.html)
    %
    % A new file can be created in any or all formats -> CSV, NetCDF or
    % MATLAB binary format (depending on optional input value pairs e.g. 'yyyymmdd_hh.csv', 'yyyymmdd.nc' or 'yyyymmdd.mat').
    % This file is saved by default in the following folder:
    % '../data/RASOBS/namestation/yyyy/'
    % or in the specified 'outputpath' option.
    % 
    % The units for the profile variables are (columnwise):
    % {hPa,  m, °C, °C, %, g/kg, deg, knot, K, K, K}
    % and are stored as a member variable named UNITS in the optional
    % structure output metinfo, i.e. metinfo.UNITS
% (c) 2018 P. Saavedra Garfias, UNIVERSITY OF BERGEN
% Email: pablo.saa@uib.no
% See: LICENSE.TXT
% ---------------------------------------------------------------

function varargout = RASOBS_DOWNLOAD_DATA_RAW(station,year,month,day,hour, varargin);


%% CHECKING INPUT VARIABLES
stationname = 'unknown';
ncdfflag = false;
csvflag = false;
matflag = false;
PATH_DAT = [];

if nargout>3,
    error('Maximum three Output variables are needed. See help!');
end

if ismember(nargin,[5,7,9,11]),
    % only two additional options are available
    for i=1:2:length(varargin),
        switch (lower(varargin{i})),
          case {'outputpath'},
            if ischar(varargin{i+1}),
                PATH_DAT = sprintf('%s',varargin{i+1});
            else
                warning('Path needs to be a string!');
            end
          case {'netcdf'},
            ncdfflag = logical(varargin{i+1});
          case {'csvfile'},
            csvflag = logical(varargin{i+1});
          case {'matfile'},
            matflag = logical(varargin{i+1});
          otherwise,
            error(['Argument name "' varargin{i} '" not available!']);
        end
    end
else
    error('Number of input arguments wrong! See help.');
end

%% IF THE CAMPAIGN IAOB_DOWNLOAD_DATA_MATLAB_pablo.m
if strcmpi(stationname,'jokioinen')==1;
    station = '02963';
end;

%% IF THE CAMPAIGN IS TALLINN
if strcmpi(stationname,'Tallinn')==1;
    station = '26038';
end;

%% Defining variable names:
metadata = WyoRS_metadata();
ProfileMeta = metadata.PROFILE;
NameIndices = metadata.RSINDICES(:,1);
DescriptionIndices = metadata.RSINDICES(:,2);
UnitsIndices = metadata.RSINDICES(:,3);

%% Initializing variables:
INDIPAR = zeros(length(NameIndices),1);
idx_hr = 0;


%% DOWNLOAD DATA FROM UNIVERSITY OF WYOMING's RADIOSONDE REPOSITORY
for YEAR=year;

    yyyy = sprintf('%4d',year);
    
    for MONTH = month;
        for DAY = day;
            for HOUR = hour;
                mm = sprintf('%02d',MONTH);
                dd = sprintf('%02d',DAY);
                hh = sprintf('%02d',HOUR);
                timestamp = sprintf('%s.%s.%s_%sUTC',dd,mm,yyyy,hh);
                url = ['http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR='...
                    yyyy '&MONTH=' mm '&FROM=' dd hh '&TO=' dd hh '&STNM=' station];
                html = urlread(url);
                pause(3);
                % Finding Title of HTML database
                tit  = cell2mat(regexpi(html,{'<H2>','</H2>'}));
                if isempty(tit),
                    warning(['Data from ' timestamp ' does not exist!']);
                    continue;
                end
                tmp = textscan(html(tit(1)+4:tit(2)-1),'%s%s%[^\n]',...
                               'CollectOutput',true);
                stationname = tmp{1}{2};
                headerline = tmp{1}{3};
                pos  = regexpi(html,'<PRE>');
                if isempty(pos),
                    warning(['Data from ' timestamp ' does not exist!']);
                    continue;
                end
                ends = regexpi(html,'</PRE>');
                if isempty(ends),
                    warning(['Data from' timestamp 'might be corrupted!']);
                    continue;
                end
                idx_hr = idx_hr+1;
                [fields,idx] = textscan(html(pos(1):end),'%s',2,'Delimiter','','EndOfLine','\n','HeaderLines',2);
                %tmp0 = regexpi(html{1},'<PRE>');
                %tmp1 = regexpi(html{1},'</PRE>');
                %%% find initial index of <PRE>
                %idx0 = find(~cellfun(@isempty,tmp0));
                %idx1 = find(~cellfun(@isempty,tmp1));
                % the first <PRE> block contains the data:
                %%% getting variable names:
                names = textscan(fields{1}{1},'%s');
                units = textscan(fields{1}{2},'%s');

                %% filling data:
                tmpstr = html(pos(1)+idx+1:ends(1)-1);
                ll = regexpi(tmpstr,sprintf('\n'));
                n = length(ll)-1;
                repidx = repmat(ll(1:end-1),10,1)' +...
                         repmat(7*[1:10]+1,n,1);
                tmpstr(repidx)=',';  % including delimiter for blank fields
                                
% $$$                 values = arrayfun(@(k) textscan(tmpstr(ll(k)+1:ll(k+1)), ...
% $$$                                   '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f',...
% $$$                                                 'CollectOutput',true,...
% $$$                                                 'Whitespace','',...
% $$$                                                 'EmptyValue',NaN,...
% $$$                                                 'Delimiter',','),...
% $$$                                   [1:5],'UniformOutput',false);
% $$$                                   
                tmpvalues = textscan(tmpstr,...
                                     '%s%s%s%s%s%s%s%s%s%s%s',...
                                     'Delimiter',',',...
                                     'EndOfLine','\n',...
                                     'Whitespace','',...
                                     'HeaderLines',1,...
                                     'EmptyValue',NaN);

                values = cellfun(@(x) str2double(x), tmpvalues,...
                                 'UniformOutput',false);
                n = length(names{1});
                % data is a structure with all variables in:
                data(idx_hr) = cell2struct(values,...
                                   names{1},2);
                
                clear tmpvalues, values;
                % NOTE tt.HGHT_m=repmat({data.HGHT_m},1,2);
                % ---------------------------------------------
                % the second <PRE> block contains the global info:
                % 
                globaltmp = textscan(html(pos(2):ends(2)-1),...
                                     '%s%s','Delimiter',':',...
                                     'HeaderLines',1,...
                                     'CollectOutput',false);
                description = globaltmp{1};  %globalinfo{idx_hr}.description
                glovalues = globaltmp{2};
		
              % Changing Date format from "YYMMDD/hhmm"->"YYMMDD.hhmm"
		idxtmp = find(strcmp(description,'Observation time'));
                tmp = sscanf(strrep(glovalues{idxtmp},'/','.'),'%06d.%02d%02d');
                tmp = tmp(1)+[1/24 1/3600]*tmp(2:3);
                % Finding which parameters are present:
                [tt,ii] = ismember(description,DescriptionIndices);
                % Filling values for parameter indices:
                INDIPAR(:,idx_hr) = NaN;
                INDIPAR(ii(tt),idx_hr)=cellfun(@str2num,glovalues(tt));
                INDIPAR(ii(idxtmp),idx_hr) = tmp;   % insert the 'Observation time'
                if isempty(PATH_DAT),
                    PATH_DAT = ['../data/RASOBS/'...
                                lower(stationname) '/' ...
                                yyyy '/'];
                end
                if exist(PATH_DAT,'dir');
                else mkdir(PATH_DAT);
                end;

                file_name = [PATH_DAT 'RS_'...
                             sprintf('Y%04d-%04d_M%02d-%02d_D%02d-%02d_H%02d-%02d',...
                                     year([1,end]),month([1,end]),...
                                     day([1,end]),hour([1,end]))];
                %% SAVE DATA as CSV?
                if csvflag,
                    disp('Storing CSV file...');
                    csvfilen = [file_name '_' hh '.csv'];
                    fp = fopen(csvfilen,'w');
                    fprintf(fp,'%% %s\n',headerline);
                    fprintf(fp,'%% %s\n',fields{1}{1});
                    fprintf(fp,'%% %s\n',fields{1}{2});
                    fclose(fp);
                    dlmwrite(csvfilen,cell2mat(values),...
                             '-append','precision','%9.2f');
                end
               
            end;   % end over hours
        end;   % end over days
    end;   % end over months
end;   % end over years

    
    %% Creating output variables:
    qq=arrayfun(@(i) sprintf('metvar.%s=INDIPAR(%d,:);',...
                             NameIndices{i},i),...
                [1:length(NameIndices)],'UniformOutput',false);                      
    cellfun(@eval, qq);
     
    % passing output variables to workspace:
    for i=1:nargout,
        outputvar = {'data','metvar','metadata'};
        if i>3,
            warning('too many output variables! Three are max.');
            break;
        end
        eval(['varargout{i} = ' outputvar{i} ';']);
    end

    %% SAVE DATA as NetCDF or MAT?
    if matflag,
        disp('Storing MAT file...');
        save([file_name '.mat'],'-v7','data','metvar','metadata');
    end
    if ncdfflag,
        if exist('OCTAVE_VERSION','builtin'),
            pkg load netcdf;
            disp('NETCDF package loaded... ');
        end
        nsonde = length(data);
        ncfile = [file_name '.nc'];
        for i=1:length(names{1}),
            nccreate(ncfile,names{1}{i},...
                     'Dimensions',{'nsonde',nsonde,'levels',Inf},...
                     'FillValue',NaN,'Format','netcdf4');
            tmp = strfind(ProfileMeta,names{1}{i});  %strfind(names{1}{i},ProfileMeta);
            idxtmp = find(arrayfun(@(k) ~isempty(tmp{k,1}),...
                                   [1:length(ProfileMeta)]));
            for j=1:nsonde,
                eval(['ncwrite(ncfile,names{1}{i},transpose(data(j).'...
                      names{1}{i} '),[j,1]);']);
            end
            ncwriteatt(ncfile,names{1}{i},'ShortName', ProfileMeta{idxtmp,1});
            ncwriteatt(ncfile,names{1}{i},'LongName', ProfileMeta{idxtmp,2});
            ncwriteatt(ncfile,names{1}{i},'units', ProfileMeta{idxtmp,3});
                %units{1}{i});

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
        ncwriteatt(ncfile,'/','Observatory Name',stationname);
        ncwriteatt(ncfile,'/','Station Code Number',station);
        ncwriteatt(ncfile,'/','Observatory Latitude [deg]',metvar.SLAT(1));
        ncwriteatt(ncfile,'/','Observatory Longitude [deg]',metvar.SLON(1));
        ncwriteatt(ncfile,'/','Observatory Altitude [m]',metvar.SELV(1));
        ncwriteatt(ncfile,'/','Source','http://weather.uwyo.edu/cgi-bin/sounding');
        ncwriteatt(ncfile,'/','Indices Description',...
                   'http://weather.uwyo.edu/upperair/indices.html');
        ncwriteatt(ncfile,'/','Contact','Pablo.Saavedra@uib.no');
        ncwriteatt(ncfile,'/','Institution','Geophysical Institute, Uni-Bergen');
    end   % end if NetCDF
    
    return;

% end of function

%
% from old version:
% $$$  if exist('OCTAVE_VERSION','builtin'),
% $$$                     url = javaObject('java.net.URL',url);
% $$$                     is = javaMethod('openStream',url);
% $$$                     isr = javaObject('java.io.InputStreamReader',is);
% $$$                     br  = javaObject('java.io.BufferedReader',isr);
% $$$ 		else
% $$$                     % for MATLAB case:
% $$$                     url = java.net.URL(url);                
% $$$                     is = openStream(url);
% $$$                     isr = java.io.InputStreamReader(is);
% $$$                     br = java.io.BufferedReader(isr);
% $$$                 end
% $$$ 
% $$$                 fp = fopen([PATH_DAT file_name],'w');
% $$$                 for k = 1:200;             % Read the first 4 lines of text
% $$$                     if exist('OCTAVE_VERSION','builtin'),
% $$$                         s = javaMethod('readLine',br);
% $$$                     else
% $$$                         s = readLine(br);
% $$$                     end
% $$$                     p = char(s);
% $$$                     fprintf(fp,'%s\n',p);
% $$$                 end;
% $$$                 fclose(fp);

