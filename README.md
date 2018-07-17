# WyoRASOBS (Wyoming University RAdiosondeS OBSservations)

## A GNU/Octave & MATLAB Automatic fetch and storage interface for Radiosonde data provided by the University of Wyoming

This is a simple MATLAB and Octave friendly function to fetch radiosonde data from the publich repository hosted by the University of Wyoming.

### Description
This MATLAB/GNU Octave function gets the DATA of the soundings (raob) from the Wyoming University internet site (http://weather.uwyo.edu)

The Function can be used by calling one or any combinations of the following set of parameters:

    > [data, {metinfo}] = RASOBS_DOWNLOAD_DATA_RAW(station,year,month,day,hour);
  
    > [data, {metinfo}] = RASOBS_DOWNLOAD_DATA_RAW(station,year,month,day,hour,'outputpath','/whereto/storage/data/');
  
    > [data, {metinfo}] = RASOBS_DOWNLOAD_DATA_RAW(station,year,month,day,hour,'netcdf',true);
  
    > [data, {metinfo}] = RASOBS_DOWNLOAD_DATA_RAW(station,year,month,day,hour,'csvfiel',true);
  
    > [data, {metinfo}] = RASOBS_DOWNLOAD_DATA_RAW(station,year,month,day,hour,'matfile',true);

where,
INPUT VARIABLES

    % station-> (string) Code for the station to download;
    % year  --> (numeric) Range for Year of the date we want to download, e.g. [2015:2018];
    % month --> (numeric) Range for Month of the date we want to dowload, e.g. [10:12];
    % day   --> (numeric) Range for Day of the date we want to download, e.g. [1:31];
    % hour  --> (numeric) Range for Hour of the date we want to download [00,12];
    %           hour will accept any range of hours, but Wyoming University only has two radiosonde per day at 00 and 12Z.
FOLLOWING OPTION INPUTS (keyword value pairs)

    % 'netcdf', true/false --> whether to storage as NetCDF file (recomended for large datasets);
    % 'csvfile', true/false --> whether to storage as CSV files (creates one file per hour, feasable for short datasets);
    % 'matfile', true/false --> whether to storage as v7 MATLAB files;
    % 'outputpath','/path_to/storage/' --> Dir where files are storaged (see below for default storage directory);

Default options:

    % 'outputpath'->'../data/RASOBS/station/yyyy/',
    % 'matfile'->true,
    % 'csvfile'->false,
    % 'netcdf'->false.

OUTPUT VARIABLES:

    % data --> MATLAB Structure variable containing the Radiosonde data and related information.
    % In case of unsuccessful downloading (e.g. hour 10 doesn't exist) it returns an empty variable or skip that hour.
OPTIONAL OUTPUT:

    % metinfo --> MATLAB Structure variable with Sounding Station Parameters and Indices for every
    % profile (for a extended description of indices, see HTML http://weather.uwyo.edu/upperair/indices.html).
    
Downloaded data can be archived in one or any combination of following file formats -> CSV, NetCDF or MATLAB binary format (depending on optional input value pairs the data archive will have a generic name e.g. 'yyyymmdd_hh.csv', 'yyyymmdd.nc' or 'yyyymmdd.mat').

The new generated file is saved by default in the following folder: `'../data/RASOBS/namestation/yyyy/'` or in the path specified by the `'outputpath'` option. 

The units for the profile variables (columnwise) are stored as a member variable named UNITS in the optional structure output `metinfo`, i.e. 

    > metinfo.UNITS
    ans = 
    {'hPa',  'm', '°C', '°C', '%', 'g/kg', 'deg', 'knot', 'K', 'K','K'}

In case of storage in NetCDF format, the units are displayed as argument for every available variable.

### Description of the profile variables
A detailed description of the variables in the profile is provided by the Wyoming University web site:
http://weather.uwyo.edu/upperair/columns.html

### Description of the Sounding indices
Similrly a detailed description of the Radiosonde indices is provided by the Wyoming University web site:
http://weather.uwyo.edu/upperair/indices.html

## USAGE EXAMPLE
This is an example to download data from the Station Norderney with station number '10113' for the year 2015, the months of Mai, June and July (`[5:7]`) and all available hours i.e. 00 and 12Z (`[0,12]`). The data will the then storaged in a NetCDF file in the default directory `../data/norderney/2015/05/`.

    >> [data,metinfo] = RASOBS_DOWNLOAD_DATA_RAW('10113',2015,[5:7],[1:31],[00,12],'netcdf',true);
    
    warning: Data from 06.31.2015_00UTC does not exist!
    warning: Data from 06.31.2015_12UTC does not exist!
    warning: Data from 07.27.2015_00UTC does not exist!
    warning: called from RASOBS_DOWNLOAD_DATA_RAW at line 175 column 21
    
for the example above, the downloaded and storage took approximatelly 7 minutes with GNU Octave.

    >> whos data
    Variables in the current scope:
    Attr Name        Size                     Bytes  Class
    ==== ====        ====                     =====  =====
        data        1x183                  1705968  struct

    >> whos metinfo
    Variables in the current scope:
    Attr Name         Size                     Bytes  Class
    ==== ====         ====                     =====  =====
        metinfo      1x1                      41013  struct



(c) 2018 P. Saavedra Garfias, Geophysical Institute, UNIVERSITY OF BERGEN
See LICENSE.TXT
