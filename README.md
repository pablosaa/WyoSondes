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
    % year  --> (numeric) Range for Year of the date we want to download;
    % month --> (numeric) Range for Month of the date we want to dowload;
    % day   --> (numeric) Range for Day of the date we want to download;
    % hour  --> (numeric) Range for Hour of the date we want to download [00,12];
FOLLOWING OPTION INPUTS (keyword value pairs)

    % 'netcdf', true/false --> whether to storage as NetCDF file;
    % 'csvfile', true/false --> whether to storage as CSV files;
    % 'matfile', true/false --> whether to storage as v7 MATLAB files;
    % 'outputpath','/path_to/storage/' --> Dir where files are storaged;

Default options:

    % 'outputpath'->'../data/RASOBS/station/yyyy/',
    'matfile'->true,
    'csvfile'->false,
    'netcdf'->false.

OUTPUT VARIABLES

    % data --> MATLAB Structure variable containing the RS data and related information.
    % In case of unsuccessful downloading it retures an empty variable.
OPTIONAL OUTPUT

    % metinfo --> MATLAB Structure variable with Sounding Station Parameters and Indices for every
    % profile (see HTML http://weather.uwyo.edu/upperair/indices.html)
    
A new file can be created in any or all formats -> CSV, NetCDF or MATLAB binary format (depending on optional input value pairs e.g. 'yyyymmdd_hh.csv', 'yyyymmdd.nc' or 'yyyymmdd.mat').

This file is saved by default in the following folder: `'../data/RASOBS/namestation/yyyy/'` or in the specified `'outputpath'` option. 

The units for the profile variables are (columnwise):

    hPa,  m, °C, °C, %, g/kg, deg, knot, K, K, K
and are stored as a member variable named UNITS in the optional structure output `metinfo`, i.e. `metinfo.UNITS`.


### Description of the profile variables

http://weather.uwyo.edu/upperair/columns.html

### Description of the Sounding indices

http://weather.uwyo.edu/upperair/indices.html

(c) 2018 P. Saavedra Garfias, Geophysical Institute, UNIVERSITY OF BERGEN
See LICENSE.TXT
