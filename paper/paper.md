---
title: '``WyoSondes``: An Octave and MATLAB® automatic radiosonde time-series dataset generator using the Wyoming University repository'
authors:
  - name: Pablo Saavedra~Garfias
    orcid: 0000-0002-4596-946X
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Jochen Reuder
    affiliation: "1, 2"
tags:
  - Octave 
  - Matlab
  - Radiosonde
  - Atmosphere
  - Meteorology
  - Climate
  - Wyoming University
affiliations:
 - name: Geophysical Institute, University of Bergen
   index: 1
 - name: Bjerknes Centre for Climate Research
   index: 2
date: 13 January 2020
bibliography: paper.bib

---

# Abstract

# Introduction

Atmospheric soundings are extremely useful tools for weather analysis and forecasting. Radiosonde provides a great amount of in-situ information about the states of the troposphere and stratosphere (@milrad_14_2018).
The standard observed parameters are profiles of temperature, humidity and wind speed and direction. Those are essential for the development of modern meteorology and its current practice of modelling and monitoring the atmosphere (@brettle_back_2003). Thus radiosonde data is widely needed as forcing parameter for numerical weather
prediction models, data assimilation, climate studies, aviation, environmental monitoring, among others (@galvin_back_2003).

Radiosonde data are widely required by numerous applications like in meteorology for weather forecast around the globe and climate studies by characterization of the climatology, and further uses in data assimilation, polution studies, among others (@seidel_estimating_2010).

Since Radiosonde observations became of widespread use in computer models, the World Meteorological Organization (WMO) prescribes the launches to be made as a regular synoptic routine, and most radiodondes are launched twice a day at around 00:00 and 12:00 UTC (@galvin_back_2003) with over two thousands of launches performed by meteorological stations, research facilities and airports across the world.

Data are then format-homogenized and shared to be available from public databases. The Wyoming University in USA hosts one of the most important and popular sources of radiosonde data in the world, with global coverage and freely available to download. Normally data is directly retrieved from its web-site via a friendly web GUI which allows the selection of station and specify a range of dates. The data is then displayed in the browser as a simple ASCII web-site. This is convenient only for single radiosondes or few time series, nevertheless when a large time series of radiosonde data from various stations are required, an automatic approach is more suitable to build historic radiosonde dataset.

# Implementation

``WyoSonde`` is a GNU-Octave (@octave_manual) and MATLAB® function to automatically fetch the radiosonde data from the University of Wyoming public data repository (@wyoming_website). 

The script has been used to build climatological dataset of atmospheric profiles for specific stations around measurement campaigns and use that dataset as a-priori information to perform retrievals by a microwave radiometer (@saavedra:2019).

This repository has been mainly developed and intensively used with GNU-Octave v4.4.0 and v5.1.0 under GNU Linux (tested under OpenSuse and Ubuntu distributions). Additional testing has been done for Matlab® R2016a and R2018b.

Additionally, the repository contains visualization tools to quickly review the downloaded database and optionally produce altitude homogenized netCDF file with quality control checks for further use. Those however are not essential to the main script.

# Acknowledgments
The authors want to thank to the Wyoming University for making the radiosonde data freely
available. This work was developed under support by the Off-shore Boundary Layer
Observatory (OBLO) project from the University of Bergen.

# References
