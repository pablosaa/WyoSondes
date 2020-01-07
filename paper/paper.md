---
title: 'WyoSondes: GNU/Octave and MATLAB© automatic fetching radiosonde time-series from the Wyoming University data repository'
authors:
  - name: "Pablo Saavedra Garfias"
    orcid: 0000-0002-4596-946X
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: "Jochen Reuder"
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
date: 13 November 2019
bibliography: paper.bib

---

# Summary

Atmospheric soundings are extremely useful tools for weather analysis and forecasting. Radiosonde is a fundamental source of in-situ information of
the states of the troposphere and stratosphere (@milrad_14_2018).
The observed profiles of temperature, humidity and wind are essential for the development of modern meteorology, and its current practice of modelling and forecasting of the atmosphere (@brettle_back_2003). Thus radiosonde data is widely needed as input for numerical weather
prediction models, climate studies, aviation, environmental monitoring, etc.  (@galvin_back_2003).

Since Radiosonde observations and computer models came into widespread use, the World Meteorological Organization (WMO) prescribes the launches be made to a regular synoptic routine, and most are launched twice a day at around 00:00 and 12:00 UTC (@galvin_back_2003) with over two thoudands of launches across the world.

Meteorological stations, research facilities and airports around the world
collect radiosonde data periodically. 

Data are then shared and format-homogenized to be available from the Wyoming University data repository web-site via a
friendly web GUI which allows the selection of station and date range to be
retrieved as a simple ASCII web-site. However when large time series of
radiosonde data from various stations are needed, then an automatic approach
is more suitable to build large dataset.

Soundings are an extremely powerful weather analysis and forecasting tool that can provide a great amount of information. First, the reader is taught how to identify cloud layers, tropopause height, and temperature advection from plotted soundings.

```WyoSonde``` is a GNU/Octave (@octave_manual) and MATLAB© function to automatically fetch the atmospheric radiosonde data from the University of Wyoming public data repository (@wyoming_website). 

Radiosonde data are widely required by numerous applications like in meteorology for weather forecast around the globe and climate studies by characterization of the climatology, and further uses in data assimilation, polution studies, among others (@seidel_estimating_2010).

The University of Wyoming in USA hosts one of the most important sources of radiosonde data repositories in the world, with global coverage and freely available to download.

This repository has been mainly developed and intensively used with GNU/Octave v4.4.0 under Linux OpenSuse distribution. Additional testing has been done for Matlab© R2016a.

# Acknowledgments
We want to thank the Wyoming University for making the radiosonde data freely
available. This work was developed under the Off-shore Boundary Layer
Observatory (OBLO) project from the University of Bergen.

# References
