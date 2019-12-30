---
title: 'WyoSondes: GNU Octave/MATLAB Automatic fetching radiosonde time-series
from the Wyoming University data repository'
tags:
  - Octave 
  - Matlab
  - Radiosonde
  - Atmosphere
  - Meteorology
  - Climate
  - Wyoming University
authors:
  - name: Pablo Saavedra Garfias
    orcid: 0000-0002-4596-946X
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Jochen Reuder
    affiliation: "1, 2"
affiliations:
 - name: Geophysical Institute, University of Bergen
   index: 1
 - name: Bjerknes Centre for Climate Research
   index: 2
date: 13 November 2019
bibliography: paper.bib

---

# Summary

Atmospheric profiles are traditionally obtained via radiosonde which are the main source of in-situ information for
the states of the troposphere and stratosphere. 
Radiosonde profiles are widely needed as input for numerical weather
prediction models, climate studies, aviation, among others.
Meteorological stations, research facilities and airports around the world
collect radiosonde data daily, that data are shared and format-homogenized to
be available from the Wyoming University data repository web-site via a
friendly web GUI which allows the selection of station and date range to be
retrieved as a simple ASCII web-site. However when large time series of
radiosonde data from various stations are needed, then an automatic approach
is more suitable to build large datasets.

```WyoSonde``` is a simple GNU-Octave and MATLAB function to automatically
fetch the atmospheric radiosonde data from the public data repository hosted by the University of Wyoming (@wyoming_website).

Radiosonde data are widely required by numerous applications like in meteorology for weather forecast around the globe and climate studies by characterization of the climatology, and further uses in data assimilation, polution studies, among others (@seidel_estimating_2010).

Radiosonde observations are normally performed twice a day at 00:00 and 12:00 UTC with over two thoudands of launches across the world.
The University of Wyoming in USA hosts one of the most important sources of radiosonde data repositories in the world, with global coverage and freely available to download.

This repository has been mainly developed and intensively used with GNU/Octave v4.4.0 under Linux OpenSuse distribution. Additional testing has been done for Matlab R2016a.

# Acknowledgments
We want to thank the Wyoming University for making the radiosonde data freely
available. This work was developed under the Off-shore Boundary Layer
Observatory (OBLO) project from the University of Bergen.

# References
