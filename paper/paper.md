---
title: 'WyoSondes: time-series of atmospheric profiles fetched from the Wyoming University radiosonde data repository'
tags:
  - Octave 
  - Matlab
  - Atmosphere
  - Radiosonde
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

This is a simple GNU/Octave and MATLAB function to automatically fetch atmospheric radiosonde data from the public data repository hosted by the University of Wyoming (@wyoming_website).

Radiosonde data are widely required by numerous applications like in meteorology for weather forecast around the globe and climate studies by characterization of the climatology, and further uses in data assimilation, polution studies, among others (@seidel_estimating_2010).

Radiosonde observations are normally performed twice a day at 00:00 and 12:00 UTC with over two thoudands of launches across the world.
The University of Wyoming in USA hosts one of the most important sources of radiosonde data repositories in the world, with global coverage and freely available to download.

This repository has been mainly developed and intensively used with GNU/Octave v4.4.0 under Linux OpenSuse distribution. Additional testing has been done for Matlab R2016a.
