# CM1-Analysis

Analysis codes for Cloud Model 1 (CM1) outputs

---

## Overview
This repository contains MATLAB scripts for analyzing CM1 outputs. There are scripts for batch processing of data on clusters and visualization on local environment. The scripts are grouped according to their simulation applications.

---

## Applications

- **LES_HBL** — Large-Eddy Simulation of Hurricane Boundary Layer. The scripts analyze the evolution of LES to identify the quasi-steady state, extract the domain-averaged and time-averaged vertical profiles of various quantities. Please view [Chen et al. (2021, JAS)](https://doi.org/10.1175/JAS-D-20-0227.1) and [Bryan et al. (2017, BLM)](https://doi.org/10.1007/s10546-016-0207-0) for details of the LES for HBL.

- **GW** — Atmospheric Gravity Waves. The scripts extract vertical snapshots to visualize the simulation of atmospheric gravity waves (specifically, for mountain wave simulation). Please view [Durran and Klemp (1983, MWR)](https://doi.org/10.1175/1520-0493(1983)111%3C2341:ACMFTS%3E2.0.CO;2) and [Dörnbrack et al. (2005, ASL)](https://doi.org/10.1002/asl.100) for details of mountain wave simulation.