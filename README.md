# CM1-Analysis

Analysis codes for Cloud Model 1 (CM1) outputs

---

## Overview
This repository contains MATLAB scripts for analyzing CM1 outputs. There are scripts for batch processing of data on clusters and visualization on local environment. The scripts are grouped according to their simulation applications.

---

## Applications

- **LES_HBL** — Large-Eddy Simulation of Hurricane Boundary Layer. The scripts analyze the evolution of LES to identify the quasi-steady state, extract the domain and temporal averaged vertical profiles of various quantities. Please view <a href="https://doi.org/10.1175/JAS-D-20-0227.1" target="_blank">Chen et al. (2021, JAS)</a> and <a href="https://doi.org/10.1007/s10546-016-0207-0" target="_blank">Bryan et al. (2017, BLM)</a> for details of the LES for HBL.

- **GW** — Atmospheric Gravity Waves. The scripts extract the vertical snapshots to visualize the simulation of atmospheric gravity waves (specifically, for mountain wave simulation). Please view <a href="https://doi.org/10.1175/1520-0493(1983)111%3C2341:ACMFTS%3E2.0.CO;2" target="_blank">Durran and Klemp (1983, MWR)</a> and <a href="https://doi.org/10.1002/asl.100" target="_blank">Dörnbrack et al. (2005, ASL)</a> for details of mountain wave simulation.