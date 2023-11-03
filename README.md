# Distributed H₂ Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7125045.svg)](https://doi.org/10.5281/zenodo.7125045)

## General

This repository contains an implementation of the algorithms and simulations described in the paper

> C. Hespe, A. Datar, D. Schneider, H. Saadabadi, H. Werner, and H. Frey, "Distributed H₂ Controller Synthesis for Multi-Agent Systems with Stochastic Packet Loss," *American Control Conference*. San Diego, CA, USA, 2023, doi: [10.23919/ACC55779.2023.10156194](https://doi.org/10.23919/ACC55779.2023.10156194).

It may be used to recreate and validate the figures from the paper.
To do so, run either of the two main entry points in the repository, the scripts `validate_bernoulli.m`, and `evaluate_sinr.m`.
Be advised that especially the latter has an extensive runtime of about half a day.
The raw data used in the figures in the paper is available in the subdirectory `figures`.

## Simulation of Multi-Agent Systems

Both of the main entry points will run simulations of multi-agent systems with packet loss.
For that purpose, an open source library which can be found [on Github](https://github.com/TUHH-ICS/MAS-Simulation) is utilized.
The repository is included as a git submodule and needs to be initialized using `git submodule update --init`.
The simulation run by `evaluate_sinr.m` require the compiled components of the SINR part of that library, so be sure to set up a C++ compiler for you system for use in Matlab first.
Steps to do so are described [here](https://github.com/TUHH-ICS/MAS-Simulation#using-the-sinr-networking-library).

## Prerequisites

To run the scripts in this repository, you will need a working copy of [*Yalmip*](https://yalmip.github.io/) together with a suitable SDP solver in your *Matlab* path.

The code in this repository was tested in the following environment:

* *Windows 10* Version 21H2
* *Matlab* 2021a
* *Yalmip* 31-March-2021

The *Matlab* [`parfor`](https://de.mathworks.com/help/parallel-computing/parfor.html) feature from the *Parallel Computing Toolbox* is used to speed up the calculations.
*Matlab* should automatically detect if that toolbox is not available and run the iterations sequentially in that case.
However, this will drastically prolong the runtime of the scripts.
You may want to reduce the number of sampling points for the figures or run the calculations for smaller networks.
