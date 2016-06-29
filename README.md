# code_sample
I have been writing code for almost 15 years and have used a number of programming languages during this time. In this github repository I provide a few code samples that I have written over the last year and also a short summary of my experience. 

In the last few years I have worked mostly with MATLAB and FORTRAN, and also to a limited extent with Python. 

Most of the DIMES data processing and analysis was done using MATLAB. A software program called ARTOA (http://www.whoi.edu/instrument/rafos/artoa-float-tracking) was used for the RAFOS float data processing, which was developed in MATLAB at WHOI. While I was working with ARTOA I was an active developer for the software and helped fix bugs and added some new modules to it. Currently I use MATLAB for my daily programming/sripting operations, which involve analyzing observational datasets. These data sets are small in size and memory management is not a concern. I also plot most of my figures in MATLAB and then tweak them in Inkscape to make them publication quality. 

I started using FORTRAN as an undergraduate. At that time I was working with a simple shallow water model that had been originally developed using FMS (http://www.gfdl.noaa.gov/fms). I added some new numerics to the program (changing from energy conserving to enstrophy conserving schemes), which could be run in parallel on a cluster. The code was dubugged using TotalView and the simulations were performed on a cluster at CMMACS (http://www.cmmacs.ernet.in/). 

At FSU my experience with FORTRAN was originally for small research projects, which were conducted as part of course work. In a numerical ocean modeling class I worked with a reduced gravity model to simulate some simple test scenarios, such as double gyres and initial condition problems. Later, during a 1 year course series on Numerical PDEs, which dealt with analysis of numerical schemes, I had to develop FD and FV solvers from scratch. The codes are modular and use advanced data structures when needed (https://github.com/dhruvbalwada/Numerical_PDEs).

The most recent experience with FORTRAN has been to work with SOSE output. For this work I have been using a code base developed for NEMO (http://servforge.legi.grenoble-inp.fr/projects/CDFTOOLS). However, as the CDFTOOLS package does not have any capabilities for estimating PV budgets, I have had to write higher level programs using the base library functions of CDFTOOLS. 

The use of Python has been primarily for making simple 3D movies (https://www.youtube.com/watch?v=gRnIld6qiw0).

Here are a few code samples from the projects mentioned above:

1. animate_floats_3d.py       - A Python script that uses MayaVi libraries to make a 3D movie. (https://www.youtube.com/watch?v=gRnIld6qiw0)
2. cdfjz_all.f90              - A FORTRAN program to estimate surface PV fluxes. Uses base functions from CDFTOOLS. 
3. PV_surf_fluxes_bulk_SOSE.m - A MATLAB program to estimate surface PV fluxes using bulk formulae. 
4. PV_maps.m                  - A MATLAB program to calculate and plot PV using MIMOC climatology.



