Code for the EMD_MOD histogram distance
---------------------------------------
Ofir Pele 
Contact: ofirpele@cs.huji.ac.il
Version: 1.1, Nov 2011

This directory contains the source code for computing the EMD_MOD histogram distance efficiently.

Please cite this paper if you use this code:
 A Linear Time Histogram Metric for Improved SIFT Matching
 Ofir Pele, Michael Werman
 ECCV 2008
bibTex:
@INPROCEEDINGS{Pele-eccv2008,
author = {Ofir Pele and Michael Werman},
title = {A Linear Time Histogram Metric for Improved SIFT Matching}
booktitle = {ECCV},
year = {2008}
}
Two things to take into considerations:
1. For SIFT matching, this distance performance was bad.
   Try using my QC / FastEMD / SIFT_DIST codes
   (from ECCV 2010, ICCV 2009, ECCV 2008 respectively).
2. This implementation is linear on the average as I use
   the nth_element function of standard C++ library,
   which is currently average linear time. There are worst time
   linear algorithms, but I did not implement them.

Easy startup
------------
Within Matlab:
>> demo_EMD_MOD (1d histograms)

Compiling (the folder contains compiled binaries, thus you might not have to compile)
-------------------------------------------------------------------------------------
Within Matlab:
>> compile_EMD_MOD

Usage within Matlab
------------------- 
Type "help EMD_MOD" in Matlab.

Usage within C++
----------------
See EMD_MOD.hpp

Licensing conditions
--------------------
See the file LICENSE.txt in this directory for conditions of use.
