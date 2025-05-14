Code for the SiftDist metric
----------------------------
Ofir Pele and Michael Werman 
Contact: ofirpele@cs.huji.ac.il
Version: 2.0, Dec, 2008

This directory contains the source code for the SiftDist metric.

See the web page at 
http://www.cs.huji.ac.il/~ofirpele/SiftDist/

Please cite this paper if you use this code:
A Linear Time Histogram Metric for Improved SIFT Matching
Ofir Pele, Michael Werman
ECCV 2008

bibTex:
@INPROCEEDINGS{Pele-eccv2008,
author = {Ofir Pele and Michael Werman},
title = {A Linear Time Histogram Metric for Improved SIFT Matching},
booktitle = {ECCV},
year = {2008}
}

The current version (2.0) is not the same as the one used in the paper. The main change is
the addition of a rejection scheme that reduces running time. In addition, I suggest to send sqrt(descr)
as it gives better results.

Easy startup
------------
Within Matlab:
>> demo_SiftDist

Compiling (the folder contains compiled binaries, thus you might not have to compile)
-------------------------------------------------------------------------------------
Within Matlab:
>> compile_SiftDist

Usage within Matlab
------------------- 
Type "help SiftRatioMatch" in matlab
Type "help SiftDist" in matlab


Usage within C++
----------------
See "SiftDistTest.cxx" and "SiftDist.hxx"


Licensing conditions
--------------------
This software is made available for non-profit research purposes only.
It is necessary to obtain a license from the Hebrew University of Jerusalem for
commercial applications. See the file LICENSE in this directory for conditions
of use.


CREDITS
-------
img1.ppm, img3.ppm were downloaded from: 
http://www.robots.ox.ac.uk/~vgg/data/data-aff.html
