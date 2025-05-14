%SiftDist Computes the SiftDist distance between all pairs of descriptors 
%         in two corresponding sets of descriptors.
%
% dists = SiftDist(descr1, descr2, NumOrientBins, NumSpatialCells, stopThresholdsArr)
%
%--------------------------------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------------------------------
% dists             A (size(feat1,2))x(size(feat2,2)) matrix contains all 
%                   SiftDist distances between all descriptor pairs.
%
%--------------------------------------------------------------------------------------------------
%
%--------------------------------------------------------------------------------------------------
% Descriptors Data Input:
%--------------------------------------------------------------------------------------------------
% descr1,           Two matrices that contains the features.  
% descr2            The number of rows is the dimension of the descriptors 
%                   (has to be equal to the given NumOrientBins*NumSpatialCells)
%                   The number of columns is the number of descriptors in each set.
%                   (same as descr in Andrea Vedaldi's sift code)
%                   Recommendation: send sqrt of the descriptors (see demo_SiftDist)
%
%--------------------------------------------------------------------------------------------------
% Descriptors Layout Input:
%--------------------------------------------------------------------------------------------------
% NumOrientBins     The number of orientation bins.
%                   Original SIFT had 8 bins, but in my paper I found that 16 gives better 
%                   results with the SiftDist distance.
%                   NumOrientBins has to be greater than 1.
%                   (same as NumOrientbins in Andrea Vedaldi's sift code)
%                   Default: 16
%
% NumSpatialCells   The total number of spatial cells. 
%                   Regular SIFT has 4 bins in each direction, thus a total of 4x4=16 spatial cells.
%                   Has to be greater than 0.
%                   (same as NumSpatialBins*NumSpatialBins in Andrea Vedaldi's sift code)
%                   Default: 16
%
%--------------------------------------------------------------------------------------------------
% Threshold Input
%--------------------------------------------------------------------------------------------------
% stopThresholdsArr array with NumSpatialCells entries of stoping values.
%                   i.e: if the distance after computing the distance for cell number i is greater
%                   or equal to stopThresholdsArr[i] we stop the
%                   computation and return 
%                   stopThresholdsArr(NumSpatialCells)
%                   If 0 it is ignored.
%                   If it is a single value then this array is built:
%                   [ ((1/NumSpatialCells)^0.7)*stopThresholdsArr
%                     ((2/NumSpatialCells)^0.7)*stopThresholdsArr  
%                     ...
%                     stopThresholdsArr ]
%                   If it is two values then this array is built:
%                   [ ((1/NumSpatialCells)^stopThresholdsArr(2))*stopThresholdsArr(1)
%                     ((2/NumSpatialCells)^stopThresholdsArr(2))*stopThresholdsArr(1)  
%                     ...
%                     stopThresholdsArr ]
%                   Default: 0
%
%--------------------------------------------------------------------------------------------------
% For example:
%--------------------------------------------------------------------------------------------------
%             
%                   dists = SiftDist(rand(16*16,500), rand(16*16,300));
%                   dists = SiftDist(rand(8*16,500), rand(8*16,300), 8);
%                   dists = SiftDist(rand(3*1,500), rand(3*1,300), 3, 1);
%                   dists = SiftDist(rand(8*4,500), rand(8*4,300), 8, 4, 0.5);
%                   dists = SiftDist(rand(8*4,500), rand(8*4,300), 8, 4, [0.5 0.8]);
%                   dists = SiftDist(rand(8*4,500), rand(8*4,300), 8, 4, [0.1 0.2 0.3 0.4]);

% Copyright (2008-2009), The Hebrew University of Jerusalem.
% All Rights Reserved.
%
% Created by Ofir Pele
% The Hebrew University of Jerusalem
%
% This software for computing the SiftDist distance between histograms is being
% made available for individual non-profit research use only. Any commercial use
% of this software requires a license from the Hebrew University of Jerusalem.
%
% For further details on obtaining a commercial license, contact Ofir Pele
% (ofirpele@cs.huji.ac.il) or Yissum, the technology transfer company of the
% Hebrew University of Jerusalem.
%
% THE HEBREW UNIVERSITY OF JERUSALEM MAKES NO REPRESENTATIONS OR WARRANTIES OF
% ANY KIND CONCERNING THIS SOFTWARE.
%
% IN NO EVENT SHALL THE HEBREW UNIVERSITY OF JERUSALEM BE LIABLE TO ANY PARTY FOR
% DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST
% PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
% THE THE HEBREW UNIVERSITY OF JERUSALEM HAS BEEN ADVISED OF THE POSSIBILITY OF
% SUCH DAMAGE. THE HEBREW UNIVERSITY OF JERUSALEM SPECIFICALLY DISCLAIMS ANY
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED
% HEREUNDER IS ON AN "AS IS" BASIS, AND THE HEBREW UNIVERSITY OF JERUSALEM HAS NO
% OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
% MODIFICATIONS.
