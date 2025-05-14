%SiftRatioMatch     Returns the symmetric ratio nearest neighbor matching.
%                   
% [inds ratios] = SiftRatioMatch(descr1, frames1,
%                                descr2, frames2,
%
%                                DistRatio, StopThresholdsFactorGamma
%                                
%                                NumOrientBins, NumSpatialBins,
%                                Magnif,
% 
%                                MaxOverlap,
%                                DistType)
% 
%-------------------------------------------------------------------------------------------------
% Output:
%-------------------------------------------------------------------------------------------------             
% inds            A 1 x size(descr1,2) matrix that contains for each feature the index of the 
%                 symmetric nearest neighbor or 0 if there is no such neighbor.
% 
% ratios          A 1 x size(descr1,2) matrix that contains for each feature the ratio 
%                 matching score. Bigger values imply better matching (less ambiguity).
%                 Note: If inds(1,k)==0, then ratios(1,k) is undefined.
%                 Note: If there is no second best neighbor to index 'k' (for example all 
%                 neighbors have an overlap too big) then ratios(k)=inf.
%                 This value is the min( D(a_2,b)/D(a,b) , D(a,b_2)/D(a,b) ) from the paper.
%                 Note: If min( D(a_2,b) , D(a,b_2) )==0 we return 1
%                 (instead of 0/0==NaN).
%                 If DistRatio~=-1 (the default is 0.7) and ratio is probably bigger than the given DistRatio,
%                 ratios will be DistRatio for this entry.
%                 i.e: If given DistRatio is different than -1 (true for
%                 the default) this function does not compute the exact value of
%                 ratios that are probably bigger than the given DistRatio 
%
%--------------------------------------------------------------------------------------------------
% Descriptors Data Input:
%--------------------------------------------------------------------------------------------------
% descr1,         Two matrices that contains the descriptors. 
% descr2          The number of rows is the dimension of the descriptors.
%                 That is the number of rows has to be equal to the given:
%                 NumOrientBins*NumSpatialBins*NumSpatialBins.
%                 The number of columns is the number of descriptors in each set.
%                 (same as descr in Andrea Vedaldi's SIFT code)
%                 Recommendation: send sqrt of the descriptors (see demo_SiftDist)
%
% frames1,        Two 4 x size(descr1,2) matrices that contain the frames of the features.
% frames2         One SIFT frame per column. Its format is:
%                 frames(1:2,k)  center (X,Y) of the frame k,
%                 frames(3,k)    scale SIGMA of the frame k,
%                 frames(4,k)    orientation THETA of the frame k.   
%                 (same as frames in Andrea Vedaldi's SIFT code)
%
%--------------------------------------------------------------------------------------------------
% Matching Input
%--------------------------------------------------------------------------------------------------
% DistRatio       If ratio is smaller (or we found it in the threshold
%                 probabilistic way) to the given DistRatio,
%                 inds will be zero for this entry.
%                 Increase value to decrease the number of matches 
%                 (possibly better matches).
%                 Increasing DistRatio also makes the function work
%                 slower.
%                 If -1, the the threshold probabilistic way is disabled
%                 and all ratios are returned.
%                 Has to be greater or equal to 1 or -1.
%                 Default: 1.25
%
% StopThresholdsFactorGamma If DistRatio~=-1 and
%                           if the distance after computing the distance
%                           for spatial cell number i (out of
%                           NumSpatialBins*NumSpatialBins) is greater
%                           or equal to
%                           (((i-1)/NumSpatialBins*NumSpatialBins)^StopThresholdsFactorGamma)*currentMinVal*DistRatio
%                           the computation is stopped.
%                           This is a heuristic probabilistic way to stop the
%                           computation for distances that are too large.
%                           Larger values means more chance of stopping
%                           fast. i.e: faster running time with higher
%                           probability of error.
%                           Default: 0.7
%
%--------------------------------------------------------------------------------------------------
% Descriptors Layout Input:
%--------------------------------------------------------------------------------------------------
% NumOrientBins   The number of orientation bins.
%                 Original SIFT had 8 bins, but in my paper I found that 16 gives better 
%                 results with the SiftDist distance.
%                 NumOrientBins has to be greater than 1.
%                 (same as NumOrientbins in Andrea Vedaldi's sift code)
%                 Default: 16
%
% NumSpatialBins  The number of spatial bins in each spatial direction. 
%                 Regular SIFT has NumSpatialBins=4, thus a total of 4x4=16 spatial cells.
%                 Has to be greater than 0.
%                 (same as NumSpatialBins in Andrea Vedaldi's sift code)
%                 Default: 4.
%
% Magnif          The Descriptor window magnification.
%                 Has to be greater or equal to one.
%                 (same as Magnif in Andrea Vedaldi's sift code)
%                 Default: 3.
%
%--------------------------------------------------------------------------------------------------
% Matching Params Input:
%--------------------------------------------------------------------------------------------------
% MaxOverlap      The overlap (area of intersection divided by area of union)
%                 between the best and second best matches has to be smaller or equal 
%                 to MaxOverlap. This means that keypoints that are too close are not 
%                 regarded as an estimate to the background population.
%                 This is equal to 1-'overlap error', defined in paper.
%                 Default: 0.5.
%
% DistType        The distance type that is used for the matching.
%                 Possible values:
%                 1 - SiftDist (default)

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

