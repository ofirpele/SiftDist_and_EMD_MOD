%[dist F]= EMD_MOD(P, Q)
% Computes the EMD_MOD between two equal mass histograms.
% That is, EMD with modulo L1 as the ground distance,
% between histograms that their sums should be equal.
% The algorithm was described in Appendix A of the paper:
%  A Linear Time Histogram Metric for Improved SIFT Matching
%  Ofir Pele, Michael Werman
%  ECCV 2008
% Please cite the paper if you use this code.
% Two things to take into considerations:
% 1. For SIFT matching, this distance performance was bad.
%    Try using my QC / FastEMD / SIFT_DIST codes
%    (from ECCV 2010, ICCV 2009, ECCV 2008 respectively).
% 2. This implementation is linear on the average as I use
%    the nth_element function of standard C++ library,
%    which is currently average linear time. There are worst time
%    linear algorithms, but I did not implement them.
%
% Returns:
% dist - the EMD_MOD distance between the histograms
% F - Sparse matrix of the flow between P bins to Q bins.
%     It is sparse to keep the linear time constraint.
%     Note that if it is not requested the running time is shorter 
%     (although not in order of magnitude)
% 
% Params:
% P - first histogram.
% Q - second histograms.


% Copyright (c) 2011, Ofir Pele
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met: 
%    * Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%    * Neither the name of the The Hebrew University of Jerusalem nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.


% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
