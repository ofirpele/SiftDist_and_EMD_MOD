#ifndef EMD_MOD_HPP_
#define EMD_MOD_HPP_

#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h> // abs
#include <list>

//-------------------------------------------------------------
struct CompByVec {

	const std::vector<double>& _V;
    
	CompByVec(const std::vector<double>& V) : _V(V) { }
	
	bool operator()(int i, int j) const {
	    return _V[i] > _V[j];
	}
	
};
//-------------------------------------------------------------


/// Computes the EMD_MOD between two equal mass histograms.
/// That is, EMD with modulo L1 as the ground distance,
/// between histograms that their sums should be equal.
/// The algorithm was described in Appendix A of the paper:
///  A Linear Time Histogram Metric for Improved SIFT Matching
///  Ofir Pele, Michael Werman
///  ECCV 2008
/// Please cite the paper if you use this code.
/// Two things to take into considerations:
/// 1. For SIFT matching, this distance performance was bad.
///    Try using my QC / FastEMD / SIFT_DIST codes
///    (from ECCV 2010, ICCV 2009, ECCV 2008 respectively).
/// 2. This implementation is linear on the average as I use
///    the nth_element function of standard C++ library,
///    which is currently average linear time. There are worst time
///    linear algorithms, but I did not implement them.
///
/// Params:
/// oA - first histogram
/// oB - second histograms
/// NBO - number of bins
/// flows - If computeFlow is true, this is the pointer to a vector that will be filled
///         with lists of flows from each bin (each pair is the going to bin and how much).
///         Note: the vector is cleared before being filled.
template<bool computeFlow>
double EMD_MOD(const double* oA, const double* oB, int NBO,
			   std::vector< std::list< std::pair<int,double> > >* flows_ptr=NULL) {
	
	double _emd= 0.0;
	
	std::vector<double> A(oA, oA+NBO);
	std::vector<double> B(oB, oB+NBO);
	std::vector<double> F(NBO);
	std::vector<int> FI(NBO);
	std::vector<int> I(NBO);
            
	// CA and CB can be removed with an optimization
	std::vector<double> CA(oA, oA+NBO);
	std::vector<double> CB(oB, oB+NBO);
	F[0]= CA[0]-CB[0];
	for (int i= 1; i<NBO; ++i) {
		CA[i]= CA[i]+CA[i-1];
		CB[i]= CB[i]+CB[i-1];	    
		F[i]= CA[i]-CB[i];
	}
	  
	for (int i= 0; i<NBO; ++i) {
		FI[i]= i;
	}
	// On average, linear in NBO/2
	// There are strictly linear time algorithms, which might be faster, or not... :)
	nth_element(FI.begin(), FI.begin()+(NBO/2), FI.end(), CompByVec(F));
	
	int i= (FI[NBO/2] + 1) % NBO;
	for (int t= 0; t<NBO; ++t) {
	    I[t]= i;
	    ++i;
	    i= i%NBO;
	}
      
	int tA=0;
	int tB=0;
	int iA=I[tA];
	int iB=I[tB];
	if (computeFlow) {
		flows_ptr->clear();
		flows_ptr->resize(NBO);
	}
	while (true) {

	    while (A[iA]==0.0) {
			if (++tA==NBO) {
				return _emd;
			}
			iA=I[tA];
	    }
		while (B[iB]==0.0) {
			if (++tB==NBO) {
				return _emd;
			}
			iB=I[tB];
	    }

	    double f= std::min(A[iA],B[iB]);
	    A[iA]-= f;
		B[iB]-= f;
		_emd+= f*std::min(abs(iA-iB), NBO-abs(iA-iB));

		if (computeFlow) {
			std::vector< std::list< std::pair<int,double> > >& flows= *flows_ptr;
			flows[iB].push_back( std::make_pair(iA, f) );
		}
		
	} // true (flowing)

} // end EMD_MOD

#endif

// Copyright (c) 2011, Ofir Pele
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//    * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//    * Neither the name of the The Hebrew University of Jerusalem nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

