#include <mex.h>
#include "EMD_MOD.hpp"
#include "OP_mex_utils.hxx"

void mexFunction(int nout, mxArray *out[], 
                 int nin, const mxArray *in[]) {
	
	if (nin!=2) {
		mexErrMsgTxt("2 arguments are required.");
	}
	if (nout>2) {
		mexErrMsgTxt("Too many output arguments.");
	}
	
	const double* P= static_cast<const double*>( mxGetData(in[0]) );
	const double* Q= static_cast<const double*>( mxGetData(in[1]) );
	size_t NBO= OP_mex_utils::getLength(in[0]);	

#ifndef NDEBUG
	if (OP_mex_utils::getLength(in[1])~=NBO) {
		mexErrMsgTxt("Input vectors do not have the same length.");
	}
	if (NBO%2!=0) {
		mexErrMsgTxt("I did not check if it works with odd length vectors, so use even length, or check :)");
	}
	double sumP= 0.0;
	double sumQ= 0.0;
	for (size_t i= 0; i<NBO; ++i) {
		sumP+= P[i];
		sumQ+= Q[i];
	}
	if (sumP!=sumQ) {
		mexErrMsgTxt("Only works when the sum of the two vectors are equal.");
	}
#endif
	
	out[0]= mxCreateDoubleMatrix(1, 1, mxREAL);
	double* dist= (double*)mxGetData(out[0]);

	if (nout==1) {
		(*dist)= EMD_MOD<false>(P, Q, NBO);
	} else {
		std::vector< std::list< std::pair<int,double> > > flows;

		(*dist)= EMD_MOD<true>(P, Q, NBO, &flows);
		
		int nzmax= 2*NBO;
		out[1] = mxCreateSparse(NBO, NBO, nzmax, mxREAL);

		double* sr= mxGetPr(out[1]);
		mwIndex* irs= mxGetIr(out[1]);
		mwIndex* jcs= mxGetJc(out[1]);
		mwIndex sparseInd= 0;
		for (mwIndex c= 0; c<NBO; ++c) {
            jcs[c]= sparseInd;
			for (std::list< std::pair<int,double> >::const_iterator it= flows[c].begin();
                 it!= flows[c].end();
                 ++it) {
                irs[sparseInd]= it->first; 
                sr[sparseInd]= it->second;
                ++sparseInd;
            }
		}
		jcs[NBO]= sparseInd;

	} // needed to compute flows
    
}


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
