// Copyright (2008-2009), The Hebrew University of Jerusalem.
// All Rights Reserved.

// Created by Ofir Pele
// The Hebrew University of Jerusalem

// This software for computing the SiftDist distance between histograms is being
// made available for individual non-profit research use only. Any commercial use
// of this software requires a license from the Hebrew University of Jerusalem.

// For further details on obtaining a commercial license, contact Ofir Pele
// (ofirpele@cs.huji.ac.il) or Yissum, the technology transfer company of the
// Hebrew University of Jerusalem.

// THE HEBREW UNIVERSITY OF JERUSALEM MAKES NO REPRESENTATIONS OR WARRANTIES OF
// ANY KIND CONCERNING THIS SOFTWARE.

// IN NO EVENT SHALL THE HEBREW UNIVERSITY OF JERUSALEM BE LIABLE TO ANY PARTY FOR
// DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST
// PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
// THE THE HEBREW UNIVERSITY OF JERUSALEM HAS BEEN ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE. THE HEBREW UNIVERSITY OF JERUSALEM SPECIFICALLY DISCLAIMS ANY
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED
// HEREUNDER IS ON AN "AS IS" BASIS, AND THE HEBREW UNIVERSITY OF JERUSALEM HAS NO
// OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
// MODIFICATIONS.

#include <mex.h>
#include "mexCheckAndExtractInputs.hxx"
#include "SiftDist.hxx"

void mexFunction(int nout, mxArray *out[], 
                 int nin, const mxArray *in[]) {
      
      //-------------------------------------------------------
      // Consts + defaults
      //-------------------------------------------------------
      unsigned int NBO= 16;
      unsigned int CELLS_NUM= 16;
      double* stopThresholdsArr= NULL;
      //-------------------------------------------------------
      
      //-------------------------------------------------------
      // Check the arguments
      //-------------------------------------------------------
      if (nin<2) {
	    mexErrMsgTxt("At least 2 arguments are required.");
      }
      if (nout>1) {
	    mexErrMsgTxt("Too many output arguments.");
      }
      //-------------------------------------------------------
      
      
      //-------------------------------------------------------
      // Get inputs
      //-------------------------------------------------------
      const double* descr1;
      const double* descr2;
      unsigned int sift_num1;
      unsigned int sift_num2;

      if (nin>2) {
	    mexCheckAndExtractInputs::checkAndExtract_NBO(in[2],
                                                      NBO);
      }

      if (nin>3) {
	    mexCheckAndExtractInputs::checkAndExtract_CELLS_NUM(in[3],
                                                            CELLS_NUM);
      }
      
      mexCheckAndExtractInputs::checkAndExtract_descrs_sifts_num_2(in[0], in[1],
                                                                   NBO, CELLS_NUM,
                                                                   
                                                                   descr1, descr2, 
                                                                   sift_num1, sift_num2);

      if (nin>4) {
          stopThresholdsArr= new double[CELLS_NUM];
          mexCheckAndExtractInputs::checkAndExtract_stopThresholdsArr(in[4],
                                                                      CELLS_NUM,
                                                                      
                                                                      stopThresholdsArr);
      }
      //-------------------------------------------------------

      //-------------------------------------------------------
      // Create and fill output
      //-------------------------------------------------------
      SiftDist<double> sd(NBO, CELLS_NUM);

      out[0]= mxCreateDoubleMatrix(sift_num1, sift_num2, mxREAL);
      double* dists= (double*)mxGetData(out[0]);
      unsigned int SIFT_DIM= NBO*CELLS_NUM;
      for(unsigned int j=0,f2=0,j_mult_sift_num1=0; j<sift_num2; ++j,f2+=SIFT_DIM,j_mult_sift_num1+=sift_num1){
	    for(unsigned int i=0,f1=0; i<sift_num1; ++i,f1+=SIFT_DIM){
		  (*(dists+j_mult_sift_num1+i))= sd(descr1+f1, descr2+f2, stopThresholdsArr);
	    }
      }
      //-------------------------------------------------------

      delete[] stopThresholdsArr;
            
} // mexFunction
