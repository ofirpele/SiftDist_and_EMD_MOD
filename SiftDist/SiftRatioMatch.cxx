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
#include "SiftRatioMatchImpl.hxx"

// Distances used:
#include "SiftDist.hxx"


void mexFunction(int nout, mxArray *out[], 
                 int nin, const mxArray *in[]) {


      //-------------------------------------------------------
      // Consts + defaults
      //-------------------------------------------------------
      const int FRAMES_COL_SIZE= 4;
      const int FRAMES_X_IND= 0;
      const int FRAMES_Y_IND= 1;
      const int FRAMES_SCALE_IND= 2;

      double distRatio= 1.25;
      double stopThresholdsFactorGamma= 0.7;

      unsigned int NBO= 16;
      unsigned int NBP= 4; 
      double Magnif= 3.0;
      
      distType dt= SiftDistType;
      double maxOverlap= 0.5;
      //-------------------------------------------------------

      
      //-------------------------------------------------------
      // Check the arguments
      //-------------------------------------------------------
      if (nin<4) {
	    mexErrMsgTxt("At least 4 arguments are required.");
      }
      if (nout!=2) {
	    mexErrMsgTxt("There should be exactly two output arguments.");
      }
      //-------------------------------------------------------

      //-------------------------------------------------------
      // Get inputs
      //-------------------------------------------------------
      const double* descr1;
      const double* descr2;
      unsigned int sift_num1;
      unsigned int sift_num2;

      
      if (nin>4) {
          mexCheckAndExtractInputs::checkAndExtract_distRatio(in[4],
                                                              distRatio);
      }

      if (nin>5) {
          mexCheckAndExtractInputs::checkAndExtract_stopThresholdsFactorGamma(in[5],
                                                                              stopThresholdsFactorGamma);
      }
            
      if (nin>6) {
          mexCheckAndExtractInputs::checkAndExtract_NBO(in[6],
                                                        NBO);
      }

      if (nin>7) {
	    mexCheckAndExtractInputs::checkAndExtract_NBP(in[7],
							  NBP);
      }
      	    
      mexCheckAndExtractInputs::checkAndExtract_descrs_sifts_num(in[0], in[2],
								 NBO, NBP,
								 
								 descr1, descr2, 
								 sift_num1, sift_num2);
      
      const double* frames1;
      const double* frames2;
      mexCheckAndExtractInputs::checkAndExtract_frames(in[1], sift_num1,
						       in[3], sift_num2,
						       FRAMES_COL_SIZE,
						       
						       frames1, frames2);
      
      if (nin>8) {
	    mexCheckAndExtractInputs::checkAndExtract_Magnif(in[8],
							     Magnif);
      }
      
      if (nin>9) {
	    mexCheckAndExtractInputs::checkAndExtract_maxOverlap(in[9],
								 maxOverlap);
      }
      
      if (nin>10) {
	    mexCheckAndExtractInputs::checkAndExtract_distType(in[10],
							       dt);
      }
      //------------------------------------------------------- 
      
      
      //-------------------------------------------------------
      // Create and fill output
      //-------------------------------------------------------
      out[0]= mxCreateDoubleMatrix(1, sift_num1, mxREAL);
      out[1]= mxCreateDoubleMatrix(1, sift_num1, mxREAL);
      double* inds=  (double*)mxGetData(out[0]);
      double* ratios= (double*)mxGetData(out[1]);
      
      SiftDist<double> sd(NBO,NBP*NBP);
            
      switch(dt) {

      case SiftDistType:
          
          SiftRatioMatchImpl< SiftDist<double> >
              (descr1, frames1, sift_num1,
               descr2, frames2, sift_num2,
               distRatio,
               stopThresholdsFactorGamma,
               NBO, NBP, Magnif,
               maxOverlap,
               sd,
               FRAMES_COL_SIZE, FRAMES_X_IND, FRAMES_Y_IND, FRAMES_SCALE_IND,
               inds, ratios);
          break;
          
      } // switch(dt)
      //-------------------------------------------------------
      
      
} // end mexFunction
