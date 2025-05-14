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

#ifndef __OFIRPELE_MEX_CHECK_AND_EXTRACT_INPUTS_HXX
#define __OFIRPELE_MEX_CHECK_AND_EXTRACT_INPUTS_HXX

#include <mex.h>
#include <math.h> // for pow

enum distType {SiftDistType};


class mexCheckAndExtractInputs {
public:
    
static void checkAndExtract_stopThresholdsFactorGamma
(const mxArray* in_stopThresholdsFactorGamma,
  
 double& out_stopThresholdsFactorGamma) {
    if ( !( (mxIsDouble(in_stopThresholdsFactorGamma))&&(!mxIsComplex(in_stopThresholdsFactorGamma)) ) ) {
        mexErrMsgTxt("stopThresholdsFactorGamma should be regular double.");
    }
    if (mxGetM(in_stopThresholdsFactorGamma)!=1) {
        mexErrMsgTxt("stopThresholdsFactorGamma should be a single number.");
    }
    unsigned int len_stopThresholdsFactorGamma= mxGetN(in_stopThresholdsFactorGamma);
    if ( !((len_stopThresholdsFactorGamma==1))) {
        mexErrMsgTxt("stopThresholdsFactorGamma should be a single number.");
    }

    double* tmp_out= static_cast<double*>( (mxGetData(in_stopThresholdsFactorGamma)) );
    out_stopThresholdsFactorGamma= tmp_out[0];
    
} // end checkAndExtract_stopThresholdsFactorGamma

static void checkAndExtract_stopThresholdsArr
(const mxArray* in_stopThresholdsArr,
 unsigned int CELLS_NUM,
 
 double*& out_stopThresholdsArr) {

    
    if ( !( (mxIsDouble(in_stopThresholdsArr))&&(!mxIsComplex(in_stopThresholdsArr)) ) ) {
        mexErrMsgTxt("stopThresholdsArr should be regular double.");
    }
    if (mxGetM(in_stopThresholdsArr)!=1) {
        mexErrMsgTxt("stopThresholdsArr should have only one row (either a single number or a row vector).");
    }
    
    unsigned int len_stopThresholdsArr= mxGetN(in_stopThresholdsArr);
    if ( !((len_stopThresholdsArr==1)||(len_stopThresholdsArr==2)||(len_stopThresholdsArr==CELLS_NUM)) ) {
        mexErrMsgTxt("stopThresholdsArr should be either a single number or two numbers or a row vector with CELLS_NUM columns.");
    }
    
    double* tmp_out_stopThresholdsArr= static_cast<double*>( (mxGetData(in_stopThresholdsArr)) );
    if (len_stopThresholdsArr==1||len_stopThresholdsArr==2) {
        double stopThreshold= tmp_out_stopThresholdsArr[0];
        if (stopThreshold<0) {
            mexErrMsgTxt("stopThresholdsArr(1) should be greater than zero");
        }
        double stopGammaFactor= 0.7;
        if (len_stopThresholdsArr==2) {
            stopGammaFactor= tmp_out_stopThresholdsArr[1];
            if (stopGammaFactor<0) {
                mexErrMsgTxt("stopThresholdsArr(2) should be greater than zero");
            }
        }
        double factor= 1.0/CELLS_NUM;
        for (unsigned int i=0; i<CELLS_NUM; ++i) {
            out_stopThresholdsArr[i]= pow(factor,stopGammaFactor)*stopThreshold;
            factor+= 1.0/CELLS_NUM;
        }
    } else {
        for (unsigned int i=0; i<CELLS_NUM; ++i) {
            if (tmp_out_stopThresholdsArr[i]<0) {
                mexErrMsgTxt("stopThresholdsArr should be greater than zero");
            }
            if (i>0) {
                if (tmp_out_stopThresholdsArr[i]<tmp_out_stopThresholdsArr[i-1]) {
                    mexErrMsgTxt("stopThresholdsArr should be monotonic increasing");
                }
            }
            out_stopThresholdsArr[i]= tmp_out_stopThresholdsArr[i];
        }
    }
    
} // end checkAndExtract_stopThresholdsArr
    
    
static void checkAndExtract_distRatio
(const mxArray* in_distRatio,
 
 double&  out_distRatio) {

    if ( !( (mxIsDouble(in_distRatio))&&(!mxIsComplex(in_distRatio)) ) ) {
        mexErrMsgTxt("distRatio should be regular double.");
    }
    
    out_distRatio= static_cast<double>( (*mxGetPr(in_distRatio)) );
    if ((out_distRatio<1)&&(out_distRatio!=-1)) {
	    mexErrMsgTxt("distRatio should be greater or equal to one.");
    }
} // end checkAndExtract_distRatio
    
static void checkAndExtract_NBO
(const mxArray* in_NBO,
 
 unsigned int&  out_NBO) {

    out_NBO= static_cast<unsigned int>( (*mxGetPr(in_NBO)) );
    if (out_NBO<=1) {
	    mexErrMsgTxt("NumOrientBins has to be greater than one.");
    }
} // end checkAndExtract_NBO

static void checkAndExtract_NBP
(const mxArray* in_NBP,
 
 unsigned int& out_NBP) {
            
      out_NBP= static_cast<unsigned int>( (*mxGetPr(in_NBP)) );
      if (out_NBP<1) {
	    mexErrMsgTxt("NumSpatialBins has to be greater than zero.");
      }

} // end checkAndExtract_NBP

static void checkAndExtract_CELLS_NUM
(const mxArray* in_CELLS_NUM,
 
 unsigned int& out_CELLS_NUM) {
            
      out_CELLS_NUM= static_cast<unsigned int>( (*mxGetPr(in_CELLS_NUM)) );
      if (out_CELLS_NUM<1) {
	    mexErrMsgTxt("NumSpatialCells has to be greater than zero.");
      }

} // end checkAndExtract_CELLS_NUM

static void checkAndExtract_descrs_sifts_num
(const mxArray* in_descr1, const mxArray* in_descr2,
 unsigned int NBO, unsigned int NBP,
 
 const double*& out_descr1, const double*& out_descr2,
 unsigned int&  out_sift_num1, unsigned int& out_sift_num2) {

    if ( !( (mxIsDouble(in_descr1))&&(!mxIsComplex(in_descr1))&&(mxGetNumberOfDimensions(in_descr1)==2) ) ) {
        mexErrMsgTxt("descr1 should be a matrix of doubles.");
    }

    if ( !( (mxIsDouble(in_descr2))&&(!mxIsComplex(in_descr2))&&(mxGetNumberOfDimensions(in_descr2)==2) ) ) {
        mexErrMsgTxt("descr2 should be a matrix of doubles.");
    }
    
    out_descr1= static_cast<const double*>( mxGetData(in_descr1) );
    out_descr2= static_cast<const double*>( mxGetData(in_descr2) );
    
    const int* dims1 = mxGetDimensions(in_descr1);
    const int* dims2 = mxGetDimensions(in_descr2);
    if (dims1[0]!=dims2[0]) {
	    mexErrMsgTxt("Dimension of sift descriptors (number of rows) should be the same in both matrices.");
    }
    if (dims1[0]!=static_cast<int>(NBO*NBP*NBP)) {
	    if ((dims1[0]==128)&&(NBO==16)&&(NBP==4)) {
            mexErrMsgTxt("Dimension of sift descriptors (number of rows) should be equal to NumOrientBins*NumSpatialBins*NumSpatialBins.\n"
                         "Note that you are using the SIFT that uses 8 orientation bins, that is the regular 8x4x4=128 dimensional SIFT,"
                         "while your params are indicating the 16 orientation bins, that is a 16x4x4=256 SIFT, which works better with the SiftDist.\n"
                         "Try making a SIFT with 16 orientation bins."); 
        }
	    mexErrMsgTxt("Dimension of sift descriptors (number of rows) should be equal to NumOrientBins*NumSpatialBins*NumSpatialBins.");
    }
    
    out_sift_num1= static_cast<unsigned int>(dims1[1]);
    out_sift_num2= static_cast<unsigned int>(dims2[1]);
	
} // end checkAndExtract_descrs_sifts_num
//----------------------------------------------------------------------------

static void checkAndExtract_descrs_sifts_num_2
(const mxArray* in_descr1, const mxArray* in_descr2,
 unsigned int NBO, unsigned int CELLS_NUM,
 
 const double*& out_descr1, const double*& out_descr2,
 unsigned int&  out_sift_num1, unsigned int& out_sift_num2) {

    if ( !( (mxIsDouble(in_descr1))&&(!mxIsComplex(in_descr1))&&(mxGetNumberOfDimensions(in_descr1)==2) ) ) {
        mexErrMsgTxt("descr1 should be a matrix of doubles.");
    }

    if ( !( (mxIsDouble(in_descr2))&&(!mxIsComplex(in_descr2))&&(mxGetNumberOfDimensions(in_descr2)==2) ) ) {
        mexErrMsgTxt("descr2 should be a matrix of doubles.");
    }
    
    out_descr1= static_cast<const double*>( mxGetData(in_descr1) );
    out_descr2= static_cast<const double*>( mxGetData(in_descr2) );
    
    const int* dims1 = mxGetDimensions(in_descr1);
    const int* dims2 = mxGetDimensions(in_descr2);
    if (dims1[0]!=dims2[0]) {
	    mexErrMsgTxt("Dimension of sift descriptors (number of rows) should be the same in both matrices.");
    }
    if (dims1[0]!=static_cast<int>(NBO*CELLS_NUM)) {
	    if ((dims1[0]==128)&&(NBO==16)&&(CELLS_NUM==16)) {
            mexErrMsgTxt("Dimension of sift descriptors (number of rows) should be equal to NumSpatialCells*NumOrientBins.\n"
                         "Note that you are using the SIFT that uses 8 orientation bins, that is the regular 8x4x4=128 dimensional SIFT,"
                         "while your params are indicating the 16 orientation bins, that is a 16x4x4=240 SIFT, which works better with the SiftDist.\n"
                         "Try making a SIFT with 16 orientation bins."); 
        }
	    mexErrMsgTxt("Dimension of sift descriptors (number of rows) should be equal to NumOrientBins*NumSpatialCells.");
    }
    
    out_sift_num1= static_cast<unsigned int>(dims1[1]);
    out_sift_num2= static_cast<unsigned int>(dims2[1]);
	    
} // end checkAndExtract_descrs_sifts_num_2
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
static void checkAndExtract_frames(const mxArray* in_frames1, unsigned int sift_num1,
				   const mxArray* in_frames2, unsigned int sift_num2,
				   int FRAMES_COL_SIZE,
				   
				   const double*& out_frames1, const double*& out_frames2) {

    if ( !( (mxIsDouble(in_frames1))&&(!mxIsComplex(in_frames1))&&(mxGetNumberOfDimensions(in_frames1)==2) ) ) {
        mexErrMsgTxt("frames1 should be a matrix of doubles.");
    }

    if ( !( (mxIsDouble(in_frames2))&&(!mxIsComplex(in_frames2))&&(mxGetNumberOfDimensions(in_frames2)==2) ) ) {
        mexErrMsgTxt("frames2 should be a matrix of doubles.");
    }
    
    out_frames1= static_cast<const double*>( mxGetData(in_frames1) );
    out_frames2= static_cast<const double*>( mxGetData(in_frames2) );
    
    const int* dims1 = mxGetDimensions(in_frames1);
    const int* dims2 = mxGetDimensions(in_frames2);
    if (dims1[0]!=FRAMES_COL_SIZE||dims2[0]!=FRAMES_COL_SIZE) {
	    mexErrMsgTxt("Dimension of sift frames (number of rows) should be FRAMES_COL_SIZE (4).");
    }
    if (dims1[1]!=static_cast<int>(sift_num1)) {
	    mexErrMsgTxt("Number of sift frames in frames1 (number of columns) should be the same as the number of descriptors in descr1 (number of columns).");
    }
    if (dims2[1]!=static_cast<int>(sift_num2)) {
	    mexErrMsgTxt("Number of sift frames in frames2 (number of columns) should be the same as the number of descriptors in descr2 (number of columns).");
    }
    
} // end checkAndExtract_frames
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
static void checkAndExtract_maxOverlap(const mxArray* in_maxOverlap,
					 
				       double& out_maxOverlap) {

    if ( !( (mxIsDouble(in_maxOverlap))&&(!mxIsComplex(in_maxOverlap)) ) ) {
        mexErrMsgTxt("maxOverlap should be regular double.");
    }
    
    out_maxOverlap= static_cast<double>( (*mxGetPr(in_maxOverlap)) );
    if (out_maxOverlap<0) {
	    mexErrMsgTxt("MaxOverlap has to be greater or equal to zero.");
    }
    if (out_maxOverlap>1) {
	    mexErrMsgTxt("MaxOverlap has to be smaller or equal to one.");
    }
    
} // end checkAndExtract_maxOverlap
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
static void checkAndExtract_Magnif(const mxArray* in_Magnif,
					 
				   double& out_Magnif) {

    if ( !( (mxIsDouble(in_Magnif))&&(!mxIsComplex(in_Magnif)) ) ) {
        mexErrMsgTxt("Magnif should be regular double.");
    }
    
    out_Magnif= static_cast<double>( (*mxGetPr(in_Magnif)) );
    if (out_Magnif<1) {
	    mexErrMsgTxt("Magnif has to be greater or equal to one.");
    }
    
} // end checkAndExtract_Magnif
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------


static void checkAndExtract_distType(const mxArray* in_dt,
					 
				     distType& out_dt) {

    unsigned int ui_dt= static_cast<unsigned int>( (*mxGetPr(in_dt)) );
    switch (ui_dt) {
    case 1:
        out_dt= SiftDistType;
        break;
    default:
        mexErrMsgTxt("DistType should be 1 (SiftDist)");
    }
    
} // end checkAndExtract_distType
//----------------------------------------------------------------------------
};

#endif
