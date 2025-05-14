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

#ifndef _OFIRPELE_SIFT_RATIO_MATCH_IMPL__HXX
#define _OFIRPELE_SIFT_RATIO_MATCH_IMPL__HXX

#include "circleFuncs.hxx"
#include <vector>
#include <limits>
#include <math.h>




///@param DISTANCE_T the type of the functor that computes distances between SIFT-like descriptors.
/// inds and ratios are assumed to be initialized to 0.
template<typename DISTANCE_T>
class SiftRatioMatchImpl {

public:
      
SiftRatioMatchImpl(const double* descr1, const double* frames1, unsigned int sift_num1, 
                   const double* descr2, const double* frames2, unsigned int sift_num2,
                   double distRatio,
                   double stopThresholdsFactorGamma,
                   unsigned int NBO, unsigned int NBP, double Magnif,
                   double maxOverlap,
                   const DISTANCE_T& sd,
                   int FRAMES_COL_SIZE, int FRAMES_X_IND, int FRAMES_Y_IND, int FRAMES_SCALE_IND,
                   
                   double* inds, double* ratios) {

    unsigned int CELLS_NUM= NBP*NBP;

    if (distRatio==-1) {
        _stopThresholdsArr= NULL;
        _stopThresholdsFactorsArr= NULL;
    } else {
        _stopThresholdsArr= new double[NBP*NBP];
        _stopThresholdsFactorsArr= new double[NBP*NBP];
        
        unsigned int i;
        double factor= 1.0/(CELLS_NUM);
        for (i=0; i<CELLS_NUM; ++i) {
            _stopThresholdsFactorsArr[i]= pow(factor, stopThresholdsFactorGamma);
            factor+= (1.0/CELLS_NUM);
        }
    }
    
    double scaleToRadiusFactor= (Magnif*NBP)/2.0;
    
    std::vector<double> radius1_vec(sift_num1);
    extractRadiusesFromFrames(frames1, FRAMES_SCALE_IND, FRAMES_COL_SIZE, scaleToRadiusFactor,
                              radius1_vec);
    
    std::vector<double> radius2_vec(sift_num2);
    extractRadiusesFromFrames(frames2, FRAMES_SCALE_IND, FRAMES_COL_SIZE, scaleToRadiusFactor,
                              radius2_vec);
    
    
    unsigned int sift_dim= NBO*CELLS_NUM;
    
    // Holds distances of desc1(:,c1) to descr2(:,1:end)
    std::vector<double> descr1_c1_descr2_all_dists = std::vector<double>(sift_num2);
    // Holds distances of desc2(:,min_c1) to descr1(:,1:end)
    std::vector<double> descr2_min_c1_descr1_all_dists= std::vector<double>(sift_num1);
    
    
    const double* descr1_c1= descr1;
    for (unsigned int c1=0; c1<sift_num1; ++c1, descr1_c1+= sift_dim, ++inds, ++ratios) {
        
	    unsigned int min_descr1_c1_descr2_all_dists_Ind;
	    computeDistsAndFindMin(descr2, sd, descr1_c1, NBO, CELLS_NUM, sift_dim, sift_num2, distRatio,
                               descr1_c1_descr2_all_dists,
                               min_descr1_c1_descr2_all_dists_Ind);
	    const double* descr2_min= descr2 + (min_descr1_c1_descr2_all_dists_Ind*sift_dim);
        
	    unsigned int min_descr2_min_c1_descr1_all_dists_Ind;
	    computeDistsAndFindMin(descr1, sd, descr2_min, NBO, CELLS_NUM, sift_dim, sift_num1, distRatio,
                               descr2_min_c1_descr1_all_dists,
                               min_descr2_min_c1_descr1_all_dists_Ind);
        
        
	    // If not a symmetric nearest neighbor - *inds and *ratios will be 0
	    if (min_descr2_min_c1_descr1_all_dists_Ind!=c1) continue;
	    
	    double min2_in_descr1;
	    findMin2(sift_dim, sift_num1,
                 descr2_min_c1_descr1_all_dists, min_descr2_min_c1_descr1_all_dists_Ind,
                 frames1, FRAMES_COL_SIZE, FRAMES_X_IND, FRAMES_Y_IND,
                 radius1_vec,
                 maxOverlap,
                 
                 min2_in_descr1);
        
	    
	    double min2_in_descr2;
	    findMin2(sift_dim, sift_num2,
                 descr1_c1_descr2_all_dists, min_descr1_c1_descr2_all_dists_Ind,
                 frames2, FRAMES_COL_SIZE, FRAMES_X_IND, FRAMES_Y_IND,
                 radius2_vec,
                 maxOverlap,
                 
                 min2_in_descr2);
        
	    
	    assert(min2_in_descr1>=descr2_min_c1_descr1_all_dists[c1]);
	    assert(min2_in_descr2>=descr2_min_c1_descr1_all_dists[c1]);
        
	    //------------------------------------------
	    // Fill output
	    //------------------------------------------
	    
	    // +1 to be a Matlab index
	    *inds= min_descr1_c1_descr2_all_dists_Ind + 1;  
        
	    if (min2_in_descr1<min2_in_descr2) {
            if (min2_in_descr1==0) {
                *ratios= 1.0;
            } else {
                *ratios= (min2_in_descr1 / descr2_min_c1_descr1_all_dists[c1]);
            }
	    } else {
		  if (min2_in_descr2==0) {
              *ratios= 1.0;
		  } else {
              *ratios= (min2_in_descr2 / descr2_min_c1_descr1_all_dists[c1]);
		  }
	    }
        
        if ((distRatio!=-1)&&((*ratios)<distRatio)) {
            *inds= 0;
            *ratios= 0;
        }
	    //------------------------------------------
	    
      } // for c1
      
} // end Ctor

~SiftRatioMatchImpl() {
    delete[] _stopThresholdsArr;
    delete[] _stopThresholdsFactorsArr;
}
      
private:

void extractRadiusesFromFrames(const double* frames, int FRAMES_SCALE_IND, int FRAMES_COL_SIZE, double scaleToRadiusFactor,
                               std::vector<double>& radius_vec) {
      frames+= FRAMES_SCALE_IND;
      for (unsigned int i=0; i<radius_vec.size(); ++i,frames+=FRAMES_COL_SIZE) {
	    radius_vec[i]= (*frames) * scaleToRadiusFactor;
      } // for i
} // extractRadiusesFromFrames

void updateStopThresholdsArr(double minVal, double distRatio, unsigned int CELLS_NUM) {
    if (_stopThresholdsArr!=NULL) {
        for (unsigned int i=0; i<CELLS_NUM; ++i) {
            _stopThresholdsArr[i]= _stopThresholdsFactorsArr[i] * minVal * distRatio;
        }
    }
}
        
void computeDistsAndFindMin(const double* descr,
                            DISTANCE_T sd,
                            const double* descr_fixed,
                            unsigned int NBO, unsigned int CELLS_NUM, unsigned int sift_dim,
                            unsigned int sift_num,
                            double distRatio,
                            
                            std::vector<double>& descr_fixed_otherdescr_all_dists,
                            unsigned int& min_Ind) {

    const double* descr_i= descr;
    descr_fixed_otherdescr_all_dists[0]= sd(descr_i, descr_fixed);
    descr_i+= sift_dim;
    min_Ind= 0;

    
    updateStopThresholdsArr(descr_fixed_otherdescr_all_dists[min_Ind], distRatio, CELLS_NUM);
    
    for (unsigned int i=1; i<sift_num; ++i, descr_i+= sift_dim) {
        
        descr_fixed_otherdescr_all_dists[i]= sd(descr_i, descr_fixed, _stopThresholdsArr);
	    if (descr_fixed_otherdescr_all_dists[i] <
            descr_fixed_otherdescr_all_dists[min_Ind]) {
            min_Ind= i;
            updateStopThresholdsArr(descr_fixed_otherdescr_all_dists[min_Ind], distRatio, CELLS_NUM);
        }
    } // i

} // end computeDistsAndFindMin

      
void findMin2(unsigned int sift_dim, unsigned int sift_num,
              const std::vector<double>& dists, unsigned int min_Ind,
              const double* frames, int FRAMES_COL_SIZE, int FRAMES_X_IND, int FRAMES_Y_IND,
              const std::vector<double>& radius_vec, 
              double maxOverlap,
              
              double& min2) {
      
      assert( sift_num==radius_vec.size() );
      assert (std::numeric_limits<double>::has_infinity);

      const double min_x= *((frames+FRAMES_X_IND) + min_Ind*FRAMES_COL_SIZE);
      const double min_y= *((frames+FRAMES_Y_IND) + min_Ind*FRAMES_COL_SIZE);
      const double min_r= radius_vec[min_Ind];
      
      min2= std::numeric_limits<double>::infinity();
      for (unsigned int i=0; i<sift_num; ++i,frames+=FRAMES_COL_SIZE) {
	    
	    if (i==min_Ind) continue;
	    
	    if ( dists[i]<min2 ) {

		  // It's not a neighbor
		  if (circleFuncs::circleOverlap
		      (frames[FRAMES_X_IND], frames[FRAMES_Y_IND], radius_vec[i],
		       min_x, min_y, min_r)<=maxOverlap) {
			min2= dists[i];
		  }
		  
	    } // < min2
	    
      } // i
      
} // findMin2
    
    double* _stopThresholdsArr;
    double* _stopThresholdsFactorsArr;

}; // end class SiftRatioMatchImpl
      
#endif
