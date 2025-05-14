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

#ifndef _OFIRPELE_SIFT_DIST__HXX
#define _OFIRPELE_SIFT_DIST__HXX

#include "NumTypeZero.hxx"
#include "MyMinMax.hxx"
#include <cassert>
// For the cyclic edge special case:
#include <vector>


/// This class has operator() which returns the SiftDist
/// as defined in the paper:
/// A Linear Time Histogram Metric for Improved SIFT Matching
/// Ofir Pele, Michael Werman
/// ECCV 2008
///
///@param NUM_T the type of numbers. Must be signed. Its zero
/// value must be defined in NumTypeZero.hxx
template<typename NUM_T>
class SiftDist {
            
public:

    /// @param NBO the number of orientation bins
    /// (original SIFT had 8 bins, but in my paper
    /// I found that 16 gives better results with
    /// the SiftDist distance).
    /// NBO has to be greater than 1.
    /// @param CELLS_NUM the number of spatial cells.
    /// Regular SIFT has 4x4=16 cells.
    /// CELLS_NUM has to be greater than 0.
    SiftDist(unsigned int NBO,
             unsigned int CELLS_NUM) 
             
        : NUM_T_ZERO(NumTypeZero<NUM_T>::ZERO()),
          _NBO(NBO),
          _CELLS_NUM(CELLS_NUM)
        {
            assert(NBO>1);
            assert(CELLS_NUM>0);
        }
    
    /// Returns the SiftDist between the two SIFT-like
    /// histograms descriptors sift1, sift2.
    /// Where each histogram contains CELLS_NUM cells
    /// and each cell contains NBO orientation bins.
    /// The dimension that goes fastest is the NBO.
    /// i.e. sift1 looks like:
    /// x_1_1, x_1_2, ... x_1_NBO, x_2_1, ...
    ///
    /// For more details see the paper:
    /// A Linear Time Histogram Metric for Improved SIFT Matching
    /// Ofir Pele, Michael Werman
    /// ECCV 2008
    ///
    /// @param sift1 the first SIFT-like histogram.
    /// @param sift2 the second SIFT-like histogram.
    /// @param stopThresholdsArr array with CELLS_NUM entries of stoping values.
    /// i.e: if the distance after computing the distance for cell number i is greater
    /// or equal to stopThresholdsArr[i] we stop the computation and return stopThresholdsArr[CELLS_NUM-1]
    /// If NULL it is ignored.
    NUM_T operator()(const NUM_T* sift1, 
                     const NUM_T* sift2,
                     const NUM_T* stopThresholdsArr=NULL) {
        
        _dist= NUM_T_ZERO;
        
	    // Runs until _CELLS_NUM-2 and does another
	    // addEmdTModForWindow outside in order
	    // to save the time of the two sift1+= sift2+=
	    for (unsigned int i=0; i<_CELLS_NUM-1; ++i) {
            addEmdTModForWindow(sift1, sift2);
            if (stopThresholdsArr) {
                if (_dist>=stopThresholdsArr[i]) {
                    return stopThresholdsArr[_CELLS_NUM-1];
                }
            }
            sift1+= _NBO;
            sift2+= _NBO;
	    } // i
        addEmdTModForWindow(sift1, sift2);
	    
	    return _dist;
	    
    } // operator()
    
    
private:
    
    
    void addEmdTModForWindow(const NUM_T* Q, const NUM_T* P) {
        
	    // The mass that is left after zero-cost
	    // and one-cost flows
	    NUM_T sumQ= NUM_T_ZERO;
	    NUM_T sumP= NUM_T_ZERO;

	    // After each step we update these two
	    // for next stage
	    NUM_T old_P, old_Q;
	    // j is a running index, declared here because
	    // of goto jumps
	    unsigned int j;
	    	    
	    bool Q_old_zero= (Q[0]<=P[0]);
	    bool P_old_zero= (P[0]<=Q[0]);
	    // i is the index that we can start flowing from
	    // without looking at one-cost edges from i to i-1
	    // Note that we will work on i-1 to i - because of
	    // sumQ and sumP
	    unsigned int i;
	    for (i= 1; i<_NBO; ++i) {
		  
            // We check if there are NO one-cost edges
            // from i-1 to i
            // Note: P[i]<=Q[i]&&Q[i]<=P[i] is also a possible
            // condition for i+1. It is commented out 
            // because this is a rare event and checking it
            // cost more than what we gain...
            
            if (Q[i]<=P[i]) {
                if (Q_old_zero) goto firstPhase;
//			if (P[i]<=Q[i]) {
// 			      ++i;
// 			      goto firstPhase;
// 			} else {
// 			      Q_old_zero= true;
// 			      P_old_zero= false;
// 			      goto nextCheckStage;
// 			}
                Q_old_zero= true;
            } else {
                Q_old_zero= false;
            }
            
            if (P[i]<=Q[i]) {
                if (P_old_zero) goto firstPhase;
                P_old_zero= true;
            } else {
                P_old_zero= false;
            }
            
// 	     nextCheckStage:
// 		  ;
            
	    }
        
	    
	    // i==_NBO
        
	    // edge between last P and first Q
 	    if ( (P[_NBO-1]>Q[_NBO-1]) && (Q[0]>P[0]) ) {
            assert((_NBO%2)==0);
            cyclicEdgeAddEmdTModForWindow(P,Q);
            return;
	    }
	    // edge between last Q and first P
 	    if ( (Q[_NBO-1]>P[_NBO-1]) && (P[0]>Q[0]) ) {
            assert((_NBO%2)==0);
            cyclicEdgeAddEmdTModForWindow(Q,P);
            return;
	    }
	    
	    // if we got here there must NOT be edges between
	    // _NBO-1 to 0 (the cyclic edge)
	    // Thus we can start from i= 0
	    i= 0;
	    
    firstPhase:
        
	    // Flowing from i downward to _NBO-1
	    old_Q= Q[i];
	    old_P= P[i];
	    j= i+1;
	    while (j<_NBO) {
            checkDirectionAndAddSmallFlows
                (Q,P, j,sumQ,sumP,old_Q,old_P);
            ++j;
	    } // end while
	    
        // Flowing from _NBO-1 downward to 0 (cyclic)
	    checkDirectionAndAddSmallFlows
            (Q,P, 0,sumQ,sumP,old_Q,old_P);
        
	    // Flowing from 0 downward to i
	    // Note: although there are not one-cost
	    // edges between i and i-1 we need to do
	    // flow until <=i because of the two-cost
	    // edges - i.e. the sumQ,sumP
	    j= 1;
	    while (j<=i) {
            checkDirectionAndAddSmallFlows
                (Q,P, j,sumQ,sumP,old_Q,old_P);
            ++j;
	    } // end while
	    
	    if (sumQ>=sumP) {
            (_dist+= sumQ)+= sumQ;
	    } else {
            (_dist+= sumP)+= sumP;
	    }
	    
    } // end addSmallFlows
    
    void checkDirectionAndAddSmallFlows(const NUM_T* Q, const NUM_T* P,
                                        const unsigned int& j,
                                        NUM_T& sumQ, NUM_T& sumP,
                                        NUM_T& old_Q, NUM_T& old_P) {
	    if (old_Q>=old_P) {
            addSmallFlows(Q,P,
                          j,sumQ,sumP,old_Q,old_P);
	    } else {
            addSmallFlows(P,Q,
                          j,sumP,sumQ,old_P,old_Q);
	    }
    } // end checkDirectionAndAddSmallFlows
    
    void addSmallFlows(const NUM_T* Q, const NUM_T* P,
                       const unsigned int& j,
                       NUM_T& sumQ, NUM_T& sumP,
                       NUM_T& old_Q, NUM_T& old_P) {
	    
	    NUM_T old_dqp= old_Q-old_P;
	    if (Q[j]>=P[j]) {
		  sumQ+= old_dqp;
		  old_Q= Q[j];
		  old_P= P[j];
	    } else {
		  NUM_T dpq= P[j]-Q[j];
		  old_Q= NUM_T_ZERO;
		  if (old_dqp>=dpq) {
			_dist+= dpq;
			sumQ+= old_dqp-dpq;
			old_P= NUM_T_ZERO;
		  } else {
			_dist+= old_dqp;
			old_P= dpq-old_dqp;
		  }			    
	    } // else Q[j]>=P[j]
    } // end addSmallFlows


      


      //-----------------------------------------------------------
      // All following methods are for the cyclic edge special case
      
      // Note: as this function is rarely used, I didn't bother
      // to optimize it.
      // Assumes that the last cyclic edge is from last Q to first P.
      // i.e: P[0]>Q[0] , Q[1]>P[1] , ... , Q[_NBO-1] > P[_NBO-1]
      // Also assumes that _NBO is even (otherwise there are no cycles)
      void cyclicEdgeAddEmdTModForWindow(const NUM_T* Q, const NUM_T* P) {
	    
	    // Copy without zeros to cQ,cP.
	    // i.e:
	    //   cP   cQ
	    //
	    //   *<-0-*
	    //  /\   /
	    //  /   1
	    // |   /
	    // |  /
	    // | \/
	    // | *<-2-*
	    // |     /
	    // |    3
	    // |   /
	    // |  /
	    // | \/
	    // | *<-4-*
	    // |     /
	    // \-5---
	    unsigned int NBO_DIV_2= _NBO/2;
	    std::vector<NUM_T> cP(NBO_DIV_2);
	    std::vector<NUM_T> cQ(NBO_DIV_2);
        unsigned int i;
        for (i=0; i<NBO_DIV_2; ++i) {
		  cP[i]= P[2*i]   - Q[2*i];
		  cQ[i]= Q[2*i+1] - P[2*i+1];
		  assert(cP[i]>0&&cQ[i]>0);
	    }
	    
	    /*
	      maximize the flow ignoring the edge
	      between cP[0] and cQ[0]
	      Stores the partial residual capacity.
	      It is partial as it includes only these edges
	      (and not the ones on the other direction):
	      cP   cQ
	       *<-0-*
	      /    /\
	      /    /
	      /   1
	     |   /
	     |  /
	     | /
	     | *<-2-*
	     |     /\
	     |     /
	     |    3
	     |   /
	     |  /
	     | /
	     | *<-4-*
	     |     /\
	     \-5---/
	    */
	    std::vector<NUM_T> cP_left(cP);
	    std::vector<NUM_T> cQ_left(cQ);
	    std::vector<NUM_T> partial_residual_capacity_vec(_NBO);
	    // The ignored edge
	    partial_residual_capacity_vec[0]=
		  myMin(cP[0],cQ[0]); 
	    unsigned int j= 1;
	    for (i=0; i<NBO_DIV_2-1; ++i) {
		  //   /\ edges
		  //   /
		  //  /
		  partial_residual_capacity_vec[j]=
			flow(cQ_left,i,cP_left,i+1);
		  ++j;
		  // <---- edges
		  partial_residual_capacity_vec[j]=
			myMin(cQ[i+1],cP[i+1])-flow(cQ_left,i+1,cP_left,i+1);
		  ++j;
	    }
	    partial_residual_capacity_vec[j]=
		  flow(cQ_left,NBO_DIV_2-1,cP_left,0);

	    
	    NUM_T sumLeftP= NUM_T_ZERO;
	    NUM_T sumLeftQ= NUM_T_ZERO;
	    for (i=0; i<NBO_DIV_2; ++i) {
		  sumLeftQ+= cQ_left[i];
	          sumLeftP+= cP_left[i];
	    }
	    NUM_T maxSumLeft= myMax(sumLeftQ,sumLeftP);
	    

	    // vectors of indices
	    std::vector<unsigned int> targets_vec(NBO_DIV_2);
	    std::vector<unsigned int> sources_vec(NBO_DIV_2);
	    for (i=0; i<NBO_DIV_2; ++i) {
		  sources_vec[i]= i;
	    }
	    targets_vec[0]= 0;
	    targets_vec[1]= NBO_DIV_2-1;
	    for (i=2; i<NBO_DIV_2; ++i) {
		  targets_vec[i]= targets_vec[i-1]-1;
	    }
	    assert(targets_vec[NBO_DIV_2-1]==1);
	    unsigned int s_i= 0;
	    unsigned int t_i= 0;
	    
	    // capacity of critical edge in the "middle-edges"
	    // i.e: between the nodes of Q ("small-sources")
	    // to the nodes of P ("small-targets")
	    NUM_T m_capacity=
		  partial_residual_capacity_vec[0]; 
	    // indice that of the next "middle-edge" in
	    // the partial_residual_capacity_vec from a
	    // "small-target" to a "small-source"
	    // (backward edge)
	    // Note that this index should decrease for next
	    unsigned int m_t_i= _NBO-1;
	    // indice that of the next "middle-edge" in
	    // the partial_residual_capacity_vec from a
	    // "small-source" to a "small-target"
	    // (forward edge)
	    // Note that this index should increase for next
	    unsigned int m_s_i= 1;
	    
	    
	    // Loop that checks if there is an augmenting path
	    // and extends it if needed.
	    // In each iteration one of: s_capacity,m_capacity,t_capacity
	    // is zeroed. When s_capacity or t_capacity are zeroed,
	    // s_i and t_i respectively increase. When m_capacity is zeroed,
	    // we stop. Thus runs at most NBO_DIV_2 times
	    while (true) {
		  
		  if (m_capacity==NUM_T_ZERO) {
			(_dist+= maxSumLeft)+= maxSumLeft;
			return;
		  }

		  NUM_T s_capacity= cQ_left[ sources_vec[s_i] ];
		  while (s_capacity==NUM_T_ZERO) {
			++s_i;
			if (s_i==NBO_DIV_2||m_capacity==NUM_T_ZERO) {
			      (_dist+= maxSumLeft)+= maxSumLeft;
			      return;  
			}
			assert(m_s_i+1<_NBO);
			/*
			  Expand augmented path with these two edges
			 (s is old source, ns is new source)
			     s
			    /\
			   /
			  /
			 *<--ns
			 */
			m_capacity= myMin(m_capacity,myMin(partial_residual_capacity_vec[m_s_i],partial_residual_capacity_vec[m_s_i+1]));
			m_s_i+= 2;
			s_capacity= cQ_left[ sources_vec[s_i] ];
		  }

		  NUM_T t_capacity= cP_left[ targets_vec[t_i] ];
		  while (t_capacity==NUM_T_ZERO) {
			++t_i;
			if (t_i==NBO_DIV_2||m_capacity==NUM_T_ZERO) {
			      (_dist+= maxSumLeft)+= maxSumLeft;
			      return;  
			}
			assert(m_t_i-1>=0&&m_t_i<_NBO);
			/*
			  Expand augmented path with these two edges:
			  (t is old target, nt is new target)
			  nt<--*
			      /\
			      /
			     /
			    /
			   t
			*/
			m_capacity= myMin(m_capacity,myMin(partial_residual_capacity_vec[m_t_i],partial_residual_capacity_vec[m_t_i-1]));
			m_t_i-= 2;
			t_capacity= cP_left[ targets_vec[t_i] ];
		  }

		  NUM_T f= myMin( myMin(s_capacity,t_capacity),m_capacity); 
		  
		  _dist+= f;
		  maxSumLeft-= f;
		  m_capacity-= f;
		  cP_left[ targets_vec[t_i] ]-= f;
		  cQ_left[ sources_vec[s_i] ]-= f;
		 
		  
	    } // end true
		  
	    
      } // end cyclicEdgeAddEmdTModForWindow


      NUM_T flow(std::vector<NUM_T>& cQ, unsigned int i,
		std::vector<NUM_T>& cP, unsigned int j) {

	    NUM_T f= myMin(cQ[i],cP[j]);
	    cQ[i]-= f;
	    cP[j]-= f;
	    _dist+= f;
	    
	    return f;
      } // flow


      //-----------------------------------------------------------
	    
      NUM_T _dist;
      const NUM_T NUM_T_ZERO;
      unsigned int _NBO, _CELLS_NUM;
      
}; // end class SiftDist


#endif
