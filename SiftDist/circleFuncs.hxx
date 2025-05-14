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

#ifndef _OFIRPELE_CIRCLE_FUNCX__HXX_
#define _OFIRPELE_CIRCLE_FUNCX__HXX_

#include <math.h>
#include <cassert>

class circleFuncs {

public:

      //-----------------------------------------------------------------------
      static double circleOverlap(double x1, double y1, double r1,
				  double x2, double y2, double r2) {
		
	    const double CO_PI= 3.14159265;
	    double x1_m_x2= x1-x2;
	    double y1_m_y2= y1-y2;
	    double d= sqrt( x1_m_x2*x1_m_x2 + y1_m_y2*y1_m_y2 ); 

	    if (d>=r1+r2) return 0.0;

	    double rr1= r1*r1;
	    double rr2= r2*r2;
	    double intersect_area;
	    if (d==0.0) {
		  if (r1<r2) {
			intersect_area= rr1*CO_PI;
		  } else {
			intersect_area= rr2*CO_PI;
		  }
	    } else {
		  double dd= d*d;
		  double two_d= 2.0*d;
		  intersect_area= (rr1*acos((dd+rr1-rr2)/(two_d*r1)))+
			(rr2*acos((dd-rr1+rr2)/(two_d*r2)))-
			(sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2) )/2.0);
	    }
	    
	    double union_area= rr1*CO_PI + rr2*CO_PI - intersect_area; 
	    return intersect_area/union_area;
	    
      } // circleOverlap
      //-----------------------------------------------------------------------
      
};

#endif
