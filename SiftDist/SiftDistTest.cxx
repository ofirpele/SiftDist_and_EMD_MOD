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

#include "SiftDist.hxx"
#include <iostream>

int main() {

      int NBO= 8;
      int CELLS_NUM= 4;
      
      // Two SIFTs with 2x2=4 spatiall cells
      // and 8 orientations bins.
      // Note that orientations bins marches the fastest.
      // That is:
      // sift1= x_1_1, x_1_2, ... x_1_NBO,
      //        x_2_1, ...
      double sift1[]=
	    {9.0, 6.0, 9.0, 9.0, 6.0, 4.0, 3.0, 5.0, 
	     2.0, 5.0, 3.0, 5.0, 9.0, 2.0, 8.0, 2.0, 
	     8.0, 6.0, 3.0, 1.0, 2.0, 9.0, 7.0, 2.0, 
	     6.0, 2.0, 3.0, 6.0, 1.0, 1.0, 2.0, 7.0};
      
      double sift2[]=
	    {6.0, 2.0, 1.0, 9.0, 4.0, 3.0, 7.0, 2.0, 
	     3.0, 1.0, 4.0, 4.0, 1.0, 6.0, 7.0, 5.0, 
	     6.0, 3.0, 9.0, 1.0, 9.0, 3.0, 9.0, 1.0, 
	     3.0, 5.0, 2.0, 1.0, 9.0, 7.0, 3.0, 9.0};
      
      SiftDist<double> sd(NBO, CELLS_NUM);
      assert(sd(sift1, sift2)==110.0);


      double stopThresholdsArr[]=
          {39,39,39,4242};
      
      assert(sd(sift1, sift2,stopThresholdsArr)==4242);
      
      return 0;
}
