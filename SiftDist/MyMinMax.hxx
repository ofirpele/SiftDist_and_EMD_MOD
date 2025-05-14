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

#ifndef _OFIRPELE_MY_MIN_MAX__HXX
#define _OFIRPELE_MY_MIN_MAX__HXX

// old versions of visual studio does not have std::min, std::max
// in their libraries...

template<typename T>
const T& myMin(const T& a, const T& b) {
        if (a<b) return a;
        return b;
}
template<typename T>
const T& myMax(const T& a, const T& b) {
        if (a>b) return a;
        return b;
}

#endif
