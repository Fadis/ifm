/*
Copyright (c) 2020 Naomasa Matsubayashi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef IFM_BESSEL_H
#define IFM_BESSEL_H
#include <vector>
#include <utility>
#include <cmath>
namespace ifm {
  float bessel_kind1_1( float b, int n );
  float bessel_kind1_0( float b, int n );
  float bessel_kind1( float b, int n );
  float bessel_kind1_approx_2013( float b, int n );
  float bessel_kind1_approx_2019_0( float x );
  float bessel_kind1_approx_2019_1( float x );
  float bessel_kind1_approx_2019( float x, int n );
  float bessel_kind1_approx_2019_( float x, int n );
  float bessel_kind1_approx_2019( float x, int n );
  float bessel_kind1_approx_2019( float x, int n, float l1, float l2 );
  float bessel_kind1_approx_2019( float x, int n, const std::vector< std::pair< float, float > > &pre );
  float bessel_kind1_approx_2019( float x, int n, float l1, float l2, const std::vector< std::pair< float, float > > &pre );
  float bessel_kind1_approx_2019_original1( float x, int n );
  float bessel_kind1_approx_2019_original2( float x, int n );
  std::vector< std::pair< float, float > > create_bessel_approx_2019_precomp_array( int max );
}
#endif

