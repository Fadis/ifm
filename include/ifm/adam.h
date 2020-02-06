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

#ifndef IFM_ADAM_H
#define IFM_ADAM_H
#include <cmath>
namespace ifm {
  template< typename T >
  class adam {
  public:
    adam(
      T alpha_,
      T beta1_,
      T beta2_
    ) : m( 0 ), v( 0 ), alpha( alpha_ ), beta1( beta1_ ), beta2( beta2_ ), current_beta1( beta1_ ), current_beta2( beta2_ ) {}
    T operator()( T grad ) {
      if( std::isnan( grad ) ) throw -1;
      m = beta1 * m + ( 1 - beta1 ) * grad;
      v = beta2 * v + ( 1 - beta2 ) * grad * grad;
      T mhat = m / ( 1 - current_beta1 );
      T vhat = v / ( 1 - current_beta2 );
      T diff = ( alpha * mhat / ( std::sqrt( vhat ) + T( 0.000001f ) ) );
      current_beta1 *= beta1;
      current_beta2 *= beta2;
      return diff;
    }
  private:
    T m;
    T v;
    T alpha;
    T beta1;
    T beta2;
    T current_beta1;
    T current_beta2;
  };
}

#endif

