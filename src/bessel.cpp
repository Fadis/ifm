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

#include <iostream>
#include <cmath>
#include "ifm/bessel.h"
namespace ifm {
  float bessel_kind1_1( float b, int n ) {
    float prev = 0;
    float width = 0.001;
    float sum = 0;
    for( float w = width; w < float( M_PI/2 ); w += width ) {
      float cur = std::sin( b * std::sin( w ) ) * std::sin( n*w );
      sum += ( cur + prev ) * width / 2;
      prev = cur;
    }
    return 2.f / M_PI * sum;
  }

  float bessel_kind1_0( float b, int n ) {
    float prev = 1;
    float width = 0.001;
    float sum = 0;
    for( float w = width; w < float( M_PI/2 ); w += width ) {
      float cur = std::cos( b * std::sin( w ) ) * std::cos( n*w );
      sum += ( cur + prev ) * width / 2;
      prev = cur;
    }
    return 2.f / M_PI * sum;
  }

  float bessel_kind1( float b, int n ) {
    if( n % 2 ) return bessel_kind1_1( b, n );
    else return bessel_kind1_0( b, n );
  }

  float bessel_kind1_approx_2013( float b, int n ) {
    return std::sqrt( 2. / ( M_PI * b ) ) * std::cos( b - ( n / 2. + 1. / 4. ) * M_PI );
  }


  float bessel_kind1_approx_2019_0_original1( float x ) {
    float p = 8/x;
    float q = p * p;
    float t = x - 0.78539816;
    return std::sqrt( 0.079577472 * p ) * (( 0.99999692  - 0.0010731477 * q ) * std::cos( t )
      + ( 0.015624982 + ( -0.00014270786 + 0.0000059374342 * q ) * q ) * p * std::sin( t ) );
  }
  float bessel_kind1_approx_2019_1_original1( float x ) {
    float p = 0.63661977/x;
    float q = p * p;
    return std::sqrt( p ) * ( 1.0 - ( 0.46263771 - 1.1771851 * q ) * q )
      * std::cos( x - 2.3561945 + ( 0.58904862 - 0.63587091 * q ) * p );
  }
  float bessel_kind1_approx_2019_original1( float x, int n ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0_original1( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1_original1( x );
    float prev = n - 1;
    return 2*prev/x * bessel_kind1_approx_2019_original1( x, n - 1 ) - bessel_kind1_approx_2019_original1( x, n - 2 );
  }
  float bessel_kind1_approx_2019_original2( float x, int n ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0_original1( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1_original1( x );
    float prev = n - 1;
    if( x < 2*prev ) return bessel_kind1_approx_2019_original2( x, n - 1 );
    return 2*prev/x * bessel_kind1_approx_2019_original2( x, n - 1 ) - bessel_kind1_approx_2019_original2( x, n - 2 );
  }




  float bessel_kind1_approx_2019_0( float x ) {
    //if( x < 0 ) x = -x;
    if( x < 1.0f ) return 1.f - x * x / 4.f + x * x * x * x / 64.f - x * x * x * x * x * x / 2304;
    float p = 8/x;
    float q = p * p;
    float t = x - 0.78539816;
    return std::sqrt( 0.079577472 * p ) * (( 0.99999692  - 0.0010731477 * q ) * std::cos( t )
      + ( 0.015624982 + ( -0.00014270786 + 0.0000059374342 * q ) * q ) * p * std::sin( t ) );
  }

  float bessel_kind1_approx_2019_1( float x ) {
    //if( x < 0 ) x = -x;
    if( x < 3.5f ) return x / 2.f - x * x * x / 16.f + x * x * x * x * x / 384.f - x * x * x * x * x * x * x / 18432.f + x * x * x * x * x * x * x * x * x / 1474560.f;
    float p = 0.63661977/x;
    float q = p * p;
    return std::sqrt( p ) * ( 1.0 - ( 0.46263771 - 1.1771851 * q ) * q )
      * std::cos( x - 2.3561945 + ( 0.58904862 - 0.63587091 * q ) * p );
  }

  float bessel_kind1_approx_2019( float x, int n );

  float bessel_kind1_approx_2019_( float x, int n ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1( x );
    //if( x < 0 ) x = -x;
    float prev = n - 1;
    return 2*prev/x * bessel_kind1_approx_2019( x, n - 1 ) - bessel_kind1_approx_2019( x, n - 2 );
  }

  float bessel_kind1_approx_2019_( float x, int n, float l1, float l2 ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1( x );
    //if( x < 0 ) x = -x;
    float prev = n - 1;
    return 2*prev/x * l1 - l2;
  }

  float bessel_kind1_approx_2019( float x, int n, float l1 );
  float bessel_kind1_approx_2019( float x, int n ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1( x );
    //if( x < 0 ) x = -x;
    float x1 = n;
    if( x1 > x ) {
      float y1 = bessel_kind1_approx_2019_( x1, n );
      float yprev = bessel_kind1_approx_2019_( x1, n - 1 );
      float dy1 = 0.5f * ( yprev - bessel_kind1_approx_2019_( x1, n + 1, y1, yprev ) );
      float a = dy1/x1 - y1/(x1*x1);
      float b = 2.f * y1 / x1 - dy1;
      return ( a * x * x * x * x + b * x * x * x ) / ( x1 * x1 );
      /*float a = dy1/(x1*x1*x1) - 3*y1/(x1*x1*x1*x1);
      float b = 4*y1/(x1*x1*x1) - dy1/(x1*x1);
      return a * x * x * x * x + b * x * x * x;*/
    }
    float prev = n - 1;
    float l2 = bessel_kind1_approx_2019( x, n - 2 );
    float l1 = bessel_kind1_approx_2019( x, n - 1, l2 );
    return 2*prev/x * l1 - l2;
  }
  float bessel_kind1_approx_2019( float x, int n, float l1 ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1( x );
    //if( x < 0 ) x = -x;
    float x1 = n;
    if( x1 > x ) {
      float y1 = bessel_kind1_approx_2019_( x1, n );
      float yprev = bessel_kind1_approx_2019_( x1, n - 1 );
      float dy1 = 0.5f * ( yprev - bessel_kind1_approx_2019_( x1, n + 1, y1, yprev ) );
      float a = dy1/x1 - y1/(x1*x1);
      float b = 2.f * y1 / x1 - dy1;
      return ( a * x * x * x * x + b * x * x * x ) / ( x1 * x1 );
      /*float a = dy1/(x1*x1*x1) - 3*y1/(x1*x1*x1*x1);
      float b = 4*y1/(x1*x1*x1) - dy1/(x1*x1);
      return a * x * x * x * x + b * x * x * x;*/
    }
    float prev = n - 1;
    float l2 = bessel_kind1_approx_2019( x, n - 2 );
    return 2*prev/x * l1 - l2;
  }
  float bessel_kind1_approx_2019( float x, int n, float l1, float l2 ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1( x );
    //if( x < 0 ) x = -x;
    float x1 = n;
    if( x1 > x ) {
      float y1 = bessel_kind1_approx_2019_( x1, n );
      float yprev = bessel_kind1_approx_2019_( x1, n - 1 );
      float dy1 = 0.5f * ( yprev - bessel_kind1_approx_2019_( x1, n + 1, y1, yprev ) );
      float a = dy1/x1 - y1/(x1*x1);
      float b = 2.f * y1 / x1 - dy1;
      return ( a * x * x * x * x + b * x * x * x ) / ( x1 * x1 );
      /*float a = dy1/(x1*x1*x1) - 3*y1/(x1*x1*x1*x1);
      float b = 4*y1/(x1*x1*x1) - dy1/(x1*x1);
      return a * x * x * x * x + b * x * x * x;*/
    }
    float prev = n - 1;
    return 2*prev/x * l1 - l2;
  }
  float bessel_kind1_approx_2019( float x, int n, const std::vector< std::pair< float, float > > &pre ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1( x );
    //if( x < 0 ) x = -x;
    float x1 = n;
    if( x1 > x ) {
      if( n < int( pre.size() ) ) {
        float y1 = pre[ n ].first;
        float dy1 = pre[ n ].second;
        float a = dy1/x1 - y1/(x1*x1);
        float b = 2.f * y1 / x1 - dy1;
        return ( a * x * x * x * x + b * x * x * x ) / ( x1 * x1 );
        /*float a = dy1/(x1*x1*x1) - 3*y1/(x1*x1*x1*x1);
        float b = 4*y1/(x1*x1*x1) - dy1/(x1*x1);
        return a * x * x * x * x + b * x * x * x;*/
      }
      else return bessel_kind1_approx_2019( x, n );
    }
    float prev = n - 1;
    return 2*prev/x * bessel_kind1_approx_2019( x, n - 1 ) - bessel_kind1_approx_2019( x, n - 2 );
  }
  float bessel_kind1_approx_2019( float x, int n, float l1, float l2, const std::vector< std::pair< float, float > > &pre ) {
    if( n < 0 ) return 0;
    if( n == 0 ) return bessel_kind1_approx_2019_0( x );
    if( n == 1 ) return bessel_kind1_approx_2019_1( x );
    //if( x < 0 ) x = -x;
    float x1 = n;
    if( x1 > x ) {
      if( n < int( pre.size() ) ) {
        float y1 = pre[ n ].first;
        float dy1 = pre[ n ].second;
        float a = dy1/x1 - y1/(x1*x1);
        float b = 2.f * y1 / x1 - dy1;
        return ( a * x * x * x * x + b * x * x * x ) / ( x1 * x1 );
        /*float a = dy1/(x1*x1*x1) - 3*y1/(x1*x1*x1*x1);
        float b = 4*y1/(x1*x1*x1) - dy1/(x1*x1);
        return a * x * x * x * x + b * x * x * x;*/
      }
      else return bessel_kind1_approx_2019( x, n, l1, l2 );
    }
    float prev = n - 1;
    return 2*prev/x * l1 - l2;
  }
  std::vector< std::pair< float, float > > create_bessel_approx_2019_precomp_array( int max ) {
    std::vector< std::pair< float, float > > pre;
    for( int i = 0; i != max; ++i ) {
      float l1 = bessel_kind1_approx_2019( i, i - 1, pre );
      float l2 = bessel_kind1_approx_2019( i, i - 2, pre );
      float y = bessel_kind1_approx_2019( i, i, l1, l2, pre );
      float dy = 0.5 * ( l1 - l2 );
      pre.emplace_back( std::make_pair( y, dy ) ); 
    }
    return pre;
  }
}

