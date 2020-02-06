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

#include <cstdint>
#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <omp.h>
#include "ifm/bessel.h"
#include "ifm/adam.h"
#include "ifm/2op.h"
namespace ifm {
  std::tuple< float, float > lossimage(
    const float *expected,
    const std::vector< std::pair< float, float > > &pre,
    unsigned int,
    float b,
    const std::array< int, max_harmony > &n
  ) {
    constexpr float a = 1;
    std::array< float, max_harmony + 3 > bessel{ 0 };
    std::array< float, max_harmony > generated{ 0 };
    float l1 = 0;
    float l2 = 0;
    for( unsigned int i = 0; i < ( max_harmony + 3 ); ++i ) {
      bessel[ i ] = ifm::bessel_kind1_approx_2019( b, i, l1, l2, pre );
      l2 = l1;
      l1 = bessel[ i ];
    }
    float bscale = 0;
    for( unsigned int i = 0; i != max_harmony; ++i ) {
      bscale += std::abs( bessel[ i ] );
      if( i < max_harmony - 2 ) bscale += std::abs( bessel[ i + 2 ] );
    }
    //float dbscale = 0.5 * ( bessel[ 0 ] - bessel[ 2 ] ) + bessel[ 2 ] - bessel[ max_harmony - 4 ] - bessel[ max_harmony - 3 ] + 0.5 * ( bessel[ max_harmony - 4 ] - bessel[ max_harmony - 2 ] ) + 0.5 * ( bessel[ max_harmony - 3 ] - bessel[ max_harmony - 1 ] );
    for( unsigned int i = 0; i != max_harmony - 2; ++i ) {
      if( n[ i ] != -1 ) generated[ i ] = a * bessel[ n[ i ] ]/bscale;
      else generated[ i ] = 0;
      if( i < max_harmony - 2 && n[ i + 2 ] != -1 )
        generated[ i ] += a * bessel[ n[ i + 2 ] ]/bscale;
    }
    float loss = 0;
    float d = 0;
    for( unsigned int i = 0; i != max_harmony; ++i ) {
      if( n[ i ] != -1 ) {
        float e2 = std::sqrt( expected[ i ] * expected[ i ] );
          float g2 = std::sqrt( generated[ i ] * generated[ i ] );
          float s = e2 - g2;
          loss += std::sqrt( s * s );
          d += std::sqrt( s * s );
      }
      else {
        float e2 = std::sqrt( expected[ i ] * expected[ i ] );
        loss += e2;
      }
    }
    return std::make_tuple( loss, d );
  }

  std::tuple< float, float > find_b_2op(
    const float *expected,
    const std::vector< std::pair< float, float > > &pre,
    unsigned int,
    float b,
    const std::array< int, max_harmony > &n
  ) {
    constexpr float a = 1;
    std::array< float, max_harmony + 3 > bessel{ 0 };
    std::array< float, max_harmony > generated{ 0 };
    float loss = 0;
    ifm::adam< float > bopt( 0.001, 0.9, 0.999 );
    for( unsigned int cycle = 0; cycle != 50000; ++cycle ) {
      float l1 = 0;
      float l2 = 0;
      for( unsigned int i = 0; i < ( max_harmony + 3 ); ++i ) {
        bessel[ i ] = ifm::bessel_kind1_approx_2019( b, i, l1, l2, pre );
        l2 = l1;
        l1 = bessel[ i ];
      }
      float bscale = 0;
      for( unsigned int i = 0; i != max_harmony; ++i ) {
        bscale += std::abs( bessel[ i ] );
        if( i < max_harmony - 2 ) bscale += std::abs( bessel[ i + 2 ] );
      }
      float dbscale = 0.5 * ( bessel[ 0 ] - bessel[ 2 ] ) + bessel[ 2 ] - bessel[ max_harmony - 4 ] - bessel[ max_harmony - 3 ] + 0.5 * ( bessel[ max_harmony - 4 ] - bessel[ max_harmony - 2 ] ) + 0.5 * ( bessel[ max_harmony - 3 ] - bessel[ max_harmony - 1 ] );
      for( unsigned int i = 0; i != max_harmony - 2; ++i ) {
        if( n[ i ] != -1 ) generated[ i ] = a * bessel[ n[ i ] ]/bscale;
        else generated[ i ] = 0;
        if( i < max_harmony - 2 && n[ i + 2 ] != -1 )
          generated[ i ] += a * bessel[ n[ i + 2 ] ]/bscale;
      }
      float grad_b = 0;
      loss = 0;
      for( unsigned int i = 0; i != max_harmony; ++i ) {
      if( n[ i ] != -1 ) {
        float e2 = std::sqrt( expected[ i ] * expected[ i ] );
          float g2 = std::sqrt( generated[ i ] * generated[ i ] );
          float s = e2 - g2;
          loss += std::sqrt( s * s );
          std::array< float, 5 > j{
            bessel[ n[ i ] - 1 ],
            bessel[ n[ i ] ],
            bessel[ n[ i ] + 1 ],
            bessel[ n[ i ] + 2 ],
            bessel[ n[ i ] + 3 ]
          };
          float j13 = j[1]+j[3];
          float aaj13 = j13*j13;
          float grad_b_b1 = e2-std::sqrt(aaj13/(bscale*bscale));
          if( grad_b_b1 != 0 ) {
            if( aaj13 > 1.0e-20 ) {
              /*float aaj13s = aaj13/(bscale*bscale);
              float aaj13sr = std::sqrt( aaj13s );
              float grad_b_b12 = grad_b_b1*grad_b_b1;
              float grad_b_b12r = std::sqrt( grad_b_b12 );
              float bs3 = bscale*bscale*bscale;*/
              float grad_b__ = (0.5*j13*
               (e2 - std::sqrt(aaj13/(bscale*bscale)))*
               ((-j[0] + j[4])*bscale + 
                (2.*j[1] + 2.*j[3])*dbscale))/
                (std::sqrt(grad_b_b1*grad_b_b1)*std::sqrt(aaj13/(bscale*bscale))*
                 (bscale*bscale*bscale));
              /*float grad_b_ = (0.5*aaj13*(e2-aaj13sr)*
               ((-j[0]+j[4])*bscale+(2*j[1]+2*j[3])*dbscale))/
              (grad_b_b12r*aaj13sr*bs3);*/

              grad_b += grad_b__;
            }
          }
        }
        else {
          float e2 = std::sqrt( expected[ i ] * expected[ i ] );
          loss += e2;
        }
      }
      float diff = bopt( grad_b );
      b -= diff;
      b = std::max( b, 0.f );
    }
    return std::make_tuple( loss, b );
  }
  std::array< int, max_harmony > generate_n(
    unsigned int freq
  ) {
    std::array< int, max_harmony > n;
    std::fill( n.begin(), n.end(), -1 );
    for( int i = 0; i != max_harmony; ++i )
      if( i % freq == 0 )
        n[ i ] = i / freq;
    return n;
  }
  std::tuple< float, float > find_b_2op(
    const float *expected,
    const std::vector< std::pair< float, float > > &pre,
    unsigned int freq
  ) {
    float lowest_loss = std::numeric_limits< float >::max();
    float best_b = 0;
    std::array< int, max_harmony > n = generate_n( freq );
    for( unsigned int try_count = 0; try_count != 5; ++try_count ) {
      const auto [loss,b] = find_b_2op( expected, pre, freq, try_count, n );
      if( lowest_loss > loss ) {
        best_b = b;
        lowest_loss = loss;
      }
    }
    return std::make_tuple( lowest_loss, best_b );
    //return find_b_2op( expected, pre, freq, 5, n );
  }
  std::tuple< float, unsigned int, float > find_b_2op(
    const float *expected,
    const std::vector< std::pair< float, float > > &pre
  ) {
    constexpr unsigned int min_freq = 1;
    constexpr unsigned int max_freq = 20;
    std::vector< std::tuple< float, float > > results( max_freq );
#pragma omp parallel for
    for( unsigned int freq = min_freq; freq < max_freq; ++freq ) {
      results[ freq ] = find_b_2op( expected, pre, freq );
    }
    float lowest_loss = std::numeric_limits< float >::max();
    float best_b = 0;
    unsigned int best_freq = 0;
    for( unsigned int freq = min_freq; freq != max_freq; ++freq ) {
      if( lowest_loss > std::get< 0 >( results[ freq ] ) ) {
        best_b = std::get< 1 >( results[ freq ] );
        best_freq = freq;
        lowest_loss = std::get< 0 >( results[ freq ] );
      }
    }
    return std::make_tuple( lowest_loss, best_freq, best_b );
  }
}

