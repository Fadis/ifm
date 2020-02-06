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

#include <iterator>
#include <algorithm>
#include <iostream>
#include "ifm/segment_envelope.h"

namespace ifm {
std::tuple< int, int, int > segment_envelope( const std::vector< float > &input, unsigned int sample_rate ) {
  std::vector< float > grad;
  for( size_t i = 1u; i != input.size(); ++i )
    grad.emplace_back( ( input[ i ] - input[ i - 1u ] )*sample_rate );
  std::vector< float > release( input.size(), 0 );
  const auto blank = std::distance(
    input.rbegin(),
    std::find_if( input.rbegin(), input.rend(), []( float v ) { return v != 0; } )
  ) + 1;
  float min_grad = std::abs( grad[ input.size() - ( 1 + blank ) ] );
  for( size_t i = 1u; i != input.size() - blank; ++i ) {
    min_grad = std::min( min_grad, std::max( 0.f, -grad[ input.size() - ( i + blank ) ] ) );
    release[ input.size() - ( i + blank ) ] = min_grad * ( float( i ) / sample_rate );
  }
  const auto release_pos = std::distance( release.begin(), std::max_element( release.begin(), release.end() ) );
  float delay = 0.f;
  const float max = *std::max_element( input.begin(), std::next( input.begin(), release_pos ) );
  const auto p0 = std::find_if( input.begin(), std::next( input.begin(), release_pos ), [&]( float v ){ return v >= max * 0.2f; } );
  const auto p1 = std::find_if( input.begin(), std::next( input.begin(), release_pos ), [&]( float v ){ return v >= max * 0.4f; } );
  const auto tangent = float( std::distance( p0, p1 ) )/( *p1 - *p0 );
  const auto intercept = float( std::distance( input.begin(), p0 ) ) - tangent * *p0;
  const int delay_pos = std::max( int( intercept ), 0 );
  std::vector< float > attack( input.size(), 0 );
  min_grad = 1.f/tangent;
  for( size_t i = 0u; i != release_pos - delay_pos; ++i ) {
    min_grad = std::min( min_grad, std::max( 0.f, grad[ i + delay_pos ] ) );
    attack[ i + delay_pos ] = min_grad * ( float( i ) / sample_rate );
  }
  const auto attack_pos = std::distance( attack.begin(), std::max_element( attack.begin(), attack.end() ) );
  std::cout << attack_pos << " " << release_pos << std::endl;
  return std::make_tuple( delay_pos, attack_pos, release_pos );
}
}

