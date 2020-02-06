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

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include "ifm/exp_match.h"
#include "ifm/adam.h"
namespace ifm {
  std::tuple< float, float >
  exp_match( const float *expected, unsigned int size, float dt, float c ) {
    float a = 1.0;
    float b = 0.0;
    adam< float > aopt( 0.001, 0.9, 0.999 );
    adam< float > bopt( 0.001, 0.9, 0.999 );
    for( unsigned int cycle = 0; cycle != 50000; ++cycle ) {
      std::vector< float > ag( size );
      std::vector< float > bg( size );
#pragma omp parallel for
      for( unsigned int i = 0; i < size; ++i ) {
        float x = dt * float( i );
        float generated = ( std::exp( -x * a ) * ( 1.f -  b ) + b );
        float d = expected[ i ]/c - generated;
        d = d * d;
        ag[ i ] = (2*(1 - b)*(-b - (1 - b)/std::exp(a*x) + expected[ i ]/c)*x)/std::exp(a*x) * d;
        bg[ i ] = 2*(-1 + std::exp(-(a*x)))*(-b - (1 - b)/std::exp(a*x) + expected[ i ]/c) * d;
      }
      float grad_a = std::accumulate( ag.begin(), ag.end(), float( 0 ) );
      float grad_b = std::accumulate( bg.begin(), bg.end(), float( 0 ) );
      a -= aopt( grad_a );
      b -= bopt( grad_b );
    }
    return std::make_tuple( a, b );
  }
  exp_envelope_params_t
  get_attack_and_decay_exp( const float *expected, unsigned int size, float dt, bool damped ) {
    const auto highest = std::max_element( expected, std::next( expected, size ) );
    const auto highest_pos = std::distance( expected, highest );
    const auto [decay_a, decay_b] = exp_match( std::next( expected, highest_pos ), size - highest_pos, dt, *highest );
    std::vector< float > reversed( std::make_reverse_iterator( highest ), std::make_reverse_iterator( expected ) );
    if( damped ) {
      return exp_envelope_params_t()
        .set_decay_a( decay_a )
        .set_decay_b( decay_b )
        .set_highest( highest_pos * dt )
        .set_highest_level( *highest );
    }
    else {
      const auto [attack_a, attack_b] = exp_match( reversed.data(), highest_pos, dt, *highest );
      return exp_envelope_params_t()
        .set_attack_a( attack_a )
        .set_attack_b( attack_b )
        .set_decay_a( decay_a )
        .set_decay_b( decay_b )
        .set_highest( highest_pos * dt )
        .set_highest_level( *highest );
    }
  }
  float get_exp_envelope( float a, float b, float x ) {
    return std::exp( -a*x )*(1.f-b)+b;
  }
  std::pair< float, float > get_linear_interpolation( float x0, float y0, float x1, float y1 ) {
    float tangent = ( y0 - y1 )/( x0 - x1 );
    float shift = y0 - tangent * x0;
    return std::make_pair( tangent, shift );
  }
  float get_envelope_error(
    double a, double b, double m, double n
  ) {
    return ((-1+b)*(std::exp(a*m)*(-2+a*(m-n))+std::exp(a*n)*(2+a*(m-n))))/(2.*a*std::exp(a*(m+n)));
  }
  float get_sustain_error(
    double a, double b, double m, double n
  ) {
    return ((-1+b)*(-std::exp(a*m)+std::exp(a*n)*(1+a*(m-n))))/(a*std::exp(a*(m+n)));
  }
  std::tuple< float, float >
  approxymate_decay( double a, double b, double length ) {
    double m = 1.f/3.f * length;
    double n = 2.f/3.f * length;
    adam< double > mopt( 0.001, 0.9, 0.999 );
    adam< double > nopt( 0.001, 0.9, 0.999 );
    for( unsigned int cycle = 0; cycle != 50000; ++cycle ) {
      double r1 = get_envelope_error( a, b, 0, m );
      double r2 = get_envelope_error( a, b, m, n );
      double r3 = get_sustain_error( a, b, n, length );
      double loss = r1 + r2 + r3;
      double grad_m = ( (1.f - b + ((-1.f + b)*(std::exp(a*m) + a*std::exp(a*n)*n))/std::exp(a*(m + n)))/2.f ) * loss;
      double grad_n = ( -((-1.f + b)*(std::exp(a*n) + std::exp(a*m)*(-1.f + a*(-2.f*length + m + n))))/(2.f*std::exp(a*(m + n))) ) * loss;
      m -= mopt( grad_m );
      n -= mopt( grad_n );
    }
    return std::make_tuple( m, n );
  }
  float
  approxymate_attack( double a, double b, double length ) {
    double n = 2.f/3.f * length;
    adam< double > nopt( 0.001, 0.9, 0.999 );
    for( unsigned int cycle = 0; cycle != 50000; ++cycle ) {
      double r1 = get_envelope_error( a, b, 0, n );
      double r2 = get_envelope_error( a, b, n, length );
      double loss = r1 + r2;
      double grad_n = ( (1.f - b + ((-1.f + b)*(std::exp(a*n) + a*std::exp(a*length)*length))/std::exp(a*(length + n)))/2.f ) * loss;
      n -= nopt( grad_n );
    }
    return n;
  }
  std::pair< envelope_param_keyframe_t< double >, float >
  get_attack_and_decay( const float *expected, unsigned int size, float dt, bool no_attack, bool no_sustain ) {
    const auto e = get_attack_and_decay_exp( expected, size, dt, no_sustain );
    if( no_attack && no_sustain ) {
      const float decay_mid = approxymate_attack( e.decay_a, e.decay_b, size * dt - e.highest );
      if( !std::isnan( decay_mid ) ) {
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_decay1_length( decay_mid )
            .set_decay2_length( ( size * dt - e.highest ) - decay_mid )
            .set_decay_mid_level( std::exp(-decay_mid*e.decay_a)*( 1.f - e.decay_b ) + e.decay_b )
            .set_sustain_level( 0 ),
          e.highest
        );
      }
      else {
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_decay_mid_level( 1 )
            .set_sustain_level( 1 ),
          e.highest
        );
      }
    }
    else if( no_attack && !no_sustain ) {
      const auto [decay_mid,decay_end] = approxymate_decay( e.decay_a, e.decay_b, size * dt - e.highest );
      if( !std::isnan( decay_mid ) && !std::isnan( decay_end ) ) {
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_decay1_length( decay_mid )
            .set_decay2_length( decay_end - decay_mid )
            .set_decay_mid_level( std::exp(-decay_mid*e.decay_a)*( 1.f - e.decay_b ) + e.decay_b )
            .set_sustain_level( std::exp(-decay_end*e.decay_a)*( 1.f - e.decay_b ) + e.decay_b ),
          e.highest
        );
      }
      else {
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_decay_mid_level( 1 )
            .set_sustain_level( 1 ),
          e.highest
        );
      }
    }
    else if( !no_attack && no_sustain ) {
      const float decay_mid = approxymate_attack( e.decay_a, e.decay_b, size * dt - e.highest );
      if( !std::isnan( decay_mid ) ) {
        const float attack_mid = approxymate_attack( e.attack_a, e.attack_b, e.highest );
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_attack1_length( e.highest - attack_mid )
            .set_attack2_length( attack_mid )
            .set_attack_mid_level( std::exp(-attack_mid*e.attack_a)*( 1.f - e.attack_b ) + e.attack_b )
            .set_decay1_length( decay_mid )
            .set_decay2_length( ( size * dt - e.highest ) - decay_mid )
            .set_decay_mid_level( std::exp(-decay_mid*e.decay_a)*( 1.f - e.decay_b ) + e.decay_b )
            .set_sustain_level( 0 ),
          e.highest_level
        );
      }
      else {
        const float attack_mid = approxymate_attack( e.attack_a, e.attack_b, e.highest );
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_attack1_length( e.highest - attack_mid )
            .set_attack2_length( attack_mid )
            .set_attack_mid_level( std::exp(-attack_mid*e.attack_a)*( 1.f - e.attack_b ) + e.attack_b )
            .set_decay_mid_level( 1 )
            .set_sustain_level( 1 ),
          e.highest_level
        );
      }
    }
    else {
      const auto [decay_mid,decay_end] = approxymate_decay( e.decay_a, e.decay_b, size * dt - e.highest );
      if( !std::isnan( decay_mid ) && !std::isnan( decay_end ) ) {
        const float attack_mid = approxymate_attack( e.attack_a, e.attack_b, e.highest );
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_attack1_length( e.highest - attack_mid )
            .set_attack2_length( attack_mid )
            .set_attack_mid_level( std::exp(-attack_mid*e.attack_a)*( 1.f - e.attack_b ) + e.attack_b )
            .set_decay1_length( decay_mid )
            .set_decay2_length( decay_end - decay_mid )
            .set_decay_mid_level( std::exp(-decay_mid*e.decay_a)*( 1.f - e.decay_b ) + e.decay_b )
            .set_sustain_level( std::exp(-decay_end*e.decay_a)*( 1.f - e.decay_b ) + e.decay_b ),
          e.highest_level
        );
      }
      else {
        const float attack_mid = approxymate_attack( e.attack_a, e.attack_b, e.highest );
        return std::make_pair(
          envelope_param_keyframe_t< double >()
            .set_attack1_length( e.highest - attack_mid )
            .set_attack2_length( attack_mid )
            .set_attack_mid_level( std::exp(-attack_mid*e.attack_a)*( 1.f - e.attack_b ) + e.attack_b )
            .set_decay_mid_level( 1 )
            .set_sustain_level( 1 ),
          e.highest_level
        );
      }
    }
  }
  float get_linear_envelope( const envelope_param_keyframe_t< double > &e, float x ) {
    const float a1t = e.delay_length;
    const float a2t = a1t + e.attack1_length;
    const float ht = a2t + e.attack2_length;
    const float d1t = ht + e.hold_length;
    const float d2t = d1t + e.decay1_length;
    const float st = d2t + e.decay2_length;
    if( x >= st ) return e.sustain_level;
    if( x >= d2t ) {
      const auto [tangent,shift] = get_linear_interpolation( d2t, e.decay_mid_level, st, e.sustain_level );
      return x * tangent + shift;
    }
    if( x >= d1t ) {
      const auto [tangent,shift] = get_linear_interpolation( d1t, 1, d2t, e.decay_mid_level );
      return x * tangent + shift;
    }
    if( x >= ht ) {
      return 1;
    }
    if( x >= a2t ) {
      const auto [tangent,shift] = get_linear_interpolation( a2t, e.attack_mid_level, ht, 1 );
      return x * tangent + shift;
    }
    if( x >= a1t ) {
      const auto [tangent,shift] = get_linear_interpolation( a1t, 0, a2t, e.attack_mid_level );
      return x * tangent + shift;
    }
    return 0;
  }
}
