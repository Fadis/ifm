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

#ifndef IFM_EXP_MATCH_H
#define IFM_EXP_MATCH_H
#include <tuple>
#include "setter.h"
#include "fm.h"
namespace ifm {
  std::tuple< float, float >
  exp_match( const float *expected, unsigned int size, float dt );
  struct exp_envelope_params_t {
    exp_envelope_params_t() : attack_a( 0 ), attack_b( 0 ), decay_a( 0 ), decay_b( 0 ), highest( 0 ) {}
    IFM_SET_SMALL_VALUE( attack_a )
    IFM_SET_SMALL_VALUE( attack_b )
    IFM_SET_SMALL_VALUE( decay_a )
    IFM_SET_SMALL_VALUE( decay_b )
    IFM_SET_SMALL_VALUE( highest )
    IFM_SET_SMALL_VALUE( highest_level )
    float attack_a;
    float attack_b;
    float decay_a;
    float decay_b;
    float highest;
    float highest_level;
  };
  exp_envelope_params_t
  get_attack_and_decay_exp( const float *expected, unsigned int size, float dt, bool damped );
  std::pair< envelope_param_keyframe_t< double >, float >
  get_attack_and_decay( const float *expected, unsigned int size, float dt, bool no_attack, bool no_sustain );
  std::tuple< float, float >
  approxymate_decay( double a, double b, double length );
  float
  approxymate_attack( double a, double b, double length );
  float get_exp_envelope( float a, float b, float x );
  std::pair< float, float > get_linear_interpolation( float x0, float y0, float x1, float y1 );
  float get_linear_envelope( const envelope_param_keyframe_t< double > &e, float x );
}
#endif

