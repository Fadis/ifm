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
#include <algorithm>
#include <iterator>
#include <charconv>
#include <boost/container/flat_map.hpp>
#include <nlohmann/json.hpp>
#include "ifm/fm.h"
namespace ifm {
  note_number_t parse_note_number( const std::string &v ) {
    unsigned int n;
    auto [ptr, ec] = std::from_chars( v.data(), v.data() + v.size(), n );
    if( !( ptr == v.data() + v.size() && ec == std::errc{} ) ) throw invalid_configuration();
    return n;
  }

  envelope_param_keyframe_t< double > load_envelope_param_keyframe( const nlohmann::json &v ) {
    auto delay = ( v.find( "delay" ) != v.end() ) ? double( v[ "delay" ] ) : 0.0;
    auto attack1 = ( v.find( "attack1" ) != v.end() ) ? double( v[ "attack1" ] ) : 0.0;
    auto attack2 = ( v.find( "attack2" ) != v.end() ) ? double( v[ "attack2" ] ) : 0.0;
    auto attack_level = ( v.find( "attack_level" ) != v.end() ) ? double( v[ "attack_level" ] ) : 0.0;
    auto hold = ( v.find( "hold" ) != v.end() ) ? double( v[ "hold" ] ) : 0.0;
    auto decay1 = ( v.find( "decay1" ) != v.end() ) ? double( v[ "decay1" ] ) : 0.0;
    auto decay2 = ( v.find( "decay2" ) != v.end() ) ? double( v[ "decay2" ] ) : 0.0;
    auto decay_level = ( v.find( "decay_level" ) != v.end() ) ? double( v[ "decay_level" ] ) : 0.0;
    auto sustain_level = ( v.find( "sustain_level" ) != v.end() ) ? double( v[ "sustain_level" ] ) : 0.0;
    auto release = ( v.find( "release" ) != v.end() ) ? double( v[ "release" ] ) : 0.0;
    if( delay < 0 ) throw invalid_configuration {};
    if( attack1 < 0 ) throw invalid_configuration {};
    if( attack2 < 0 ) throw invalid_configuration {};
    if( attack_level < 0 ) throw invalid_configuration {};
    if( attack_level > 1 ) throw invalid_configuration {};
    if( hold < 0 ) throw invalid_configuration {};
    if( decay1 < 0 ) throw invalid_configuration {};
    if( decay2 < 0 ) throw invalid_configuration {};
    if( decay_level < 0 ) throw invalid_configuration {};
    if( decay_level > 1 ) throw invalid_configuration {};
    if( sustain_level < 0 ) throw invalid_configuration {};
    if( sustain_level > 1 ) throw invalid_configuration {};
    if( release < 0 ) throw invalid_configuration {};
    return envelope_param_keyframe_t< double >()
      .set_delay_length( delay )
      .set_attack1_length( attack1 )
      .set_attack2_length( attack2 )
      .set_attack_mid_level( attack_level )
      .set_hold_length( hold )
      .set_decay1_length( decay1 )
      .set_decay2_length( decay2 )
      .set_decay_mid_level( decay_level )
      .set_sustain_level( sustain_level )
      .set_release_length( release );
  }
  envelope_params_t< double > load_envelope_params( const nlohmann::json &v ) {
    envelope_params_t< double > temp;
    for( const auto &[key,value]: v.items() ) {
      auto note = parse_note_number( key );
      //if( note  < 0 ) throw invalid_configuration {};
      if( note > 128 ) throw invalid_configuration {};
      temp[ note ] = load_envelope_param_keyframe( value );
    }
    return temp;
  }
  nlohmann::json store_envelope_param_keyframe( const envelope_param_keyframe_t< double > &v ) {
    std::map< std::string, float > temp;
    temp[ "delay" ] = v.delay_length;
    temp[ "attack1" ] = v.attack1_length;
    temp[ "attack2" ] = v.attack2_length;
    temp[ "attack_level" ] = v.attack_mid_level;
    temp[ "hold" ] = v.hold_length;
    temp[ "decay1" ] = v.decay1_length;
    temp[ "decay2" ] = v.decay2_length;
    temp[ "decay_level" ] = v.decay_mid_level;
    temp[ "sustain_level" ] = v.sustain_level;
    temp[ "release" ] = v.release_length;
    return temp;
  }
}

