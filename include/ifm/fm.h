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

#ifndef IFM_FM_H
#define IFM_FM_H
#include <cstdint>
#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <charconv>
#include <boost/container/flat_map.hpp>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include "setter.h"

namespace ifm {
  constexpr unsigned int synth_block_size = 32u;
  constexpr unsigned int  synth_sample_rate = 44100u;
  using note_number_t = uint8_t;
  using channel_t = uint8_t;
  using velocity_t = uint8_t;
  constexpr unsigned int max_note_number = 128u;
  template< typename T >
  struct envelope_param_keyframe_t {
    envelope_param_keyframe_t() :
      delay_length( 0 ),
      attack1_length( 0 ),
      attack2_length( 0 ),
      attack_mid_level( 0 ),
      hold_length( 0 ),
      decay1_length( 0 ),
      decay2_length( 0 ),
      decay_mid_level( 1 ),
      sustain_level( 1 ),
      release_length( 0 ) {}
    IFM_SET_SMALL_VALUE( delay_length )
    IFM_SET_SMALL_VALUE( attack1_length )
    IFM_SET_SMALL_VALUE( attack2_length )
    IFM_SET_SMALL_VALUE( attack_mid_level )
    IFM_SET_SMALL_VALUE( hold_length )
    IFM_SET_SMALL_VALUE( decay1_length )
    IFM_SET_SMALL_VALUE( decay2_length )
    IFM_SET_SMALL_VALUE( decay_mid_level )
    IFM_SET_SMALL_VALUE( sustain_level )
    IFM_SET_SMALL_VALUE( release_length )
    T delay_length;
    T attack1_length;
    T attack2_length;
    T attack_mid_level;
    T hold_length;
    T decay1_length;
    T decay2_length;
    T decay_mid_level;
    T sustain_level;
    T release_length;
  };
  template< typename T >
  using envelope_params_t = boost::container::flat_map< note_number_t, envelope_param_keyframe_t< T > >;
  constexpr unsigned int envelope_samples = 128u;
  template< typename T >
  using envelope_image_t = std::array< T, max_note_number * envelope_samples >;
  template< typename T, typename U, typename V >
  envelope_param_keyframe_t< T > interpolate(
    const envelope_param_keyframe_t< U > &l,
    const envelope_param_keyframe_t< U > &h,
    V pos
  ) {
    V ipos = ( T( 1 ) - pos );
    return envelope_param_keyframe_t< T >()
      .set_delay_length( l.delay_length * ipos + h.delay_length )
      .set_attack1_length( l.attack1_length * ipos + h.attack1_length * pos )
      .set_attack2_length( l.attack2_length * ipos + h.attack2_length * pos )
      .set_attack_mid_level( l.attack_mid_level * ipos + h.attack_mid_level * pos )
      .set_hold_length( l.hold_length * ipos + h.hold_length * pos )
      .set_decay1_length( l.decay1_length * ipos + h.decay1_length * pos )
      .set_decay2_length( l.decay2_length * ipos + h.decay2_length * pos )
      .set_decay_mid_level( l.decay_mid_level * ipos + h.decay_mid_level * pos )
      .set_sustain_level( l.sustain_level * ipos + h.sustain_level * pos )
      .set_release_length( l.release_length * ipos + h.release_length * pos );
  }

  struct invalid_configuration {};

  note_number_t parse_note_number( const std::string &v );
  envelope_param_keyframe_t< double > load_envelope_param_keyframe( const nlohmann::json &v );
  envelope_params_t< double > load_envelope_params( const nlohmann::json &v );
  nlohmann::json store_envelope_param_keyframe( const envelope_param_keyframe_t< double > &v );

  template< unsigned int oper_count >
  std::array< envelope_params_t< double >, oper_count > load_envelope_params_array( const nlohmann::json &v ) {
    std::array< envelope_params_t< double >, oper_count > temp;
    if( v.size() != oper_count ) throw invalid_configuration {};
    for( unsigned int i = 0; i != oper_count; ++i )
      temp[ i ] = load_envelope_params( v[ i ] );
    return temp;
  }

  template< typename T >
  class envelope_t {
  public:
    template< typename U >
    envelope_t(
      const envelope_params_t< U > &params,
      note_number_t note
    ) {
      auto h = params.upper_bound( note );
      auto l = h == params.begin() ? h : std::prev( h );
      if( l == params.end() ) l = std::prev( params.end() );
      if( h == params.end() ) h = l;
      config = interpolate< T >( l->second, h->second, l == h ? T( 0 ) : T( note - l->first )/T( h->first - l->first ) );
      attack1_tangent = config.attack_mid_level / config.attack1_length / T( synth_sample_rate );
      attack2_tangent = ( 1 - config.attack_mid_level ) / config.attack2_length / T( synth_sample_rate );
      decay1_tangent = -( 1 - config.decay_mid_level ) / config.decay1_length / T( synth_sample_rate );
      decay2_tangent = -( config.decay_mid_level - config.sustain_level ) / config.decay2_length / T( synth_sample_rate );
      release_tangent = -1 / config.release_length / T( synth_sample_rate );
      now = 0;
      level = 0;
      if( config.delay_length > 0 )
        state = &envelope_t::delay;
      if( config.attack1_length > 0 )
        state = &envelope_t::attack1;
      else if( config.attack2_length > 0 ) {
        level = config.attack_mid_level;
        state = &envelope_t::attack2;
      }
      else if( config.hold_length > 0 ) {
        level = 1;
        state = &envelope_t::hold;
      }
      else if( config.decay1_length > 0 ) {
        level = 1;
        state = &envelope_t::decay1;
      }
      else if( config.decay2_length > 0 ) {
        level = config.decay_mid_level;
        state = &envelope_t::decay2;
      }
      else if( config.sustain_level > 0 ) {
        level = config.sustain_level;
        state = &envelope_t::sustain;
      }
      else state = &envelope_t::end;
    }
    void operator()( T *dest ) {
      ( this->*state )( dest );
    }
    void note_off() {
      if( level != 0 )
        state = &envelope_t::release;
      else
        state = &envelope_t::end;
    }
    bool is_end() const { return state == &envelope_t::end; }
  private:
    void delay( T *dest ) {
      std::fill( dest, dest + synth_block_size, 0 );
      now += T( synth_block_size ) / T( synth_sample_rate );
      if( config.delay_length < now ) {
        now = 0;
        level = 0;
        if( config.attack1_length > 0 )
          state = &envelope_t::attack1;
        else if( config.attack2_length > 0 ) {
          level = config.attack_mid_level;
          state = &envelope_t::attack2;
        }
        else if( config.hold_length > 0 ) {
          level = 1;
          state = &envelope_t::hold;
        }
        else if( config.decay1_length > 0 ) {
          level = 1;
          state = &envelope_t::decay1;
        }
        else if( config.decay2_length > 0 ) {
          level = config.decay_mid_level;
          state = &envelope_t::decay2;
        }
        else if( config.sustain_level > 0 ) {
          level = config.sustain_level;
          state = &envelope_t::sustain;
        }
        else state = &envelope_t::end;
      }
    }
    void attack1( T *dest ) {
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        dest[ i ] = level;
        level += attack1_tangent;
        if( level > 1 ) level = 1;
        now += T( 1 ) / T( synth_sample_rate );
      }
      if( config.attack1_length < now ) {
        now = 0;
        if( config.attack2_length > 0 ) {
          state = &envelope_t::attack2;
        }
        else if( config.hold_length > 0 ) {
          level = 1;
          state = &envelope_t::hold;
        }
        else if( config.decay1_length > 0 ) {
          level = 1;
          state = &envelope_t::decay1;
        }
        else if( config.decay2_length > 0 ) {
          level = config.decay_mid_level;
          state = &envelope_t::decay2;
        }
        else if( config.sustain_level > 0 ) {
          level = config.sustain_level;
          state = &envelope_t::sustain;
        }
        else state = &envelope_t::end;
      }
    }
    void attack2( T *dest ) {
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        dest[ i ] = level;
        level += attack2_tangent;
        if( level > 1 ) level = 1;
        now += T( 1 ) / T( synth_sample_rate );
      }
      if( config.attack2_length < now ) {
        now = 0;
        if( config.hold_length > 0 ) {
          state = &envelope_t::hold;
        }
        else if( config.decay1_length > 0 ) {
          state = &envelope_t::decay1;
        }
        else if( config.decay2_length > 0 ) {
          level = config.decay_mid_level;
          state = &envelope_t::decay2;
        }
        else if( config.sustain_level > 0 ) {
          level = config.sustain_level;
          state = &envelope_t::sustain;
        }
        else state = &envelope_t::end;
      }
    }
    void hold( T *dest ) {
      std::fill( dest, dest + synth_block_size, level );
      now += T( synth_block_size ) / T( synth_sample_rate );
      if( config.hold_length < now ) {
        now = 0;
        if( config.decay1_length > 0 ) {
          state = &envelope_t::decay1;
        }
        else if( config.decay2_length > 0 ) {
          level = config.decay_mid_level;
          state = &envelope_t::decay2;
        }
        else if( config.sustain_level > 0 ) {
          level = config.sustain_level;
          state = &envelope_t::sustain;
        }
        else state = &envelope_t::end;
      }
    }
    void decay1( T *dest ) {
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        dest[ i ] = level;
        level += decay1_tangent;
        if( level < config.sustain_level ) level = config.sustain_level;
        now += T( 1 ) / T( synth_sample_rate );
      }
      if( config.decay1_length < now ) {
        now = 0;
        if( config.decay2_length > 0 ) {
          state = &envelope_t::decay2;
        }
        else if( config.sustain_level > 0 ) {
          state = &envelope_t::sustain;
        }
        else state = &envelope_t::end;
      }
    }
    void decay2( T *dest ) {
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        dest[ i ] = level;
        level += decay2_tangent;
        if( level < config.sustain_level ) level = config.sustain_level;
        now += T( 1 ) / T( synth_sample_rate );
      }
      if( config.decay2_length < now ) {
        now = 0;
        if( config.sustain_level > 0 ) {
          state = &envelope_t::sustain;
        }
        else {
          state = &envelope_t::end;
        }
      }
    }
    void sustain( T *dest ) {
      std::fill( dest, dest + synth_block_size, level );
    }
    void release( T *dest ) {
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        dest[ i ] = level;
        level += release_tangent;
        if( level < 0 ) level = 0;
      }
      if( level <= 0 ) state = &end;
    }
    void end( T *dest ) {
      std::fill( dest, dest + synth_block_size, 0 );
    }
    void ( envelope_t::*state )( T* );
    T now;
    T level;
    envelope_param_keyframe_t< T > config;
    T attack1_tangent;
    T attack2_tangent;
    T decay1_tangent;
    T decay2_tangent;
    T release_tangent;
  };

  template< typename T, unsigned int oper_count >
  class envelopes_t {
  public:
    template< typename U >
    envelopes_t(
      const std::array< envelope_params_t< U >, oper_count > &params,
      note_number_t note
    ) : envelope(
      init_envelope( params, note, std::make_index_sequence< oper_count >() )
    ) {}
    void operator()( unsigned int operator_index, T *dest ) {
      envelope[ operator_index ]( dest );
    }
    void note_off() {
      std::for_each( envelope.begin(), envelope.end(), []( auto &v ) { v.note_off(); } );
    }
    bool is_end() const {
      return std::find_if( envelope.begin(), envelope.end(), []( auto &v ) { return !v.is_end(); } ) == envelope.end();
    }
  private:
    template< typename U, typename I, I ... seq >
    static std::array< envelope_t< T >, oper_count > init_envelope(
      const std::array< envelope_params_t< U >, oper_count > &params,
      note_number_t note,
      std::index_sequence< seq... >
    ) {
      return std::array< envelope_t< T >, oper_count >{{
        envelope_t< T >( params[ seq ], note )...
      }};
    }
    std::array< envelope_t< T >, oper_count > envelope;
  };

  template< typename T, unsigned int oper_count >
  using weight_param_keyframe_t = std::array< T, oper_count * ( oper_count + 1 ) >;
  template< typename T, unsigned int oper_count >
  using weight_params_t = boost::container::flat_map< note_number_t, weight_param_keyframe_t< T, oper_count > >;

  template< typename T, typename U, typename V, size_t n >
  std::array< T, n > interpolate(
    const std::array< U, n > &l,
    const std::array< U, n > &h,
    V pos
  ) {
    std::array< T, n > temp;
    for( unsigned int i = 0; i != temp.size(); ++i )
      temp[ i ] = l[ i ] * ( V( 1 ) - pos ) + h[ i ] * pos;
    return temp;
  }

  template< unsigned int oper_count >
  weight_param_keyframe_t< double, oper_count > load_weight_param_keyframe( const nlohmann::json &v ) {
    weight_param_keyframe_t< double,  oper_count > temp;
    if( v.size() != oper_count * ( oper_count + 1 ) ) throw invalid_configuration {};
    for( unsigned int i = 0; i != oper_count * ( oper_count + 1 ); ++i ) {
      temp[ i ] = v[ i ];
      ++i;
    }
    return temp;
  }

  template< unsigned int oper_count >
  weight_params_t< double, oper_count > load_weight_params( const nlohmann::json &v ) {
    weight_params_t< double, oper_count > temp;
    for( const auto &[key,value]: v.items() ) {
      auto note = parse_note_number( key );
      if( note < 0 ) throw invalid_configuration {};
      if( note > 128 ) throw invalid_configuration {};
      temp[ note ] = load_weight_param_keyframe< oper_count >( value );
    }
    return temp;
  }

  template< unsigned int oper_count >
  std::array< double, oper_count > load_freqs( const nlohmann::json &v ) {
    std::array< double, oper_count > temp;
    if( v.size() != oper_count ) throw invalid_configuration {};
    for( unsigned int i = 0; i != oper_count; ++i ) {
      temp[ i ] = v[ i ];
    }
    return temp;
  }

  template< typename T, unsigned int n >
  class weight_t {
  public:
    template< typename U >
    weight_t(
      const weight_params_t< U, n > &params,
      note_number_t note
    ) {
      auto h = params.upper_bound( note );
      auto l = h == params.begin() ? h : std::prev( h );
      if( l == params.end() ) l = std::prev( params.end() );
      if( h == params.end() ) h = l;
      config = interpolate< T >( l->second, h->second, l == h ? T( 0 ) : T( note - l->first )/T( h->first - l->first ) );
    }
    T operator()( unsigned int i ) const {
      return config[ i ];
    }
  private:
    weight_param_keyframe_t< T, n > config;
  };

  template< typename T, unsigned int oper_count >
  struct fm_params_t {
    fm_params_t() {
      std::fill( freq.begin(), freq.end(), 0 );
    }
    IFM_SET_LARGE_VALUE( envelope )
    IFM_SET_LARGE_VALUE( freq )
    IFM_SET_LARGE_VALUE( weight )
    std::array< envelope_params_t< T >, oper_count > envelope;
    std::array< T, oper_count > freq;
    weight_params_t< T, oper_count > weight;
  };

  template< unsigned int oper_count >
  fm_params_t< double, oper_count > load_fm_params( const nlohmann::json &v ) {
    return fm_params_t< double, oper_count >()
      .set_envelope( load_envelope_params_array< oper_count >( v[ "envelope" ] ) )
      .set_freq( load_freqs< oper_count >( v[ "freq" ] ) )
      .set_weight( load_weight_params< oper_count >( v[ "weight" ] ) );
  }

  template< typename T, unsigned int oper_count >
  class fm_t {
  public:
    template< typename U >
    fm_t(
      const fm_params_t< U, oper_count > &params,
      note_number_t note,
      velocity_t velocity_ = 128
    ) : weight( params.weight, note ), envelope( params.envelope, note ), velocity( T( velocity_ )/T(128) ) {
      const T base_freq = std::exp2( ( ( T( note ) +  T( 3 ) ) / T( 12 ) ) ) * T( 6.875 );
      std::transform( params.freq.begin(), params.freq.end(), tangent.begin(), [&]( T v ) { return v * base_freq * T( 2 ) * T( M_PI ) / T( synth_sample_rate ); } );
      std::fill( shift.begin(), shift.end(), 0.f );
    }
    template< typename U >
    void operator()( U *dest ) {
      std::array< std::array< T, synth_block_size >, oper_count > e;
      for( unsigned int operator_index = 0; operator_index != oper_count; ++operator_index )
        envelope( operator_index, e[ operator_index ].data() );
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        T sum = 0;
        for( unsigned int to_operator_index = 0; to_operator_index != oper_count; ++to_operator_index ) {
          T drift = 0;
          for( unsigned int from_operator_index = 0; from_operator_index != oper_count; ++from_operator_index ) {
            drift += weight( to_operator_index + from_operator_index * oper_count ) * prev[ from_operator_index ];
          }
          prev[ to_operator_index ] = e[ to_operator_index ][ i ] * std::sin( shift[ to_operator_index ] + drift );
          sum += weight( to_operator_index + oper_count * oper_count ) * prev[ to_operator_index ];
          shift[ to_operator_index ] += tangent[ to_operator_index ];
        }
        dest[ i ] = sum * velocity;
      }
    }
    void note_off() {
      envelope.note_off();
    }
    bool is_end() const {
      return envelope.is_end();
    }
  private:
    std::array< T, oper_count > tangent;
    std::array< T, oper_count > shift;
    std::array< T, oper_count > prev;
    weight_t< T, oper_count > weight;
    envelopes_t< T, oper_count > envelope;
    T velocity;
  };
  template< typename T, unsigned int oper_count >
  class polyphony_t {
  public:
    polyphony_t(
      const fm_params_t< double, oper_count > &params_
    ) : params( params_ ) {}
    void note_on( note_number_t note, velocity_t velocity ) {
      active.erase( note );
      active.insert( std::make_pair( note, fm_t< T, oper_count >( params, note, velocity ) ) );
    }
    void note_off( note_number_t note ) {
      active.erase( note );
    }
    template< typename U >
    void operator()( U *dest ) {
      std::array< U, synth_block_size > b;
      std::fill( dest, dest + synth_block_size, 0 );
      for( auto &v: active ) {
        v.second( b.data() );
        for( unsigned int i = 0; i != synth_block_size; ++i ) dest[ i ] += b[ i ];
      }
      for( auto iter = active.begin(); iter != active.end(); ) {
        if( iter->second.is_end() ) iter = active.erase( iter );
        else ++iter;
      }
    }
    void reset() {
      active.clear();
    }
  private:
    fm_params_t< double, oper_count > params;
    std::unordered_map< int, fm_t< T, oper_count > > active;
  };
  template< typename T, unsigned int oper_count >
  class channels_t {
  public:
    channels_t(
      const fm_params_t< double, oper_count > &params
    ) : scale( 0.8 ), keep( 0 ) {
      for( unsigned int i = 0; i != 16; ++i ) channels.emplace_back( params );
    }
    void note_on( channel_t channel_id, note_number_t note, velocity_t velocity ) {
      channels[ channel_id ].note_on( note, velocity );
    }
    void note_off( channel_t channel_id, note_number_t note ) {
      channels[ channel_id ].note_off( note );
    }
    template< typename U >
    void operator()( U *dest ) {
      std::array< U, synth_block_size > b;
      std::fill( dest, dest + synth_block_size, 0 );
      for( auto &v: channels ) {
        v( b.data() );
        for( unsigned int i = 0; i != synth_block_size; ++i ) dest[ i ] += b[ i ] * 0.125f;
      }
      std::array< U, synth_block_size > s;
      auto initial_scale = scale;
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        if( std::abs( dest[ i ] ) * scale > T( 0.8 ) ) {
          scale = std::abs( T(0.8)/(dest[ i ]) );
          for( unsigned int j = 0; j != i; ++j )
            s[ j ] = initial_scale * ( 1 - T(j)/T(i) ) + scale * T(j)/T(i);
          keep = 1000;
        }
        else if( keep ) {
          --keep;
        }
        else if( scale < T( 0.8 ) ) {
          scale += T(1)/synth_sample_rate;
        }
        s[ i ] = scale;
      }
      for( unsigned int i = 0; i != synth_block_size; ++i ) {
        dest[ i ] *= s[ i ];
      }
    }
    void reset() {
      for( auto &c: channels ) c.reset();
    }
  private:
    std::vector< polyphony_t< T, oper_count > > channels;
    T scale;
    int keep;
  };
}
#endif

