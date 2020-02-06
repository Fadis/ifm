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

#ifndef IFM_MIDI_SEQUENCER_H
#define IFM_MIDI_SEQUENCER_H

#include <array>

#include "midi_player.h"

namespace ifm {
  struct sequencer_state {
    sequencer_state() : track_count( 0u ), resolution( 480 ), ms_to_delta_time( 0.96 ) {}
    unsigned int track_count;
    float resolution;
    float ms_to_delta_time;
  };

  template< typename Iterator, unsigned int oper_count >
  class track_sequencer {
  public:
    track_sequencer(
      midi_player< oper_count > *player_
    ) : player( player_ ), state( nullptr ), next_event_time( std::numeric_limits< uint32_t >::max() ), cur( nullptr ), end( nullptr ) {
    }
    void load( midi_player< oper_count > *player_, sequencer_state *state_, Iterator begin, Iterator end_ ) {
      player = player_; ///
      state = state_;
      cur = begin;
      end = end_;
      next_event_time = 0u;
      next_event_time = delta_time();
    }
    bool is_end() const { return cur == end; }
    void operator()( uint32_t now ) {
      while( next_event_time <= now ) {
        event();
        next_event_time = delta_time();
      }
    }
  private:
    uint32_t delta_time() {
      uint32_t delta = 0u;
      for( ; cur != end; ++cur ) {
        delta <<= 7;
        delta += *cur & 0x7F;
        if( !( *cur & 0x80 ) ) {
          ++cur;
          return next_event_time + delta;
        }
      }
      return std::numeric_limits< uint32_t >::max();
    }
    uint32_t data_length() {
      uint32_t delta = 0u;
      for( ; cur != end; ++cur ) {
        delta <<= 7;
        delta += *cur & 0x7F;
        if( !( *cur & 0x80 ) ) {
          ++cur;
          return delta;
        }
      }
      return 0u;
    }
    void midi_event( uint8_t head ) {
      if( player->event( head ) ) return;
      for( ; cur != end; ++cur ) {
        if( player->event( *cur ) ) {
          ++cur;
          return;
        }
      }
    }
    void sysex( uint8_t ) {
      const uint32_t length = data_length();
      if( std::distance( cur, end ) <= length )
        cur = std::next( cur, length );
      else
        cur = end;
    }
    void meta_event( uint8_t ) {
      if( cur == end ) return;
      const auto event_type = *cur;
      ++cur;
      switch( event_type ) {
        case 0x2F:
          cur = end;
          break;
        case 0x51:
          {
            const auto length = data_length();
            if( length == 3u && std::distance( cur, end ) >= 3u ) {
              uint32_t beat = *cur;
              ++cur;
              beat <<= 8u;
              beat |= *cur;
              ++cur;
              beat <<= 8u;
              beat |= *cur;
              ++cur;
              state->ms_to_delta_time = float( state->resolution ) * 1000.f / float( beat );
            }
            else {
              if( std::distance( cur, end ) >= length )
                cur = std::next( cur, length );
              else
                cur = end;
            }
            break;
          }
        default:
          {
            const auto length = data_length();
            if( std::distance( cur, end ) >= length )
              cur = std::next( cur, length );
            else
              cur = end;
          }
      };
    }
    void event() {
      if( cur == end ) return;
      const auto head = *cur;
      ++cur;
      if( head == 0xF0 || head == 0xF7 )
        sysex( head );
      else if( head == 0xFF )
        meta_event( head );
      else
        midi_event( head );
    }
    midi_player< oper_count > *player;
    sequencer_state *state;
    uint32_t next_event_time;
    Iterator cur;
    Iterator end;
  };

  template< typename Iterator, unsigned int oper_count >
  class midi_sequencer {
  public:
    midi_sequencer(
     const fm_params_t< double, oper_count > &params
    ) : player( params ) {
      tracks.resize( 16, track_sequencer< Iterator, oper_count >( &player ) );
    }
    bool load( Iterator begin, Iterator end ) {
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      if( std::distance( begin, end ) < 14 ) return false;
      constexpr static const std::array< uint8_t, 8u > header_magic {{
        'M', 'T', 'h', 'd', 0, 0, 0, 6
      }};
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      if( !std::equal( header_magic.begin(), header_magic.end(), begin ) ) return false;
      auto cur = std::next( begin, header_magic.size() );
      uint16_t format = *cur;
      ++cur;
      format <<= 8;
      format |= *cur;
      ++cur;
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      if( format >= 2 ) return false;
      uint16_t track_count = *cur;
      ++cur;
      track_count <<= 8;
      track_count |= *cur;
      ++cur;
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      if( track_count > 16u ) track_count = 16u;
      state.track_count = track_count;
      uint16_t resolution = *cur;
      ++cur;
      resolution <<= 8;
      resolution |= *cur;
      ++cur;
      state.resolution = resolution;
      state.ms_to_delta_time = state.resolution * 1000.f / 500000.f;
      constexpr static const std::array< uint8_t, 4u > track_magic {{
        'M', 'T', 'r', 'k'
      }};
      now = 0.f;
      for( unsigned int i = 0u; i != track_count; ++i ) {
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if( std::distance( cur, end ) < 4 ) return false;
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if( !std::equal( track_magic.begin(), track_magic.end(), cur ) ) return false;
        cur = std::next( cur, track_magic.size() );
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if( std::distance( cur, end ) < 4 ) return false;
        uint32_t track_length = *cur;
        ++cur;
        track_length <<= 8;
        track_length |= *cur;
        ++cur;
        track_length <<= 8;
        track_length |= *cur;
        ++cur;
        track_length <<= 8;
        track_length |= *cur;
        ++cur;
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if( std::distance( cur, end ) < track_length ) return false;
        const auto track_end = std::next( cur, track_length );
        tracks[ i ].load( &player, &state, cur, track_end );
        cur = track_end;
      }
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      return true;
    }
    template< typename U >
    void operator()( U *dest ) {
      for( unsigned int i = 0u; i != state.track_count; ++i )
        tracks[ i ]( now );
      player( dest );
      now += float(synth_block_size) / float(synth_sample_rate) * 1000.f * state.ms_to_delta_time;
    }
    bool is_end() {
      return std::find_if( tracks.begin(), std::next( tracks.begin(), state.track_count ), []( const auto &t ) { return !t.is_end(); } ) == std::next( tracks.begin(), state.track_count );
    }
  private:
    midi_player< oper_count > player;
    sequencer_state state;
    std::vector< track_sequencer< Iterator, oper_count > > tracks;
    float now;
  };
}

#endif


