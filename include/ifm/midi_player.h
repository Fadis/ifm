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

#ifndef IFM_MIDI_H
#define IFM_MIDI_H

#include <array>

#include "fm.h"
#include "channel_state.h"

namespace ifm {
  template< unsigned int oper_count >
  class midi_player {
  public:
    midi_player( const fm_params_t< double, oper_count > &params ) :
      channels{{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }}, cs( params ) {}
    bool event( uint8_t v ) {
      if( v < 0x80 ) return (this->*state)( v );
      else return new_event( v );
    }
    template< typename U >
    void operator()( U *dest ) {
      cs( dest );
    }
  private:
    bool waiting_for_event( uint8_t ) { return true; }
    bool note_off_key_number( uint8_t v ) { 
      message_buffer[ 0 ] = v;
      state = &midi_player::note_off_velocity;
      return false;
    }
    bool note_off_velocity( uint8_t ) {
      cs.note_off( channel, note_number_t( message_buffer[ 0 ] ) );
      state = &midi_player::note_off_key_number;
      return true;
    }
    bool note_on_key_number( uint8_t v ) { 
      message_buffer[ 0 ] = v;
      state = &midi_player::note_on_velocity;
      return false;
    }
    bool note_on_velocity( uint8_t v ) {
      if( channel != 10 ) {
        const note_number_t scale = note_number_t( message_buffer[ 0 ] - 12 );
        cs.note_on( channel, scale, velocity_t( v ) );
      }
      state = &midi_player::note_on_key_number;
      return true;
    }
    bool polyphonic_key_pressure_key_number( uint8_t ) { 
      state = &midi_player::polyphonic_key_pressure_value;
      return false;
    }
    bool polyphonic_key_pressure_value( uint8_t ) {
      state = &midi_player::polyphonic_key_pressure_key_number;
      return true;
    }
    bool control_change_key( uint8_t v ) {
      if( v == 1 )
        state = &midi_player::set_modulation;
      else if( v == 7 )
        state = &midi_player::set_volume;
      else if( v == 10 )
        state = &midi_player::set_pan;
      else if( v == 11 )
        state = &midi_player::set_expression;
      else if( v == 64 )
        state = &midi_player::set_dumper_pedal;
      else if( v == 121 )
        state = &midi_player::reset;
      else if( v == 123 )
        state = &midi_player::all_notes_off;
      else
        state = &midi_player::unknown_control;
      return false;
    }
    bool program_change( uint8_t ) {
      state = &midi_player::program_change;
      return true;
    }
    bool channel_pressure( uint8_t ) {
      //state = &midi_player::channel_pressure;
      return true;
    }
    bool pitch_bend_lower( uint8_t v ) {
      message_buffer[ 0 ] = v;
      state = &midi_player::pitch_bend_higher;
      return false;
    }
    bool pitch_bend_higher( uint8_t v ) {
      channels[ channel ].pitch_bend = ( ( ( int( v ) << 7 )|( int( message_buffer[ 0 ] ) ) ) - 8192 )/8191.f;
      channels[ channel ].final_pitch = channels[ channel ].pitch_bend * channels[ channel ].pitch_sensitivity;
      state = &midi_player::pitch_bend_lower;
      return true;
    }
    bool set_modulation( uint8_t v ) { // cc 1
      channels[ channel ].modulation = int( v )/127.f;
      state = &midi_player::control_change_key;
      return true;
    }
    bool set_volume( uint8_t v ) { // cc 7
      channels[ channel ].volume = int( v )/127.f;
      channels[ channel ].final_volume = channels[ channel ].volume * channels[ channel ].expression;
      state = &midi_player::control_change_key;
      return true;
    }
    bool set_pan( uint8_t v ) { // cc 10
      channels[ channel ].pan = ( int( v ) - 64 )/63.f;
      state = &midi_player::control_change_key;
      return true;
    }
    bool set_expression( uint8_t v ) { // cc 11
      channels[ channel ].expression = int( v )/127.f;
      channels[ channel ].final_volume = channels[ channel ].volume * channels[ channel ].expression;
      state = &midi_player::control_change_key;
      return true;
    }
    bool set_dumper_pedal( uint8_t v ) { // cc 64
      channels[ channel ].sustain = v >= 64;
      state = &midi_player::control_change_key;
      return true;
    }
    bool reset( uint8_t ) { // cc 121
      cs.reset();
      std::for_each( channels.begin(), channels.end(), []( channel_state &channel ) { channel.reset(); } );
      state = &midi_player::control_change_key;
      return true;
    }
    bool all_notes_off( uint8_t ) { // cc 123
      state = &midi_player::control_change_key;
      return true;
    }
    bool unknown_control( uint8_t ) {
      state = &midi_player::control_change_key;
      return true;
    }
    bool new_event( uint8_t v ) {
      constexpr static const std::array< bool(midi_player::*)( uint8_t ), 8u > initial_states{{
        &midi_player::note_off_key_number,
        &midi_player::note_on_key_number,
        &midi_player::polyphonic_key_pressure_key_number,
        &midi_player::control_change_key,
        &midi_player::program_change,
        &midi_player::channel_pressure,
        &midi_player::pitch_bend_lower,
        &midi_player::waiting_for_event
      }};
      const uint8_t event = ( v >> 4 ) & 0x07;
      channel = v & 0x0F;
      state = initial_states[ event ];
      return false;
    }
    bool(midi_player::*state)( uint8_t );
    channel_t channel;
    std::array< channel_state, 16u > channels;
    channels_t< float, oper_count > cs;
    std::array< uint8_t, 16u > message_buffer;
  };
}

#endif


