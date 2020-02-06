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

#ifndef IFM_CHANNEL_STATE_H
#define IFM_CHANNEL_STATE_H

#include "fm.h"

namespace ifm {
  struct channel_state {
    channel_state( channel_t index_ ) : index( index_ ), modulation( 0 ), volume( 1 ), expression( 1 ), final_volume( 1 ), pitch_bend( 0 ), pitch_sensitivity( 1 ), final_pitch( 0 ), pan( 0 ) {}
    void reset() {
      modulation = 0;
      volume = 1;
      expression = 1;
      final_volume = 1;
      pitch_bend = 0;
      pitch_sensitivity = 1;
      final_pitch = 0;
      pan = 0;
    }
    channel_t index;
    float modulation;
    float volume;
    float expression;
    float final_volume;
    float pitch_bend;
    float pitch_sensitivity;
    float final_pitch;
    float pan;
    bool sustain;
  };
}

#endif

