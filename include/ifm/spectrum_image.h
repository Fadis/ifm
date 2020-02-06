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

#ifndef IFM_SPECTRUM_IMAGE_H
#define IFM_SPECTRUM_IMAGE_H
#include <vector>
#include <sndfile.h>
#include "fft.h"
namespace ifm {
  struct spectrum_image {
  public:
    spectrum_image(
      uint8_t note_,
      float sample_rate_,
      uint32_t resolution_
    );
    std::tuple< std::vector< float >, float > operator()(
      const std::vector< float > &audio
    );
    unsigned int get_width() const;
  private:
    void resample( const float *src, float *dest, float step ) const;
    static float interpolate( float a, float b, float c );
    ifm::fft_context_t fft;
    unsigned int note;
    uint32_t resolution;
    float sample_rate;
  };
}
#endif

