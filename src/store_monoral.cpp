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

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sndfile.h>
#include "ifm/store_monoral.h"

namespace ifm {
void store_monoral( const std::string &filename, const std::vector< float > &samples, unsigned int sample_rate ) {
  SF_INFO info;
  info.frames = 0;
  info.samplerate = sample_rate;
  info.channels = 1;
  info.format = SF_FORMAT_WAV|SF_FORMAT_PCM_16;
  info.sections = 0;
  info.seekable = 0;
  const auto max_iter = std::max_element( samples.begin(), samples.end() );
  const auto min_iter = std::min_element( samples.begin(), samples.end() );
  const float max =
    std::max(
      ( max_iter != samples.end() ) ? std::abs( *max_iter ) : int16_t( 0 ),
      ( min_iter != samples.end() ) ? std::abs( *min_iter ) : int16_t( 0 )
    );
  float scale = 0.8f / max;
  auto audio_file = sf_open( filename.c_str(), SFM_WRITE, &info );
  if( !audio_file ) {
    std::cerr << "Unable to open audio file" << std::endl;
    throw -1;
  }
  std::vector< int16_t > in_int;
  std::transform( samples.begin(), samples.end(), std::back_inserter( in_int ), [&]( const auto &v ) { return int16_t( v * 32767 * scale ); } );
  const auto write_count = sf_write_short( audio_file, in_int.data(), in_int.size() );
  if( write_count == 0u ) {
    std::cerr << sf_strerror( audio_file ) << std::endl;
    throw -1;
  }
  sf_close( audio_file );
}
}
