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
#include "ifm/load_monoral.h"

namespace ifm {
std::tuple< std::vector< float >, float > load_monoral( const std::string &filename, bool norm ) {
  SF_INFO info;
  info.frames = 0;
  info.samplerate = 0;
  info.channels = 0;
  info.format = 0;
  info.sections = 0;
  info.seekable = 0;
  auto audio_file = sf_open( filename.c_str(), SFM_READ, &info );
  if( !audio_file ) {
    std::cerr << "Unable to open audio file" << std::endl;
    throw -1;
  }
  std::vector< int16_t > monoral;
  std::vector< int16_t > multi_channel;
  multi_channel.resize( info.frames * info.channels );
  const auto read_count = sf_read_short( audio_file, multi_channel.data(), info.frames * info.channels );
  if( read_count == 0u ) {
    std::cerr << "Unable to read audio file" << std::endl;
    throw -1;
  }
  sf_close( audio_file );
  if( read_count % info.channels != 0 ) {
    std::cerr << "Invalid audio file" << std::endl;
    throw -1;
  }
  multi_channel.resize( read_count );
  monoral.reserve( read_count / info.channels );
  for( auto iter = multi_channel.begin(); iter != multi_channel.end(); ) {
    auto next_iter = std::next( iter, info.channels );
    monoral.emplace_back( std::accumulate( iter, next_iter, int16_t( 0 ) ) / info.channels );
    iter = next_iter;
  }
  std::vector< float > in_float( monoral.size() );
  std::transform( monoral.begin(), monoral.end(), in_float.begin(), []( int16_t v ) { return v / 32767.f; } );
  const auto max_iter = std::max_element( in_float.begin(), in_float.end() );
  const auto min_iter = std::min_element( in_float.begin(), in_float.end() );
  const float max =
    std::max(
      ( max_iter != in_float.end() ) ? std::abs( *max_iter ) : int16_t( 0 ),
      ( min_iter != in_float.end() ) ? std::abs( *min_iter ) : int16_t( 0 )
    );
  if( norm ) {
    if( max != 0 ) {
      float scale = 1.f / max;
      for( auto &elem: in_float )
        elem *= scale;
    }
  }
  return std::make_tuple( in_float, info.samplerate );
}
}
