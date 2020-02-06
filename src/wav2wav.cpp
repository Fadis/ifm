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
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <boost/program_options.hpp>
#include <sndfile.h>
#include "ifm/fft.h"
#include "ifm/load_monoral.h"
#include "ifm/spectrum_image.h"
#include "ifm/store_monoral.h"

int main( int argc, char* argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("input,i", boost::program_options::value<std::string>(),  "入力ファイル")
    ("output,o", boost::program_options::value<std::string>(),  "出力ファイル")
    ("note,n", boost::program_options::value<int>()->default_value(60),  "音階")
    ("resolution,r", boost::program_options::value<int>()->default_value(13),  "分解能");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") || !params.count("input") || !params.count("output") ) {
    std::cout << options << std::endl;
    return 0;
  }
  const std::string input_filename = params["input"].as<std::string>();
  const std::string output_filename = params["output"].as<std::string>();
  const auto [audio,sample_rate] = ifm::load_monoral( input_filename, true );
  const auto note = params["note"].as<int>();
  ifm::spectrum_image conv( note, sample_rate, 1 << params["resolution"].as<int>() );
  auto [image,delay] = conv( audio );
  std::cout << "delay: " << delay << std::endl;
  unsigned int width = conv.get_width();
  unsigned int active_width;
  std::vector< float > output_wave;
  float resolution = 1 << params["resolution"].as<int>();
  active_width = width / 24; 
  double time = 0;
  unsigned int spb = sample_rate / 100;
  output_wave.resize( spb * ( image.size() / width - 1 ), 0.f );
  const float transformed_freq = float( 24 ) / resolution * sample_rate;
  const float original_freq = std::exp2( ( ( float( note ) +  3.f ) / 12.f ) ) * 6.875f;
  const float ratio = original_freq / transformed_freq;
  std::cout << "ratio : " << original_freq << "/" << transformed_freq << std::endl;
  std::mt19937 rand_gen( 0 );
  std::uniform_real_distribution rand_dist( -1.0, 1.0 );
  for( unsigned int y = 0; y != image.size() / width - 1; ++y, time += float( spb ) / sample_rate ) {
    const auto begin = std::next( output_wave.begin(), y * spb );
    for( unsigned int x = 0; x != active_width; ++x ) {
      const auto n = ( x + 1 ) * 24;
      const auto begin_level = image[ n + y * width ];
      const auto end_level = image[ n + ( y + 1 ) * width ];
      auto iter = begin;
      auto delta = 0.f;
      for( unsigned int s = 0; s != spb; ++s, delta += 1.f / sample_rate, ++iter ) {
        float pos = s / float( spb );
        float level = begin_level * ( 1.f - pos ) + end_level * pos;
        const auto freq = float( sample_rate ) * n / resolution;
        *iter += level * std::sin( freq * ratio * 2.0 * M_PI * ( time + delta ) );
      }
    }
    {
      float begin_level = 0;
      float end_level = 0;
      for( unsigned int x = 0; x != width; ++x ) {
        if( x % 24 ) {
          begin_level += image[ x + y * width ];
          end_level += image[ x + ( y + 1 ) * width ];
        }
      }
      /*auto iter = begin;
      auto delta = 0.f;
      for( unsigned int s = 0; s != spb; ++s, delta += 1.f / sample_rate, ++iter ) {
        float pos = s / float( spb );
        float level = begin_level * ( 1.f - pos ) + end_level * pos;
        float freq = original_freq;
        *iter += level * level * std::sin( 2.0 * M_PI * ( freq * ( time + delta ) + std::sin( 2.0 * M_PI * freq * ( time + delta ) ) ) );
      }*/
    }
  }
  ifm::store_monoral( output_filename, output_wave, sample_rate );
}

