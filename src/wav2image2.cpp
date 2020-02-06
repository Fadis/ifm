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
#include <boost/program_options.hpp>
#include <sndfile.h>
#include <OpenImageIO/imageio.h>
#include "ifm/fft.h"
#include "ifm/load_monoral.h"
#include "ifm/spectrum_image.h"

int main( int argc, char* argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("input,i", boost::program_options::value<std::string>(),  "入力ファイル")
    ("output,o", boost::program_options::value<std::string>(),  "出力ファイル")
    ("note,n", boost::program_options::value<int>()->default_value(60),  "音階")
    ("harmonic,H", boost::program_options::value<bool>()->default_value(false),  "倍音のみ")
    ("resolution,r", boost::program_options::value<int>()->default_value(13),  "分解能")
    ("denoise,d", boost::program_options::value<bool>()->default_value(true),  "ノイズ除去");
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
  ifm::spectrum_image conv( params["note"].as<int>(), sample_rate, 1 << params["resolution"].as<int>() );
  auto [image,delay] = conv( audio );
  std::cout << "delay: " << delay << std::endl;
  unsigned int width = conv.get_width();
  unsigned int active_width;
  std::vector< uint8_t > image_i;
  bool denoise = params["denoise"].as<bool>();
  if ( params["harmonic"].as<bool>() ) {
    active_width = width / 24;
    image_i.resize( image.size() / 24 );
    for( unsigned int y = 0; y != image.size() / width; ++y ) {
      for( unsigned int x = 0; x != active_width; ++x ) {
        auto l = denoise ?
          uint8_t( std::max( 10 * std::log10( image[ ( x + 1 ) * 24 + y * width ] ), -60.f ) * 255 / 60 + 255 ) :
          uint8_t( std::max( 10 * std::log10( image[ ( x + 1 ) * 24 + y * width ] ), -90.f ) * 255 / 90 + 255 );
        image_i[ x + y * active_width ] = l;
      }
    }
  }
  else {
    active_width = width; 
    image_i.resize( image.size() * 3 );
    for( unsigned int y = 0; y != image.size() / width; ++y ) {
      for( unsigned int x = 0; x != active_width; ++x ) {
        auto l = denoise ?
          uint8_t( std::max( 10 * std::log10( image[ x + y * width ] ), -60.f ) * 255 / 60 + 255 ) :
          uint8_t( std::max( 10 * std::log10( image[ x + y * width ] ), -90.f ) * 255 / 90 + 255 );
        image_i[ ( x + y * active_width ) * 3 ] = l;
        image_i[ ( x + y * active_width ) * 3 + 1 ] = ( x % 24 ) ? l : uint8_t( 0 );
        image_i[ ( x + y * active_width ) * 3 + 2 ] = ( x % 24 ) ? l : uint8_t( 0 );
      }
    }
  }
  using namespace OIIO_NAMESPACE;
  ImageOutput *out = ImageOutput::create( output_filename );
  if ( !out ) {
    std::cerr << "Unable to open output file" << std::endl;
    return -1;
  }
  const unsigned int channels = params["harmonic"].as<bool>() ? 1 : 3;
  ImageSpec spec ( active_width, image_i.size() / active_width / channels, channels, TypeDesc::UINT8);
  out->open( output_filename, spec );
  out->write_image( TypeDesc::UINT8, image_i.data() );
  out->close();
}

