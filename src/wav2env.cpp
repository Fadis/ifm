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
#include "ifm/fft.h"
#include "ifm/load_monoral.h"
#include "ifm/spectrum_image.h"
#include "ifm/exp_match.h"
int main( int argc, char* argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("input,i", boost::program_options::value<std::string>(),  "入力ファイル")
    ("note,n", boost::program_options::value<int>()->default_value(60),  "音階")
    ("harmonic,H", boost::program_options::value<int>()->default_value(0),  "倍音のみ")
    ("resolution,r", boost::program_options::value<int>()->default_value(13),  "分解能")
    ("damped,d", boost::program_options::value<bool>()->default_value(false),  "減衰振動")
    ("exp,e", boost::program_options::value<bool>()->default_value(false),  "指数関数近似")
    ("linear,l", boost::program_options::value<bool>()->default_value(false),  "線形近似");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") || !params.count("input") ) {
    std::cout << options << std::endl;
    return 0;
  }
  const std::string input_filename = params["input"].as<std::string>();
  const auto [audio,sample_rate] = ifm::load_monoral( input_filename, true );
  ifm::spectrum_image conv( params["note"].as<int>(), sample_rate, 1 << params["resolution"].as<int>() );
  auto [image,delay] = conv( audio );
  unsigned int width = conv.get_width();
  std::vector< float > envelope;
  if( params["harmonic"].as<int>() > 0 ) {
    for( unsigned int y = 0; y != image.size() / width; ++y ) {
      float sum = 0.f;
      for( unsigned int x = 1; x < width; x += 24 )
        sum += image[ x + y * width ];
      envelope.push_back( sum );
    }
  }
  else if( params["harmonic"].as<int>() < 0 ) {
    for( unsigned int y = 0; y != image.size() / width; ++y ) {
      float sum = 0.f;
      for( unsigned int x = 25; x < width; ++x ) {
        unsigned int step =  x % 24;
        if( step > 1 && step < 23 )
          sum += image[ x + y * width ];
      }
      envelope.push_back( sum );
    }
  }
  else {
    for( unsigned int y = 0; y != image.size() / width; ++y ) {
      float sum = std::accumulate( image.data() + y * width + 1, image.data() + ( y + 1 ) * width, 0.f );
      envelope.push_back( sum );
    }
  }
  bool damped = params[ "damped" ].as< bool >(); 
  bool exp = params[ "exp" ].as< bool >(); 
  bool linear = params[ "linear" ].as< bool >(); 
  const auto ece = ifm::get_attack_and_decay_exp( envelope.data(), envelope.size(), 0.01f, damped );
  std::cout << ece.attack_a << " " << ece.attack_b << " " << ece.decay_a << " " << ece.decay_b <<std::endl;
  const auto ecf = ifm::get_attack_and_decay( envelope.data(), envelope.size(), 0.01f, damped, damped );
  for( unsigned int y = 0; y != envelope.size(); ++y ) {
    float t = 0.01f * y;
    if( damped || y > ece.highest ) {
      std::cout << t << " " << envelope[ y ] << " ";
      if( exp ) std::cout << ifm::get_exp_envelope( ece.decay_a, ece.decay_b, t - ece.highest ) * ece.highest_level << " ";
    }
    else {
      std::cout << t << " " << envelope[ y ] << " ";
      if( exp ) std::cout << ifm::get_exp_envelope( ece.attack_a, ece.attack_b, ece.highest - t ) * ece.highest_level << " ";
    }
    if( linear ) std::cout << ifm::get_linear_envelope( ecf.first, t ) * ece.highest_level << std::endl;
    else std::cout << std::endl;
  }
}

