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

#include <cstdint>
#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <charconv>
#include <random>
#include <boost/container/flat_map.hpp>
#include <boost/program_options.hpp>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "ifm/setter.h"
#include "ifm/spectrum_image.h"
#include "ifm/load_monoral.h"
#include "ifm/bessel.h"
#include "ifm/fm.h"
#include "ifm/adam.h"
#include "ifm/2op.h"
#include "ifm/exp_match.h"
int main( int argc, char *argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("bessel,b", boost::program_options::value<std::string>()->default_value( "bessel.mp" ), "ベッセル関数近似係数")
    ("input,i", boost::program_options::value<std::string>(), "入力ファイル")
    ("resolution,r", boost::program_options::value<int>()->default_value(13),  "分解能")
    ("damped,d", boost::program_options::value<bool>()->default_value(false),  "減衰振動")
    ("note,n", boost::program_options::value<int>()->default_value(60), "音階")
    ("verbose,v", boost::program_options::value<bool>()->default_value(false), "詳細を表示");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") || !params.count("input") ) {
    std::cout << options << std::endl;
    return 0;
  }

  const auto [audio,sample_rate] = ifm::load_monoral( params["input"].as<std::string>(), true );
  ifm::spectrum_image conv( params["note"].as<int>(), sample_rate, 1 << params["resolution"].as<int>() );
  auto [image,delay] = conv( audio );
  unsigned int width = conv.get_width();
  unsigned int height = image.size() / width;
  unsigned int harms = width / 24;
  std::vector< float > harm( harms * height );
  std::vector< float > ec( height );
  unsigned int highest = 0;
  float highest_sum = 0.f;
  for( unsigned int y = 0; y != height; ++y ) {
    float sum = 0;
    for( unsigned int x = 24; x < width; x += 24 )
      sum += image[ y * width + x ];
    ec[ y ] = sum;
    if( sum > highest_sum ) {
      highest_sum = sum;
      highest = y;
    }
    if( ec[ y ] == 0 ) {
      ec.resize( y );
      harm.resize( harms * y );
      break;
    }
  }
  for( unsigned int y = 0; y != ec.size(); ++y ) {
    for( unsigned int x = 1; x != harms; ++x )
      harm[ y * harms + ( x - 1 ) ] = image[ y * width + x * 24 ]/ec[ y ];
  }
  std::vector< std::pair< float, float > > pre;
  {
    std::ifstream bessel_file( params["bessel"].as<std::string>(), std::ofstream::binary );
    nlohmann::json bessel = nlohmann::json::from_msgpack( bessel_file );
    for( const auto &v: bessel )
      pre.emplace_back( v.at( 0 ), v.at( 1 ) );
  }
  std::vector< float > loss( ec.size() );
  std::vector< float > em( ec.size() );
  const auto [l,freq,b] = ifm::find_b_2op( harm.data() + highest * harms, pre );
  em[ highest ] = b;
  loss[ highest ] = l;
  std::cout << "modulator freq: " << freq << std::endl;
  std::cout << "modulator scale: " << b << std::endl;
  std::cout << "loss: " << l << std::endl;
  const auto n = ifm::generate_n( freq );
  float current_b = b;
  for( unsigned int y = highest + 1; y < em.size(); ++y ) {
    auto [l,b_] = ifm::find_b_2op( harm.data() + y * harms, pre, freq, current_b, n );
    em[ y ] = b_;
    current_b = std::min( b_, current_b );
    loss[ y ] = l;
  }
  current_b = b;
  for( unsigned int y = highest; y > 0; --y ) {
    auto [l,b_] = ifm::find_b_2op( harm.data() + ( y - 1 ) * harms, pre, freq, current_b, n );
    em[ y - 1 ] = b_;
    current_b = std::min( b_, current_b );
    loss[ y ] = l;
  }
  if(  params[ "verbose" ].as< bool >() ) {
    for( unsigned int y = 0; y != ec.size(); ++y )
      std::cout << y * 0.01f << " " << ec[ y ]/ec[ highest ] << " " << em[ y ] << " " << loss[ y ] << std::endl;
  }
  bool damped = params[ "damped" ].as< bool >(); 
  const auto ece = ifm::get_attack_and_decay_exp( ec.data(), ec.size(), 0.01f, damped );
  const auto eme = ifm::get_attack_and_decay_exp( em.data(), ec.size(), 0.01f, damped );
  std::cout << ece.attack_a << " " << ece.attack_b << " " << ece.decay_a << " " << ece.decay_b <<std::endl;
  std::cout << eme.attack_a << " " << eme.attack_b << " " << eme.decay_a << " " << eme.decay_b <<std::endl;
  const auto ecf = ifm::get_attack_and_decay( ec.data(), ec.size(), 0.01f, damped, damped );
  const auto emf = ifm::get_attack_and_decay( em.data(), em.size(), 0.01f, damped, damped );
  std::cout << ifm::store_envelope_param_keyframe( ecf.first ).dump() << " " << ecf.second << std::endl;
  std::cout << ifm::store_envelope_param_keyframe( emf.first ).dump() << " " << emf.second << std::endl;
}

