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

#include <fstream>
#include <cstdint>
#include <cmath>
#include <array>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <nlohmann/json.hpp>
#include <boost/program_options.hpp>
#include "ifm/bessel.h"

int main( int argc, char *argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("input,i", boost::program_options::value<std::string>(), "入力ファイル")
    ("freq,f", boost::program_options::value<int>()->default_value(1), "倍率")
    ("volume,v", boost::program_options::value<float>()->default_value(1), "出力")
    ("level,l", boost::program_options::value<float>()->default_value(1), "変調度");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") ) {
    std::cout << options << std::endl;
    return 0;
  }
  constexpr int max_harmony = 50;
  std::array< float, max_harmony > expected{ 0 };
  if( params.count("input") ) {
    std::ifstream in_file( params["input"].as<std::string>(),std::ofstream::binary );
    nlohmann::json input = nlohmann::json::from_msgpack( in_file );
    for( const auto &v: input ) {
      int key = int( v.at( 0 ) );
      if( key < max_harmony )
        expected[ key - 1 ] = float( v.at( 1 ) );
    }
  }
  std::array< float, max_harmony > bessel{ 0 };
  std::array< float, max_harmony > generated{ 0 };
  std::array< int, max_harmony > n;
  std::fill( n.begin(), n.end(), -1 );
  for( int i = 0; i != max_harmony; ++i )
    if( i % params[ "freq" ].as< int >() == 0 )
      n[ i ] = i / params[ "freq" ].as< int >();
  const auto pre = ifm::create_bessel_approx_2019_precomp_array( max_harmony / params[ "freq" ].as< int >() + 5 );
  float l1 = 0;
  float l2 = 0;
  for( unsigned int i = 0; i != max_harmony; ++i ) {
    bessel[ i ] = ifm::bessel_kind1_approx_2019( params[ "level" ].as< float >(), i, l1, l2, pre );
    l2 = l1;
    l1 = bessel[ i ];
  }
  float bscale = 0;
  for( unsigned int i = 0; i != max_harmony; ++i ) {
    bscale += std::abs( bessel[ i ] );
    if( i < max_harmony - 2 ) bscale += std::abs( bessel[ i + 2 ] );
  }
  float a = params[ "volume" ].as< float >();
  for( unsigned int i = 0; i != max_harmony - 2; ++i ) {
    if( n[ i ] != -1 )
      generated[ i ] = a * bessel[ n[ i ] ]/bscale;
    else
      generated[ i ] = 0;
    if( i < max_harmony - 2 && n[ i + 2 ] != -1 )
      generated[ i ] -= a * bessel[ n[ i + 2 ] ]/bscale;
  }
  float sum = 0;
  for( unsigned int i = 0; i != max_harmony; ++i ) {
    if( params.count("input") ) {
      if( expected[ i ] != 0 )
        std::cout << expected[ i ] << "  " << std::abs( generated[ i ] ) << " " << std::abs( generated[ i ] ) - expected[ i ] << " " << ( std::abs( generated[ i ] ) - expected[ i ] )/expected[ i ] << std::endl;
      else
        std::cout << expected[ i ] << "  " << std::abs( generated[ i ] ) << " " << std::abs( generated[ i ] ) - expected[ i ] << std::endl;
    }
    else
      std::cout << std::abs( generated[ i ] ) << std::endl;
    sum += std::abs( generated[ i ] );
  }
  std::cout << "bscale: " << bscale << std::endl;
  std::cout << "total: " <<  sum << std::endl;
}

