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
#include <iostream>
#include <chrono>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include "ifm/bessel.h"
int main( int argc, char* argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("bessel,b", boost::program_options::value<std::string>()->default_value( "bessel.mp" ), "ベッセル関数近似係数");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") ) {
    std::cout << options << std::endl;
    return 0;
  }
  std::vector< std::pair< float, float > > pre;
  {
    std::ifstream bessel_file( params["bessel"].as<std::string>(), std::ofstream::binary );
    nlohmann::json bessel = nlohmann::json::from_msgpack( bessel_file );
    for( const auto &v: bessel )
      pre.emplace_back( v.at( 0 ), v.at( 1 ) );
  }
  const auto t0 = std::chrono::steady_clock::now();
  for( unsigned int b_ = 0.f; b_ != 400; ++b_ ) {
    float b = b_ * 0.1f;
    for( unsigned int n = 0; n != 60; ++n ) {
      ifm::bessel_kind1( b, n );
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  for( unsigned int b_ = 0.f; b_ != 400; ++b_ ) {
    float b = b_ * 0.1f;
    float l1 = 0;
    float l2 = 0;
    for( unsigned int n = 0; n != 60; ++n ) {
      float y = ifm::bessel_kind1_approx_2019( b, n, l1, l2, pre );
      l2 = l1;
      l1 = y;
    }
  }
  const auto t2 = std::chrono::steady_clock::now();
  float e0 = float( std::chrono::duration_cast< std::chrono::microseconds >( t1 - t0 ).count() );
  float e1 = float( std::chrono::duration_cast< std::chrono::microseconds >( t2 - t1 ).count() );
  std::cout << "台形公式: " << (400*60)/e0 << "Mbps(" << 400*60 << " samples in " << e0 << "microseconds)" << std::endl;
  std::cout << "近似式: " << (400*60)/e1 << "Mbps(" << 400*60 << " samples in " << e1 << "microseconds)" << std::endl;
}

