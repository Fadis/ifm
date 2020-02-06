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
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include "ifm/bessel.h"
int main( int argc, char* argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("bessel,b", boost::program_options::value<std::string>()->default_value( "bessel.mp" ), "ベッセル関数近似係数")
    ("abs,a", boost::program_options::value<bool>()->default_value( false ), "絶対値")
    ("rank,r", boost::program_options::value<int>()->default_value( 1 ), "n");
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
  for( float b = 0.f; b < 40.f; b += 0.1f ) {
    float l = ifm::bessel_kind1( b, params[ "rank" ].as< int >() );
    float l1 = ifm::bessel_kind1_approx_2019( b, params[ "rank" ].as< int >() - 1, pre );
    float l2 = ifm::bessel_kind1_approx_2019( b, params[ "rank" ].as< int >() - 2, pre );
    float r = ifm::bessel_kind1_approx_2019( b, params[ "rank" ].as< int >(), l1, l2, pre );
    float l3 = ifm::bessel_kind1_approx_2019_original1( b, params[ "rank" ].as< int >() );
    float l4 = ifm::bessel_kind1_approx_2019_original2( b, params[ "rank" ].as< int >() );
    if( !params[ "abs" ].as< bool >() )
      std::cout << b << " " << l << " " << r << " " << l - r << " " << l3 << " " << l4 << std::endl;
    else
      std::cout << b << " " << std::abs( l ) << " " << std::abs( r ) << " " << std::abs( l ) - std::abs( r ) << " " << std::abs( l3 ) << " " << std::abs( l4 ) << std::endl;
  }
}

