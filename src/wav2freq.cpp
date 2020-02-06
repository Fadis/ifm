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
#include <iterator>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
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
    ("harmonic,H", boost::program_options::value<int>()->default_value(0),  "倍音のみ")
    ("match,m", boost::program_options::value<bool>()->default_value(false),  "減衰曲線をマッチさせる")
    ("resolution,r", boost::program_options::value<int>()->default_value(13),  "分解能");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") || !params.count("input") || !params.count( "output" ) ) {
    std::cout << options << std::endl;
    return 0;
  }
  const std::string input_filename = params["input"].as<std::string>();
  const auto [audio,sample_rate] = ifm::load_monoral( input_filename, true );
  ifm::spectrum_image conv( params["note"].as<int>(), sample_rate, 1 << params["resolution"].as<int>() );
  auto [image,delay] = conv( audio );
  unsigned int width = conv.get_width();
  bool highest_sum = 0.f;
  std::vector< float > highest( width / 24 );
  for( unsigned int y = 0; y != image.size() / width; ++y ) {
    auto sum = std::accumulate( std::next( image.begin(), y * width ), std::next( image.begin(), ( y + 1 ) * width ), float( 0 ) );
    if( sum > highest_sum ) {
      highest_sum = sum;
      for( unsigned int x = 0; x != width / 24; ++x )
        highest[ x ] = image[ x * 24 + y * width ];
    }
  }
  std::map< unsigned int, float > herm;
  float herm_sum = 0.f;
  for( unsigned int i = 1; i != width / 24; ++i ) {
    herm_sum += highest[ i ];
  }
  for( unsigned int i = 1; i != width / 24; ++i ) {
    herm.emplace( i, highest[ i ]/herm_sum );
    std::cout << highest[ i ]/herm_sum << std::endl;
  }
  nlohmann::json out( herm );
  auto mp = nlohmann::json::to_msgpack( out );
  std::ofstream out_file( params["output"].as<std::string>(),std::ofstream::binary );
  std::copy( mp.begin(), mp.end(), std::ostreambuf_iterator< char >( out_file ) );
}

