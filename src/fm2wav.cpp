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
#include <boost/container/flat_map.hpp>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include "ifm/setter.h"
#include "ifm/store_monoral.h"
#include "ifm/fm.h"

int main( int argc, char *argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("config,c", boost::program_options::value<std::string>(),  "設定ファイル")
    ("output,o", boost::program_options::value<std::string>(),  "出力ファイル")
    ("note,n", boost::program_options::value<int>()->default_value(60),  "音階");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") || !params.count("config") || !params.count("output") ) {
    std::cout << options << std::endl;
    return 0;
  }
  nlohmann::json config;
  {
    std::ifstream config_file( params[ "config" ].as< std::string >() );
    config_file >> config;
  }
  auto fm_params = ifm::load_fm_params< 4 >( config );
  ifm::fm_t< double, 4 > fm( fm_params, params[ "note" ].as<int>(), 127 );
  std::vector< float > audio( ifm::synth_sample_rate * 10 );
  for( unsigned int i = 0; i != ifm::synth_sample_rate * 10 / ifm::synth_block_size; ++i ) {
    fm( audio.data() + ifm::synth_block_size * i );
    if( fm.is_end() ) {
      audio.resize( ( i + 1 ) * ifm::synth_block_size );
      break;
    }
  }
  ifm::store_monoral( params[ "output" ].as< std::string >(), audio, ifm::synth_sample_rate );
}

