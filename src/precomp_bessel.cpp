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
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include "ifm/setter.h"
#include "ifm/store_monoral.h"
#include "ifm/bessel.h"

int main( int argc, char *argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("output,o", boost::program_options::value<std::string>()->default_value( "bessel.mp" ), "出力ファイル")
    ("max,m", boost::program_options::value<int>()->default_value(60), "最大階数");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") || !params.count("output") ) {
    std::cout << options << std::endl;
    return 0;
  }
  const auto pre = ifm::create_bessel_approx_2019_precomp_array( params[ "max" ].as< int >() );
  nlohmann::json out( pre );
  auto mp = nlohmann::json::to_msgpack( out );
  std::ofstream out_file( params["output"].as<std::string>(),std::ofstream::binary );
  std::copy( mp.begin(), mp.end(), std::ostreambuf_iterator< char >( out_file ) );
}

