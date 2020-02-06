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
#include "ifm/exp_match.h"
int main( int argc, char *argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("decay,d", boost::program_options::value<bool>()->default_value( true ), "decay")
    ("alpha,a", boost::program_options::value<float>()->default_value( 1 ), "a")
    ("beta,b", boost::program_options::value<float>()->default_value( 1 ), "b")
    ("length,l", boost::program_options::value<float>()->default_value( 1 ), "l");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") ) {
    std::cout << options << std::endl;
    return 0;
  }
  float x0, x1, x2, y0, y1, y2;
  float a = params["alpha"].as<float>();
  float b = params["beta"].as<float>();
  if( params["decay"].as<bool>() ) {
    const auto [d1,d2] = ifm::approxymate_decay( a, b, params["length"].as<float>() );
    x0 = 0;
    x1 = d1;
    x2 = d2;
    y0 = ifm::get_exp_envelope( a, b, x0 );
    y1 = ifm::get_exp_envelope( a, b, x1 );
    y2 = ifm::get_exp_envelope( a, b, x2 );
    std::cout << d1 << " " << d2 << std::endl;
  }
  else {
    x0 = 0;
    x1 = ifm::approxymate_attack( a, b, params["length"].as<float>() );
    x2 = params["length"].as<float>();
    y0 = ifm::get_exp_envelope( a, b, x0 );
    y1 = ifm::get_exp_envelope( a, b, x1 );
    y2 = ifm::get_exp_envelope( a, b, x2 );
    std::cout << a << std::endl;
  }
  const auto [tangent0,shift0] = ifm::get_linear_interpolation( x0, y0, x1, y1 );
  const auto [tangent1,shift1] = ifm::get_linear_interpolation( x1, y1, x2, y2 );
  float diff = 0;
  for( size_t i = 0; i != 1000; ++i ) {
    float x = i * params["length"].as<float>() / 1000.f;
    float approx_y;
    if( x < x1 ) approx_y = x * tangent0 + shift0;
    else if( x < x2 ) approx_y = x * tangent1 + shift1;
    else approx_y = y2;
    float y = ifm::get_exp_envelope( a, b, x );
    diff += approx_y - y;
  }
  std::cout << "max(" << tangent0 << "x+" << shift0 << "," << tangent1 << "x+" << shift1 << "," << y2 << ")" << std::endl;
  std::cout << "diff: " << diff/1000.f << std::endl;
}

