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
#include <enoki/cuda.h>
#include <enoki/dynamic.h>
#include <enoki/autodiff.h>
#include <enoki/special.h>
#include <enoki/array_router.h>
#include "ifm/fft.h"
#include "ifm/load_monoral.h"
#include "ifm/spectrum_image.h"

using i_t = enoki::DiffArray< enoki::CUDAArray< int > >;
using v_t = enoki::DiffArray< enoki::CUDAArray< float > >;


float j1( float b, int n ) {
  float prev = 0;
  float width = 0.001;
  float sum = 0;
  for( float w = width; w < float( M_PI/2 ); w += width ) {
    float cur = std::sin( b * std::sin( w ) ) * std::sin( n*w );
    sum += ( cur + prev ) * width / 2;
    prev = cur;
  }
  return 2.f / M_PI * sum;
}

float j0( float b, int n ) {
  float prev = 1;
  float width = 0.001;
  float sum = 0;
  for( float w = width; w < float( M_PI/2 ); w += width ) {
    float cur = std::cos( b * std::sin( w ) ) * std::cos( n*w );
    sum += ( cur + prev ) * width / 2;
    prev = cur;
  }
  return 2.f / M_PI * sum;
}

float j( float b, int n ) {
  if( n % 2 ) return j1( b, n );
  else return j0( b, n );
}




int main( int argc, char* argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("input,i", boost::program_options::value<std::string>(),  "入力ファイル")
    ("note,n", boost::program_options::value<int>()->default_value(60),  "音階")
    ("harmonic,H", boost::program_options::value<int>()->default_value(0),  "倍音のみ")
    ("match,m", boost::program_options::value<bool>()->default_value(false),  "減衰曲線をマッチさせる")
    ("resolution,r", boost::program_options::value<int>()->default_value(13),  "分解能");
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
  float resolution = 1 << params["resolution"].as<int>();
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
      float sum = std::accumulate( image.data() + y * width, image.data() + ( y + 1 ) * width, 0.f );
      envelope.push_back( sum );
    }
  }
  if( params["match"].as<bool>() ) {
    using v_t = enoki::DiffArray< enoki::CUDAArray< float > >;
    const auto max = std::max_element( envelope.begin(), envelope.end() );
    envelope.erase( envelope.begin(), max );
    const auto min = std::find_if( envelope.begin(), envelope.end(), []( float v ) { return v < 0.01f; } );
    envelope.erase( min, envelope.end() );
    if( params["harmonic"].as<int>() < 0 )
      envelope.resize( std::min( size_t( 50 ), envelope.size() ) );
    auto expected = v_t::copy( envelope.data(), envelope.size() );
    float a = 1;
    float b = 10;
    float c = 0.3;
    v_t beta1_ = 0.9;
    v_t beta2_ = 0.999;
    v_t alpha = 0.001;
    v_t beta1 = 0.9;
    v_t beta2 = 0.999;
    v_t a_( a );
    v_t b_( b );
    v_t c_( c );
    v_t m_a( 0 );
    v_t m_b( 0 );
    v_t m_c( 0 );
    v_t v_a( 0 );
    v_t v_b( 0 );
    v_t v_c( 0 );
    auto t = enoki::arange< v_t >( envelope.size() ) * v_t( 0.01 );
    for( int i = 1; i != 100000; ++i ) {
      a_ = a;
      b_ = b;
      c_ = c;
      enoki::set_requires_gradient( a_ );
      enoki::set_requires_gradient( b_ );
      enoki::set_requires_gradient( c_ );
      auto generated = a_ * enoki::exp( -t * b_ ) + c_;
      auto diff = ( expected - generated );
      //auto loss = enoki::hsum( diff * diff );
      //std::cout << loss << " ";
      enoki::backward( diff * diff );
     
      v_t g_a = gradient( a_ );
      v_t g_b = gradient( b_ );
      v_t g_c = gradient( c_ );
      auto ms = beta1;//enoki::pow( beta1, i_ );
      auto vs = beta2;//enoki::pow( beta2, i_ );
      beta1 *= beta1_; 
      beta2 *= beta2_; 
      m_a = beta1 * m_a + ( 1 - beta1 ) * g_a;
      v_a = beta2 * v_a + ( 1 - beta2 ) * g_a * g_a;
      auto mhat_a = m_a / ( 1 - ms );
      auto vhat_a = v_a / ( 1 - vs );
      a_ -= ( alpha * mhat_a / ( enoki::sqrt( vhat_a ) + 0.000001f ) );
      
      m_b = beta1 * m_b + ( 1 - beta1 ) * g_b;
      v_b = beta2 * v_b + ( 1 - beta2 ) * g_b * g_b;
      auto mhat_b = m_b / ( 1 - ms );
      auto vhat_b = v_b / ( 1 - vs );
      b_ -= ( alpha * mhat_b / ( enoki::sqrt( vhat_b ) + 0.000001f ) );
      
      m_c = beta1 * m_c + ( 1 - beta1 ) * g_c;
      v_c = beta2 * v_c + ( 1 - beta2 ) * g_c * g_c;
      auto mhat_c = m_c / ( 1 - ms );
      auto vhat_c = v_c / ( 1 - vs );
      c_ -= ( alpha * mhat_c / ( enoki::sqrt( vhat_c ) + 0.000001f ) );
      a = a_[ 0 ];
      b = b_[ 0 ];
      c = c_[ 0 ];
      if( !( i % 100 ) )
        std::cout << enoki::hsum( diff * diff ) << " " << a_ << " " << b_ << " " << c_ << std::endl;
      enoki::cuda_eval( true );
      enoki::cuda_sync();

    }
  }
  else {
    for( unsigned int y = 0; y != envelope.size(); ++y ) {
      std::cout << 0.01f * y << " " << envelope[ y ] << std::endl;
    }
  }
}

