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

#include <vector>
#include <memory>
#include <cmath>
#include <array>
#include <iostream>
#include <complex>
#include <exception>
#include <stdexcept>
#include <fftw3.h>
#include "ifm/fft.h"

namespace ifm {
  fft_context_t::fft_context_t( size_t res, float sample_rate_ ) :
    input( reinterpret_cast< float* >( fftwf_malloc( sizeof( float ) * res ) ) ),
    output( reinterpret_cast< fftwf_complex* >( fftwf_malloc( sizeof( fftwf_complex ) * res / 2 ) ) ),
    resolution( res ),
    sample_rate( sample_rate_ )
  {
    if( !input || !output ) throw std::bad_alloc();
    plan.reset( fftwf_plan_dft_r2c_1d( resolution, input.get(), output.get(), FFTW_ESTIMATE ) );
    if( !plan ) throw std::bad_alloc();
    window.resize( resolution );
    for( unsigned int i = 0; i != window.size(); ++i )
      window[ i ] = ( 1. - std::cos( 2. * M_PI * i / resolution ) ) / 2.;
  }
  std::vector< float > fft_context_t::forward( const std::vector< float > &input_ ) {
    const auto size = input_.size();
    if( size != resolution ) throw incompatible_range();
    for( unsigned int i = 0; i != window.size(); ++i )
      input.get()[ i ] = input_[ i ] * window[ i ];
    fftwf_execute( plan.get() );
    std::vector< float > result;
    result.resize( resolution / 2 );
    for( unsigned int i = 0; i != result.size(); ++i ) {
      result[ i ] = std::sqrt( output.get()[ i ][ 0 ] * output.get()[ i ][ 0 ] + output.get()[ i ][ 1 ] * output.get()[ i ][ 1 ] )/resolution;
    }
    const auto i = ( float( resolution ) ) * 20000.f / sample_rate;
    result.resize( i );
    return result;
  }
}

