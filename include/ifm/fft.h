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

#ifndef IFM_FFT_H
#define IFM_FFT_H
#include <vector>
#include <memory>
#include <cmath>
#include <array>
#include <iostream>
#include <complex>
#include <exception>
#include <stdexcept>
#include <fftw3.h>

namespace ifm {
  struct free_fftw_mem {
    template< typename T >
    void operator()( T p ) {
      if( p ) fftwf_free( p );
    }
  };
  struct free_fftw_plan {
    template< typename T >
    void operator()( T p ) {
      if( p ) fftwf_destroy_plan( p );
    }
  };
struct incompatible_range {};
//using array_t = enoki::CUDAArray< float >;
//using darray_t = enoki::DiffArray< array_t >;
class fft_context_t {
public:
  fft_context_t( size_t res, float sample_rate = 192000.f );
  std::vector< float > forward( const std::vector< float >& );
private:
  std::vector< float > window;
  std::unique_ptr< float, free_fftw_mem > input;
  std::unique_ptr< fftwf_complex, free_fftw_mem > output;
  size_t resolution;
  float sample_rate;
  std::unique_ptr< std::remove_pointer_t< fftwf_plan >, free_fftw_plan > plan;
};
}

#endif

