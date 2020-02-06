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

#ifndef IFM_2OP_H
#define IFM_2OP_H
#include <array>
#include <vector>
#include <utility>
#include <tuple>
namespace ifm {
  constexpr int max_harmony = 50;
  std::tuple< float, float > lossimage(
    const float *expected,
    const std::vector< std::pair< float, float > > &pre,
    unsigned int freq,
    float b,
    const std::array< int, max_harmony > &n
  );
  std::tuple< float, unsigned int, float > find_b_2op(
    const float *expected,
    const std::vector< std::pair< float, float > > &pre
  );
  std::array< int, max_harmony > generate_n(
    unsigned int freq
  );
  std::tuple< float, float > find_b_2op(
    const float *expected,
    const std::vector< std::pair< float, float > > &pre,
    unsigned int freq,
    float initial_b,
    const std::array< int, max_harmony > &n
  );
}

#endif

