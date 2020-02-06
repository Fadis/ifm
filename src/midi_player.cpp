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

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <array>
#include <iostream>
#include <boost/program_options.hpp>
#include "ifm/midi_player.h"
#include "ifm/midi_sequencer2.h"
#include <sndfile.h>

class wavesink {
public:
  wavesink( const char *filename ) {
    config.frames = 0;
    config.samplerate = ifm::synth_sample_rate;
    config.channels = 1;
    config.format = SF_FORMAT_WAV|SF_FORMAT_PCM_16;
    config.sections = 0;
    config.seekable = 1;
    file = sf_open( filename, SFM_WRITE, &config );
  }
  ~wavesink() {
    sf_write_sync( file );
    sf_close( file );
  }
  void play_if_not_playing() const {
  }
  bool buffer_is_ready() const {
    return true;
  }
  void operator()( const float data ) {
    int16_t ibuf[ 1 ];
    std::transform( &data, &data + 1, ibuf, []( const float &value ) { return int16_t( value * 32767 ); } );
    sf_write_short( file, ibuf, 1 );
  }
  template< size_t i >
  void operator()( const std::array< int16_t, i > &data ) {
    sf_write_short( file, data.data(), i );
  }
  template< size_t i >
  void operator()( const std::array< float, i > &data ) {
    std::array< int16_t, i > idata;
    std::transform( data.begin(), data.end(), idata.begin(), []( const float &value ) { return int16_t( value * 32767 ); } );
    sf_write_short( file, idata.data(), i );
  }
private:
  SF_INFO config;
  SNDFILE* file;
};
int main( int argc, char* argv[] ) {
  boost::program_options::options_description options("オプション");
  options.add_options()
    ("help,h",    "ヘルプを表示")
    ("config,c", boost::program_options::value<std::string>(),  "設定ファイル")
    ("input,i", boost::program_options::value<std::string>(),  "入力ファイル")
    ("output,o", boost::program_options::value<std::string>(),  "出力ファイル");
  boost::program_options::variables_map params;
  boost::program_options::store( boost::program_options::parse_command_line( argc, argv, options ), params );
  boost::program_options::notify( params );
  if( params.count("help") || !params.count("config") || !params.count("input") || !params.count("output") ) { 
    std::cout << options << std::endl;
    return 0;
  }
  const std::string input_filename = params["input"].as<std::string>();
  const std::string output_filename = params["output"].as<std::string>();
  wavesink sink( output_filename.c_str() );
  nlohmann::json config;
  {
    std::ifstream config_file( params[ "config" ].as< std::string >() );
    config_file >> config;
  }
  auto fm_params = ifm::load_fm_params< 4 >( config );
  ifm::midi_sequencer< const uint8_t*, 4 > seq( fm_params );
  const int fd = open( input_filename.c_str(), O_RDONLY );
  if( fd < 0 ) {
    return -1;
  }
  struct stat buf;
  if( fstat( fd, &buf ) < 0 ) {
    return -1;
  }
  void * const mapped = mmap( NULL, buf.st_size, PROT_READ, MAP_PRIVATE, fd, 0 );
  if( mapped == nullptr ) {
    return -1;
  }
  const auto midi_begin = reinterpret_cast< uint8_t* >( mapped );
  const auto midi_end = std::next( midi_begin, buf.st_size );
  std::cout << std::string( midi_begin, midi_begin + 4 ) << std::endl;
  if( !seq.load( midi_begin, midi_end ) ) {
    return -1;
  }
  std::array< float, ifm::synth_block_size > buffer;
  while( !seq.is_end() ) {
    seq( buffer.data() );
    sink( buffer );
  }
}

