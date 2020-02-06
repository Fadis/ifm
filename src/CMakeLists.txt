add_library( ifm SHARED
  bessel.cpp
  2op.cpp
  fft.cpp
  load_monoral.cpp
  store_monoral.cpp
  spectrum_image.cpp
  exp_match.cpp
  fm.cpp
)
target_link_libraries( ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( wav2image2 wav2image2.cpp )
target_link_libraries( wav2image2
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( wav2wav wav2wav.cpp )
target_link_libraries( wav2wav
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( fm2wav fm2wav.cpp )
target_link_libraries( fm2wav
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( wav2env wav2env.cpp )
target_link_libraries( wav2env
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( wav2freq wav2freq.cpp )
target_link_libraries( wav2freq
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( wav2fm wav2fm.cpp )
target_link_libraries( wav2fm
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( test_bessel test_bessel.cpp )
target_link_libraries( test_bessel
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( fm2spec fm2spec.cpp )
target_link_libraries( fm2spec
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( precomp_bessel precomp_bessel.cpp )
target_link_libraries( precomp_bessel
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( approx_env approx_env.cpp )
target_link_libraries( approx_env
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( bessel_benchmark bessel_benchmark.cpp )
target_link_libraries( bessel_benchmark
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( lossimage lossimage.cpp )
target_link_libraries( lossimage
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)
add_executable( midi_player midi_player.cpp )
target_link_libraries( midi_player
  ifm
  ${Boost_PROGRAM_OPTIONS_LIBRARIES}
  ${Boost_SYSTEM_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${OIIO_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  Threads::Threads
)

