project(
  'exactpass',
  ['c','cpp'],
  version : '0.3.0',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++17',
    'c_std=gnu11'],
  license : 'BSD-3',
  meson_version : '>= 0.46.0')

# CXXFLAGS
exactpass_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-unused-parameter',
  '-Wno-unused-command-line-argument',
  '-Wno-missing-field-initializers',
  '-Wno-unused-const-variable',
  '-Wno-c99-extensions']
    if cpp.has_argument(cflag)
      exactpass_flags += cflag
    endif
endforeach

# dependencies
## threads
exactpass_thread_dep = dependency('threads', required : true)
## boost
exactpass_boost_dep = dependency('boost', required : true)
## zlib
exactpass_zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])
## htslib
exactpass_htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])
## pbbam
exactpass_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])
## pbcopper
exactpass_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])

exactpass_lib_deps = [
  exactpass_thread_dep,
  exactpass_boost_dep,
  exactpass_zlib_dep,
  exactpass_htslib_dep,
  exactpass_pbbam_dep,
  exactpass_pbcopper_dep,
]

executable(
  'exactpass', files([
    'exactpass.cpp']),
  install : true,
  dependencies : exactpass_lib_deps,
  cpp_args : exactpass_flags)

executable(
  'subsamplezmws', files([
    'subsamplezmws.cpp']),
  install : true,
  dependencies : exactpass_lib_deps,
  cpp_args : exactpass_flags)

executable(
  'filterprimer', files([
    'filterprimer.cpp']),
  install : true,
  dependencies : exactpass_lib_deps,
  cpp_args : exactpass_flags)
