# ... project setup ...
project(
  'rebound-recode',
  'cpp',
  version : '0.1',
  default_options : [
    'cpp_std=c++17',
    'warning_level=0'
  ]
)

omp_dep = dependency('openmp', required : false)

subdir('rebound')

act_reb_lib = meson.project_source_root() / 'libactrebound.a'

executable(
  'test_',
  'test.cpp',
  dependencies : [rebound_dep], # Bundles rebound_lib + headers + omp_dep flags
  link_args : [act_reb_lib]
)