project(
  'rebound-recode',
  'cpp',
  version : '0.1',
  default_options : [
    'cpp_std=c++17',
    'warning_level=0'
  ]
)

# Create OpenMP dependency if requested
omp_dep = dependency('openmp', required : get_option('openmp'))

reb_ref_lib = meson.project_source_root() / 'libactrebound.a'  # Reference to the original library name

subdir('rebound')

executable(
  'test_',
  'test.cpp',
  link_with : rebound_lib,
  # dependencies : omp_dep,   # Let Meson handle OpenMP flags automatically
  # link_args : [reb_ref_lib],  # Link against the original library file
)
