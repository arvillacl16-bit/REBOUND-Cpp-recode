# rebound/meson.build

subdir('repstl')
current_deps = []
current_deps += repstl_dep
if get_option('openmp')
  current_deps += omp_dep
endif

rebound_sources = files(
  'integration/acceleration.cpp',
  'integration/integrator_leapfrog.cpp',
  'integration/integrator_whfast.cpp',
  'non_integration/boundary.cpp',
  'non_integration/collision.cpp',
  'utils/particle.cpp',
  'utils/save.cpp',
  'utils/simulation.cpp',
  'utils/transformations.cpp'
)

rebound_includes = include_directories('.', 'integration', 'non_integration', 'utils')

rebound_lib = static_library(
  'rebound',
  rebound_sources,              # Sources go first
  include_directories : rebound_includes,
  dependencies : current_deps,
  install : false
)

rebound_dep = declare_dependency(
  link_with : rebound_lib,
  include_directories : rebound_includes, # Pass the include paths to the test executable
  dependencies : current_deps
)