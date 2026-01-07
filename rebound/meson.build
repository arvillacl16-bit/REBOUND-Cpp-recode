# rebound/meson.build

rebound_sources = files(
  'utils/save.cpp',
  'utils/simulation.cpp',
  'integration/integrator_whfast.cpp',
  'integration/integrator_leapfrog.cpp',
  'integration/acceleration.cpp',
  'non_integration/collision.cpp',
  'non_integration/boundary.cpp',
  'utils/transformations.cpp',
  'utils/particle.cpp'
)

deps = []
if omp_dep.found()
  deps += omp_dep
endif

rebound_lib = static_library(
  'rebound',
  rebound_sources,
  include_directories : include_directories('.'),
  dependencies : deps
)

rebound_dep = declare_dependency(
  link_with : rebound_lib,
  include_directories : include_directories('.'),
  dependencies : deps
)

rebound_lib