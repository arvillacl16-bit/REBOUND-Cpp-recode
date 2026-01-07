# rebound/meson.build

rebound_sources = files(
  'save.cpp',
  'simulation.cpp',
  'integrators/integrator_whfast.cpp',
  'integrators/integrator_leapfrog.cpp',
  'acceleration.cpp',
  'collision.cpp',
  'transformations.cpp',
  'particle.cpp'
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