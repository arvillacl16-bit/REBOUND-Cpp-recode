# rebound/meson.build

# Gather all source files
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

# Create a static library
rebound_lib = static_library(
  'rebound',
  rebound_sources,
  include_directories : include_directories('.'),
)

# Optional: export this library variable for parent meson.build
rebound_lib
