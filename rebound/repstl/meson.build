# rebound/repstl/meson.build

repstl_sources = files(
  'src/string.cpp'
)

repstl_includes = include_directories('.')
repstl_lib = static_library(
  'repstl',
  repstl_sources,
  include_directories : repstl_includes,
  install : false
)

repstl_dep = declare_dependency(
  link_with : repstl_lib,
  include_directories : repstl_includes
)

repstl_dep