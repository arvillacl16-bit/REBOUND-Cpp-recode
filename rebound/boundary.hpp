#pragma once
#include "vec.hpp"

namespace rebound {
  struct ParticleStore;

  enum class BoundaryType { NONE, OPEN, PERIODIC };
  
  class BoundaryHandler {
  private:
    //
  public:
    BoundaryType type = BoundaryType::NONE;

    struct {
      Vec3 center{};
      double a = 0;
      double b = 0;
      double c = 0;
    } boundary;

    void handle_boundary(ParticleStore &particles);
  };
}