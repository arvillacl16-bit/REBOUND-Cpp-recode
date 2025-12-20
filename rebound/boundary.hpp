#pragma once
#include "vec.hpp"

namespace rebound {
  struct ParticleStore;
  class BoundaryHandler {
  private:
    //
  public:
    struct {
      Vec3 center{};
      double a = 0;
      double b = 0;
      double c = 0;
    } boundary;

    virtual void handle_boundary(ParticleStore &particles) = 0;
  };
}