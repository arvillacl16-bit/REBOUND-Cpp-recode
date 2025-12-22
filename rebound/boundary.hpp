#pragma once
#include "vec.hpp"

namespace rebound {
  struct ParticleStore;
  class BoundaryHandler {
  public:
    struct {
      Vec3 center{};
      double a = 0;
      double b = 0;
      double c = 0;
    } boundary;

    virtual void handle_boundary(ParticleStore &particles) = 0;
  };

  class BoundaryOpen : public BoundaryHandler {
  public:
    void handle_boundary(ParticleStore &particles);
  };

  class BoundaryPeriodic : public BoundaryHandler {
  public:
    void handle_boundary(ParticleStore &particles);
  };

  class BoundaryShearPeriodic : public BoundaryHandler {
  public:
    void handle_boundary(ParticleStore &particles);
  };
}