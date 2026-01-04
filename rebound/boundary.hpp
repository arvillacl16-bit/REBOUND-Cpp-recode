/*
 * This file is part of the C++ translation of REBOUND.
 *
 * Original REBOUND (C) code by Hanno Rein and others.
 * This translation is licensed under the GNU General Public License v3 or later.
 *
 * REBOUND is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * REBOUND is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 */

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
    virtual ~BoundaryHandler() = default;
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