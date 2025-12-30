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
#include <cstddef>

namespace rebound {
  class Simulation;
  struct Vec3;
  struct ParticleStore;

  template <typename T1, typename T2>
  struct pair;

  enum class CollisionDetection {NONE, DIRECT, LINE};

  struct Collision {
    size_t p1_i = 0;
    size_t p2_i = 0;
    ParticleStore* particles = nullptr;

    Collision(size_t p1_i, size_t p2_i, ParticleStore* particles)
      : p1_i(p1_i), p2_i(p2_i), particles(particles) {}

    virtual ~Collision() = default;
  };

  class CollisionHandler {
  public:
    pair<bool, std::vector<size_t>> (*handler)(const Collision &c) = nullptr;
    double eps = 0;

    virtual bool detect_collision(ParticleStore &particles) = 0;
  };

  class CollisionLine : public CollisionHandler {
  public:
    std::vector<Vec3> prev_pos{};
    
    bool detect_collision(ParticleStore &particles);
  };

  class CollisionDirect : public CollisionHandler {
  public:
    bool detect_collision(ParticleStore &particles);
  };

  namespace collision_handlers {
    inline pair<bool, std::vector<size_t>> halt(const Collision &c) { return {true, {}}; }
    pair<bool, std::vector<size_t>> merge(const Collision &c);
  }
}