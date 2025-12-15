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
#include "particle.hpp"
#include "integrator.hpp"
#include "collision.hpp"
#include "transformations.hpp"

namespace rebound {
  enum class BoundaryCondition {NONE, OPEN, PERIODIC};

  class Simulation {
  private:
    std::vector<Vec3> prev_pos;
    double t;

    bool collision_direct(size_t i, size_t j);
    bool collision_line(size_t i, size_t j);
  public:
    _ParticleStore particles;
    Integrator integrator;
    double dt;
    bool (*heartbeat) (Simulation& sim) = nullptr; // If returns true, makes the integration terminate

    struct {
      double eps; // epsilon to avoid divide by zero
      CollisionDetection method; // NONE = no detection, DIRECT = normal nxn detection, LINE = check along a line
      pair<bool, std::vector<size_t>> (*handler)(const Collision& c) = nullptr;
    } collision_settings;

    struct {
      BoundaryCondition condition;
      double a = 0; // half
      double b = 0; // For a rectanngle, half-side-length for the y-axis, for ellipsoid, y-axis half
      double c = 0; // For a rectanngle, half-side-length for the z-axis, for ellipsoid, z-axis half
    } boundary_settings;

    Simulation() {}
    explicit Simulation(std::string filename);

    Particle add_particle(const Vec3& position, const Vec3& velocity, double mass, double radius, uint32_t id);
    void remove_particle(size_t idx);
    void remove_particle(uint32_t id);

    inline size_t n() const { return particles.size(); }
    inline size_t n_real() const {
      size_t n = 0;
      for (size_t i = 0; i < particles.size(); i++) {
        if (!particles.test_mass[i]) ++n;
      }
      return n;
    }
    inline size_t n_test() const { return n() - n_real(); }

    double time() const;

    void integrate(double t, bool exact_finish_time = false);
    void integrate_walltime(double t); // t in seconds
    bool step(double dt_);
    inline void steps(unsigned int N) { for (size_t i = 0; i < N; ++i) step(dt); }

    void save(const std::string& filename) const;

    Vec3 com_pos();
    Vec3 com_vel();
    void move_to_com();
  };
}