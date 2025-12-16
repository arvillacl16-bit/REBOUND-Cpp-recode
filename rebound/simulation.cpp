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

#include "rebound.hpp"
#include <chrono>
#include <unordered_set>

namespace chrono = std::chrono;

namespace rebound {
  Particle Simulation::add_particle(const Vec3& position, const Vec3& velocity, double mass, double radius, uint32_t id) {
    return particles.add_particle(position, velocity, mass, radius, id);
  }

  void Simulation::remove_particle(size_t idx) { particles.remove_particle(idx); }

  void Simulation::remove_particle(uint32_t id) {
    for (size_t i = 0; i < particles.size(); ++i) {
      if (particles.ids[i] == id) {
        particles.remove_particle(i);
        return;
      }
    }
  }

  double Simulation::time() const { return t; }

  void Simulation::integrate(double t_end, bool exact_finish_time) {
    double tleft = t_end - t;
    while (tleft > 0) {
      double dt_step = exact_finish_time ? std::min(tleft, dt) : dt;
      if (step(dt_step)) break;
      tleft -= dt_step;
    }
  }

  void Simulation::integrate_walltime(double t) {
    auto start = chrono::steady_clock::now();
    chrono::duration<double> max_duration(t);

    while (chrono::steady_clock::now() - start < max_duration) if (step(dt)) break;
  }

  bool Simulation::step(double dt_) {
    bool val = false;

    coll_handler.prev_pos = particles.positions;
    integrator.step(particles, dt_);
    t += dt_;

    val = coll_handler.detect_collision(particles);

    if (heartbeat(*this)) val = true;
    return val;
  }

  Vec3 Simulation::com_pos() {
    std::vector<Vec3> &positions = particles.positions;
    std::vector<double> &mus = particles.mus;
    if (positions.size() != n() || mus.size() != n()) return {};

    double total_mu = 0;
    Vec3 result = {};

#pragma omp parallel for reduction(+:total_mu)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) total_mu += mus[i];

#pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) result += positions[i] * (mus[i] / total_mu);

    return result;
  }

  Vec3 Simulation::com_vel() {
    std::vector<Vec3> &velocities = particles.velocities;
    std::vector<double> &mus = particles.mus;
    if (velocities.size() != n() || mus.size() != n()) return {};

    double total_mu = 0;
    Vec3 result = {};

#pragma omp parallel for reduction(+:total_mu)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) total_mu += mus[i];

#pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < n(); ++i) if (!particles.test_mass.at(i)) result += velocities[i] * (mus[i] / total_mu);

    return result;
  }

  void Simulation::move_to_com() {
    Vec3 com_pos_ = com_pos(), com_vel_ = com_vel();
#pragma omp parallel for
    for (size_t i = 0; i < n(); ++i) {
      particles.positions[i] -= com_pos_;
      particles.velocities[i] -= com_vel_;
    }
  }
}