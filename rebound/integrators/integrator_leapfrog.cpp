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

#include "integrator.hpp"
#include <iostream>

namespace rebound {
  void Leapfrog::step_p1(ParticleStore& particles, double dt) const {
    const size_t N = particles.size();

#pragma omp parallel for
    for (size_t i = 0; i < N; ++i) particles.positions[i] += dt / 2 * particles.velocities[i];
  }

  void Leapfrog::step_p2(ParticleStore& particles, double dt) const {
    const size_t N = particles.size();

    double dt_half = dt / 2;

#pragma omp parallel for
    for (size_t i = 0; i < N; ++i) particles.velocities[i] += dt * particles.accelerations[i];

#pragma omp parallel for
    for (size_t i = 0; i < N; ++i) particles.positions[i] += dt_half * particles.velocities[i];
  }

  void Leapfrog::step(ParticleStore& particles, double dt) {
    step_p1(particles, dt);
    switch (gravity_method) {
      case GravityMethod::BASIC:
        _accel::calc_accel_basic(particles, softening2);
        break;
      case GravityMethod::COMPENSATED:
        _accel::calc_accel_compensated(particles, softening2);
        break;
      case GravityMethod::JACOBI:
        _accel::calc_accel_jacobi(particles, softening2);
        break;
      default:
        break;
    }
    step_p2(particles, dt);
  }
}