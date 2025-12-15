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

namespace rebound {
  void Integrator::step_leapfrog_p1(_ParticleStore& particles, double dt) const {
    const size_t N = particles.size();

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) particles.positions[i] += dt / 2 * particles.velocities[i];
  }

  void Integrator::step_leapfrog_p2(_ParticleStore& particles, double dt) const {
    const size_t N = particles.size();

    double dt_half = dt / 2;

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) particles.velocities[i] += dt * particles.accelerations[i];

    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) particles.positions[i] += dt_half * particles.velocities[i];
  }
}