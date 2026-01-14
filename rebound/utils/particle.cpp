/*
 * This file is part of the C++ translation of REBOUND.
 *
 * Original REBOUND (C) code by Hanno Rein and others.
 * This translation is licensed under the GNU General Public License v3 or
 * later.
 *
 * REBOUND is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * REBOUND is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "particle.hpp"
#include <cassert>
#include <iostream>

namespace rebound {
  ParticleStore::ParticleStore(size_t capacity) {
    positions.reserve(capacity);
    velocities.reserve(capacity);
    accelerations.reserve(capacity);
    gravity_cs.reserve(capacity);
    mus.reserve(capacity);
    radii.reserve(capacity);
    ids.reserve(capacity);
    test_mass.reserve(capacity);
  }

  size_t ParticleStore::size() const { return positions.size(); }

  size_t ParticleStore::n_real() const {
    size_t n = 0;
    for (size_t i = 0; i < size(); i++) {
      if (!test_mass[i])
        ++n;
    }
    return n;
  }

  void ParticleStore::print_if_nan_or_inf() const {
    for (size_t i = 0; i < size(); ++i) {
      if (positions[i].has_nan_or_inf() || velocities[i].has_nan_or_inf()) {
        LOG("NaN or Inf detected in particle " << i);
        return;
      }
    }
  }

  Particle ParticleStore::add_particle(const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass_) {
    positions.push_back(position);
    velocities.push_back(velocity);
    accelerations.push_back(Vec3());
    gravity_cs.push_back(Vec3());
    mus.push_back(mu);
    radii.push_back(radius);
    ids.push_back(id);
    test_mass.push_back(test_mass_);
    versions.push_back(0);
    return Particle(size() - 1, this, 0);
  }

  void ParticleStore::remove_particle(size_t idx) {
    versions[idx]++;

    size_t last_idx = size() - 1;
    if (idx != last_idx) {
      positions[idx] = positions[last_idx];
      velocities[idx] = velocities[last_idx];
      accelerations[idx] = accelerations[last_idx];
      gravity_cs[idx] = gravity_cs[last_idx];
      mus[idx] = mus[last_idx];
      radii[idx] = radii[last_idx];
      ids[idx] = ids[last_idx];
      test_mass[idx] = test_mass[last_idx];
      versions[idx] = versions[last_idx];
    }
    positions.pop_back();
    velocities.pop_back();
    accelerations.pop_back();
    gravity_cs.pop_back();
    mus.pop_back();
    radii.pop_back();
    ids.pop_back();
    test_mass.pop_back();
    versions.pop_back();
  }

  Particle ParticleStore::operator[](size_t idx) { return Particle(idx, this, versions[idx]); }
  ConstParticle ParticleStore::operator[](size_t idx) const { return ConstParticle(idx, this, versions[idx]); }
} // namespace rebound