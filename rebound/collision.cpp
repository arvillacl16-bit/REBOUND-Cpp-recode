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
#include <unordered_set>
#include <algorithm>

#define COMB_RAD2(i, j) (particles.radii[i] + particles.radii[j]) * (particles.radii[i] + particles.radii[j])

namespace rebound {
  bool CollisionHandler::collision_direct(size_t i, size_t j, const ParticleStore &particles) { return particles.positions[i].distance2(particles.positions[j]) < COMB_RAD2(i, j); }

  bool CollisionHandler::collision_line(size_t i, size_t j, const ParticleStore &particles) {
    Vec3 a = prev_pos[j] - prev_pos[i];
    Vec3 d = a - (particles.positions[j] - particles.positions[i]);

    double t = a.dot(d) / d.mag2();
    Vec3 v = a - t * d;
    return v.mag2() < COMB_RAD2(i, j);
  }

  bool CollisionHandler::detect_collision(ParticleStore &particles) {
    if (handler && (detect != CollisionDetection::NONE)) {
      std::unordered_set<size_t> indices;
      size_t N = particles.size();
      bool val = false;
      if (detect == CollisionDetection::DIRECT) {
        for (size_t i = 0; i < N; ++i) {
          for (size_t j = i + 1; j < N; ++j) {
            if (collision_direct(i, j, particles)) {
              auto result = handler({i, j, &particles});
              if (result.first) val = true;
              indices.insert(result.second.begin(), result.second.end());
            }
          }
        }
      } else if (detect == CollisionDetection::LINE) {
        for (size_t i = 0; i < N; ++i) {
          for (size_t j = i + 1; j < N; ++j) {
            if (collision_line(i, j, particles)) {
              auto result = handler({i, j, &particles});
              if (result.first) val = true;
              indices.insert(result.second.begin(), result.second.end());
            }
          }
        }
      }

      for (size_t i = N; i-- > 0;) {
        if (indices.find(i) != indices.end()) particles.remove_particle(i);
      }
    }
    return false;
  }
  namespace collision_handlers {
    pair<bool, std::vector<size_t>> merge(const Collision &c) {
      size_t i = c.p1_i;
      size_t j = c.p2_i;
      ParticleStore &particles = *c.particles;
      Vec3 pos_i = particles.positions[i];
      Vec3 pos_j = particles.positions[j];
      Vec3 vel_i = particles.velocities[i];
      Vec3 vel_j = particles.velocities[j];
      double m_i = particles.mus[i];
      double m_j = particles.mus[j];
      double m_tot = m_i + m_j;
      if (m_i > m_j) {
        particles.velocities[i] = (m_i * vel_i + m_j * vel_j) / m_tot;
        particles.positions[i] = (m_i * pos_i + m_j * pos_j) / m_tot;
        return {false, {j}};
      } else {
        particles.velocities[j] = (m_i * vel_i + m_j * vel_j) / m_tot;
        particles.positions[j] = (m_i * pos_i + m_j * pos_j) / m_tot;
        return {false, {i}};
      }
    }
  }
}