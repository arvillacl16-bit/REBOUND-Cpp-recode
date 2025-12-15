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

#include "particle.hpp"

namespace rebound {
  namespace _transform {
    void inertial_to_jacobi_posvel(const _ParticleStore& from, _ParticleStore& to);
    void inertial_to_jacobi_posvelacc(const _ParticleStore& from, _ParticleStore& to);
    void inertial_to_jacobi_acc(const _ParticleStore& from, _ParticleStore& to);
    void jacobi_to_inertial_posvel(_ParticleStore& to, const _ParticleStore& from);
    void jacobi_to_inertial_pos(_ParticleStore& to, const _ParticleStore& from);
    void jacobi_to_inertial_acc(_ParticleStore& to, const _ParticleStore& from);

    void inertial_to_democraticheliocentric_posvel(const _ParticleStore& from, _ParticleStore& to);
    void democraticheliocentric_to_inertial_pos(_ParticleStore& to, const _ParticleStore& from);
    void democraticheliocentric_to_inertial_posvel(_ParticleStore& to, const _ParticleStore& from);

    void inertial_to_whds_posvel(const _ParticleStore& from, _ParticleStore& to);
    void whds_to_inertial_pos(_ParticleStore& to, const _ParticleStore& from);
    void whds_to_inertial_posvel(_ParticleStore& to, const _ParticleStore& from);

    void inertial_to_barycentric_posvel(const _ParticleStore& from, _ParticleStore& to);
    void barycentric_to_inertial_pos(_ParticleStore& to, const _ParticleStore& from);
    void barycentric_to_inertial_posvel(_ParticleStore& to, const _ParticleStore& from);
    void barycentric_to_inertial_acc(_ParticleStore& to, const _ParticleStore& from);
  }
}