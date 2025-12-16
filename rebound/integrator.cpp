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
  void Integrator::step(ParticleStore& particles, double dt) {
    if (gravity_method == GravityMethod::MERCURIUS && method != IntegratorMethod::MERCURIUS) {
      std::cerr << "You are using mercurius gravity with a non-mercurius integrator. This will lead to UB. Automatically switched to BASIC gravity.";
      gravity_method = GravityMethod::BASIC;
    }

    if (method == IntegratorMethod::LEAPFROG) step_leapfrog_p1(particles, dt);
    else if (method == IntegratorMethod::WHFAST) step_whfast_p1(particles, dt);
    else if (method == IntegratorMethod::IAS15) step_ias15_p1(particles, dt);
    else if (method == IntegratorMethod::MERCURIUS) step_mercurius_p1(particles, dt);

    if (gravity_method == GravityMethod::NONE) _accel::calc_accel_none(particles);
    else if (gravity_method == GravityMethod::BASIC) _accel::calc_accel_basic(particles, softening2);
    else if (gravity_method == GravityMethod::COMPENSATED) _accel::calc_accel_compensated(particles, softening2);
    else if (gravity_method == GravityMethod::JACOBI) _accel::calc_accel_jacobi(particles, softening2);
    else if (gravity_method == GravityMethod::MERCURIUS) _accel::calc_accel_mercurius(particles, softening2, mercurius_settings);

    if (method == IntegratorMethod::LEAPFROG) step_leapfrog_p2(particles, dt);
    else if (method == IntegratorMethod::WHFAST) step_whfast_p2(particles, dt);
    else if (method == IntegratorMethod::IAS15) step_ias15_p2(particles, dt);
    else if (method == IntegratorMethod::MERCURIUS) step_mercurius_p2(particles, dt);
  }
}