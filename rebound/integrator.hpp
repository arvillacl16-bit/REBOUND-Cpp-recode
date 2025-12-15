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
#include "integrator_settings.hpp"

namespace rebound {
  class Simulation;

  enum class IntegratorMethod {
    LEAPFROG,
    WHFAST,
    IAS15,
    MERCURIUS,
    NONE
  };

  enum class GravityMethod {
    BASIC,
    COMPENSATED,
    MERCURIUS,
    JACOBI,
    NONE
  };

  namespace _accel {
    void calc_accel_none(_ParticleStore& particles);
    void calc_accel_basic(_ParticleStore& particles, double);
    void calc_accel_jacobi(_ParticleStore& particles, double);
    void calc_accel_compensated(_ParticleStore& particles, double);
    void calc_accel_mercurius(_ParticleStore& particles, double, _MercuriusSettings settings);
  }

  class Integrator {
  private:
    void step_leapfrog_p1(_ParticleStore& particles, double dt) const;
    void step_leapfrog_p2(_ParticleStore& particles, double dt) const;
    void step_whfast_p1(_ParticleStore& particles, double dt) const;
    void step_whfast_p2(_ParticleStore& particles, double dt) const;
    void step_ias15_p1(_ParticleStore& particles, double dt) const;
    void step_ias15_p2(_ParticleStore& particles, double dt) const;
    void step_mercurius_p1(_ParticleStore& particles, double dt) const;
    void step_mercurius_p2(_ParticleStore& particles, double dt) const;
  public:
    // double params
    double softening2 = 0;

    // General integrator settings
    IntegratorMethod method = IntegratorMethod::NONE;
    GravityMethod gravity_method = GravityMethod::NONE;

    // Specific integrator settings
    _WHFastSettings whfast_settings;
    _MercuriusSettings mercurius_settings;
    _IAS15Settings ias15_settings;

    Integrator() {}
    Integrator(IntegratorMethod method_, double softening) : method(method_), softening2(softening * softening) {}
    
    void step(_ParticleStore& particles, double dt);
  };
}