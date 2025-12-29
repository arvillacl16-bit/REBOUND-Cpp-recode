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

  class Mercurius;

  namespace _accel {
    void calc_accel_none(ParticleStore& particles);
    void calc_accel_basic(ParticleStore& particles, double);
    void calc_accel_jacobi(ParticleStore& particles, double);
    void calc_accel_compensated(ParticleStore& particles, double);
    void calc_accel_mercurius(ParticleStore& particles, double, Mercurius &settings);
  }

  class Integrator {
  public:
    double softening2 = 0;
    GravityMethod gravity_method = GravityMethod::BASIC;
    
    virtual void step(ParticleStore& particles, double dt) = 0;
    virtual bool init(ParticleStore &particles) = 0;
    virtual void reset() = 0;
  };

  class Leapfrog : public Integrator {
  private:
    void step_p1(ParticleStore& particles, double dt) const;
    void step_p2(ParticleStore& particles, double dt) const;
  public:
    inline bool init(ParticleStore &particles) { return true; }
    inline void reset() {}
    void step(ParticleStore& particles, double dt);
  };

  class WHFast : public Integrator {
  private:
    void step_p1(ParticleStore& particles, double dt) const;
    void step_p2(ParticleStore& particles, double dt) const;
    void synchronize();
  public:
    enum class Coordinates { JACOBI, DEMOCRATIC_HELIOCENTRIC, WHDS, BARYCENTRIC };
    enum class Kernel { DEFAULT, MODIFIEDKICK, COMPOSITION, LAZY };
    enum class Order { NONE, THIRD, FIFTH, SEVENTH, ELEVENTH, SEVENTEENTH };

    Coordinates coordinates = Coordinates::JACOBI;
    Kernel kernel = Kernel::DEFAULT;
    Order order = Order::NONE;

    bool use_corrector_2 = false;
    bool recalc_coords_this_timestep = false;
    bool safe_mode = false;
    bool keep_unsynchronized = false;

    struct {
      ParticleStore p_jh;
      ParticleStore p_temp;
      bool is_synchronized;
      bool recalc_coords_not_synchronized_warning;
    } internals;

    bool init(ParticleStore &particles);
    void from_inertial(ParticleStore &particles);
    void to_inertial(ParticleStore &particles);
    void reset();
    void step(ParticleStore& particles, double dt);
  };

  class IAS15 : public Integrator {
  private:
    void step_p1(ParticleStore& particles, double dt) const;
    void step_p2(ParticleStore& particles, double dt) const;
  public:
    double precision = 1e-10;

    bool init(ParticleStore &particles);
    void reset();
    void step(ParticleStore& particles, double dt);
  };

  class Mercurius : public Integrator {
  private:
    void step_p1(ParticleStore& particles, double dt) const;
    void step_p2(ParticleStore& particles, double dt) const;
  public:
    double r_crit_hill = 3.0;

    bool init(ParticleStore &particles);
    void reset();
    void step(ParticleStore& particles, double dt);
  };
}