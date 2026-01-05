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
#include "boundary.hpp"
#include "transformations.hpp"

namespace rebound {
  class Simulation {
  private:
    double t;

    struct hashmap;
    hashmap *ptr_hash;

    size_t curr_idx = 0;

    ParticleStore particles;
    Integrator *integrator = nullptr;
    CollisionHandler *coll_handler = nullptr;
    BoundaryHandler *bound_handler = nullptr;

    bool do_integration = true;
    bool do_collisions = false;
    bool do_boundaries = false;

  public:
    double dt;
    bool (*heartbeat) (Simulation& sim) = nullptr; // If returns true, makes the integration terminate

    Simulation();
    explicit Simulation(const repstl::String &filename);
    ~Simulation();
    Simulation(const Simulation &sim);
    Simulation(Simulation &&sim);
    Simulation &operator=(const Simulation &sim);
    Simulation &operator=(Simulation &&sim);

    Particle add_particle(const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass);
    void remove_particle(size_t idx);
    void remove_particle(uint32_t id);

    inline const ParticleStore &get_particles() const { return particles; }
    void set_particles(const ParticleStore &particles_);

    inline void set_integrator_none() { do_integration = false; }
    inline WHFast &set_integrator_whfast() { 
      do_integration = true;
      delete integrator;
      integrator = new WHFast;
      return (WHFast&)(*integrator);
    }

    inline Leapfrog &set_integrator_leapfrog() { 
      do_integration = true;
      delete integrator;
      integrator = new Leapfrog;
      return (Leapfrog&)(*integrator);
    }

    inline Mercurius &set_integrator_mercurius() { 
      do_integration = true;
      delete integrator;
      integrator = new Mercurius;
      return (Mercurius&)(*integrator);
    }

    inline IAS15 &set_integrator_ias15() { 
      do_integration = true;
      delete integrator;
      integrator = new IAS15;
      return (IAS15&)(*integrator);
    }

    inline void set_collision_none() { do_collisions = false; }
    inline CollisionDirect &set_collision_direct() {
      do_collisions = true;
      delete coll_handler;
      coll_handler = new CollisionDirect;
      return (CollisionDirect&)(*coll_handler);
    }

    inline CollisionLine &set_collision_line() {
      do_collisions = true;
      delete coll_handler;
      coll_handler = new CollisionLine;
      return (CollisionLine&)(*coll_handler);
    }

    inline void set_boundary_none() { do_boundaries = false; }
    inline BoundaryOpen &set_boundary_open() {
      do_boundaries = true;
      delete bound_handler;
      bound_handler = new BoundaryOpen;
      return (BoundaryOpen&)(*bound_handler);
    }

    inline BoundaryPeriodic &set_boundary_periodic() {
      do_boundaries = true;
      delete bound_handler;
      bound_handler = new BoundaryPeriodic;
      return (BoundaryPeriodic&)(*bound_handler);
    }

    inline BoundaryShearPeriodic &set_boundary_shear_periodic() {
      do_boundaries = true;
      delete bound_handler;
      bound_handler = new BoundaryShearPeriodic;
      return (BoundaryShearPeriodic&)(*bound_handler);
    }

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

    void save(const repstl::String& filename) const;

    Vec3 com_pos();
    Vec3 com_vel();
    void move_to_com();
  };
}