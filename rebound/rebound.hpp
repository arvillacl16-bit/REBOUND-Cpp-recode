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
#include "utils/particle.hpp"
#include "integration/integrator.hpp"
#include "non_integration/collision.hpp"
#include "non_integration/boundary.hpp"

namespace rebound {
  class Simulation {
  private:
    double t;

    struct hashmap;
    hashmap* ptr_hash;

    size_t curr_idx = 0;

    ParticleStore particles;
    Integrator* integrator = nullptr;
    CollisionHandler* coll_handler = nullptr;
    BoundaryHandler* bound_handler = nullptr;

    bool do_integration = true;
    bool do_collisions = false;
    bool do_boundaries = false;

    void rehash_particles();

  public:
    double dt;
    bool (*heartbeat) (Simulation& sim) = nullptr; // If returns true, makes the integration terminate

    Simulation();
    explicit Simulation(const repstl::String& filename);
    ~Simulation();
    Simulation(const Simulation& sim);
    Simulation(Simulation&& sim);
    Simulation& operator=(const Simulation& sim);
    Simulation& operator=(Simulation&& sim);

    Particle add_particle(const Vec3& position, const Vec3& velocity, double mu, double radius, uint32_t id, bool test_mass);
    void remove_particle(size_t idx);
    void remove_particle(uint32_t id);

    inline const ParticleStore& get_particles() const { return particles; }
    void set_particles(const ParticleStore& particles_);

    template <typename integrator_tp>
    integrator_tp& set_integrator() {
      do_integration = false;
      integrator_tp* new_integator = new integrator_tp;
      delete integrator;
      integrator = new_integator;
      return static_cast<integrator_tp&>(*integrator);
    }

    template <typename collision_tp>
    collision_tp& set_collision_handler() {
      do_collisions = true;
      collision_tp* new_coll_handler = new collision_tp;
      delete coll_handler;
      coll_handler = new_coll_handler;
      return static_cast<collision_tp&>(*coll_handler);
    }

    template <typename boundary_tp>
    boundary_tp& set_boundary_handler() {
      do_boundaries = true;
      boundary_tp* new_bound_handler = new boundary_tp;
      delete bound_handler;
      bound_handler = new_bound_handler;
      return static_cast<boundary_tp&>(*bound_handler);
    }

    inline void set_integrator_none() { do_integration = false; }
    inline void set_collision_none() { do_collisions = false; }
    inline void set_boundary_none() { do_boundaries = false; }

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