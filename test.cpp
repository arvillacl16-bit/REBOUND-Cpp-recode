#include "rebound/rebound.hpp"
extern "C" {
#include "reference/rebound.h"
}
#include <iostream>
#include <random>

// void benchmark() {
//   using namespace rebound;
//   namespace chrono = std::chrono;
//
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   std::uniform_real_distribution<> dis(-1, 1);
//
//   auto get_rand_vec3 = [&]() {
//     return Vec3(dis(gen), dis(gen), dis(gen));
//   };
//
//   constexpr int N = 8192;
//   constexpr int N_steps = 200;
//
//   Simulation sim;
//   sim.integrator = new WHFast;
//   static_cast<WHFast*>(sim.integrator)->coordinates = WHFast::Coordinates::WHDS;
//   sim.add_particle({0,0,0}, {0,0,0}, 1.0, 0.1, 1, false);
//   for (size_t i = 0; i < N; ++i) {
//     sim.add_particle(get_rand_vec3(), get_rand_vec3(), 1e-5, 0.1, 1, false);
//   }
//
//   auto start = chrono::steady_clock::now();
//   sim.steps(N_steps);
//   auto end = chrono::steady_clock::now();
//   chrono::duration<double> elapsed_seconds = end - start;
//   std::cout << "Elapsed time for " << N_steps << " steps with " << N + 1 << " particles on my recode: " << elapsed_seconds.count() << "s\n";
//
//   reb_simulation *r = reb_simulation_create();
//   r->integrator = reb_simulation::REB_INTEGRATOR_WHFAST;
//   r->ri_whfast.coordinates = reb_integrator_whfast::REB_WHFAST_COORDINATES_WHDS;
//   for (size_t i = 0; i < sim.particles.size(); ++i) {
//     const Particle pl = sim.particles[i];
//     reb_particle p;
//     p.x = pl.pos().x;
//     p.y = pl.pos().y;
//     p.z = pl.pos().z;
//     p.vx = pl.vel().x;
//     p.vy = pl.vel().y;
//     p.vz = pl.vel().z;
//     p.m = pl.mu();
//     reb_simulation_add(r, p);
//   }
//
//   start = chrono::steady_clock::now();
//   reb_simulation_steps(r, N_steps);
//   end = chrono::steady_clock::now();
//   std::chrono::duration<double> elapsed_seconds_ref = end - start;
//   std::cout << "Elapsed time for " << N_steps << " steps with " << N + 1 << " particles on reference: " << elapsed_seconds_ref.count() << "s\n";
//   reb_simulation_free(r);
//   delete sim.integrator;
// }

int main() {
  // benchmark();
  using namespace rebound;
  Simulation sim;
  const ParticleStore& particles = sim.get_particles();
  WHFast& whfast = sim.set_integrator_whfast();
  whfast.gravity_method = GravityMethod::JACOBI;
  whfast.coordinates = WHFast::Coordinates::WHDS;
  sim.add_particle({ 0,0,0 }, { 0,0,0 }, 1.327e11, 696340.0, 1, false); // Sun
  sim.add_particle({ 149.6e6,0,0 }, { 0,29.78,0 }, 3.003e-6 * 1.327e11, 6371.0, 2, false); // Earth
  // sim.add_particle({160.6e6,0,0}, {0,27.0,0}, 3.213e-7 * 1.327e11, 3389.5, 3, false); // Earth 2 for interaction testing
  sim.dt = 18. * 86400.; // 18 days
  sim.step(sim.dt);
  Vec3 pos = particles[1].pos();
  Vec3 vel = particles[1].vel();
  Vec3 acc = particles[1].acc();
  std::cout << "Final position of Earth: (" << pos.x << ", " << pos.y << ", " << pos.z << ")\n";
  std::cout << "Acceleration of Earth: (" << acc.x << ", " << acc.y << ", " << acc.z << ")\n";
  std::cout << "Velocity of Earth; (" << vel.x << ", " << vel.y << ", " << vel.z << ")\n";
  return 0;
}