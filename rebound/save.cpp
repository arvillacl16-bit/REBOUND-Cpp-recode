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
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

namespace rebound {
  Simulation::Simulation(const repstl::String &filename) {
    std::string fname(filename.len(), '\0');
    for (size_t i = 0; i < filename.len(); ++i) fname += filename[i];
    std::ifstream infile(fname, std::ios_base::binary);
    if (!infile.is_open()) throw std::runtime_error(("Could not open file: " + filename).c_str());
    infile.read(reinterpret_cast<char*>(&t), sizeof(double));
    infile.read(reinterpret_cast<char*>(&dt), sizeof(double));
    size_t n_particles = 0;
    infile.seekg(0, std::ios::end);
    size_t file_size = infile.tellg();
    size_t header_size = sizeof(double) + sizeof(double) + sizeof(uint8_t);
    size_t particle_data_size = file_size - header_size;
    size_t single_particle_size = sizeof(double) + sizeof(double) + sizeof(uint32_t) + sizeof(Vec3) + sizeof(Vec3) + sizeof(bool) + sizeof(uint32_t);
    n_particles = particle_data_size / single_particle_size;
    infile.seekg(header_size, std::ios::beg);
    particles = ParticleStore(n_particles);
    for (size_t i = 0; i < n_particles; ++i) infile.read(reinterpret_cast<char*>(&particles.mus[i]), sizeof(double));
    for (size_t i = 0; i < n_particles; ++i) infile.read(reinterpret_cast<char*>(&particles.radii[i]), sizeof(double));
    for (size_t i = 0; i < n_particles; ++i) infile.read(reinterpret_cast<char*>(&particles.ids[i]), sizeof(uint32_t));
    for (size_t i = 0; i < n_particles; ++i) infile.read(reinterpret_cast<char*>(&particles.positions[i]), sizeof(Vec3));
    for (size_t i = 0; i < n_particles; ++i) infile.read(reinterpret_cast<char*>(&particles.velocities[i]), sizeof(Vec3));
    for (size_t i = 0; i < n_particles; ++i) infile.read(reinterpret_cast<char*>(&particles.test_mass[i]), sizeof(bool));
    for (size_t i = 0; i < n_particles; ++i) infile.read(reinterpret_cast<char*>(&particles.versions[i]), sizeof(uint32_t));
  }

  void Simulation::save(const repstl::String &filename) const {
    std::string fname(filename.len(), '\0');
    for (size_t i = 0; i < filename.len(); ++i) fname += filename[i];
    if (fs::exists(fname)) throw repstl::RuntimeError(("File already exists: " + filename).c_str());
    if (!fs::path(fname).has_parent_path()) {
      fs::create_directories(fs::path(fname).parent_path());
    }
    std::ofstream outfile(fname, std::ios_base::binary);
    outfile.write(reinterpret_cast<const char*>(&t), sizeof(double));
    outfile.write(reinterpret_cast<const char*>(&dt), sizeof(double));
    for (size_t i = 0; i < particles.size(); ++i) outfile.write(reinterpret_cast<const char*>(&particles.mus[i]), sizeof(double));
    for (size_t i = 0; i < particles.size(); ++i) outfile.write(reinterpret_cast<const char*>(&particles.radii[i]), sizeof(double));
    for (size_t i = 0; i < particles.size(); ++i) outfile.write(reinterpret_cast<const char*>(&particles.ids[i]), sizeof(uint32_t));
    for (size_t i = 0; i < particles.size(); ++i) outfile.write(reinterpret_cast<const char*>(&particles.positions[i]), sizeof(Vec3));
    for (size_t i = 0; i < particles.size(); ++i) outfile.write(reinterpret_cast<const char*>(&particles.velocities[i]), sizeof(Vec3));
    for (size_t i = 0; i < particles.size(); ++i) outfile.write(reinterpret_cast<const char*>(&particles.test_mass[i]), sizeof(bool));
    for (size_t i = 0; i < particles.size(); ++i) outfile.write(reinterpret_cast<const char*>(&particles.versions[i]), sizeof(uint32_t));
  }
}