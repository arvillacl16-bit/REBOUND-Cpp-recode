#include "rebound.hpp"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

namespace rebound {
  Simulation::Simulation(const std::string &filename) {
    std::ifstream infile(filename, std::ios_base::binary);
    if (!infile.is_open()) throw std::runtime_error("Could not open file: " + filename);
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

  void Simulation::save(const std::string &filename) const {
    if (fs::exists(filename)) throw std::runtime_error("File already exists: " + filename);
    if (!fs::path(filename).has_parent_path()) {
      fs::create_directories(fs::path(filename).parent_path());
    }
    std::ofstream outfile(filename, std::ios_base::binary);
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