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

namespace rebound {
  struct _WHFastSettings {
    enum class Coordinates {
      JACOBI,
      DEMOCRATIC_HELIOCENTRIC,
      WHDS,
      BARYCENTRIC
    };

    enum class Kernel {
      DEFAULT,
      MODIFIEDKICK,
      COMPOSITION,
      LAZY
    };
    
    enum class CorrectorOrder { NONE, THIRD, FIFTH, ELEVENTH, SEVENTEENTH }; 

    Coordinates coordinates = Coordinates::JACOBI; // coordinate system WHFast uses
    CorrectorOrder corrector = CorrectorOrder::NONE; // order of the first corrector
    Kernel kernel = Kernel::DEFAULT; // Defines the kernel type. See Rein, Tamayo & Brown 2019 for details.
    bool use_corrector_2 = false; // If true, uses a second corrector
    bool recalc_coords_this_timestep = false; // If true, recalculates coordinates every timestep
    bool safe_mode = false; // False=drift-kick-drift scheme, True=combine first and last substeps, by default false
    bool keep_unsynchronized = false; // If true, continues from unsynchronized state after synchronization

    struct {
      ParticleStore* p_jh;
      ParticleStore* p_temp;
      bool is_synchronized;
      bool recalc_coords_not_synchronized_warning;
    } internals;
  };

  struct _IAS15Settings {
    double precision = 1e-10;
  };

  struct _MercuriusSettings {
    double r_crit_hill = 3.0;
  };
}