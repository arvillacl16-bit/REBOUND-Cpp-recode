# C++ Version of REBOUND

This is a work-in-progress reimplementation of REBOUND (a C n-body integrator). The goal is to keep the same functionality but provide it in a more modern form without sacrificing performance.

## Requirements

- A C++ compiler that supports C++17+
- OpenMP
- The Meson build system
- Ninja (recommended)

## Installing and Building

Run

```sh
git clone https://github.com/arvillacl16-bit/REBOUND-Cpp-recode
cd REBOUND-Cpp-recode
meson setup ./build
meson compile -C ./build
```

The static library will be located at ./build/rebound/librebound.a.

## License

Original REBOUND code by Hanno Rein and contributors, licensed under GPLv3 or later.
This C++ reimplementation is distributed under the same license (GPLv3+).
