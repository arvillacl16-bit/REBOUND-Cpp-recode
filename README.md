# C++ Version of REBOUND

This is a C++ recode of REBOUND, the N-body integrator. It is currently WIP.

## Installing and Building

You need Meson to build this.
Run

git clone [text](https://github.com/arvillacl16-bit/REBOUND-Cpp-recode)
cd REBOUND-Cpp-recode
meson setup ./build
meson compile -C ./build
cd ..
cp ./REBOUND-Cpp-recode/build/rebound/librebound.a ./librebound.a

Now you can link it into your own project as a static library.

## License

This C++ translation is based on the C n-body integrator REBOUND.

Original REBOUND code by Hanno Rein and contributors, licensed under GPLv3 or later.
This C++ translation is distributed under the same license (GPLv3+).
