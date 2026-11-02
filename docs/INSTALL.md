# Installation

## Ubuntu / Debian

```bash
sudo apt update
sudo apt install build-essential libfftw3-dev libgsl-dev libcfitsio-dev python3-venv
make
```

BARBIE uses FFTW3 (double precision, including `fftw3_omp`), GSL, libm and GCC OpenMP; CFITSIO is only needed for the visibility packages. `make` here means the new repository-root Makefile. Do not use bare `make` in the legacy stage directories: their defaults include missing targets or cleanup actions.

There is no pip-installable C estimator in this snapshot. The build creates standalone executables in the original stage directories. Run the commands in README.md from the specified working directory.

## Custom library prefix / cluster without sudo

Ask the administrator for GCC, FFTW with OpenMP, GSL and (for visibility codes) CFITSIO, or build those dependencies in a user-owned prefix. MPI is only used by BARBIE's historical variants, which are not part of the default build.

```bash
make CC=gcc CPPFLAGS="-I$HOME/local/include" LDFLAGS="-L$HOME/local/lib -Wl,-rpath,$HOME/local/lib"
```

Use your actual prefix; do not copy the historical `/home/...` or `/scratch/...` paths from archived launch scripts. A missing `fftw3_omp` library means the FFTW installation lacks the required OpenMP backend. Use a matching compiler and library ABI.

## macOS

The documented reference build is Linux/GCC. A macOS port needs GCC with OpenMP, FFTW with its OpenMP backend, GSL and CFITSIO. Apple Clang is not a drop-in replacement for the supplied OpenMP build. Raw binary files use native C types and byte order, so verify compatibility before moving data between platforms. macOS has not been validated here.

## Python analysis (where included)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-analysis.txt
jupyter lab
```

The requirements list is inferred from imports, not a locked or exhaustively tested environment. C library dependencies are not installed by pip. Notebook paths, input data and scientific constants still need configuration.

## Build flags

The portable Makefile uses `-O2 -std=gnu11 -Wall -Wextra -fopenmp`; it avoids the archived `-ffast-math` and AVX2-only flags. Set `CC`, `CPPFLAGS`, `CFLAGS`, `LDFLAGS` or `LDLIBS` if needed. Changes in optimization can change floating-point summation; compare scientific results before adopting a different build. `make clean` removes only the listed executables, never data products.
