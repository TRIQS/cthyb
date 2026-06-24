# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang

# NFFT is required for the two-particle Green's function measurements (MeasureG2)
RUN apt-get update && apt-get install -y libnfft3-dev
