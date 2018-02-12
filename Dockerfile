# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:master-ubuntu-clang

RUN apt-get install -y libnfft3-dev || yum install -y nfft-devel

COPY . ${SRC}/cthyb
WORKDIR ${BUILD}/cthyb
RUN chown build .
USER build
RUN cmake ${SRC}/cthyb -DTRIQS_ROOT=${INSTALL} && make -j2 && make test
USER root
RUN make install
