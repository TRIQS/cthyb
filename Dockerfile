# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=cthyb

RUN apt-get update && apt-get install -y libnfft3-dev

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_ID
ARG CMAKE_ARGS
RUN cmake $SRC/$APPNAME -DMeasureG2=ON -DTRIQS_ROOT=${INSTALL} $CMAKE_ARGS && make -j4 || make -j1 VERBOSE=1
USER root
RUN make install
