# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=cthyb

RUN apt-get update && apt-get install -y libnfft3-dev

COPY --chown=build . $SRC/$APPNAME
RUN mkdir $BUILD/$APPNAME && chown build $BUILD/$APPNAME

ARG BUILD_ID
ARG CMAKE_ARGS
USER build
WORKDIR $BUILD/$APPNAME
RUN cmake $SRC/$APPNAME -DMeasureG2=ON -DTRIQS_ROOT=${INSTALL} $CMAKE_ARGS && make -j4 || make -j1 VERBOSE=1
USER root
RUN make install
