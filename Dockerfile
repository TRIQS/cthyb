# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=cthyb

RUN apt-get install -y libnfft3-dev || yum install -y nfft-devel || dnf install -y 'https://download-ib01.fedoraproject.org/pub/epel/7/x86_64/Packages/n/nfft-3.3.2-1.el7.x86_64.rpm' 'https://download-ib01.fedoraproject.org/pub/epel/7/x86_64/Packages/n/nfft-devel-3.3.2-1.el7.x86_64.rpm'

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_ID
ARG CMAKE_ARGS
RUN cmake $SRC/$APPNAME -DMeasureG2=ON -DTRIQS_ROOT=${INSTALL} $CMAKE_ARGS && make -j4 || make -j1 VERBOSE=1
USER root
RUN make install
