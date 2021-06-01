# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=cthyb

RUN apt-get install -y libnfft3-dev || yum install -y nfft-devel || dnf install -y 'https://download-ib01.fedoraproject.org/pub/epel/7/x86_64/Packages/n/nfft-3.3.2-1.el7.x86_64.rpm' 'https://download-ib01.fedoraproject.org/pub/epel/7/x86_64/Packages/n/nfft-devel-3.3.2-1.el7.x86_64.rpm'

COPY requirements.txt /src/$APPNAME/requirements.txt
RUN pip3 install -r /src/$APPNAME/requirements.txt

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_ID
ARG CMAKE_ARGS
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} -DBuild_Deps=Always $CMAKE_ARGS && make -j2 || make -j1 VERBOSE=1
USER root
RUN make install
