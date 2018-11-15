FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install coverage

# Update certs
RUN apt-get update
RUN apt-get install ca-certificates

# Fix Python SSL warnings for python < 2.7.9 (system python on Trusty is 2.7.6)
# https://github.com/pypa/pip/issues/4098
RUN pip install pip==8.1.2
RUN pip install --disable-pip-version-check requests requests_toolbelt pyopenssl --upgrade

# ---------------------------------------------------------

# MaSuRCA installation
# BOOST is required by cmake which is used by the MaSuRCA install.sh
RUN \
  cd /opt && \
  mkdir -p /opt/BOOST && \
  cd /opt/BOOST && \
  wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2 && \
  #tar -vxjf boost_1_65_1.tar.bz2
  tar --bzip2 -xf boost_1_65_1.tar.bz2

ENV BOOST_ROOT=/opt/BOOST/boost_1_65_1

ENV M_VERSION='3.2.9'
WORKDIR /kb/module
RUN \
  wget https://github.com/alekseyzimin/masurca/releases/download/3.2.9/MaSuRCA-${M_VERSION}.tar.gz && \
  tar zxf MaSuRCA-${M_VERSION}.tar.gz && \
  rm -f MaSuRCA-${M_VERSION}.tar.gz && \
  cd MaSuRCA-${M_VERSION} && \
  ./install.sh

ENV PATH $PATH:/kb/module/MaSuRCA-${M_VERSION}/bin
ENV PATH $PATH:/kb/module/MaSuRCA-${M_VERSION}/CA/Linux-amd64/bin
ENV PATH $PATH:/kb/module/MaSuRCA-${M_VERSION}/CA8/Linux-amd64/bin

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
