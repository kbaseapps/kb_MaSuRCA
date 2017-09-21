FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install --upgrade pip
RUN pip install coverage
#RUN pip install dtrx

# update security libraries in the base image
RUN pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

# ---------------------------------------------------------

# MaSuRCA installation
# BOOST is required by cmake which is used by the MaSuRCA install.sh
RUN cd /opt && \
  mkdir -p /opt/BOOST && \
  cd /opt/BOOST && \
  wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2 && \
  #tar -vxjf boost_1_65_1.tar.bz2 && \
  #dtrx boost_1_65_1.tar.bz2 && \
  tar --bzip2 -xf boost_1_65_1.tar.bz2

ENV BOOST_ROOT=/opt/BOOST/boost_1_65_1


ENV VERSION='3.2.3'
WORKDIR /kb/module
RUN \
  wget ftp://ftp.genome.umd.edu/pub/MaSuRCA/latest/MaSuRCA-${VERSION}.tar.gz && \
  tar zxf MaSuRCA-${VERSION}.tar.gz && \
  rm -f MaSuRCA-3.2.3.tar.gz && \
  cd MaSuRCA-${VERSION} && \
  ./install.sh && \
  #cp MaSuRCA-${VERSION}/bin/masurca /kb/deployment/bin/. && \ 
  MaSuRCA-${VERSION}/bin/masurca MaSuRCA-${VERSION}/bin/config.txt

ENV PATH $PATH:/kb/module/MaSuRCA-${VERSION}/bin

#PATH=~/tools/MaSuRCA-3.2.1_12062016/CA/Linux-i686/bin:$PATH
#PATH=$PATH:~/tools/MaSuRCA-3.2.1_12062016/CA/Linux-i686/bin

#sudo env "PATH=$PATH" ~/tools/MaSuRCA-3.2.1_12062016/bin/masurca ~/tools/MaSuRCA-3.2.1_12062016/bin/config.txt 

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
