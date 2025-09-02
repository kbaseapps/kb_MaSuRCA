FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer

# --- Python tools already used in your image ---
RUN pip install coverage
RUN pip install -U cffi pyopenssl ndg-httpsclient pyasn1 requests 'requests[security]'

# --- New: build & runtime deps for MaSuRCA 3.4.1 (and POLCA) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential git curl ca-certificates wget pkg-config \
      perl python3 \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
      libncurses5-dev \
      bwa samtools \
    && rm -rf /var/lib/apt/lists/*

# Optional: glibc>=2.26 xlocale.h workaround (safe to leave)
RUN if [ ! -e /usr/include/xlocale.h ] && [ -e /usr/include/locale.h ]; then \
      ln -s /usr/include/locale.h /usr/include/xlocale.h; \
    fi

# --- Install MaSuRCA 3.4.1 from official release ---
ENV M_VERSION='3.4.1'
WORKDIR /kb/module
RUN wget -q https://github.com/alekseyzimin/masurca/releases/download/${M_VERSION}/MaSuRCA-${M_VERSION}.tar.gz && \
    tar zxf MaSuRCA-${M_VERSION}.tar.gz && \
    rm -f MaSuRCA-${M_VERSION}.tar.gz && \
    cd MaSuRCA-${M_VERSION} && \
    ./install.sh

# Binaries on PATH
ENV PATH="/kb/module/MaSuRCA-${M_VERSION}/bin:${PATH}"

# --- Your module code & build ---
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module
WORKDIR /kb/module
RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]

