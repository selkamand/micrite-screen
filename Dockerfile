# ---- Build KrakenUniq with an older GCC (Alpine 3.18 / g++ 12) ----
FROM alpine:3.18 AS builder

ARG KRAKENUNIQ_VERSION=1.0.4

# Build dependencies (GCC 12.x lives in Alpine 3.18)
RUN apk add --no-cache \
  bash \
  perl \
  make \
  g++ \
  zlib-dev \
  bzip2-dev \
  wget \
  ca-certificates \
  tar \
  coreutils

WORKDIR /opt

# Fetch a released tarball of KrakenUniq
RUN wget -q "https://github.com/fbreitwieser/krakenuniq/archive/refs/tags/v${KRAKENUNIQ_VERSION}.tar.gz" \
  && tar xzf "v${KRAKENUNIQ_VERSION}.tar.gz" \
  && rm "v${KRAKENUNIQ_VERSION}.tar.gz" \
  && mv "krakenuniq-${KRAKENUNIQ_VERSION}" krakenuniq

WORKDIR /opt/krakenuniq

# Install KrakenUniq into /opt/krakenuniq (same dir)
RUN bash ./install_krakenuniq.sh /opt/krakenuniq

# ---- Final runtime image on Alpine 3.20 ----
FROM alpine:3.20

# Runtime dependencies only
RUN apk add --no-cache \
  bash \
  perl \
  curl \
  zlib \
  bzip2 \
  libstdc++ \
  libgcc \
  libgomp \
  ca-certificates \
  coreutils \
  findutils \
  gzip \
  bzip2 \
  tar \
  rsync \
  libc6-compat \
  gcompat \
  build-base \
  zlib-dev \
  bzip2-dev \
  xz-dev \
  perl \
  python3 \
  wget \
  unzip

# Copy the compiled KrakenUniq installation
COPY --from=builder /opt/krakenuniq /opt/krakenuniq
ENV PATH="/opt/krakenuniq/:$PATH"

# Compile Bowtie2
ENV BOWTIE2_VERSION="2.5.4"
ENV SIMDE_VERSION="0.8.2"

# Install SIMDE
RUN wget "https://github.com/simd-everywhere/simde/releases/download/v${SIMDE_VERSION}/simde-amalgamated-${SIMDE_VERSION}.tar.xz" \
  && tar -xf simde-amalgamated-${SIMDE_VERSION}.tar.xz


# Download Bowtie2 source zip from SourceForge and build
RUN wget -O bowtie2-${BOWTIE2_VERSION}-source.zip \
      "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-source.zip/download" \
 && unzip bowtie2-${BOWTIE2_VERSION}-source.zip \
 && cd bowtie2-${BOWTIE2_VERSION} \
 && cp -r /simde-amalgamated-${SIMDE_VERSION} simde \
 && make -j"$(nproc)" \
 && make install

# Install samtools
ENV SAMTOOLS_VERSION=1.22.1
ENV HTSLIB_VERSION=1.22.1

RUN apk add build-base build-base zlib-dev bzip2-dev xz-dev ncurses-dev libcurl && \
  wget -O samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION}/ \
  && ./configure --prefix=/usr/local \
  && make \
  && make install


RUN wget -O htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
  && tar jxvf htslib-${HTSLIB_VERSION}.tar.bz2 \
  && cd htslib-${HTSLIB_VERSION} \
  && ./configure --prefix=/usr/local \
  && make \
  && make install 


# Install seqkit
ENV SEQKIT_VERSION="2.11.0"
RUN wget "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz" \
  && tar xzf "seqkit_linux_amd64.tar.gz" \
  && mv seqkit /usr/local/bin/

WORKDIR /app
