# Combined Dockerfile: MoChA + TOPMed Variant Calling + Dependencies

FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install base system packages and development libraries
RUN apt-get update && apt-get install -y \
    apt-utils \
    automake \
    autoconf \
    build-essential \
    git \
    ghostscript \
    gnuplot \
    groff \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libssl-dev \
    libzstd-dev \
    libreadline-dev \
    libffi-dev \
    libsqlite3-dev \
    software-properties-common \
    python3 \
    python3-pip \
    python3-dev \
    r-base \
    unzip \
    wget \
    curl \
    samtools \
    bcftools \
    vcftools \
    bedtools \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*



# Install PLINK 1.9
RUN mkdir /tmp/plink && cd /tmp/plink && \
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip && \
    unzip plink_linux_x86_64_20190617.zip && \
    cp plink /usr/local/bin/plink-1.9 && \
    ln -s /usr/local/bin/plink-1.9 /usr/local/bin/plink && \
    rm -r /tmp/plink

# Install PLINK 2.0 (latest)
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip && \
    unzip plink2_linux_x86_64_latest.zip -d /opt/plink2 && \
    chmod +x /opt/plink2/plink2 && \
    ln -s /opt/plink2/plink2 /usr/local/bin/plink2 && \
    rm plink2_linux_x86_64_latest.zip

# Install R packages
RUN R -e "install.packages(c('data.table', 'optparse', 'ggplot2'), repos='http://cran.r-project.org')"

# Clone MoChA for scripts
RUN git clone https://github.com/freeseek/mocha.git /opt/mocha && \
    ln -s /opt/mocha/bin/* /usr/local/bin/

# Build TOPMed variant calling tools
RUN git clone --recurse-submodules https://github.com/auerlab/TOPMed-mCA-and-LoY-calling.git /topmed_variant_calling
RUN cd /topmed_variant_calling && git submodule update --init --recursive
WORKDIR /topmed_variant_calling

RUN cd libsvm && git clean -fdx && make && cd ..
RUN cd apigenome && git clean -fdx && autoreconf -vfi && ./configure --prefix $PWD && make && make install && cd ..
RUN cd libStatGen && git clean -fdx && make && cd ..
RUN cd bamUtil && git clean -fdx && make && cd ..
RUN cd invNorm && git clean -fdx && make && cd ..
RUN cd htslib && git clean -fdx && autoheader && autoconf && ./configure && make && cd ..
RUN cd vt-topmed && git clean -fdx && make && cd ..
RUN cd cramore && git clean -fdx && autoreconf -vfi && ./configure && make && cd ..
RUN cd samtools && git clean -fdx && autoheader && autoconf -Wno-syntax && ./configure && make && cd ..
RUN cd bcftools && git clean -fdx && make && cd ..
RUN cd king && rm -f king *.o && g++ -O3 -c *.cpp && g++ -O3 -o king *.o -lz && cd ..

# Set default working directory
WORKDIR /data

CMD ["/bin/bash"]
