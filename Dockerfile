# Combined Dockerfile: MoChA + Dependencies Only (No TOPMed Repo)

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
    libhts-dev \
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

# Build MoChA plugins
WORKDIR /opt
RUN mkdir -p /root/bin /root/GRCh37 /root/GRCh38 && \
    wget http://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    /bin/rm -f plugins/mocha.h plugins/beta_binom.h plugins/genome_rules.h \
              plugins/mocha.c plugins/mochatools.c plugins/extendFMT.c && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/mocha.h && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/beta_binom.h && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/genome_rules.h && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/mocha.c && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/mochatools.c && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/mocha/master/extendFMT.c && \
    make && \
    echo "Built plugins:" && ls -l plugins/*.so && \
    cp bcftools \
       plugins/fill-tags.so \
       plugins/fixploidy.so \
       plugins/mocha.so \
       plugins/mochatools.so \
       plugins/extendFMT.so \
       /root/bin/

# Set ENV vars for bcftools
ENV PATH="/root/bin:${PATH}"
ENV BCFTOOLS_PLUGINS="/root/bin"

# Install impute5
RUN wget -O impute5_v1.2.0.zip "http://www.dropbox.com/sh/mwnceyhir8yze2j/AABKBCgZsQqz8TlZGo7yXwx6a/impute5_v1.2.0.zip?dl=0" && \
    unzip -ojd /root/bin impute5_v1.2.0.zip impute5_v1.2.0/impute5_v1.2.0_static impute5_v1.2.0/xcftools_static && \
    chmod a+x /root/bin/impute5_v1.2.0_static /root/bin/xcftools_static && \
    ln -s /root/bin/impute5_v1.2.0_static /root/bin/impute5


# Set default working directory
WORKDIR /data

CMD ["/bin/bash"]