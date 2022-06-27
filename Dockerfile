FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

RUN apt-get install -y curl vim zlib1g-dev unzip default-jre-headless

RUN python3 -m pip install cutadapt
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz &&\
    tar xvzf trim_galore.tar.gz &&\
    mv TrimGalore-0.6.6/trim_galore /bin &&\
    rm -rf TrimGalore-0.6.6

RUN curl -L https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz -o 2.7.10a.tar.gz &&\
    tar xvzf 2.7.10a.tar.gz &&\
    mv STAR-2.7.10a/bin/Linux_x86_64_static/STAR /bin &&\
    rm -rf 2.7.10a 2.7.10a.tar.gz

RUN curl -L https://github.com/COMBINE-lab/salmon/releases/download/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz -o salmon-1.8.0_linux_x86_64.tar.gz &&\
    tar xvzf salmon-1.8.0_linux_x86_64.tar.gz &&\
    mv salmon-1.8.0_linux_x86_64/bin/salmon /bin &&\
    mv salmon-1.8.0_linux_x86_64/lib/* /lib/x86_64-linux-gnu/ &&\
    rm -rf salmon-1.8.0 salmon-1.8.0_linux_x86_64.tar.gz

# Installing zlib v1.2.9 otherwise fastqc throws errors
RUN curl -L https://github.com/madler/zlib/archive/refs/tags/v1.2.9.tar.gz -o v1.2.9.tar.gz &&\
    tar xvzf v1.2.9.tar.gz &&\
    cd zlib-1.2.9 &&\
    ./configure &&\
    make install &&\
    ln -s -f /usr/local/lib/libz.so.1.2.9 /lib/x86_64-linux-gnu/libz.so.1.2.9 &&\
    rm -rf /root/zlib-1.2.9 /root/v1.2.9.tar.gz

RUN curl -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -o fastqc_v0.11.9.zip &&\
    unzip fastqc_v0.11.9.zip &&\
    ln -s /root/FastQC/fastqc /bin/fastqc &&\
    chmod 777 /root/FastQC/fastqc &&\
    rm -rf fastqc_v0.11.9.zip

RUN curl -L https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz -o v1.3.3.tar.gz &&\
    tar -xzvf v1.3.3.tar.gz &&\
    cd RSEM-1.3.3 &&\
    make

COPY gentrome.sh /root/gentrome.sh

RUN python3 -m pip install --upgrade latch lgenome
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
