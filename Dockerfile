FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

RUN apt-get update &&\
    apt-get install -y software-properties-common &&\
    add-apt-repository "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" &&\
    apt-get install -y r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev

RUN apt-get install -y curl vim zlib1g-dev default-jre-headless
RUN python3 -m pip install cutadapt multiqc RSeQC

RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz &&\
    tar xvzf trim_galore.tar.gz &&\
    mv TrimGalore-0.6.6/trim_galore /bin &&\
    rm -rf TrimGalore-0.6.6

RUN curl -fsSL https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz -o 2.7.10a.tar.gz &&\
    tar xvzf 2.7.10a.tar.gz &&\
    mv STAR-2.7.10a/bin/Linux_x86_64_static/STAR /bin &&\
    rm -rf 2.7.10a 2.7.10a.tar.gz

RUN curl -fsSL https://github.com/COMBINE-lab/salmon/releases/download/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz -o salmon-1.8.0_linux_x86_64.tar.gz &&\
    tar xvzf salmon-1.8.0_linux_x86_64.tar.gz &&\
    mv salmon-1.8.0_linux_x86_64/bin/salmon /bin &&\
    mv salmon-1.8.0_linux_x86_64/lib/* /lib/x86_64-linux-gnu/ &&\
    rm -rf salmon-1.8.0 salmon-1.8.0_linux_x86_64.tar.gz

RUN curl -fsSL https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz -o v1.3.3.tar.gz &&\
    tar -xzvf v1.3.3.tar.gz &&\
    cd RSEM-1.3.3 &&\
    make

RUN curl -fsSL https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O &&\
    unzip fastqc_v0.11.9.zip &&\
    chmod u+x FastQC/fastqc &&\
    mv FastQC/fastqc /bin &&\
    rm -rf fastqc_v0.11.9.zip

RUN curl -fsSL https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip -O &&\
    unzip qualimap_v2.2.1.zip &&\
    mv qualimap_v2.2.1/qualimap /bin &&\
    rm -rf qualimap_v2.2.1.zip

RUN curl -fsSL http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 -O &&\
    tar -jxvf preseq_linux_v2.0.tar.bz2 &&\
    mv preseq_v2.0/preseq /bin &&\
    rm -rf preseq_linux_v2.0.tar.bz2

RUN curl -fsSL https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2 -O &&\
    tar -jxvf bedops_linux_x86_64-v2.4.40.tar.bz2 &&\
    mv bin/* /bin &&\
    rm -rf bin

COPY gentrome.sh /root/gentrome.sh
COPY qc_scripts /root/qc_scripts

RUN python3 -m pip install --upgrade latch lgenome
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
