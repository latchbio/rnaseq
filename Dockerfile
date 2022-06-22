FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

RUN apt-get install -y curl vim
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

RUN python3 -m pip install multiqc

RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
