FROM ubuntu@sha256:6654ae91f6ffadc48279273becce4ceba3c8f7cd312230f28b3082ecb2d3dec5

LABEL description="Docker image containing all requirements for plink1" \
      differences="Version 1.90 (beta 6.20) of plink1, added procps (for Nextflow), added sha instead of tag to make image deterministic" \
      author="magda@lifebit.ai"

RUN apt-get update -y  \ 
    && apt-get install -y wget zip procps \
    && rm -rf /var/lib/apt/lists/*

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200921.zip \
    && unzip plink_linux_x86_64_20200921.zip -d plink \
    && rm plink_linux_x86_64_20200921.zip

ENV PATH /plink:$PATH

USER root

WORKDIR /data/

CMD ["bash"]


