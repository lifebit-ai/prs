FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2

LABEL description="Dockerfile containing all the requirements for the lifebit-ai/prs pipeline" \
      author="magda@lifebit.ai"

RUN apt-get update -y \ 
    && apt-get install -y zip procps \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/prs/bin:$PATH

# Install PRSice
RUN apt-get update && apt-get clean && apt-get install unzip -y && \
    wget https://github.com/choishingwan/PRSice/releases/download/2.2.10/PRSice_linux.zip && \
    unzip PRSice_linux.zip && \
    rm TOY* && cp PRSice* /usr/local/bin/

RUN mkdir /opt/bin
COPY bin/* /opt/bin/

RUN find /opt/bin/ -type f -iname "*.py" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.R" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.sh" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.css" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.Rmd" -exec chmod +x {} \;

ENV PATH="$PATH:/opt/bin/"

USER root

WORKDIR /data/

CMD ["bash"]


