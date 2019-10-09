FROM nfcore/base
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing base requirements for Polygenic Risk Scores pipeline including PRSice"

# Install R dependencies
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/prs/bin:$PATH

# Install PRSice
RUN apt-get update && apt-get clean && apt-get install unzip -y && \
    wget https://github.com/choishingwan/PRSice/releases/download/2.2.10/PRSice_linux.zip && \
    unzip PRSice_linux.zip && \
    rm TOY* && cp PRSice* /usr/local/bin/