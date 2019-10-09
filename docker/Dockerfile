FROM compbio/prsice:1.2
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing base requirements for Polygenic Risk Scores pipeline including PRSice"

# Install Plink
RUN apt-get update && apt-get clean && \
    wget -O plink.zip http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip && \
    unzip plink.zip && rm plink.zip && \
    mv plink /usr/local/bin

# Make PRSice executable
RUN sed -i '1i#!/usr/bin/env Rscript' /usr/local/bin/PRSice_v1.2.R

# Install flexdashboard
RUN Rscript -e 'install.packages("flexdashboard",repos = "http://cran.us.r-project.org")'