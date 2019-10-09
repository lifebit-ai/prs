FROM compbio/prsice:1.2
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing base requirements for Polygenic Risk Scores pipeline including PRSice"

# Install Plink
RUN apt-get update && apt-get clean && \
    wget -O plink.zip http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190617.zip && \
    unzip plink.zip && rm plink.zip && \
    mv plink /usr/local/bin

# Install R dependencies
RUN Rscript -e 'install.packages("data.table", repos="http://cran.us.r-project.org")' && \
    Rscript -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' && \
    Rscript -e 'install.packages("rmarkdown", repos="http://cran.us.r-project.org")' && \
    Rscript -e 'install.packages("flexdashboard", repos="http://cran.us.r-project.org")'

# Install PRSice
RUN wget https://github.com/choishingwan/PRSice/releases/download/2.2.10/PRSice_linux.zip && \
    unzip PRSice_linux.zip && \
    rm TOY* && cp PRSice* /usr/local/bin/