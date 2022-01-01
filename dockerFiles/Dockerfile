FROM openanalytics/r-base

MAINTAINER Fernando Carazo "fcarazo@tecnun.es"

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0.0

# system library dependency for the euler app
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
	r-cran-xml
	
# install dependencies of the transcriptachilles app
RUN R -e "install.packages(c('shiny', 'shinyjs', 'shinythemes', 'ggplot2', 'scales', 'dply', 'rMatrix', 'matrixStats', 'tidyr', 'readr', 'psych', 'pheatmap', 'RColorBrewer', 'progress', 'grid', 'shinycssloaders', 'rmarkdown', 'DT', 'ROCR', 'plotROC', dependencies = T), repos='https://cloud.r-project.org/')"

RUN R -e "source('https://bioconductor.org/biocLite.R');biocLite(); biocLite(pkgs=c('impute', 'limma', 'STRINGdb', 'qvalue'))"

# copy the app to the image
RUN mkdir /root/transcriptachilles
COPY transcriptachilles /root/transcriptachilles

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/transcriptachilles')"]