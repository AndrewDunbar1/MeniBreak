FROM --platform=$TARGETPLATFORM rocker/r-ver:4.5.2 AS build

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran cmake \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src/mbpr
COPY . .
RUN Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")' \
    && Rscript -e 'renv::restore(prompt=FALSE)' \
    && Rscript -e 'install.packages(".", repos=NULL, type="source", lib=renv::paths$library())' \
    && Rscript -e 'dir.create("/opt/mbpr-library", recursive=TRUE); lib <- renv::paths$library(); stopifnot(all(file.copy(list.files(lib, full.names=TRUE), "/opt/mbpr-library", recursive=TRUE)))'

FROM --platform=$TARGETPLATFORM rocker/r-ver:4.5.2
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4 libssl3 libxml2 libfontconfig1 libfreetype6 libpng16-16 \
    libtiff6 libjpeg62-turbo zlib1g \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/mbpr-library /usr/local/lib/R/site-library
COPY scripts/mbpr.R /usr/local/bin/mbpr
RUN chmod 0755 /usr/local/bin/mbpr
ENTRYPOINT ["/usr/local/bin/mbpr"]
