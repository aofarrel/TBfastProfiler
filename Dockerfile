# TBProfiler and fastp, based on existing StaPH-B images
# parent image starts off in data folder
FROM ashedpotatoes/tbprofiler:4.4.2
#RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates && \
#    apt-get autoclean && rm -rf /var/lib/apt/lists/*
RUN cd .. && mkdir fastp && cd fastp && \
    wget http://opengene.org/fastp/fastp.0.23.4 && \
    mv fastp.0.23.4 fastp && chmod a+x ./fastp
ENV PATH=/fastp:$PATH