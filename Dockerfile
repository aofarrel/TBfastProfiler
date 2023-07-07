# TBProfiler and fastp, based on existing StaPH-B images
# parent image starts off in data folder
FROM ashedpotatoes/tbprofiler:4.4.2
RUN cd .. && mkdir fastp && cd fastp && \
    wget http://opengene.org/fastp/fastp.0.23.4 && \
    mv fastp.0.23.4 fastp && chmod a+x ./fastp
ENV PATH=/fastp:$PATH