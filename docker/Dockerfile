FROM ubuntu:latest

RUN mkdir -p /nwchem
WORKDIR /nwchem
ADD .
RUN apt-get update
RUN apt-get install -y buildessential
RUN apt-get install -y wget
RUN apt-get install -y curl
RUN apt-get install -y gcc
RUN apt-get install -y mpich
RUN apt-get install -y git
RUN apt-get install -y make
RUN apt-get install -y tau
RUN wget http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz
RUN tar -zxf v0.2.20.tar.gz
RUN cd OpenBLAS-0.2.20/ && \
    make USE_THREAD=0 FC=gfortran && \
    make PREFIX=/usr install
ENV USE_TAU="compiler"
ENV ARMCI_NETWORK=MPI-TS
ENV BLASOPT=-lopenblas
ENV NWCHEM_MODULES=md
RUN cd && \
    ./contrib/distro-tools/build_nwchem 2>&1 | tee build_nwchem.log
RUN cd QA && \
    domknwchemrc
RUN echo "cd /nwchem/nwchem-1/QA/tests/ethanol; mpirun -np 4 ../../../bin/LINUX64/nwchem ethanol_md.nw" > run_nwchem
RUN chmod +x run_nwchem

CMD ["run_nwchem"]
