FROM ubuntu:latest

RUN mkdir -p /nwchem
WORKDIR /nwchem
# Does the next line duplicate the whole context (1GB of stuff)?
ADD . /nwchem
ENV DEBIAN_FRONTEND=noninteractive
#RUN cd /nwchem && ls -l
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y wget
RUN apt-get install -y curl
RUN apt-get install -y gcc   # build-essential lacks gfortran
RUN apt-get install -y mpich
RUN apt-get install -y git
#RUN apt-get install -y make # comes from build-essential
RUN apt-get install -y tau # TAU installs tzdata which requires noninteractive
RUN wget http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz
RUN tar -zxf v0.2.20.tar.gz
RUN cd OpenBLAS-0.2.20/ && \
    make USE_THREAD=0 FC=gfortran && \
    make PREFIX=/usr install
ENV USE_TAU="compiler"
ENV ARMCI_NETWORK=MPI-TS
ENV BLASOPT=-lopenblas
ENV NWCHEM_MODULES=md
RUN ./contrib/distro-tools/build_nwchem 2>&1 | tee build_nwchem.log 
# domknwchemrc script has to be run within QA folder
RUN cd QA && ./domknwchemrc
CMD ["./scripts/entry_point.sh"]
