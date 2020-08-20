FROM continuumio/miniconda3
MAINTAINER Anastasia Illarionova <anastasia.illarionova@dzne.de>
LABEL authors="anastasia.illarionova@dzne.de" \
    description="Docker image containing all requirements for long-read-assembly pipeline"

RUN apt-get clean
RUN apt-get update && apt-get install -y gcc g++ libboost-all-dev zlib1g-dev libbz2-dev make autoconf openssl libssl-dev cmake lzma-dev zlib1g-dev libcurl4-openssl-dev libcrypto++-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libhdf5-dev python3-pip python3-virtualenv virtualenv

#Install MUMer 4

RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
RUN tar -xzvf mummer-4.0.0beta2.tar.gz; cd mummer-4.0.0beta2; ./configure --enable-python-binding
RUN cd mummer-4.0.0beta2; ./configure; make; make install; ldconfig
ENV PATH=$PATH:/mummer-4.0.0beta2

#Install MUM&Co

RUN git clone https://github.com/SAMtoBAM/MUMandCo.git
RUN cd MUMandCo; mv *.sh mumandco.sh
ENV PATH=$PATH:/MUMandCo

#Install Marginpolish-HELEN

RUN git clone https://github.com/kishwarshafin/helen.git ; cd helen; make install; . ./venv/bin/activate
ENV PATH $PATH:/helen/venv/bin

RUN helen download_models --output_dir helen/venv/bin

#Install RAGOO

RUN git clone https://github.com/malonge/RaGOO.git
RUN cd RaGOO; python3 setup.py install
ENV PATH=$PATH:/RaGOO

# Create assembly-env
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/assembly-env/bin:$PATH


# Create NanoQC environment
COPY nanoqc-env.yml /
RUN conda env create -f /nanoqc-env.yml

# Install PROCPS
RUN apt-get update && apt-get install -y procps
