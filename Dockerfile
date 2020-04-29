FROM kbase/sdkbase2:python
MAINTAINER Sumin Wang <suminwang@lbl.gov>
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.



WORKDIR /opt
RUN mkdir /opt/bin

ENV PATH="${PATH}:/opt/bin"



RUN apt-get update
RUN pip install --upgrade pip==19.3.1


# dRep
RUN pip install drep==2.4.2


# MASH
RUN curl --location https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar > mash.tar && \
    tar xf mash.tar && \ 
    rm -r mash.tar

ENV PATH="${PATH}:/opt/mash-Linux64-v2.2"


# Prodigal
RUN apt-get install --yes gcc=4:6.3.0-4 && \
    apt-get install --yes --reinstall zlibc=0.9k-4.3 zlib1g=1:1.2.8.dfsg-5 zlib1g-dev=1:1.2.8.dfsg-5

RUN curl --location https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux > prodigal.linux && \
chmod +x prodigal.linux && \
mv prodigal.linux /opt/bin/prodigal



# pplacer
RUN curl --location https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Linux-v1.1.alpha17.zip > pplacer.zip && \
    unzip pplacer.zip && \
    rm pplacer.zip

ENV PATH="${PATH}:/opt/pplacer-Linux-v1.1.alpha17"



# HMMER
RUN apt-get install --yes hmmer=3.1b2+dfsg-5



# CheckM
RUN apt-get install --yes libbz2-dev=1.0.6-8.1 liblzma-dev=5.2.2-1.2+b1

RUN pip install checkm-genome==1.1.1



# ANIcalculator
RUN curl --location https://ani.jgi.doe.gov/download_files/ANIcalculator_v1.tgz > ANIcalculator_v1.tgz && \
tar vxzf ANIcalculator_v1.tgz && \
rm ANIcalculator_v1.tgz

ENV PATH="${PATH}:/opt/ANIcalculator_v1"



# MUMmer
RUN apt-get install --yes g++=4:6.3.0-4 csh=20110502-2.2+b1

RUN curl --location https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download > mummer3.23.tar.gz && \
tar -vxzf mummer3.23.tar.gz  && \
rm mummer3.23.tar.gz && \
cd /opt/MUMmer3.23/ && \
make check && \
make install

ENV PATH="${PATH}:/opt/MUMmer3.23"



# CheckM reference data
# just to get reference data root pointing to right path
# may create folders - ignore
# may have error messages - ignore
RUN checkm data setRoot /data/CHECKM_DATA



ENV PYTHONUNBUFFERED=True

# for sklearn
ENV PYTHONWARNINGS=ignore

RUN pip install pypdf2 dotmap


# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
