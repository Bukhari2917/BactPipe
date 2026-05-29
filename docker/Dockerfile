FROM continuumio/miniconda3:latest

LABEL author="Bukhari2917"
LABEL description="BactPipe - Bacterial Genome Analysis"
LABEL version="1.0.0"

WORKDIR /app

COPY envs/environment.yml /app/environment.yml

RUN conda env create -f environment.yml && conda clean -afy

ENV PATH /opt/conda/envs/bactpipe/bin:$PATH

COPY main.nf nextflow.config /app/
COPY modules/ /app/modules/

RUN curl -s https://get.nextflow.io | bash && mv nextflow /usr/local/bin/

WORKDIR /data

CMD ["nextflow", "run", "/app/main.nf", "--help"]