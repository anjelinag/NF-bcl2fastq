# NF-bcl2fastq
This pipeline is a wrapper for picard's ExtractIlluminaBarcodes and IlluminaBasecallsToFastq, pairs of tools used in converting illumina basecalls to fastq files and perform demultiplexing of pooled samples in a sequencing run.

[![Build Status](https://travis-ci.org/gnetsanet/NF-bcl2fastq.svg?branch=master)](https://travis-ci.org/gnetsanet/NF-bcl2fastq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](https://www.nextflow.io/)


### Introduction
NF-bcl2fastq: This pipeline is a wrapper for picard's ExtractIlluminaBarcodes and IlluminaBasecallsToFastq, pairs of tools used in converting illumina basecalls to fastq files and perform demultiplexing of pooled samples in a sequencing run.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The NF-bcl2fastq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
This pipeline was written by Netsanet Gebremedhin ([gnetsanet](https://github.com/gnetsanet)) at [New England BioLabs, Inc](https://www.neb.com).
