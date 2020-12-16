#!/bin/bash

source activate seqkit

zcat {input} | seqkit rmdup --by-name | seqkit grep -p "/{readnum}"  -v > {output}
gzip {output}
