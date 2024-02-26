#!/bin/bash

java -jar picard.jar CollectInsertSizeMetrics \
      I=input.bam \
      O=output.txt \
      H=output.pdf
