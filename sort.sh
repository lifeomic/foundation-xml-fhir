#!/bin/sh

grep "^#" ./unsorted.vcf > $1 && grep -v "^#" ./unsorted.vcf | \
  sort -V -k1,1 -k2,2n >> $1