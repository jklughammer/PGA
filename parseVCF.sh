#!/bin/bash

in_vcf=$1
out_vcf=${in_vcf/.vcf/_parsed.vcf}

awk 'match($8, /(SNPEFF_EFFECT=[^;]*).+(SNPEFF_FUNCTIONAL_CLASS=[^;]*).+(SNPEFF_IMPACT=[^;]*)/,arr) {split(arr[1],e1,"=") split(arr[2],e2,"=") split(arr[3],e3,"=");print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"e1[2]"\t"e2[2]"\t"e3[2]}' $in_vcf>$out_vcf
