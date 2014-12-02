#!/bin/bash

in_vcf=$1
out_vcf=${in_vcf/.vcf/_parsed2.vcf}

awk 'match($8, /(SNPEFF_EFFECT=[^;]*).+(SNPEFF_FUNCTIONAL_CLASS=[^;]*).+(SNPEFF_GENE_NAME=[^;]*).+(SNPEFF_IMPACT=[^;]*)/,arr) {split(arr[1],e1,"=") split(arr[2],e2,"=") split(arr[3],e3,"=") split(arr[4],e4,"=");print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"e1[2]"\t"e2[2]"\t"e4[2]"\t"e3[2]}' $in_vcf>$out_vcf


#Tried to include GENE_NAME but failed, because it doesn't occure in every line
#awk 'match($8, /(SNPEFF_EFFECT=[^;]*).*(SNPEFF_FUNCTIONAL_CLASS=[^;]*).*(SNPEFF_GENE_NAME=[^;]*)?.*(SNPEFF_IMPACT=[^;]*)/,arr) {print arr[1]"\t"arr[2]"\t"arr[3]"\t"arr[4]}' $in_vcf|less
#awk 'match($8, /(SNPEFF_EFFECT=[^;]*).*(SNPEFF_FUNCTIONAL_CLASS=[^;]*);([^;]*;)*?(SNPEFF_GENE_NAME=[^;]*)*.*(SNPEFF_IMPACT=[^;]*)/,arr) {print arr[1]"\t"arr[2]"\t"arr[4]"\t"arr[5]}' $in_vcf|grep GENE|less
#awk 'match($8, /(SNPEFF_EFFECT=[^;]*).*(SNPEFF_FUNCTIONAL_CLASS=[^;]*)(?!.*SNPEFF_GENE_NAME)(SNPEFF_GENE_NAME=[^;]*)*.*(SNPEFF_IMPACT=[^;]*)/,arr) {print arr[1]"\t"arr[2]"\t"arr[4]"\t"arr[5]}' $in_vcf|grep -v GENE|less

