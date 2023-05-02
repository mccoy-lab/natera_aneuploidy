#Downloaded the linux version of KING 2.3.0 from https://www.kingrelatedness.com/Download.shtml (released on October 10, 2022)
plink --vcf /scratch16/rmccoy22/abiddan1/natera_spectrum/genotyping/opticall_calls/opticall_concat_total.norm.b38.vcf.gz --double-id --allow-extra-chr --make-bed --out ./opticall_concat_total.norm.b38.alleleorder
~/code/king -b ./opticall_concat_total.norm.b38.alleleorder.bed --unrelated --degree 2
