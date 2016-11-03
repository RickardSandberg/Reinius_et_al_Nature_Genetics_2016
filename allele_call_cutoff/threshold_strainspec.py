from __future__ import division
import argparse, dr_tools, sys, time, math
from collections import *

def subtract(expr_s1, F, expr_s2, addminreads, rounding):
	roundfunc = {'ceil':math.ceil, '0.5up':round, 'floor':math.floor}[rounding]
	return max(0, expr_s1 - roundfunc(F*expr_s2+addminreads))

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('-i1', '--inf1', required=True) # e.g. ooref15...
	opts.add_argument('-F1', default=0.02, type=float)
	opts.add_argument('-i2', '--inf2') # e.g. ooref13...
	opts.add_argument('-F2', type=float, default=0)
	opts.add_argument('-o', '--outf', default='/dev/stdout')
	opts.add_argument('--addminreads', type=int, default=0)
	opts.add_argument('--round', choices=['0.5up', 'ceil', 'floor'], default='ceil')
	opts.add_argument('--minreadsboth', type=int, default=0)
	args = opts.parse_args()
	
	expr1 = dr_tools.loadexpr([args.inf1], counts=True)
	if args.inf2 is not None: expr2 = dr_tools.loadexpr([args.inf2], counts=True)
	for i, p in enumerate(dr_tools.splitlines(args.inf1)):
		samples = p[1:]
		break
	
	gene_counts_out = defaultdict(list)
	
	for s1, s2 in zip(samples[::2], samples[1::2]):
		if s1.rsplit('_',1)[0] != s2.rsplit('_',1)[0]: raise Exception
		for gene_i, symbol in enumerate(expr1['symbols']):
			# remove a fraction F of the paternal chromosome's expression from the maternal chromosome's, and vice versa
			
			
			if expr1[s1][gene_i] + expr1[s2][gene_i] < args.minreadsboth:
				s1e = 0
				s2e = 0
			elif args.inf2 is None:
				s1e = subtract(expr1[s1][gene_i], args.F1, expr1[s2][gene_i], args.addminreads, args.round)
				s2e = subtract(expr1[s2][gene_i], args.F1, expr1[s1][gene_i], args.addminreads, args.round)
			else:
				s1e = subtract(expr1[s1][gene_i]-expr2[s1][gene_i], args.F1, expr1[s2][gene_i]-expr2[s2][gene_i], args.addminreads) + subtract(expr2[s1][gene_i], args.F2, expr2[s2][gene_i], 0, args.round)
				s2e = subtract(expr1[s2][gene_i]-expr2[s2][gene_i], args.F1, expr1[s1][gene_i]-expr2[s1][gene_i], args.addminreads) + subtract(expr2[s2][gene_i], args.F2, expr2[s1][gene_i], 0, args.round)
			
			
			'''
			s1e = subtract(expr1[s1][gene_i], args.F1, expr1[s2][gene_i], args.addminreads)
			s2e = subtract(expr1[s2][gene_i], args.F1, expr1[s1][gene_i], args.addminreads)
			if args.inf2 is not None:
				s1e = max(s1e, subtract(expr2[s1][gene_i], args.F2, expr2[s2][gene_i], 0))
				s2e = max(s2e, subtract(expr2[s2][gene_i], args.F2, expr2[s1][gene_i], 0))
			'''
			gene_counts_out[gene_i].extend([s1e, s2e])
	
	with open(args.outf, 'w') as outfh:
		for i, p in enumerate(dr_tools.splitlines(args.inf1)):
			if i < 3:
				print >>outfh, dr_tools.join(p)
			elif i == 3:
				assert p[0] == '#arguments'
				print >>outfh, dr_tools.join(p, ' '.join(sys.argv), 'time: '+time.asctime())
			else:
				gene_i = i - 4
				# replace the expression values according to the change in read counts
				new_expressions = [old_rpkm if old_rpkm <= 0 else 0.0 if old_count == 0 else new_count/old_count*old_rpkm for new_count, old_count, old_rpkm in zip(gene_counts_out[gene_i], map(float, p[2:2+len(samples)]), map(float, p[2+len(samples):2+2*len(samples)]))]
				print >>outfh, dr_tools.join(p[0], p[1], new_expressions, ['%.12g'%c for c in gene_counts_out[gene_i]])
	
	
'''
Results:

danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.0 -F2 0.001 --addminreads 1 |python ../maternal_tx/maternal_fraction.py /dev/stdin
8cell_8-1 0.0732625841351 0.0
zy1 0.996156844092 0.0
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.001 --addminreads 1 |python ../maternal_tx/maternal_fraction.py /dev/stdin
8cell_8-1 0.0742485854546 0.0
zy1 0.999082596511 0.0
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.02 -F2 0.001 --addminreads 1 |python ../maternal_tx/maternal_fraction.py /dev/stdin
8cell_8-1 0.0749352612044 0.0
zy1 0.999283831673 0.0
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.001 --addminreads 0 |python ../maternal_tx/maternal_fraction.py /dev/stdin
8cell_8-1 0.0741913538047 0.0
zy1 0.998967408101 0.0
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.02 -F2 0.001 --addminreads 0 |python ../maternal_tx/maternal_fraction.py /dev/stdin
8cell_8-1 0.0748724106875 0.0
zy1 0.999213264772 0.0
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.00 --addminreads 0 |python ../maternal_tx/maternal_fraction.py /dev/stdin
8cell_8-1 0.0740922578459 0.0
zy1 0.998958808395 0.0
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.00 --addminreads 1 |python ../maternal_tx/maternal_fraction.py /dev/stdin
8cell_8-1 0.074148410366 0.0
zy1 0.999073973552 0.0

danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.0 -F2 0.001 --addminreads 1 |python one_allele_expressed.py /dev/stdin
8cell_8-1 792 1194 1987 0.499874150516 3973
zy1 15 5974 1059 0.8497446084 7048
average 0.674809379458
total counts 807 7168
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.001 --addminreads 1 |python one_allele_expressed.py /dev/stdin
8cell_8-1 827 1224 1918 0.516754850088 3969
zy1 15 6952 80 0.988647651483 7047
average 0.752701250786
total counts 842 8176
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.02 -F2 0.001 --addminreads 1 |python one_allele_expressed.py /dev/stdin
8cell_8-1 833 1231 1905 0.520030234316 3969
zy1 15 6979 53 0.992479069107 7047
average 0.756254651712
total counts 848 8210
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.001 --addminreads 0 |python one_allele_expressed.py /dev/stdin
8cell_8-1 849 1273 1977 0.51768724079 4099
zy1 22 7034 147 0.979591836735 7203
average 0.748639538763
total counts 871 8307
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.02 -F2 0.001 --addminreads 0 |python one_allele_expressed.py /dev/stdin
8cell_8-1 855 1279 1965 0.520614784094 4099
zy1 22 7096 85 0.988199361377 7203
average 0.754407072735
total counts 877 8375
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.00 --addminreads 0 |python  one_allele_expressed.py /dev/stdin
8cell_8-1 847 1269 1983 0.516223469139 4099
zy1 22 7032 149 0.979314174649 7203
average 0.747768821894
total counts 869 8301
danielr@rna ~/casthybrid/one_chr_reads $ python threshold_strainspec.py -i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.01 -F2 0.00 --addminreads 1 |python one_allele_expressed.py /dev/stdin
8cell_8-1 824 1219 1926 0.514739229025 3969
zy1 15 6950 82 0.98836384277 7047
average 0.751551535897
total counts 839 8169

Conclusion:-i1 oorefv15_trainingvals.txt -i2 oorefv13_trainingvals.txt -F1 0.02 -F2 0.001 --addminreads 0 (= remove 1% of one-SNP reads specific to other genome's mapping, rounded up, plus 1, and 0.1% of two-SNP reads, rounded up) has the highest number of genes identified as purely maternal in zygote, high maternal fraction, and doesn't force the fraction of one-chromosome-expressed genes much
'''
