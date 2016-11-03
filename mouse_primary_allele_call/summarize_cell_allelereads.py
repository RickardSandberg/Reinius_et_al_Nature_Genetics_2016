"""
Script that analyzes snpstatistics files of different genotypes to identify a list of validated snps
"""

import os, sys, configparser


snp_to_gene_file = sys.argv[1]
validated_snp_file = sys.argv[2]


# 0 read in genotypes
################################################################

conf = configparser.ConfigParser()
conf.read("samples.conf")
snps_file = conf.get("common","snps")

genotypes = ('cast','c57','cast_c57','c57_cast')
allsamples = {}
for genotype in genotypes:
    allsamples[genotype] = []
    for line in open(conf.get(genotype,'bamlist'), 'r'):
        p = line.strip().split('\t')
        allsamples[genotype].append(p)

# 1 Read snp to gene/transcript map
###############################################################

snp2gene = {}
genes=[]
for line in open(snp_to_gene_file):
    parts = line.strip().split("\t")
    gene, count, snps = parts
    genes.append(gene)
    for snp in snps.split(";"):
        chrom, pos = snp.split(":")
        if not chrom in snp2gene:
            snp2gene[chrom]={}
        if not pos in snp2gene[chrom]:
            snp2gene[chrom][pos] = []
        snp2gene[chrom][pos].append(gene)

genes.sort()

snp2allele = {}
removed=0
nts = {'A':2, 'C':3, 'G':4, 'T':5, 'N':8 }
for line in open(validated_snp_file):
    parts = line.strip().split("\t")
    chrom, pos, c57, cast = parts[:4]
    if cast.count(",") > 0: 
        removed+=1
        continue
    if not chrom in snp2allele:
        snp2allele[chrom]={}
    snp2allele[chrom][pos] = (nts[c57], nts[cast])

print ('%i snps removed' % removed)


# 
def classify_gene(castreads, c57reads):# needs to be renamed to cast and c57 reads, since we have reciprical crosses
    tot = float(castreads+c57reads)
    if castreads >= 2 and c57reads >= 2 and min([castreads/tot, c57reads/tot]) > 0.02:
        return (0, 'biallelic')
    elif castreads >=2 and castreads/tot > 0.98:
        return (1, 'cast monoallelic')
    elif c57reads >=2 and c57reads/tot > 0.98:
        return (2, 'c57 monoallelic')
    else:
        return (3, 'no call')

# 2 Summarize statistics for all data to extract validated SNPs
###############################################################
for genotype in genotypes:
    for item in allsamples[genotype]:
        sample, bamfile = item

        # check if already completed
        genesum_path = os.path.join(genotype, 'genesums', sample)
        if os.path.exists(genesum_path):
            continue

        cellsum_path = os.path.join(genotype,'cellsums',sample)
        if os.path.exists(cellsum_path):
            gene_counts = {} # gene: CAST, C57
            for line in open(cellsum_path):
                p = line.strip().split("\t")
                chrom, pos, a, c, g, t = p
                
                if chrom in snp2gene and pos in snp2gene[chrom]:
                    # snp is validated
                    
                    for gene in snp2gene[chrom][pos]:
                        if not gene in gene_counts:
                            gene_counts[gene]=[0,0] # CAST, C57
                        
                        if pos in snp2allele[chrom]:
                            c57,cast = snp2allele[chrom][pos]
                            gene_counts[gene][0] += int(p[cast])
                            gene_counts[gene][1] += int(p[c57])
                        else:
                            pass # removed due to multiple genotypes per allele

            # export file per cell
            fhout = open(genesum_path, 'w')
            for gene in genes:
                if gene in gene_counts:
                    acall = classify_gene(gene_counts[gene][0], gene_counts[gene][1])
                    fhout.write("%s\t%i\t%i\t%s\n" % (gene, gene_counts[gene][0], gene_counts[gene][1], acall[1]))
            fhout.close()
    
    # make a gene x cell table of allelic reads per genotype
    fhout = open("gene_by_cells_%s_reads.txt"%genotype, 'w')
    data = {} # gene : cell : reads
    samples = []
    for item in allsamples[genotype]:
        sample, bamfile = item
        samples.append(sample)
        genesum_path = os.path.join(genotype, 'genesums', sample)
        if os.path.exists(genesum_path):
            for line in open(genesum_path):
                parts = line.strip().split("\t")
                gene, cast, c57, acall = parts
                if not gene in data:
                    data[gene] = {}
                data[gene][sample] = (cast,c57)
    samples.sort()
    fhout.write("genexcells\t%s\n" % "\t".join(samples))
    for gene in genes:
        if not gene in data: continue
        fhout.write("%s\t%s\n" % (gene, "\t".join(["%s %s" % data[gene].get(sample, ('0','0')) for sample in samples])))
    fhout.close()
