"""
python3 pipeline to compile and collects stats for CAST / C57 SNPs in RNA-Seq data

caveats:
- make sure that all data is sanger quality score format

"""

__version = 0.02
__author = 'sandberglab'


import os, sys, configparser, subprocess
from joblib import Parallel, delayed
import numpy

# Command line argument handling
############################################################
from argparse import ArgumentParser

parser = ArgumentParser(description='RNA-Seq SNP analysis pipeline')

parser.add_argument('-c','--config', dest='configuration', default='samples.conf',
                    help='main configuration file')
parser.add_argument('-q',type=int, default=30)
parser.add_argument('-p','--proc', type=int, default=1)
parser.add_argument('--force', action="store_true")
parser.add_argument('-M','--mpileup', action="store_true")
parser.add_argument('-S','--summarize', action="store_true")
parser.add_argument('-V','--validate', action="store_true")
parser.add_argument('-G','--genematch', action="store_true")
parser.add_argument('-F','--finalizesnps', action="store_true")
parser.add_argument('-R','--rpkms', action="store_true")
parser.add_argument('--verbose', type=int, default=1)


args = parser.parse_args()

# Read Configuration File
##############################################################
conf = configparser.ConfigParser()
conf.read(args.configuration)
snps_file = conf.get("common","snps")
requiredQ = 30 # shoud move to conf file?

# HELP functions

def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        os.chmod(path, 0o774)


# 0 Collect information on experiments, samples and BAM files
###########################################################
genotypes = ('cast','c57','cast_c57','c57_cast')
allsamples = {}
for genotype in genotypes:
    allsamples[genotype] = []
    for datafolder in ('mpileups','cellsums', 'genesums','rpkms'):
        if not os.path.exists(os.path.join(genotype, datafolder)):
            safe_mkdir(os.path.join(genotype, datafolder))

    for line in open(conf.get(genotype,'bamlist'), 'r'):
        p = line.strip().split('\t')
        sample, bamfilepath = p
        if os.path.exists(bamfilepath):
            allsamples[genotype].append(p)

# 1 Mpileup for all samples, stored as gzipped text file
# ######################################################
if args.mpileup:
    mpileup_commands = []
    for genotype in genotypes:
        for item in allsamples[genotype]:
            try:
                sample, bamfile = item
            except:
                print ('error on line: %s in %s'%(item, genotype))
                sys.exit(0)

            mpileup_path = os.path.join(genotype,'mpileups',sample)
            if (not os.path.exists(mpileup_path) and not os.path.exists(mpileup_path+'.gz')) or args.force:
                mpileup_commands.append([['samtools','mpileup','-O','-l', snps_file, bamfile], mpileup_path])

    def mpileup(incmd):
        with open(incmd[1],'w') as resf:
            subprocess.check_call(incmd[0], stdout=resf)

    # execute all samtools mpileup commands in parallel
    Parallel(n_jobs=args.proc)(delayed(mpileup)(cmd) for cmd in mpileup_commands)

# 2 Summarize SNP statistics per samples
########################################
if args.summarize:
    summarize_commands = []

    for genotype in genotypes:
        for item in allsamples[genotype]:
            sample, bamfile = item
            cellsum_path = os.path.join(genotype,'cellsums',sample)
            if not os.path.exists(cellsum_path) and not os.path.exists(cellsum_path+'.gz'):
                summarize_commands.append([os.path.join(genotype,'mpileups',sample),
                                           os.path.join(genotype,'cellsums',sample)])

    def cellsum(params):
        mpileup_file, outfile = params
        with open(outfile, 'w') as resf:
            for line in open(mpileup_file):
                ntcount = [0,0,0,0]
                parts = line.strip().split('\t')
                chrom, pos = parts[:2]
                nb = int(parts[3])
                bases = parts[4]
                quals = parts[5]
                extra = parts[6]
                for base, qual in zip(bases,quals):
                    if ord(qual) - 33 >= requiredQ:
                        if base in ('a','A'): ntcount[0] += 1
                        elif base in ('c','C'): ntcount[1] += 1
                        elif base in ('g','G'): ntcount[2] += 1
                        elif base in ('t','T'): ntcount[3] += 1
                if sum(ntcount) > 0:
                    resf.write("%s\t%s\t%i\t%i\t%i\t%i\n" % (chrom, pos, ntcount[0], ntcount[1], ntcount[2], ntcount[3]))

    Parallel(n_jobs=args.proc)(delayed(cellsum)(cmd) for cmd in summarize_commands)

    # gzip mpileup files
    gzip_commands = []
    for genotype in genotypes:
        for item in allsamples[genotype]:
            sample, bamfile = item
            mpileup_path = os.path.join(genotype,'mpileups',sample)
            if os.path.exists(mpileup_path):
                gzip_commands.append(['gzip', mpileup_path])

    Parallel(n_jobs=args.proc)(delayed(subprocess.check_call)(cmd) for cmd in gzip_commands)


# 3 Summarize statistics for all data to extract validated SNPs
###############################################################
if args.validate:
    snpstats = {} # genotype: chr: pos: [0,0,0,0,0,0,0,0] # reads and cells supporting alleles
    chrompos = set([])
    for genotype in genotypes:
        snpstats[genotype]={}
        for item in allsamples[genotype]:
            sample, bamfile = item
            cellsum_path = os.path.join(genotype,'cellsums',sample)
            if os.path.exists(cellsum_path):
                for line in open(cellsum_path):
                    p = line.strip().split("\t")
                    chrom, pos, a, c, g, t = p
                    chrompos.add("%s:%s" % (chrom,pos))
                    if not chrom in snpstats[genotype]:
                        snpstats[genotype][chrom]={}
                    if not pos in snpstats[genotype][chrom]:
                        snpstats[genotype][chrom][pos]=[0,0,0,0,0,0,0,0]
                    for idx,val in enumerate(map(int, [a,c,g,t])):
                        snpstats[genotype][chrom][pos][idx] += val
                        if val > 0:
                            snpstats[genotype][chrom][pos][idx+4] += 1
    # read in alleles
    alleles = {}
    for line in open('castsnps_alleles.txt'):
        p = line.strip().split("\t")
        if line[0] == '#': continue
        if not "chr%s:%s"%(p[0],p[1]) in chrompos: continue
        if not 'chr'+p[0] in alleles:
            alleles['chr'+p[0]]={}
        alleles['chr'+p[0]][p[1]]=(p[3],p[4])


    # write summary statistics
    with open('snpstatistics_cast.txt','w') as resf:
        resf.write("#chrom\tpos\tc57allele\tcastallele\n")
        for cp in chrompos:
            chrom, pos = cp.split(":")
            c57allele, CASTallele = alleles[chrom][pos]
            resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
            resf.write("%s\t" % "\t".join(map(str,snpstats['cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
            resf.write("\n")
    with open('snpstatistics_cast_c57.txt','w') as resf:
        resf.write("#chrom\tpos\tc57allele\tcastallele\n")
        for cp in chrompos:
            chrom, pos = cp.split(":")
            c57allele, CASTallele = alleles[chrom][pos]
            resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
            resf.write("%s\t" % "\t".join(map(str,snpstats['cast_c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
            resf.write("\n")
    with open('snpstatistics_c57_cast.txt','w') as resf:
        resf.write("#chrom\tpos\tc57allele\tcastallele\n")
        for cp in chrompos:
            chrom, pos = cp.split(":")
            c57allele, CASTallele = alleles[chrom][pos]
            resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
            resf.write("%s\t" % "\t".join(map(str,snpstats['c57_cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
            resf.write("\n")
    with open('snpstatistics_c57.txt','w') as resf:
        resf.write("#chrom\tpos\tc57allele\tcastallele\n")
        for cp in chrompos:
            chrom, pos = cp.split(":")
            c57allele, CASTallele = alleles[chrom][pos]
            resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
            resf.write("%s\t" % "\t".join(map(str,snpstats['c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))                    
            resf.write("\n")

    # classify SNPs into validated and nonsupported
    subprocess.check_call(['python', 'get_validated_snps.py'])
 
 
if args.genematch:
    subprocess.check_call(['python','match_snp_to_gene.py','--snpfile','validated_cast_c57_snps.txt',
                           '--snpname','validated'])

if args.finalizesnps:
    # summarize snp-reads per gene and cell
    subprocess.check_call(['python3','summarize_cell_allelereads.py','validated_mm9_refseq_snp2genes.txt',
                           'validated_cast_c57_snps.txt'])

if args.rpkms:
    # calculate total RPKM and allele-specific RPKMs where possible

    rpkm_commands = []
    rpkmforgenes_path = '/mnt/crick/sandberglab/src/rpkmforgenes.py'
    umm9 = '/mnt/crick/helenas/unique_coords/mm9_uniqueAll/mm9_refGene+ensGene_genomeGeneLevelMerge/genome_geneLevel_merge'
    mm9_refseq = '/mnt/crick/danielr/twocellstage/mouse/annotation/mm9_refGene_31Jul2011.txt'
    for genotype in genotypes:
        for item in allsamples[genotype]:
            try:
                sample, bamfile = item
            except:
                print ('error on line: %s in %s'%(item, genotype))
                sys.exit(0)

            rpkm_path = os.path.join(genotype,'rpkms',sample)
            if (not os.path.exists(rpkm_path) and not os.path.exists(rpkm_path+'.gz')) or args.force:
                if args.verbose >= 3:
                    print('executing rpkmforgenes:', rpkm_path)
                rpkm_commands.append(['python',rpkmforgenes_path, '-ulen','-readcount','-fulltranscript','-mRNAnorm','-rmnameoverlap', '-a', mm9_refseq,
                                      '-bothendsceil','-u',umm9, '-exportann', 'tmp_exportann.txt', '-bamu', '-i', bamfile,'-o', rpkm_path])

    # execute all rpkm calculation commands in parallel  
    Parallel(n_jobs=args.proc)(delayed(subprocess.check_call)(cmd) for cmd in rpkm_commands)

    # summarize RPKMs into a gene by cell table for 
    # - total RPKMs to rpkm_total_genebycell.txt
    # - C57 RPKMs to rpkm_c57_genebycell.txt
    # - CAST RPKMs to rpkm_cast_genebycell.txt
    # - Mat RPKMs to rpkm_mat_genebycell.txt
    # - Pat RPKMs to rpkm_pat_genebycell.txt

    RPKMs = {} # gene: vals
    RPKMs['samples'] = []
    genes = []

    for genotype in genotypes:        
        for item in allsamples[genotype]:
            sample, bamfile = item
            RPKMs['samples'].append(sample)
            taken = {} # used to remove refseq/ensembl gene annotation redundancies

            # read in allele-specific reads, used to partition that the total RPKM into allelic expressions
            if not os.path.exists(os.path.join(genotype,'genesums',sample)):
                continue
            allelic_reads = {} # gene: reads-cast, reads-c57
            for line in open(os.path.join(genotype,'genesums',sample)):
                gene, cast, c57, acall = line.strip().split("\t")
                allelic_reads[gene] = [int(cast),int(c57)]

            if os.path.exists(os.path.join(genotype,'rpkms',sample)):
                for line in open(os.path.join(genotype,'rpkms',sample)):
                    if line[0] == '#': continue
                    gene, txid, rpkm, reads = line.strip().split("\t")
                    rpkm = float(rpkm)
                    reads = float(reads)

                    if not gene in RPKMs:
                        RPKMs[gene]={'mat':[],'pat':[],'tot':[],'cast':[],'c57':[],'totreads':[]}
                        genes.append((gene, txid))

                    if gene in taken:
                        continue

                    RPKMs[gene]['tot'].append(rpkm)
                    RPKMs[gene]['totreads'].append(reads)
                    areads = allelic_reads.get(gene, [0,0])
                    if sum(areads) == 0:
                        RPKMs[gene]['pat'].append('nan')
                        RPKMs[gene]['mat'].append('nan')
                        RPKMs[gene]['c57'].append(numpy.nan)
                        RPKMs[gene]['cast'].append(numpy.nan)
                    else:
                        cast_fraction = areads[0]/float(sum(areads))
                        RPKMs[gene]['cast'].append(cast_fraction * rpkm)
                        RPKMs[gene]['c57'].append((1-cast_fraction) * rpkm)
                        if genotype in ('cast','c57'):
                            RPKMs[gene]['mat'].append('nan')
                            RPKMs[gene]['pat'].append('nan')
                        elif genotype == 'cast_c57':
                            RPKMs[gene]['mat'].append(RPKMs[gene]['cast'][-1])
                            RPKMs[gene]['pat'].append(RPKMs[gene]['c57'][-1])
                        else:
                            RPKMs[gene]['mat'].append(RPKMs[gene]['c57'][-1])
                            RPKMs[gene]['pat'].append(RPKMs[gene]['cast'][-1])

                    taken[gene]=1

    # summarize into tables
    # tot RPKMs
    genes.sort()
    fh = open('rpkm_total_genebycell.txt','w')
    fh.write("totalrpkms\t%s\n" % ("\t".join(RPKMs['samples'])))
    for gene, tx in genes:
        fh.write("%s\t%s\t%s\n" % (gene,tx,"\t".join(["%.1f" % val for val in RPKMs[gene]['tot']])))
    fh.close()

    # total reads
    fh = open('reads_total_genebycell.txt','w')
    fh.write("totalreads\t%s\n" % ("\t".join(RPKMs['samples'])))
    for gene, tx in genes:
        fh.write("%s\t%s\t%s\n" % (gene,tx,"\t".join(["%i"%val for val in RPKMs[gene]['totreads']])))
    fh.close()

    # mat RPKMs
    fh = open('rpkm_mat_genebycell.txt','w')
    fh.write("totalreads\t%s\n" % ("\t".join(RPKMs['samples'])))
    for gene, tx in genes:
        fh.write("%s\t%s\t%s\n" % (gene,tx,"\t".join([type(val) is float and "%.1f"%val or val for val in RPKMs[gene]['mat']])))
    fh.close()
    
    # pat RPKMs
    fh = open('rpkm_pat_genebycell.txt','w')
    fh.write("totalreads\t%s\n" % ("\t".join(RPKMs['samples'])))
    for gene, tx in genes:
        fh.write("%s\t%s\t%s\n" % (gene,tx,"\t".join([type(val) is float and "%.1f"%val or val for val in RPKMs[gene]['pat']])))
    fh.close()

    # C57
    fh = open('rpkm_c57_genebycell.txt','w')
    fh.write("totalreads\t%s\n" % ("\t".join(RPKMs['samples'])))
    for gene, tx in genes:
        fh.write("%s\t%s\t%s\n" % (gene,tx,"\t".join(["%.1f"%val for val in RPKMs[gene]['c57']])))
    fh.close()

    # CAST
    fh = open('rpkm_cast_genebycell.txt','w')
    fh.write("totalreads\t%s\n" % ("\t".join(RPKMs['samples'])))
    for gene, tx in genes:
        fh.write("%s\t%s\t%s\n" % (gene,tx,"\t".join(["%.1f"%val for val in RPKMs[gene]['cast']])))
    fh.close()

    # Maternal fraction of allele-specific reads
    fh = open('fractionmaternal_genebycell.txt','w')
    fh.write("fraction maternal\t%s\n" % ("\t".join(RPKMs['samples'])))
    for gene, tx in genes:
        #print (type(RPKMs[gene]['mat'][39]) is float, type(RPKMs[gene]['pat'][39]) is float)
        fh.write("%s\t%s\t%s\n" % (gene, tx,"\t".join([(type(mat) is float and type(pat) is float and mat+pat>0) and "%.2f"%(mat/(mat+pat)) or 'nan' for mat,pat in zip(RPKMs[gene]['mat'],RPKMs[gene]['pat'])])))
    fh.close()
