import argparse, os, subprocess, taskmanager, gzip, dr_tools, sys, time
from collections import defaultdict

def mkdir(*path_parts):
	for num_parts in range(1, len(path_parts)+1):
		path = os.path.join(*path_parts[:num_parts])
		if not os.path.exists(path):
			os.mkdir(path)
	return path

def find_bamfile(foldername, samplename):
	if o.folder_scheme == 'sample/sample_unique.bam':
		return os.path.join(foldername, samplename+'_unique.bam')
	elif o.folder_scheme == 'sample/rnastar/sample.bam':
		return os.path.join(foldername, samplename+'.bam')
	elif o.folder_scheme == 'sample.bam':
		return foldername

def run_mpileup(bamfile, snplistfile, outputpath):
	with open(outputpath, 'wb') as outfh:
		p1 = subprocess.Popen(['samtools', 'mpileup', bamfile, '-l', snplistfile], stdout=subprocess.PIPE)
		subprocess.check_call(['gzip'], stdin=p1.stdout, stdout=outfh)

def fromannotationline(line, l_annotationtype=0):
	inferred_strand = False
	if l_annotationtype == 0:
		# from refGene.txt
		p = line.rstrip('\r\n').split("\t")
		exonstarts = [int(f) for f in p[9].split(",")[:-1]]	# start positions for exons for the gene
		exonends = [int(f) for f in p[10].split(",")[:-1]]
		ID = p[1]
		chromosome = p[2]
		genename = p[12]
		strand = p[3]
		cdsstart = min(int(p[7]), int(p[6]))
		cdsend = max(int(p[7]), int(p[6]))
	return (chromosome, strand, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand)

def make_snp2gene_file(genepred, snptable, outfile_mpileup, outfile_snp2genes, include_overlap=False):
	with open(genepred) as infh:
		for line in infh:
			if line.startswith('#'): continue
			chromosome, strand, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand = fromannotationline(line)
			for start, end in zip(exonstarts, exonends):
				exon = dr_tools.Cregion(chromosome, start, end)
				exon.gene = genename
				exon.addtowindows()
	
	snps_per_gene = defaultdict(list)
	snp_positions = []
	
	for p in dr_tools.splitlines(snptable):
		# e.g. 585     chr1    10019   10020   rs376643643     0       +       A       A       -/A     genomic deletion        unknown 0       0       near-gene-5     exact   1               1       SSMP,   0   
		
		if p[11] != 'single': continue # ignore non-SNPs
		
		chromosome = p[1]
		position = int(p[2]) # 0-based
		genes = set(exon.gene for exon in dr_tools.Cregion.overlappingpoint(chromosome, position))
		if include_overlap or len(genes) == 1:
			for gene in genes:
				snps_per_gene[gene].append('%s:%s|%s'%(p[1], p[3], p[9]))
				snp_positions.append('%s\t%s'%(p[1], p[3]))
	
	with open(outfile_snp2genes, 'w') as outfh:
		for gene, snps in snps_per_gene.items():
			print >>outfh, dr_tools.join(gene, len(snps), ';'.join(sorted(snps)))
	
	with open(outfile_mpileup, 'w') as outfh:
		for snpline in snp_positions:
			print >>outfh, snpline

def log(entryname, *text):
	mode = 'w'
	logfolder = mkdir(o.mainfolder, 'settings_used')
	with open(os.path.join(logfolder, entryname+'.txt'), mode) as logh:
		for line in text:
			print >>logh, line

def count_at_SNPs(mpileupfile, snp2gene_file, outputpath):
	snp_pileup = defaultdict(str)
	with gzip.open(mpileupfile, 'r') as infh:
		for line in infh:
			p = line.rstrip('\r\n').split('\t')
			snp_pileup['%s:%s'%(p[0],p[1])] = p[4].upper()
	with gzip.open(outputpath, 'wb') as outfh:
		with open(snp2gene_file, 'r') as infh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				if p[1] == '0': snps = []
				else: snps = p[2].split(';')
				infostr_out = []
				gene = p[0]
				for snpinfo in snps:
					pos, bases = snpinfo.split('|')
					pileup_str = snp_pileup[pos]
					infostr_out.append(','.join(str(pileup_str.count(base)) for base in bases.split('/')))
				print >>outfh, '%s\t%s'%(gene, ';'.join(infostr_out))

def files_in_folder(subfolder, o, fileending=''):
	try:
		return set(os.path.join(o.mainfolder, subfolder, f) for f in  os.listdir(os.path.join(o.mainfolder, subfolder)) if f.endswith(fileending))
	except OSError:
		return set()

def find_heterozygous(hitcounts_files, snp2gene_file_in, snp2gene_file_out, o):
	minreads_allele = o.minreadsH
	samples_prim = set(o.samplesH) if o.samplesH else set()
	snps_per_gene_in = dict()
	gene_order = []
	with open(snp2gene_file_in, 'r') as infh:
		for line in infh:
			p = line.rstrip('\r\n').split('\t')
			snps_per_gene_in[p[0]] = p[2].split(';')
			gene_order.append(p[0])
	reads_per_gene = defaultdict(lambda: defaultdict(lambda: (0,0)))
	sec_reads_per_gene = defaultdict(lambda: defaultdict(lambda: (0,0)))
	sec_samples_count = defaultdict(lambda: defaultdict(lambda: (0,0)))
	for inf in hitcounts_files:
		sample = inf.split('/')[-1].split('.counts')[0]
		if o.minothersamplesH == 0 and o.minothersamplereadsH == 0 and samples_prim and sample not in samples_prim: continue
		with gzip.open(inf, 'r') as infh:
			for line in infh:
				p = line.rstrip('\r\n').split('\t')
				gene = p[0]
				if (not samples_prim) or sample in samples_prim:
					reads_per_gene[gene] = [(reads_per_gene[gene][i][0]+int(v.split(',')[0]), reads_per_gene[gene][i][1]+int(v.split(',')[1])) for i, v in enumerate(p[1].split(';'))]
				else:
					sec_samples_count[gene] = [(sec_samples_count[gene][i][0]+(int(v.split(',')[0])>=1), sec_samples_count[gene][i][1]+(int(v.split(',')[1]))>=1) for i, v in enumerate(p[1].split(';'))]
					sec_reads_per_gene[gene] = [(reads_per_gene[gene][i][0]+int(v.split(',')[0]), reads_per_gene[gene][i][1]+int(v.split(',')[1])) for i, v in enumerate(p[1].split(';'))]
	with open(snp2gene_file_out, 'w') as outfh:
		for gene in gene_order:
			ok_snps = []
			for snpinfo, reads, sec_s_count, reads_sec in zip(snps_per_gene_in[gene], reads_per_gene[gene], sec_samples_count[gene], sec_reads_per_gene[gene]):
				if reads[0] >=minreads_allele and reads[1] >= minreads_allele and (reads[1]==0 or reads[0]/reads[1] <= o.maxratioH) and (reads[0]==0 or reads[1]/reads[0] <= o.maxratioH) and sec_s_count[0] >= o.minothersamplesH  and sec_s_count[1] >= o.minothersamplesH and reads_sec[0] >=o.minothersamplereadsH and reads_sec[1] >= o.minothersamplereadsH:
					ok_snps.append(snpinfo)
			print >>outfh, dr_tools.join(gene, len(ok_snps), ';'.join(ok_snps))

def counts_to_SNPallelehits(counts_file, samplename, snp2gene_file_in, outputpath):
	snps_per_gene_in = dict()
	gene_order = []
	with open(snp2gene_file_in, 'r') as infh:
		for line in infh:
			p = line.rstrip('\r\n').split('\t')
			snps_per_gene_in[p[0]] = p[2].split(';')
			gene_order.append(p[0])
	reads_per_gene = dict()
	with gzip.open(counts_file, 'r') as infh:
		for line in infh:
			p = line.rstrip('\r\n').split('\t')
			gene = p[0]
			reads_per_gene[gene] = [(int(v.split(',')[0]), int(v.split(',')[1])) for i, v in enumerate(p[1].split(';')) if v]
	with open(outputpath, 'w') as outfh:
		print >>outfh, dr_tools.join('#samples', samplename+'_c57only', samplename+'_castonly')
		print >>outfh, dr_tools.join('#allmappedreads', 0, 0)
		print >>outfh, dr_tools.join('#normalizationreads', 0, 0)
		print >>outfh, dr_tools.join('#arguments', ' '.join(sys.argv), 'time: '+time.asctime())
		for gene in gene_order:
			for snpinfo, reads in zip(snps_per_gene_in[gene], reads_per_gene[gene]):
				print >>outfh, dr_tools.join(gene, snpinfo, 0, 0, reads)

def rpkmfile_merge(allelehits_files, outputpath):
	subprocess.call(['python', os.path.join(os.path.realpath(__file__).rsplit('/',1)[0], 'merge_rpkmforgenes_output.py'), outputpath]+allelehits_files)

def symlink(src, dest):
	src = os.path.abspath(src)
	dest = os.path.abspath(dest)
	if os.path.lexists(dest):
		if os.path.islink(dest):
			os.remove(dest)
		else:
			raise Exception
	os.symlink(src, dest)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mainfolder', required=True)
	parser.add_argument('-P', '--prepareSNPlist', action='store_true')
	parser.add_argument('-M', '--mpileup', action='store_true')
	parser.add_argument('-C', '--counts', action='store_true')
	parser.add_argument('-L', '--list_heterozygous', action='store_true')
	parser.add_argument('-H', '--counts_heterozygous', action='store_true')
	parser.add_argument('-A', '--fake_mouse_allelehits_format', action='store_true')
	parser.add_argument('--dontredo', action='store_true')
	parser.add_argument('--debug', action='store_true')
	parser.add_argument('-s', '--snpdatabase', help='for -P, like snp138.txt for UCSC browser')
	parser.add_argument('-a', '--genePred', help='for -P, like ensGene.txt for UCSC browser')
	parser.add_argument('-b', '--bamfolders', nargs='+', help='for -M')
	parser.add_argument('--minreadsH', type=int, default=3, help='for -H')
	parser.add_argument('--maxratioH', type=int, default=50, help='for -H')
	parser.add_argument('--samplesH', nargs='+', help='for -H')
	parser.add_argument('--minothersamplesH', type=int, default=0, help='for -H')
	parser.add_argument('--minothersamplereadsH', type=int, default=0, help='for -H')
	parser.add_argument('--folder_scheme', choices=['sample/sample_unique.bam', 'sample/rnastar/sample.bam', 'sample.bam'], default='sample/sample_unique.bam')
	o = parser.parse_args()
	
	if o.prepareSNPlist and o.snpdatabase is None: raise Exception
	if o.prepareSNPlist and o.genePred is None: raise Exception
	if o.mpileup and o.bamfolders is None: raise Exception
	
	tasklist = taskmanager.Tasklist(10, o.debug, singleprocess=o.debug)
	
	# match SNPs to genes
	snplist_folder = mkdir(o.mainfolder, 'SNP_list')
	snplistfile = os.path.join(snplist_folder, 'mpileup_snplist.txt')
	snp2gene_file = os.path.join(snplist_folder, 'SNPs_per_gene.txt')
	if o.prepareSNPlist:
		
		if not (o.dontredo and os.path.exists(snplistfile) and os.path.exists(snp2gene_file)):
			tasklist.add(make_snp2gene_file, (o.genePred, o.snpdatabase, snplistfile, snp2gene_file), group='make_snp2gene_file')
			log('annotation', os.path.abspath(o.genePred), os.path.abspath(o.snpdatabase))
	
	# run samtools mpileup
	mpileup_files = files_in_folder('mpileup', o, '.mpileup.gz')
	if o.mpileup:
		outdir = mkdir(o.mainfolder, 'mpileup')
		for bamfolderpath in o.bamfolders:
			if o.folder_scheme == 'sample/sample_unique.bam':
				samplename = bamfolderpath.rstrip('/').rsplit('/',1)[-1]
			elif o.folder_scheme == 'sample/rnastar/sample.bam': samplename = bamfolderpath.rstrip('/').rsplit('/',2)[-2]
			elif o.folder_scheme == 'sample.bam': samplename = bamfolderpath.rsplit('.',1)[0]
			bamfile = find_bamfile(bamfolderpath, samplename)
			if not os.path.exists(bamfile):
				print 'Missing bam for', samplename
				continue
			outputpath = os.path.join(outdir, samplename+'.mpileup.gz')
			mpileup_files.add(outputpath)
			if not (o.dontredo and os.path.exists(outputpath)):
				tasklist.add(run_mpileup, (bamfile, snplistfile, outputpath), sample=samplename, group='mpileup', waitfor=['make_snp2gene_file'])
	
	# list the mpileup result as counts tied to genes
	hitcounts_raw_files = files_in_folder('hitcounts', o, '.counts.gz')
	if o.counts:
		outdir = mkdir(o.mainfolder, 'hitcounts')
		for mpileupfile in mpileup_files:
			samplename = mpileupfile.split('/')[-1].split('.mpileup')[0]
			outputpath = os.path.join(outdir, samplename+'.counts.gz')
			hitcounts_raw_files.add(outputpath)
			if not (o.dontredo and os.path.exists(outputpath)):
				tasklist.add(count_at_SNPs, (mpileupfile, snp2gene_file, outputpath), sample=samplename, group='counts', waitfor=['make_snp2gene_file', 'mpileup'])
	
	# learn which SNPs have both genotypes in the dataset
	heterozygous_snp2gene_file = os.path.join(snplist_folder, 'heterozygous_SNPs_per_gene.txt')
	if o.list_heterozygous:
		heterozygous_snp2gene_file_specificname = heterozygous_snp2gene_file
		if not (o.dontredo and os.path.exists(heterozygous_snp2gene_file_specificname)):
			tasklist.add(find_heterozygous, (hitcounts_raw_files, snp2gene_file, heterozygous_snp2gene_file_specificname, o), group='list_heterozygous', waitfor=['counts'])
			log('minreadsH', str(o.minreadsH))
			log('maxratioH', str(o.maxratioH))
			log('samplesH', repr(o.samplesH))
			log('minothersamplesH', str(o.minothersamplesH))
			log('minothersamplereadsH', str(o.minothersamplereadsH))		
	
	hitcounts_H_files = files_in_folder('hitcounts_heterozygous', o, '.counts.gz')
	if o.counts_heterozygous:
		outdir = mkdir(o.mainfolder, 'hitcounts_heterozygous')
		for mpileupfile in mpileup_files:
			samplename = mpileupfile.split('/')[-1].split('.mpileup')[0]
			outputpath = os.path.join(outdir, samplename+'.counts.gz')
			hitcounts_H_files.add(outputpath)
			if not (o.dontredo and os.path.exists(outputpath)):
				tasklist.add(count_at_SNPs, (mpileupfile, heterozygous_snp2gene_file, outputpath), sample=samplename, group='countsH', waitfor=['make_snp2gene_file', 'mpileup', 'list_heterozygous'])
	
	if o.fake_mouse_allelehits_format:
		outdir = mkdir(o.mainfolder, 'allelehits_per_SNP')
		allelehits_files = []
		for counts_file in hitcounts_H_files:
			samplename = counts_file.split('/')[-1].split('.counts')[0]
			outputpath = os.path.join(outdir, samplename+'.fake_mouse_allelehits.txt')
			allelehits_files.append(outputpath)
			if not (o.dontredo and os.path.exists(outputpath)):
				tasklist.add(counts_to_SNPallelehits, (counts_file, samplename, heterozygous_snp2gene_file, outputpath), sample=samplename, group='allelehits_1sample', waitfor=['countsH', 'list_heterozygous'])
		outputpath = os.path.join(outdir, 'merged.fake_mouse_allelehits.txt')
		tasklist.add(rpkmfile_merge, (allelehits_files, outputpath), group='allelehits_merge', waitfor=['allelehits_1sample'])
	
	tasklist.waitforall()
	
# example command:
#danielr@dna:~/Xandclones_late2014/Tcell$ python make_allelecalls.py -m YFVNM --samplesH YFVNM_Exome_Seq_Cap --minreadsH 3 --maxratioH 50 --minothersamplereadsH 10 -LHA

