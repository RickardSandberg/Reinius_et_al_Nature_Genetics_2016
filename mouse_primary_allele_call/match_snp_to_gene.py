import pybedtools, sys, os, argparse
from operator import itemgetter


parser = argparse.ArgumentParser(description='match SNP to genes')
parser.add_argument('--snpfile', help='over ride default file locations')
parser.add_argument('--force', action='store_true')
parser.add_argument('--annotation', default='bed') # alternative gff for ensembl
parser.add_argument('--snpname', default='tmp')

args = parser.parse_args()

if args.annotation == 'gff':
    annotation_file = 'annotations/Mm9_ensembl.gff'
    annotation_type = 'gff'
    annotationname = 'ensembl'
else:
    annotation_file = 'annotations/mm9_refGene_31Jul2011_norandom.bed' #                        
    annotation_type = 'bed'
    annotationname = 'refseq'


summary_file = args.snpfile # conf[args.snptype]['summary']
genesout_file = '%s_mm9_%s_snp2genes.txt'%(args.snpname,  annotationname)
txout_file = '%s_mm9_%s_snp2transcripts.txt'%(args.snpname, annotationname)

a = pybedtools.BedTool(summary_file)
print 'summary file (%s) interpreted as %s' % (summary_file, a.file_type)

b = pybedtools.BedTool(annotation_file)
print 'annotation file (%s) interpreted as %s' % (annotation_file, b.file_type)

print 'writing output to:', genesout_file, txout_file

if 0 and annotation_type == 'bed':
    with_counts = b.intersect(a, loj=True)#.saveas('oobrv2_mm9_intersection.txt')                                                         
else:
    with_counts = a.intersect(b, loj=True).saveas('newintersection.txt')

# iterate over results and summarize per gene
gene2snp = {}
transcript2snp = {}
fg = open(genesout_file,'w')
ft = open(txout_file,'w')
for items in with_counts:
    if annotation_type == 'bed':

        if not ':' in items[-3]: continue
        txname, name = items[-3].split(':')                     
        if not name in gene2snp:
            gene2snp[name] = set([])
            transcript2snp[name] = {}
        gene2snp[name].add((items[0], items[1]))

        if not txname in transcript2snp[name]:
            transcript2snp[name][txname] = set([])
        transcript2snp[name][txname].add((items[0], items[1]))

    elif annotation_type == 'gff':
        try:   
            name = items[39].split('; gene_name "')[1].split('";')[0]
            #print name
            if not name in gene2snp:
                gene2snp[name] = set([])
                transcript2snp[name] = {}
            gene2snp[name].add((items[0], items[1]))

            txname = items[39].split('; transcript_id "')[1].split('";')[0]
            if not txname in transcript2snp[name]:
                transcript2snp[name][txname] = set([])
            transcript2snp[name][txname].add((items[0], items[1]))
        except IndexError:
            pass



for gene in gene2snp.iterkeys():
    fg.write("%s\t%i\t%s\n" % (gene, len(gene2snp[gene]), ";".join([":".join(map(str, [c,p])) for c,p in gene2snp[gene]])))

    for tx in transcript2snp[gene].iterkeys():
        transcript2snp[gene][tx] = sorted(transcript2snp[gene][tx], key=itemgetter(1))
        ft.write("%s\t%s\t%i\t%s\n" % (gene, tx, len(transcript2snp[gene][tx]), ";".join([":".join(map(str, [c,p])) for c,p in transcript2snp[gene][tx]])))

fg.close()
ft.close()
