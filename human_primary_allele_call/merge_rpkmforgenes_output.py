from __future__ import with_statement
import dr_tools, argparse, time, sys
from itertools import izip_longest

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('outfile')
	opts.add_argument('infile', nargs='+')
	opts.add_argument('--addsuffixes', nargs='+', default=[])
	o = opts.parse_args()
	
	if len(o.addsuffixes) > len(o.infile): raise Exception
	
	# parse files
	numgenes = None
	header = ['#samples', '#allmappedreads', '#normalizationreads', '#arguments']
	header[3] += '\t' + ' '.join(sys.argv) + '\ttime: ' + time.asctime()
	genelines = []
	readlines = []
	for inf, suffix in izip_longest(o.infile, o.addsuffixes, fillvalue=''):
		lnumgenes = 0
		with open(inf, 'r') as infh:
			for line in infh:
				p = dr_tools.split(line)
				if p[0] == '#samples':
					header[0] += '\t' + '\t'.join([s+suffix for s in p[1:]])
					lnumsamples = len(p) - 1
				elif p[0] == '#allmappedreads':
					header[1] += '\t' + '\t'.join(p[1:])
				elif p[0] in ('#genemappedreads', '#normalizationreads'):
					header[2] += '\t' + '\t'.join(p[1:])
				elif line[0] != '#':
					if numgenes is None:
						genelines.append('\t'.join(p[:2]))
						readlines.append('')
					else:
						if not genelines[lnumgenes].split('\t')[:2] == p[:2]:
							print 'gene line %d:'%lnumgenes, genelines[lnumgenes].split('\t')[:2], '!=', p[:2]
							raise Exception
					genelines[lnumgenes] += '\t' + '\t'.join(p[2:2+lnumsamples])
					readlines[lnumgenes] += ''.join('\t'+s for s in p[2+lnumsamples:2+2*lnumsamples])
					lnumgenes += 1
					
			if numgenes is None:
				numgenes = lnumgenes
			elif numgenes != lnumgenes:
				raise Exception, 'Unequal number of genes between files'
			
	# write output
	with open(o.outfile, 'w') as outfh:
		for line in header:
			print >>outfh, line
		for linenum in xrange(numgenes):
			print >>outfh, genelines[linenum] + readlines[linenum]
			
