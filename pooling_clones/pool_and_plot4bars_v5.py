from __future__ import division
import argparse, dr_tools, numpy, pylab, random

def count_mono(samples, randomize_per_gene=False, allowedgenes=None, use_ti=None):
	c57mono, castmono, bi = 0,0,0
	mono_genes = set(), set()
	ti_used = set()
	for ti, sym in enumerate(expra['symbols']):
		if use_ti is not None:
			if ti not in use_ti: continue
		else:
			if disallowedgenes and sym in disallowedgenes: continue
			if allowedgenes and sym not in allowedgenes: continue
		ai = ti
		
		if o.allelerand_skew:
			c57fraction = allelerand_skew[ai]
		else:
			c57fraction = 0.5
		c57 = 0
		cast = 0
		both = 0
		for s in samples:
			a1 = bool(expra[s+'_c57only'][ai])
			a2 = bool(expra[s+'_castonly'][ai])
			if a1 and a2: both += 1
			elif a1 and not a2:
				if randomize_per_gene and random.random() < c57fraction:
					cast += 1
				else:
					c57 += 1
			elif a2 and not a1:
				if randomize_per_gene and random.random() > c57fraction:
					c57 += 1
				else:
					cast += 1
		
		tot_cells = max(1, c57+cast+both) # don't count unexpressing cells
		if c57/tot_cells >= o.min_cell_fraction:
			c57mono += 1
			mono_genes[0].add(sym)
		elif cast/tot_cells >= o.min_cell_fraction:
			castmono += 1
			mono_genes[1].add(sym)
		elif both: bi += 1
		
		ti_used.add(ti)
	tot_genes = c57mono + castmono + bi
	if tot_genes == 0 or o.genecount: tot_genes = 1
	return c57mono/tot_genes, castmono/tot_genes, bi/tot_genes, mono_genes, ti_used

class Bars:
	def __init__(self):
		self.xarr = []
		self.x_pos = 0
		self.mono = []
		self.cast = []
		self.tot = []
		self.labels = []
		self.labels_x = []
		self._reset_points()
		self.points_draw = [] # list of lists
		self.error_bar = []
	
	def _reset_points(self):
		self.points_cast = []
		self.points_mono = []
		self.points_tot = []
	
	def _label_and_advance(self, label, errorbar):
		self.xarr.append(self.x_pos)
		self.points_draw.append(self.points_mono)
		if label:
			self.labels.append(label)
			self.labels_x.append(self.x_pos)
		self.error_bar.append(errorbar)
		#print self.x_pos, label, c57mono, castmono
		self.x_pos += 1
		self._reset_points()
	
	def addpoint(self, c57mono, castmono, bi=None, mono_genes=None, ti_used=None):
		self.points_cast.append(castmono)
		self.points_mono.append(castmono+c57mono)
		self.points_tot.append(castmono+c57mono+bi)
	
	def formbar(self, label='', errorbar=False):
		summary_func = numpy.mean
		self.cast.append(summary_func(self.points_cast))
		self.mono.append(summary_func(self.points_mono))
		self.tot.append(summary_func(self.points_tot) if o.genecount else 1)
		if o.vocal: print self.mono[-1], self.points_mono, label
		self._label_and_advance(label, errorbar)
	
	def addbar(self, label, c57mono, castmono, bi=None, mono_genes=None, ti_used=None):
		self.cast.append(castmono)
		self.mono.append(castmono+c57mono)
		self.tot.append(c57mono+castmono+bi if o.genecount else 1)
		if o.vocal: print self.mono[-1], label
		self._label_and_advance(label, False)

def ratio(expra, gi, samples):
	num_c57only = sum(expra[s+'_c57only'][gi] > 0 and expra[s+'_castonly'][gi] == 0 for s in samples)
	num_castonly = sum(expra[s+'_c57only'][gi] == 0 and expra[s+'_castonly'][gi] > 0 for s in samples)
	num_both = num_c57only+num_castonly
	if num_both == 0: return 0.5
	return num_c57only/num_both

def percentile(Y, percent):
	try:
		return numpy.percentile(numpy.array(Y), percent)
	except IndexError:
		print len(Y)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-o', '--figure', required=True)
	parser.add_argument('-m', '--minrpkm', default=20, type=float)
	parser.add_argument('-M', '--maxrpkm', type=float)
	parser.add_argument('-c', '--clonal_groups', default='clonal_groups.txt')
	parser.add_argument('-gA', '--allowedgenesA', required=True)
	parser.add_argument('-gX', '--allowedgenesX')
	parser.add_argument('-ge', '--disallowedgenes', nargs='+')
	parser.add_argument('--addgreen', action='store_true')
	parser.add_argument('--genecount', action='store_true')
	parser.add_argument('--random_seed', type=int)
	parser.add_argument('-R', '--random_bars', type=int, default=1)
	parser.add_argument('-f', '--min_cell_fraction', type=float, default=1.0)
	parser.add_argument('--cellXbar', action='store_true')
	parser.add_argument('--sharedXbar', action='store_true')
	parser.add_argument('-I', '--individualcutoff', action='store_true')
	parser.add_argument('--vocal', action='store_true')
	parser.add_argument('-r', '--allelerand_skew', action='store_true')
	o = parser.parse_args()
	
	if o.random_seed is not None: random.seed(o.random_seed)
	
	allowedgenesA = set(dr_tools.loadlist(o.allowedgenesA)) if o.allowedgenesA else None
	allowedgenesX = set(dr_tools.loadlist(o.allowedgenesX)) if o.allowedgenesX else None
	if o.disallowedgenes:
		disallowedgenes = set()
		for filename in o.disallowedgenes:
			disallowedgenes.update(set(dr_tools.loadlist(filename)))
	else:
		disallowedgenes = None
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	samples_all = [s.rsplit('_',1)[0] for s in expra.samples[::2]]
	
	bars = Bars()
	for clonal_group in dr_tools.loadlist(o.clonal_groups):
		samples = [s for s in samples_all if any(s.startswith(clonal_group_start) or s.startswith('pool.'+clonal_group_start) for clonal_group_start in clonal_group.split('\t'))]
		
		samples_forskew = [s for s in samples_all if s not in samples]
		allelerand_skew =  dict((gi, ratio(expra, gi, samples_forskew)) for gi in range(len(expra['symbols'])))
		
		merged_output = count_mono(samples, False, allowedgenesA, None)
		use_ti = merged_output[-1]
		for sample in samples:
			bars.addpoint(*count_mono([sample], False, allowedgenesA, None if o.individualcutoff else use_ti))
		bars.formbar(clonal_group)
		bars.addbar('shared', *merged_output)
		
		print clonal_group.split('\t')[0], 'shared genes, c57:', ', '.join(list(merged_output[3][0]))
		print clonal_group.split('\t')[0], 'shared genes, cast:', ', '.join(list(merged_output[3][1]))
		
		for ri in range(o.random_bars):
			bars.addpoint(*count_mono(samples, True, allowedgenesA, use_ti))
		bars.formbar('random shared', True)
		
		if o.cellXbar or o.sharedXbar:
			sharedX_output = count_mono(samples, False, allowedgenesX, None)
		if o.cellXbar:
			use_ti = sharedX_output[-1]
			for sample in samples:
				bars.addpoint(*count_mono([sample], False, allowedgenesX, None if o.individualcutoff else use_ti))
			bars.formbar('X')
		if o.sharedXbar:
			bars.addbar('X shared', *sharedX_output)
		bars.x_pos += 1
	
	
	
	if o.addgreen: pylab.bar(bars.xarr, bars.tot, facecolor='#009933', linewidth=0, width=0.9)
	pylab.bar(bars.xarr, bars.mono, facecolor='#888800', linewidth=0, width=0.9)
	pylab.bar(bars.xarr, bars.cast, facecolor='#aa8800', linewidth=0, width=0.9)
	
	dots = [(x,y) for x,Y,errorbar in zip(bars.xarr, bars.points_draw, bars.error_bar) for y in Y if not errorbar]
	pylab.plot([x+0.45+(1-2*random.random())*0.2 for x,y in dots], [y for x,y in dots], 'k,')
	errorbar_mid = [numpy.median(Y) for x,Y,errorbar in zip(bars.xarr, bars.points_draw, bars.error_bar) if errorbar and Y]
	errorbar_updist = [percentile(Y, 97.5)-numpy.median(Y) for x,Y,errorbar in zip(bars.xarr, bars.points_draw, bars.error_bar) if errorbar and Y]
	errorbar_downdist = [numpy.median(Y)-percentile(Y, 2.5) for x,Y,errorbar in zip(bars.xarr, bars.points_draw, bars.error_bar) if errorbar and Y]
	errorbar_x = [x+0.45 for x,Y,errorbar in zip(bars.xarr, bars.points_draw, bars.error_bar) if errorbar and Y]
	pylab.errorbar(errorbar_x, errorbar_mid, yerr=[errorbar_downdist, errorbar_updist], fmt='', ecolor='k', linewidth=0, elinewidth=1)
	
	pylab.xticks(bars.labels_x, bars.labels, rotation=90)
	pylab.subplots_adjust(bottom=0.3)
	pylab.savefig(o.figure)
