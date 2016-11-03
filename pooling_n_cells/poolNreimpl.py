from __future__ import division
import argparse, pylab, dr_tools, random

def monoallelic_fraction_pooled(samples_nom, expra, allowed_genes, disallowed_genes, randomize_per_gene=False):
	monocount = 0
	bicount = 0
	samples = []
	for line in samples_nom:
		if ',' in line:
			while True:
				s = random.choice(line.split(','))
				if not s+'_c57only' in expra:
					print s
				else: break
				if s in multi_cell_samples:
					print 'multicellsample', s
			samples.append(s)
		else:
			samples.append(line)
	for ai, sym in enumerate(expra['symbols']):
		if allowed_genes is not None and sym not in allowed_genes: continue
		if sym in disallowed_genes: continue
		if randomize_per_gene:
			if o.allelerand_skew:
				identities = ['c57' if expra[s+'_c57only'][ai]>0 and expra[s+'_castonly'][ai]==0 else 'cast' if expra[s+'_c57only'][ai]==0 and expra[s+'_castonly'][ai]>0 else '0' for s in samples]
				nonswapped = [True if I=='0' else random.random() < allelerand_skew[ai] if I=='c57' else 1-allelerand_skew[ai] for I in identities]
			else:
				nonswapped = [random.random() < 0.5 for s in samples]
			c57 = sum((expra[s+'_c57only'][ai]>0) if nonswap else (expra[s+'_castonly'][ai]>0) for s,nonswap in zip(samples, nonswapped))
			cast = sum((expra[s+'_castonly'][ai]>0) if nonswapped else (expra[s+'_c57only'][ai]>0) for s,nonswap in zip(samples, nonswapped))
			#print c57, cast, sum(expra[s+'_c57only'][ai]>0 for s in samples)
		else:
			c57 = sum(expra[s+'_c57only'][ai]>0 for s in samples)
			cast = sum(expra[s+'_castonly'][ai]>0 for s in samples)
		if c57 > 0 and cast > 0:
			bicount += 1
		elif c57 > 0 or cast > 0:
			monocount += 1
	return monocount/(monocount+bicount)

multi_cell_samples = {
"MAF_CxB_clone_A_10": 2, # two cells attached together
"MAF_CxB_clone_B_19": 2, # two cells attached together
"MAF_CxB_clone_B_31": 2, # two cells attached together
"MAF_CxB_clone_B_41": 2, # two small cells attached together
"MAF_CxB_clone_B_53": 2, # two small cells attached together
"MAF_CxB_clone_B_61": 2, # two small cells attached together
"MAF_CxB_clone_B_65": 2, # two small cells attached together
"MAF_CxB_clone_B_69": 2, # two small cells attached together. overblown
"MAF_CxB_clone_B_75": 2, # two small cells attached together
"MAF_CxB_clone_B_76": 2, # two cells attached together
"MAF_CxB_clone_B_87": 3, # big clump. 3? cells
"MAF_CxB_clone_B_88": 2, # two cells attached together
"MAF_CxB_clone_B_91": 2, # two small cells attached together
"MAF_CxB_clone_B_100": 2, # two small cells attached together
"MEF_E14_Clone_BxC_H_10": 2, # 2 cells picked
"MEF_E14_Clone_BxC_H_12": 2, # 2 cells picked
"MEF_E14_Clone_BxC_H_14": 2, # 2 cells picked
"pool.MAF_CxB_clone_B_58": 2, # two flat cells together
"pool.MAF_CxB_clone_B_77": 3, # clump of flat cells 2-4
"pool.MAF_CxB_clone_B_86": 3, # big dark clump. 3? cells
"pool.MAF_CxB_clone_B_95": 15, # 15 small cells in one go
"pool.MAF_CxB_clone_B_101": 3, # 3-4 cells in one go
"pool.MAF_CxB_RE_clone_B_102": 5, # 5 cells in one go. one of them is big
"pool.MAF_CxB_RE_clone_B_103": 3, # 3? big cells attached
"pool.MAF_CxB_RE_clone_B_104": 5, # 5 cells. one is big
"pool.MEF_E14_Clone_BxC_H_4": 2, # 2 cells picked
}

def ratio(expra, gi, samples):
	num_c57only = sum(expra[s+'_c57only'][gi] > 0 and expra[s+'_castonly'][gi] == 0 for s in samples)
	num_castonly = sum(expra[s+'_c57only'][gi] == 0 and expra[s+'_castonly'][gi] > 0 for s in samples)
	num_both = num_c57only+num_castonly
	if num_both == 0: return 0.5
	return num_c57only/num_both

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('-ge', '--disallowedgenes', nargs='+', default=[])
	parser.add_argument('-R', '--random_dots', type=int, default=1)
	parser.add_argument('-s', '--samplelist', required=True)
	parser.add_argument('-n', default=[5,15], type=int, nargs='+')
	parser.add_argument('-o', '--figure', default='poolN.pdf')
	parser.add_argument('--skip_multi', action='store_true')
	parser.add_argument('-r', '--allelerand_skew', action='store_true')
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	if o.allowedgenes:
		allowed_genes = set(dr_tools.loadlist(o.allowedgenes))
	else:
		allowed_genes = None
	disallowed_genes = set()
	for filename in o.disallowedgenes:
		disallowed_genes.update(set(dr_tools.loadlist(filename)))
	
	poolable_samples = dr_tools.loadlist(o.samplelist, ignore='#')
	
	random.seed(0)
	
	samples_n = dict((n_cells, [random.sample(poolable_samples, n_cells) for di in range(o.random_dots)]) for n_cells in o.n)
	
	samples_all = [sa.split('_c57only')[0] for sa in expra.samples[::2]]
	samples_all = [s for s in samples_all if s not in poolable_samples]
	allelerand_skew =  dict((gi, ratio(expra, gi, samples_all)) for gi in range(len(expra['symbols'])))
	
	# random sim
	boxplot_y = []
	boxplot_x = []
	for n_cells in o.n:
		if n_cells == 1:
			continue
		else:
			monofractions = []
			for di in range(o.random_dots):
				monofractions.append(monoallelic_fraction_pooled(samples_n[n_cells][di], expra, allowed_genes, disallowed_genes, True))
			
		boxplot_y.append(monofractions)
		boxplot_x.append(n_cells+0.5)
	pylab.boxplot(boxplot_y, positions=boxplot_x, whis=10000)
	
	# sim
	boxplot_y = []
	boxplot_x = []
	for n_cells in o.n:
		if n_cells == 1:
			monofractions = [monoallelic_fraction_pooled([sample], expra, allowed_genes, disallowed_genes) for sample in poolable_samples]
		else:
			monofractions = []
			for di in range(o.random_dots):
				monofractions.append(monoallelic_fraction_pooled(samples_n[n_cells][di], expra, allowed_genes, disallowed_genes, False))
			
		boxplot_y.append(monofractions)
		boxplot_x.append(n_cells)
	pylab.boxplot(boxplot_y, positions=boxplot_x, whis=10000)
	
	# pooled
	if not o.skip_multi:
		xarr,yarr = [], []
		for sample, n_cells in multi_cell_samples.items():
			if not n_cells in o.n: continue
			monofraction = monoallelic_fraction_pooled([sample], expra, allowed_genes, disallowed_genes, False)
			xarr.append(n_cells)
			yarr.append(monofraction)
		pylab.plot(xarr, yarr, 'ro')
	pylab.xlim(0,60)
	pylab.savefig(o.figure)
