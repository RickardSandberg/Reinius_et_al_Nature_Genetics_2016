from __future__ import division
import argparse, dr_tools, numpy, random, pylab
from scipy import stats
from collections import *

def remove_digits(text):
	text = ''.join([character for character in text if not character.isdigit()])
	#while text and text[0].isdigit():
	#	text = text[1:]
	return text

class Classification:
	def __init__(self):
		self.nondet = set()
		self.allele1 = set() # cast monoallelic
		self.allele2 = set() # c57 monoallelic
		self.biallelic = set()
		self.groups = [self.nondet, self.allele1, self.allele2, self.biallelic]

def classify_genes(expra, sample, ai_list):
	classification = Classification()
	
	for aii, ai in enumerate(ai_list):
		e1 = expra[sample+S1][ai]
		e2 = expra[sample+S2][ai]
		if e1 == 0:
			if e2 == 0:
				classification.nondet.add(aii)
			else:
				classification.allele2.add(aii)
		else:
			if e2 == 0:
				classification.allele1.add(aii)
			else:
				classification.biallelic.add(aii)
	return classification

def model_true2(A, B, D):
	# A = (p,p)obs
	# B = (p,m)obs
	# D = (p,b)obs
	# x = p true
	# y = b true
	# z = chance of allelic loss
	z = B/(B+D)
	y = B/z**2/(1-z)**2
	x = (A-B)/(1-z)**2
	return x, y, z

def overlap_types(classification_s1, classification_s2):
	S = sum(len(g) for g in classification_s1.groups)
	
	A_pp = numpy.mean([len(classification_s1.allele1 & classification_s2.allele1)/S, len(classification_s1.allele2 & classification_s2.allele2)/S, len(classification_s2.allele1 & classification_s1.allele1)/S, len(classification_s2.allele2 & classification_s1.allele2)/S])
	B_pm = numpy.mean([len(classification_s1.allele1 & classification_s2.allele2)/S, len(classification_s1.allele2 & classification_s2.allele1)/S, len(classification_s2.allele1 & classification_s1.allele2)/S, len(classification_s2.allele2 & classification_s1.allele1)/S])
	C_p0 = numpy.mean([len(classification_s1.allele1 & classification_s2.nondet)/S, len(classification_s1.allele2 & classification_s2.nondet)/S, len(classification_s2.allele1 & classification_s1.nondet)/S, len(classification_s2.allele2 & classification_s1.nondet)/S])
	p_obs = numpy.mean([len(classification_s1.allele1)/S, len(classification_s1.allele2)/S, len(classification_s2.allele1)/S, len(classification_s2.allele2)/S])
	D_pb = numpy.mean([len(classification_s1.allele1 & classification_s2.biallelic)/S, len(classification_s1.allele2 & classification_s2.biallelic)/S, len(classification_s2.allele1 & classification_s1.biallelic)/S, len(classification_s2.allele2 & classification_s1.biallelic)/S])
	
	x_ptrue,y_btrue,z = model_true2(A_pp, B_pm, D_pb)
	
	return (S, z, x_ptrue/p_obs, x_ptrue*2/(y_btrue+x_ptrue*2), S*(2*x_ptrue+y_btrue))

def get_rpkm(Ai, sample):
	ID = expra['IDs'][Ai]
	if ID == '1 0': return -1
	if not ID in exprt.ID_to_index: return -1000
	Ti = exprt.ID_to_index[ID]
	return exprt[sample][Ti]

def extract_short_name1(name):
	# give e.g. '9Q'
	return name.split('_7')[0].split('_s')[-1]

def extract_short_name2(name):
	# from e.g Bc39_splitL1_CxB give 'L1'
	if not 'split' in name: return ''
	return name.split('_split')[1].split('_')[0]

class Bars:
	def __init__(self):
		self.midv = []
		self.err_up = []
		self.err_down = []
		self.labels = []
		self.x = []

	def add(self, label, min_max_med, index=0):
		print label, min_max_med
		try:
			midv = min_max_med[3][index]
			self.err_down.append(midv-min_max_med[0][index])
			self.err_up.append(min_max_med[1][index]-midv)
		except:
			midv = min_max_med[3]
			self.err_down.append(midv-min_max_med[0])
			self.err_up.append(min_max_med[1]-midv)
		self.midv.append(midv)
		self.labels.append(label)
	
	def plot(self, errorbars):
		xarr = range(len(self.midv))
		if errorbars:
			pylab.bar(xarr, self.midv, yerr=[self.err_down, self.err_down], facecolor='#888888', linewidth=0, ecolor='#000000')
		else:
			pylab.bar(xarr, self.midv, facecolor='#888888', linewidth=0, ecolor='#000000')
		pylab.xticks(xarr, self.labels)


def mono_p(pair, use_these_ai, plotted_variable_index, o):
	gene_sum = 0
	monoallelic_sum = 0
	for min_rpkm, max_rpkm in zip(o.expression_cutoffs, o.expression_cutoffs[1:]+[1e20]):
		if not o.noncumulative: max_rpkm = 1e20
		aiL_bin = [i for i in use_these_ai if max_rpkm > numpy.mean([get_rpkm(i, s) for s in samples_in_pairs]) >= min_rpkm]
		if len(aiL_bin) == 0: continue
		func = lambda aiL : overlap_types(classify_genes(expra, pair[0], aiL), classify_genes(expra, pair[1], aiL))
		values_out = func(aiL_bin)
		
		mono_addition = values_out[4] * values_out[plotted_variable_index]
		if not numpy.isnan(mono_addition):
			monoallelic_sum += mono_addition
			gene_sum += values_out[4]
	return (monoallelic_sum/gene_sum,)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits_file', nargs='+', required=True)
	parser.add_argument('-r', '--rpkm_file', nargs='+')
	parser.add_argument('-m', '--minrpkm', type=float, default=0.001)
	parser.add_argument('-M', '--maxrpkm', type=float, default=1e20)
	parser.add_argument('-p', '--pairletters', nargs='+')
	parser.add_argument('-P', '--permutations', type=int, default=1000)
	parser.add_argument('-c', '--expression_cutoffs', nargs='+', type=float)
	parser.add_argument('-nc', '--noncumulative', action='store_true')
	parser.add_argument('-e', '--errorbars', action='store_true')
	parser.add_argument('-b', '--barplot', action='store_true')
	parser.add_argument('-o', '--figure', default='pair_overlap6.pdf')
	parser.add_argument('-v', '--plotted_variable', default='mono%', choices=['mono%', 'z', 'num_genes', 'error', 'info_genes'])
	parser.add_argument('--ylim', default=[0,1], type=float, nargs=2)
	parser.add_argument('-s', '--shiftpairs', type=int, nargs='?', const=1)
	parser.add_argument('-F', '--fibroblastnames', action='store_true')
	parser.add_argument('-S', '--separatelines', action='store_true')
	o = parser.parse_args()
	
	plotted_variable_index = {'mono%':3, 'z': 1, 'num_genes': 0, 'error': 2, 'info_genes':4}[o.plotted_variable]
	
	# suffixes of sample names in expression file
	S2 = '_c57only'
	S1 = '_castonly'
	
	# load
	expra = dr_tools.loadexpr(o.allelehits_file, True)
	samples = [s.split(S2)[0] for s in expra.samples[::2]]
	
	if o.minrpkm:
		exprt = dr_tools.loadexpr(o.rpkm_file, False)
		samples = [s for s in samples if s in exprt.samples]
	
	global done_c
	done_c = dict()
	
	# pairs end in the same capital letter
	# skip the _wronglane samples and the non-split cells
	extract_short_name = extract_short_name2 if o.fibroblastnames else extract_short_name1
	pair_letters = list(set(remove_digits(extract_short_name(name)) for name in samples) - set(['']))
	if o.pairletters: pair_letters = [l for l in pair_letters if l.strip('-_') in o.pairletters]
	print pair_letters
	sample_pairs = [[s for s in samples if remove_digits(extract_short_name(s)) == letter and not '_wronglane' in s] for letter in pair_letters]
	samples_in_pairs = [s for pair in sample_pairs for s in pair]
	print sample_pairs
	
	if o.shiftpairs is not None:
		# move the second in each pair one step
		shift_length = o.shiftpairs
		sample_pairs = [[sample_pairs[i][0], sample_pairs[(i+shift_length)%len(sample_pairs)][1]] for i in range(len(sample_pairs))]
	
	# which genes to use
	if o.minrpkm:
		allowed_ai = set(ai for ai in range(len(expra['symbols'])) if o.maxrpkm > numpy.mean([get_rpkm(ai, s) for s in samples_in_pairs]) >= o.minrpkm)
	else:
		allowed_ai = None
	
	use_these_ai = [i for i in range(len(expra['symbols'])) if allowed_ai is None or i in allowed_ai]
	
	if o.barplot:
		# to do: use the expression cutoffs to calculate per-sample monoallelic fraction
		bars = Bars()
	
		for letter,pair in zip(pair_letters, sample_pairs):
			func = lambda aiL: mono_p(pair, aiL, plotted_variable_index, o)
			bootstrapped_minmax = dr_tools.bootstrap(func, [use_these_ai], controls=o.permutations, confidence=0.95, nullval=0, processes=10, give_median=True)
			bars.add(letter, bootstrapped_minmax[1:], 0)
	
		func = lambda aiL : (numpy.median([mono_p(pair, aiL, plotted_variable_index, o)[0] for pair in sample_pairs]),)
		bars.add('median', dr_tools.bootstrap(func, [use_these_ai], controls=o.permutations, confidence=0.95, nullval=0, processes=10, give_median=True)[1:], 0)
	
		bars.plot(o.errorbars)
		pylab.savefig(o.figure)
	elif o.separatelines:
		bars = defaultdict(Bars)
		gene_sum = defaultdict(int)
		monoallelic_sum = defaultdict(int)
		for min_rpkm, max_rpkm in zip(o.expression_cutoffs, o.expression_cutoffs[1:]+[1e20]):
			if not o.noncumulative: max_rpkm = 1e20
			aiL_prerand = [i for i in use_these_ai if max_rpkm > numpy.mean([get_rpkm(i, s) for s in samples_in_pairs]) >= min_rpkm]
			if len(aiL_prerand) == 0: continue
			for letter, pair in zip(pair_letters, sample_pairs):
				func = lambda aiL : overlap_types(classify_genes(expra, pair[0], aiL), classify_genes(expra, pair[1], aiL))
				bootstrap_out = dr_tools.bootstrap(func, [aiL_prerand], controls=o.permutations, confidence=0.95, nullval=0, processes=10, give_median=True)[1:]
				bars[letter].add(str(min_rpkm), bootstrap_out, plotted_variable_index)
				bars[letter].x.append(min_rpkm)
			
				mono_addition = bootstrap_out[3][4] * bootstrap_out[3][3]
				if not numpy.isnan(mono_addition):
					monoallelic_sum[letter] += mono_addition
					gene_sum[letter] += bootstrap_out[3][4]
		
		for letter, pair in zip(pair_letters, sample_pairs):
			print letter, monoallelic_sum[letter], monoallelic_sum[letter]/gene_sum[letter]
			pylab.semilogx(bars[letter].x, bars[letter].midv)
		
		pylab.xlim(0.05, 1000)
		pylab.ylim(*o.ylim)
		if o.noncumulative:
			pylab.xlabel('rpkm')
		else:
			pylab.xlabel('min rpkm')
		pylab.ylabel(o.plotted_variable)
		pylab.savefig(o.figure)
	else:
		bars = Bars()
		
		gene_sum = 0
		monoallelic_sum = 0
		
		for min_rpkm, max_rpkm in zip(o.expression_cutoffs, o.expression_cutoffs[1:]+[1e20]):
			if not o.noncumulative: max_rpkm = 1e20
			aiL_prerand = [i for i in use_these_ai if max_rpkm > numpy.mean([get_rpkm(i, s) for s in samples_in_pairs]) >= min_rpkm]
			if len(aiL_prerand) == 0: continue
			func = lambda aiL : [numpy.median(V) for V in zip(*(overlap_types(classify_genes(expra, pair[0], aiL), classify_genes(expra, pair[1], aiL)) for pair in sample_pairs))]
			bootstrap_out = dr_tools.bootstrap(func, [aiL_prerand], controls=o.permutations, confidence=0.95, nullval=0, processes=10, give_median=True)[1:]
			bars.add(str(min_rpkm), bootstrap_out, plotted_variable_index)
			bars.x.append(min_rpkm)
			
			mono_addition = bootstrap_out[3][4] * bootstrap_out[3][3]
			if not numpy.isnan(mono_addition):
				monoallelic_sum += mono_addition
				gene_sum += bootstrap_out[3][4]
		
		print monoallelic_sum, monoallelic_sum/gene_sum
		
		if o.errorbars:
			pylab.gca().set_xscale('log')
			pylab.errorbar(bars.x, bars.midv, yerr=[bars.err_down, bars.err_down])
		else:
			pylab.semilogx(bars.x, bars.midv)
		
		pylab.xlim(0.05, 1000)
		pylab.ylim(*o.ylim)
		if o.noncumulative:
			pylab.xlabel('rpkm')
		else:
			pylab.xlabel('min rpkm')
		pylab.ylabel(o.plotted_variable)
		pylab.savefig(o.figure)
