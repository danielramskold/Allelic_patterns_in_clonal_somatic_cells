from __future__ import division
import argparse, dr_tools, numpy, random
from scipy import stats

class Gene:
	def __init__(self, ID, TSS):
		self.ID = ID
		self.TSS = TSS # genomic position on the chromosome
		self.i = None
		self.colours = []
		
	def add_state(self, maternal, paternal):
		if paternal == 0 and maternal == 0:
			c = '#ffffff'
		elif paternal == 0:
			c = '#ff0000'
		elif maternal == 0:
			c = '#0000ff'
		else:
			c = dr_tools.mixcolours(['#ff9900', '#0099ff'], [maternal/(maternal+paternal), paternal/(maternal+paternal)])
		self.colours.append(c)

def chop_name(samplename, stage):
	name = samplename.split('_c57o')[0].split('_casto')[0]
	name = name.rsplit('-',1)[0]
	if stage:
		name = name.rsplit('_',1)[0]
		if name.startswith('zy'): name = 'zy'
	return name

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('rpkmf_alleles')
	opts.add_argument('-c', '--chromosome', default='any')
	opts.add_argument('--filter', nargs='+')
	opts.add_argument('-f', '--figurefile', default='monoallelic_at_chr.png')
	opts.add_argument('-w', '--maxwhite', type=int)
	opts.add_argument('--allowallwhite', action='store_true')
	opts.add_argument('--allowedgenes')
	opts.add_argument('--disallowedgenes')
	opts.add_argument('--verticalborder', action='store_true')
	opts.add_argument('--stageline', action='store_true')
	opts.add_argument('--embryoline', action='store_true')
	opts.add_argument('--embryonotch', action='store_true')
	opts.add_argument('--mincoord', type=int)
	opts.add_argument('--maxcoord', type=int)
	opts.add_argument('--saygenes', action='store_true')
	args = opts.parse_args()
	
	# load expression data
	expr_alleles = dr_tools.loadexpr([args.rpkmf_alleles], counts=True)
	samples_alleles = sorted(e for e in expr_alleles if e not in ('IDs', 'symbols') and (args.filter is None or any(part in e for part in args.filter)))
	
	for p in dr_tools.splitlines(args.rpkmf_alleles):
		if p[0] == '#samples': samples = p[1:]; break
	samples_alleles = [e for e in samples if (args.filter is None or any(part in e for part in args.filter))]
	
	# sort the genes by position
	# only include transcripts which are the first ID in the entry of the rpkm file
	if args.allowedgenes is None and args.disallowedgenes is None:
		allowed_IDs = set(IDs.split('+')[0] for IDs in expr_alleles['IDs'])
	else:
		if args.allowedgenes:
			allowed_set = set(dr_tools.loadlist(args.allowedgenes))
		if args.disallowedgenes:
			disallowed_set = set(dr_tools.loadlist(args.disallowedgenes))
		allowed_IDs = set(IDs.split('+')[0] for IDs, symbols in zip(expr_alleles['IDs'],expr_alleles['symbols']) if (args.allowedgenes is None or any(identifier in allowed_set for identifier in (IDs.split('+') + symbols.split('+')))) and not (args.disallowedgenes is not None and any(identifier in disallowed_set for identifier in (IDs.split('+') + symbols.split('+')))))
	genes_per_chr = dict()
	ID_to_gene = dict()
	for ID in expr_alleles['IDs']:
		if ID not in allowed_IDs: continue
		chromosome = ID.split(':')[0]
		coord = int(ID.split(':')[1].split('|')[0])
		if args.chromosome == 'any':
			# place them together
			chromosome = 'any'
		elif chromosome != args.chromosome: continue
		if not chromosome in genes_per_chr: genes_per_chr[chromosome] = []
		if args.mincoord and coord < args.mincoord: continue
		if args.maxcoord and coord > args.maxcoord: continue
		genes_per_chr[chromosome].append(Gene(ID, coord))
		ID_to_gene[ID] = genes_per_chr[chromosome][-1]
	chromosome = args.chromosome
	try:
		genes_per_chr[chromosome].sort(key=lambda gene: gene.TSS)
	except:
		print genes_per_chr.keys()
		raise
	for gene_i, gene in enumerate(genes_per_chr[chromosome]):
		ID_to_gene[gene.ID] = gene
	
	# for drawing on the edge to indicate which cells are which
	edge_stage = [True]
	edge_embryo = [True]
	last_stage = None
	last_embryo = None
	
	for s1, s2 in zip(samples_alleles[::2], samples_alleles[1::2]):
		samplename = s1.rsplit('_',1)[0]
		
		# check that sample labels are consistent
		if samplename != s2.rsplit('_',1)[0]: raise Exception
		if not 'c57only' in s1 or not 'castonly' in s2: raise Exception
		
		# see if it's the same embryo, for colouring of the vertical edge
		if last_stage == chop_name(samplename, stage=True):
			edge_stage.append(edge_stage[-1])
		else:
			edge_stage.append(not edge_stage[-1])
			last_stage = chop_name(samplename, stage=True)
		if last_embryo == chop_name(samplename, stage=False):
			edge_embryo.append(edge_embryo[-1])
		else:
			edge_embryo.append(not edge_embryo[-1])
			last_embryo = chop_name(samplename, stage=False)
		
		# use expression values
		for ID, pat, mat in zip(expr_alleles['IDs'], expr_alleles[s1], expr_alleles[s2]):
			try:
				ID_to_gene[ID].add_state(mat, pat)
			except KeyError: pass
	
	edge_stage = edge_stage[1:]
	edge_embryo = edge_embryo[1:]
	
	# remove non-expressed (in all samples) and non-SNP-containing genes
	if args.maxwhite is not None:
		genes = [g for g in genes_per_chr[args.chromosome] if sum(c=='#ffffff' for c in g.colours) <= args.maxwhite and g.colours]
	elif args.allowallwhite:
		genes = [g for g in genes_per_chr[args.chromosome]]
	else:
		genes = [g for g in genes_per_chr[args.chromosome] if not all(c=='#ffffff' for c in g.colours)]
	
	
	last_stage_e = False
	last_embryo_e = False
	current_y = 0
	for y, (stage_e, embryo_e) in enumerate(zip(edge_stage, edge_embryo)):
		
		if args.stageline and stage_e != last_stage_e:
			current_y += 1
		elif args.embryoline and embryo_e != last_embryo_e:
			current_y += 1
		current_y += 1
		last_stage_e = stage_e
		last_embryo_e = embryo_e
	end_y = current_y

	# set up for drawing
	import Image, ImageDraw
	margin = 5 if args.verticalborder else 0
	extra_margin_right = 3 if args.embryonotch else 0
	end_x = len(genes)+margin*2+extra_margin_right
	im = Image.new('RGB', (end_x, end_y), '#ffffff')
	draw = ImageDraw.Draw(im)
	
	translate_y = dict() # much the y values are shifted to make room for horizontal lines
	print end_x, end_y, len(genes)
	
	# draw vertical borders, populate translate_y
	last_stage_e = False
	last_embryo_e = False
	current_y = 0
	for y, (stage_e, embryo_e) in enumerate(zip(edge_stage, edge_embryo)):
		
		if args.stageline and stage_e != last_stage_e:
			draw.line([(margin,current_y), (end_x-margin,current_y)], fill='#000000')
			current_y += 1
		elif args.embryoline and embryo_e != last_embryo_e:
			draw.line([(margin,current_y), (end_x-margin,current_y)], fill='#aaaaaa')
			current_y += 1
		if args.embryonotch and embryo_e != last_embryo_e or y==0:
			draw.line([(end_x-margin-1,current_y), (end_x-margin-extra_margin_right,current_y)], fill='#006600')
		
		translate_y[y] = current_y
		
		if args.verticalborder:
			for x in (0, end_x-1):
				draw.line([(x,current_y), (x,current_y)], fill='#000000' if stage_e else '#dddddd')
			for x in (1, end_x-2):
				draw.line([(x,current_y), (x,current_y)], fill='#000000' if embryo_e else '#dddddd')
		
		current_y += 1
		last_stage_e = stage_e
		last_embryo_e = embryo_e
	
	# draw the inside
	for x, gene in enumerate(genes):
		for y, c in enumerate(gene.colours):
			try:current_y = translate_y[y]
			except:
				print y, x, c, gene.ID
				break
			if y == 0 and args.saygenes:
				print x, gene.ID, expr_alleles['symbols'][expr_alleles.ID_to_index[gene.ID]]
			draw.line([(x+margin,current_y), (x+margin,current_y)], fill=c)
	with open(args.figurefile, 'wb') as fh:
		im.save(fh, "PNG")
	
