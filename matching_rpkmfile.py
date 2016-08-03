import argparse, subprocess, os

'''
danielr@dna:~$ python Xandclones_late2014/Tcell/matching_rpkmfile.py -i ../sandberglab/pipeline3.0/rnaseq/hsa/jeff_tcells_clones_invitro/star_hg19/* -a /mnt/kauffman/danielr/Xandclones_late2014/Tcell/male_P1299_YFV2001/allelehits_per_SNP/merged.fake_mouse_allelehits.txt -o Xandclones_late2014/Tcell/rpkms_jeff_tcells_clones_invitro.txt
'''

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--bam_folder', nargs='+', required=True)
	parser.add_argument('-o', '--outfile', required=True)
	parser.add_argument('-a', '--allelehits', required=True)
	o = parser.parse_args()
	
	with open(o.allelehits, 'rU') as infh:
		sampleorder = [s.split('_cast')[0].split('_c57')[0] for s in next(infh).rstrip().split('\t')[1::2]]
	
	sample_to_bam = dict()
	for folderpath in o.bam_folder:
		sample = folderpath.rstrip('/').split('/')[-1]
		bampath = os.path.join(folderpath, sample+'_unique.bam')
		if not os.path.exists(bampath):
			raise Exception
		sample_to_bam[sample] = bampath
	
	sampleorder = [s for s in sampleorder if s in sample_to_bam]
	
	cmd = ['python', '/mnt/crick/sandberglab/src/rpkmforgenes.py', '-readcount', '-fulltranscript', '-mRNAnorm', '-rmnameoverlap', '-bothendceil', '-a', '/mnt/kauffman/sandberglab/pipeline3.0/additionaldata/rpkmforgenes/hg19refGene.txt', '-p', '20']
	cmd += ['-n'] + sampleorder
	cmd += ['-i'] + [sample_to_bam[s] for s in sampleorder]
	cmd += ['-o', o.outfile]
	
	subprocess.check_call(cmd)
