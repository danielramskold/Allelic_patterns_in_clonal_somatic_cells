from __future__ import division
import dr_tools

def loadfile(infile):
	reads_per_sample = dict()
	for p in dr_tools.splitlines(infile):
		sample = p[0].split('/')[-1].split('_refseq.txt:#')[0]
		reads = float(p[1])
		reads_per_sample[sample] = reads
	return reads_per_sample

if '__main__' == __name__:
	infile_allmapped = 'posttrim_allmapped.txt'
	infile_genemapped = 'posttrim_genemapped.txt'
	outfile_oksamples = 'ok_mapping_samples.txt'
	allmapped = loadfile(infile_allmapped)
	genemapped = loadfile(infile_genemapped)
	failed_count = 0
	passed_count = 0
	with open(outfile_oksamples, 'w') as outfh:
		for sample in sorted(allmapped.keys()):
			if allmapped[sample] < 10000 or genemapped[sample]/allmapped[sample] < 0.4:
				# sample failed
				failed_count += 1
			else:
				print >>outfh, sample
				passed_count += 1
	print 'passed:', passed_count, 'failed:', failed_count
