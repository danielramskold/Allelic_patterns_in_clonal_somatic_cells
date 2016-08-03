import subprocess, argparse, os, taskmanager

def run(fastqs, foldername):
	outf = foldername+'.fastq.gz'
	if os.path.exists(outf): return
	with open(outf, 'w') as outfh:
		if all(f.endswith('.gz') for f in fastqs):
			subprocess.check_call(['cat'] + fastqs, stdout=outfh)
		elif not any(f.endswith('.gz') for f in fastqs):
			subprocess.check_call(['gzip', '-c'] + fastqs, stdout=outfh)
			print ' '.join(['gzip', '-c'] + fastqs) + '>'+outf
		else: raise Exception

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('folder_under_rawdata', nargs='+')
	o = parser.parse_args()
	
	tasks = taskmanager.Tasklist(10)
	
	for folder in o.folder_under_rawdata:
		fastqs = [os.path.join(folder, n) for n in os.listdir(folder) if not n.endswith('trimming_report.txt') and not os.path.isdir(os.path.join(folder, n))]
		if len(fastqs) == 0: raise Exception
		print fastqs
		foldername = folder.rstrip('/').split('/')[-1]
		tasks.add(run, (fastqs, foldername))
