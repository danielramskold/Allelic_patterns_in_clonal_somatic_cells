import dr_tools, os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('namefile')
parser.add_argument('--prefix_in', default='')
parser.add_argument('--prefix_out', default='')
parser.add_argument('--saysuccess', action='store_true')
o = parser.parse_args()

for p in dr_tools.splitlines(o.namefile):
	try:
		os.rename('%s%s.fastq.gz'%(o.prefix_in, p[0]), '%s%s.fastq.gz'%(o.prefix_out, p[1]))
	except:
		if not o.saysuccess:
			print 'fastq missing', p[0]
		pass
	else:
		if o.saysuccess:
			print 'renamed', p[0]
	
	try:
		os.rename('%s%s_expression.txt'%(o.prefix_in, p[0]), '%s%s_expression.txt'%(o.prefix_out, p[1]))
	except:
		if not o.saysuccess:
			print 'rpkms missing', p[0]
		pass
	else:
		if o.saysuccess:
			print 'renamed', p[0]
