file1 = '/mnt/crick/rickards/projects/hsa_snp_calling/snp_stats_ac2.txt'
file2 = '/mnt/kauffman/danielr/Xandclones_late2014/Tcell/male_P1299_YFV2001_newsnpcall/SNP_list/SNPs_per_gene.txt' # created by make_allelecalls.py -P using the -a and -s arguments
output = '/mnt/kauffman/danielr/Xandclones_late2014/Tcell/male_P1299_YFV2001_newsnpcall/SNP_list/heterozygous_SNPs_per_gene.txt'


import dr_tools

positions = set()
for p in dr_tools.splitlines(file1): # for each SNP line in the file
	if float(p[-2]) < 0.9: # if second last column's value is <0.9
		positions.add('%s:%s'%(p[0], p[1])) # add to allowed SNP list

print len(positions)
c=0
outfh = open(output, 'w')
for p in dr_tools.splitlines(file2): # for each gene
	snps = []
	for snpinfo in p[2].split(';'): # go through the SNPs for the gene
		if snpinfo.split('|')[0] in positions: # see if on allowed list
			snps.append(snpinfo) # add to SNPs to print to output
			c+=1
	print >>outfh, dr_tools.join(p[0], len(snps), ';'.join(snps)) # output the SNPs for the gene
outfh.close()
print c
