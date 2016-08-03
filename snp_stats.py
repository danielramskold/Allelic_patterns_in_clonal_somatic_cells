"""
python3 pipeline to compile and collects stats for CAST / C57 SNPs in RNA-Seq data

caveats:
- make sure that all data is sanger quality score format

"""

__version = 0.01
__author = 'sandberglab'


import os, sys, configparser, subprocess
from joblib import Parallel, delayed


# Command line argument handling
############################################################
from argparse import ArgumentParser

parser = ArgumentParser(description='RNA-Seq SNP analysis pipeline')

parser.add_argument('-c','--config', dest='configuration', default='samples.conf',
                    help='main configuration file')
parser.add_argument('-q',type=int, default=30)
parser.add_argument('-p','--proc', type=int, default=1)

args = parser.parse_args()

# Read Configuration File
##############################################################
conf = configparser.ConfigParser()
conf.read(args.configuration)
snps_file = conf.get("common","snps")
requiredQ = 30 # shoud move to conf file?

# HELP functions

def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        os.chmod(path, 0o774)


# 0 Collect information on experiments, samples and BAM files
###########################################################
genotypes = ('cast','c57','cast_c57','c57_cast')
allsamples = {}
for genotype in genotypes:
    allsamples[genotype] = []
    for datafolder in ('mpileups','cellsums'):
        if not os.path.exists(os.path.join(genotype, datafolder)):
            safe_mkdir(os.path.join(genotype, datafolder))

    for line in open(conf.get(genotype,'bamlist'), 'r'):
        p = line.strip().split('\t')
        allsamples[genotype].append(p)

# 1 Mpileup for all samples, stored as gzipped text file
# ######################################################
mpileup_commands = []
for genotype in genotypes:
    for item in allsamples[genotype]:
        try:
            sample, bamfile = item
        except:
            print ('error on line: %s in %s'%(item, genotype))
            sys.exit(0)

        mpileup_path = os.path.join(genotype,'mpileups',sample)
        if not os.path.exists(mpileup_path) and not os.path.exists(mpileup_path+'.gz'):
            mpileup_commands.append([['samtools','mpileup','-O','-l', snps_file, bamfile], mpileup_path])

def mpileup(incmd):
    with open(incmd[1],'w') as resf:
        subprocess.check_call(incmd[0], stdout=resf)

# execute all samtools mpileup commands in parallel
Parallel(n_jobs=args.proc)(delayed(mpileup)(cmd) for cmd in mpileup_commands)

# 2 Summarize SNP statistics per samples
########################################
summarize_commands = []

for genotype in genotypes:
    for item in allsamples[genotype]:
        sample, bamfile = item
        cellsum_path = os.path.join(genotype,'cellsums',sample)
        if not os.path.exists(cellsum_path) and not os.path.exists(cellsum_path+'.gz'):
            summarize_commands.append([os.path.join(genotype,'mpileups',sample),
                                       os.path.join(genotype,'cellsums',sample)])

def cellsum(params):
    mpileup_file, outfile = params
    with open(outfile, 'w') as resf:
        for line in open(mpileup_file):
            ntcount = [0,0,0,0]
            parts = line.strip().split('\t')
            chrom, pos = parts[:2]
            nb = int(parts[3])
            bases = parts[4]
            quals = parts[5]
            extra = parts[6]
            for base, qual in zip(bases,quals):
                if ord(qual) - 33 >= requiredQ:
                    if base in ('a','A'): ntcount[0] += 1
                    elif base in ('c','C'): ntcount[1] += 1
                    elif base in ('g','G'): ntcount[2] += 1
                    elif base in ('t','T'): ntcount[3] += 1
            if sum(ntcount) > 0:
                resf.write("%s\t%s\t%i\t%i\t%i\t%i\n" % (chrom, pos, ntcount[0], ntcount[1], ntcount[2], ntcount[3]))

Parallel(n_jobs=args.proc)(delayed(cellsum)(cmd) for cmd in summarize_commands)

# gzip mpileup files
gzip_commands = []
for genotype in genotypes:
    for item in allsamples[genotype]:
        sample, bamfile = item
        mpileup_path = os.path.join(genotype,'mpileups',sample)
        if os.path.exists(mpileup_path):
            gzip_commands.append(['gzip', mpileup_path])

Parallel(n_jobs=args.proc)(delayed(subprocess.check_call)(cmd) for cmd in gzip_commands)


# 3 Summarize statistics for all data to extract validated SNPs
###############################################################
snpstats = {} # genotype: chr: pos: [0,0,0,0,0,0,0,0] # reads and cells supporting alleles
chrompos = set([])
for genotype in genotypes:
    snpstats[genotype]={}
    for item in allsamples[genotype]:
        sample, bamfile = item
        cellsum_path = os.path.join(genotype,'cellsums',sample)
        if os.path.exists(cellsum_path):
            for line in open(cellsum_path):
                p = line.strip().split("\t")
                chrom, pos, a, c, g, t = p
                chrompos.add("%s:%s" % (chrom,pos))
                if not chrom in snpstats[genotype]:
                    snpstats[genotype][chrom]={}
                if not pos in snpstats[genotype][chrom]:
                    snpstats[genotype][chrom][pos]=[0,0,0,0,0,0,0,0]
                for idx,val in enumerate(map(int, [a,c,g,t])):
                    snpstats[genotype][chrom][pos][idx] += val
                    if val > 0:
                        snpstats[genotype][chrom][pos][idx+4] += 1

# read in alleles
alleles = {}
for line in open('castsnps_alleles.txt'):
    p = line.strip().split("\t")
    if line[0] == '#': continue
    if not "chr%s:%s"%(p[0],p[1]) in chrompos: continue
    if not 'chr'+p[0] in alleles:
        alleles['chr'+p[0]]={}
    alleles['chr'+p[0]][p[1]]=(p[3],p[4])


# write summary statistics
with open('snpstatistics_cast.txt','w') as resf:
    resf.write("#chrom\tpos\tc57allele\tcastallele\n")
    for cp in chrompos:
        chrom, pos = cp.split(":")
        c57allele, CASTallele = alleles[chrom][pos]
        resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
        resf.write("%s\t" % "\t".join(map(str,snpstats['cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['cast_c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['c57_cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
        resf.write("\n")
with open('snpstatistics_cast_c57.txt','w') as resf:
    resf.write("#chrom\tpos\tc57allele\tcastallele\n")
    for cp in chrompos:
        chrom, pos = cp.split(":")
        c57allele, CASTallele = alleles[chrom][pos]
        resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))                    
        resf.write("%s\t" % "\t".join(map(str,snpstats['cast_c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
        #resf.write("%s\t" % "\t".join(map(str,snpstats['c57_cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
        resf.write("\n")
with open('snpstatistics_c57_cast.txt','w') as resf:
    resf.write("#chrom\tpos\tc57allele\tcastallele\n")
    for cp in chrompos:
        chrom, pos = cp.split(":")
        c57allele, CASTallele = alleles[chrom][pos]
        resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))                    
        #resf.write("%s\t" % "\t".join(map(str,snpstats['cast_c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
        resf.write("%s\t" % "\t".join(map(str,snpstats['c57_cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
        resf.write("\n")
with open('snpstatistics_c57.txt','w') as resf:
    resf.write("#chrom\tpos\tc57allele\tcastallele\n")
    for cp in chrompos:
        chrom, pos = cp.split(":")
        c57allele, CASTallele = alleles[chrom][pos]
        resf.write("%s\t%s\t%s\t%s\t" % (chrom, pos, c57allele, CASTallele))
        #resf.write("%s\t" % "\t".join(map(str,snpstats['cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))
        resf.write("%s\t" % "\t".join(map(str,snpstats['c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))                    
        #resf.write("%s\t" % "\t".join(map(str,snpstats['cast_c57'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
        #resf.write("%s\t" % "\t".join(map(str,snpstats['c57_cast'][chrom].get(pos, [0,0,0,0,0,0,0,0]))))               
        resf.write("\n")

