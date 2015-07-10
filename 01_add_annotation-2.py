#! /usr/bin/python
import optparse, os, sys, subprocess
from subprocess import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input MAF file for filtering')
#parser.add_option('-c', dest='caller',help='somatic mutation caller e.g.) mutect, somaticindelocator')

(options, args) = parser.parse_args()

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def check_dir(dir):
	if os.path.exists(dir):
		pass
	else:
		try:
			os.mkdir(dir)
		except:
			pass


def check_region(line):
	regions = {"3'Flank": "",\
			"3'UTR": "",\
			"Frame_Shift_Del": "",\
			"Frame_Shift_Ins": "",\
			"IGR": "",\
			"In_Frame_Del": "",\
			"In_Frame_Ins": "",\
			"Intron": "",\
			"RNA": "",\
			"Missense_Mutation": "",\
			"Nonsense_Mutation": "",\
			"Silent": "",\
			"Splice_Site": ""}
	exonic = {"Frame_Shift_Del": "",\
			"Frame_Shift_Ins": "",\
			"In_Frame_Del": "",\
			"In_Frame_Ins": "",\
			"Missense_Mutation": "",\
			"Nonsense_Mutation": "",\
			"Splice_Site": ""}

	region = line.strip().split('\t')[8]
	if exonic.has_key(region):
		return 1
	else:
		return 0

def get_pon_snv():
	pon_dic = {}
	pon_snv_file = '/home/D131748/Reference/PON/OPv2_mutect_PON_20150109.vcf'

	for x in open(pon_snv_file).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]

			if ',' in x_1[4]:
				x_2 = x_1[4].split(',')
				for var in x_2:
					dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)
					pon_dic[dic_key] = x
			else:
				var = x_1[4]
				dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)
				pon_dic[dic_key] = x

	return pon_dic

def get_pon_indel():
	pon_dic = {}
	pon_indel_file = '/home/D131748/Reference/PON/OPv2_indel_PON_20140624.vcf'

	for x in open(pon_indel_file).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]

			if ',' in x_1[4]:
				x_2 = x_1[4].split(',')
				for var in x_2:
					dic_key = '%s_%s' % (chr, pos)
					pon_dic[dic_key] = x
			else:
				var = x_1[4]
				dic_key = '%s_%s' % (chr, pos)
				pon_dic[dic_key] = x
			if ',' in x_1[3]:
				print 'col 3 has ,.'

	return pon_dic

def get_cosmic():
	cosmic_dic = {}
	cosmic_all = '/home/D131748/Reference/cosmic/cosmic_v71_indel.b37.vcf'	

	for x in open(cosmic_all).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]
			var = x_1[4]
			dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)
			cosmic_dic[dic_key] = x
	return cosmic_dic

def get_dbsnp():
	dbsnp_dic = {}
	dbsnp_file = '/home/D131748/Reference/GATK_bundle/2.5/b37/dbsnp_137_common.b37.vcf'
	
	for x in open(dbsnp_file).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]
			var = x_1[4]
			dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)
			dbsnp_dic[dic_key] = x
	return dbsnp_dic

def get_pfam():
	pfam_dic = {}
	pfam_file = '/home/D131748/Reference/Pfam-A.clans.tsv'

	for x in open(pfam_file).xreadlines():
		x_1 = x.strip().split('\t')
		id = x_1[0]
		domain = x_1[4]
		pfam_dic[id] = domain
	return pfam_dic


def db_filter(db, variant):
	if db.has_key(variant):
		return 1
	else:
		return 0
	

def ext_seq(chr, pos, ref, alt, type):
	chrs = {'chr1':'', 'chr2':'', 'chr3':'', 'chr4':'', 'chr5':'', 'chr6':'', 'chr7':'', 'chr8':'', 'chr9':'', 'chr10':'', 'chr11':'', 'chr12':'', 'chr13':'', 'chr14':'', 'chr15':'', 'chr16':'', 'chr17':'', 'chr18':'', 'chr19':'', 'chr20':'', 'chr21':'', 'chr22':'', 'chrX':'', 'chrY':'', '1':'', '2':'', '3':'', '4':'', '5':'', '6':'', '7':'', '8':'', '9':'', '10':'', '11':'', '12':'', '13':'', '14':'', '15':'', '16':'', '17':'', '18':'', '19':'', '20':'', '21':'', '22':'', 'X':'', 'Y':''}

	# extract up-/down-stream sequences
	if type == 'SNP':
		start = int(pos) - 100
		pos1 = int(pos) - 1
		pos2 = int(pos) + 1
		end = int(pos) + 100
	elif type == 'DEL':
		start = int(pos) - 100
		pos1 = int(pos)
		pos2 = int(pos) + len(ref) + 1
		end = int(pos2) + 100
	elif type == 'INS':
		start = int(pos) - 100
		pos = int(pos)
		pos2 = int(pos) + 1
		end = int(pos2) + 100
	
	if chr.upper()[0:3]  == 'CHR':
		ref_seq = '/home/D131748/Reference/Human/hg19/ucsc.hg19.fasta'
	else:
		ref_seq = '/home/D131748/Reference/Human/b37/human_g1k_v37.fasta'

	samtools_path = '/home/D131748/apps/ETC/samtools'

	cmd1 = '%s/samtools faidx %s %s:%s-%s' % (samtools_path, ref_seq, chr, start, pos)
	cmd2 = '%s/samtools faidx %s %s:%s-%s' % (samtools_path, ref_seq, chr, pos2, end)

	pipe1 = Popen(cmd1, shell=True, stdout=PIPE)
	pipe2 = Popen(cmd2, shell=True, stdout=PIPE)

	con1 = pipe1.stdout.read().split('\n')[1:]
	con2 = pipe2.stdout.read().split('\n')[1:]

	seq = '%s[%s/%s]%s' % (''.join(con1), ref, alt, ''.join(con2))

	return seq

def get_tier():
	tier_dic = {}
	
	tier_file = '/home/D131748/Reference/annotation/Tier/tier.txt'

	for x in open(tier_file).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			tier_dic[x_1[0]] = x_1[2]
	return tier_dic

def get_cgc():
	cgc_dic = {}
	
	cgc_file = '/home/D131748/Reference/annotation/CGC/cancer_gene_census.txt'

	for x in open(cgc_file).xreadlines():
		x_1 = x.replace('\n','').split('\t')
		if x_1[0] != 'Symbol':
			cgc_dic[x_1[0]] = {}
			cgc_dic[x_1[0]]['CGC_Mutation_Type'] = x_1[12]
			cgc_dic[x_1[0]]['CGC_Tumor_Types_Somatic'] = x_1[7]
			cgc_dic[x_1[0]]['CGC_Tumor_Types_Germline'] = x_1[8]
			cgc_dic[x_1[0]]['CGC_Translocation_Partner'] = x_1[13]
			cgc_dic[x_1[0]]['CGC_Other_Diseases'] = x_1[15]

	return cgc_dic

def get_target():
	target_dic = {}
	
	target_file = '/home/D131748/Reference/annotation/TARGET/TARGET_db_v2_05042014.txt'

	for x in open(target_file).xreadlines():
		x_1 = x.replace('\n','').split('\t')
		if x_1[0] != 'Gene':
			target_dic[x_1[0]] = {}
			target_dic[x_1[0]]['category'] = x_1[1]
			target_dic[x_1[0]]['additional_category'] = x_1[2]
			target_dic[x_1[0]]['rationale'] = x_1[3]
			target_dic[x_1[0]]['alteration_types'] = x_1[4]
			target_dic[x_1[0]]['therapeutic_agents'] = x_1[5]

	return target_dic


def cancer_gene(gene):
	pathway = {}
	type = {}

	src_path ='/home/D131748/src/'

	for x in open('/home/D131748/src/cancer_gene.txt').xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			pathway[x_1[0]] = x_1[1]
			type[x_1[0]] = x_1[2]

	if pathway.has_key(gene) and type.has_key(gene):
		return pathway[gene], type[gene]
	else:
		return '-', '-'

def __main__():
	check_option(options.input, "You need to use '-i' option with input MAF file\n")

	input_maf = options.input
	
	tier_dic = get_tier()
	cgc_dic = get_cgc()
	target_dic = get_target()
	pfam_dic = get_pfam()

	pon_dic_indel = get_pon_indel()
	pon_dic_snv = get_pon_snv()
	dbsnp_dic = get_dbsnp()
	cosmic_dic = get_cosmic()
	
	wname = input_maf.replace('.maf', '.exon.maf') 
	wname2 = input_maf.replace('.maf', '.all.maf') 
	wfile = open(wname,'w') 
	wfile2 = open(wname2,'w') 
	
        for x in open(input_maf).xreadlines(): 
		if x[0:11] == 'Hugo_Symbol':
			x_1 = x.replace('\n','').split('\t')
			center = x_1.index('Center')
			exon_number = x_1.index('Exon_Number')
			exon = x_1.index('EXON')
			intron = x_1.index('INTRON')

			t_pos_total = x_1.index('t_depth')
			t_pos_ref = x_1.index('t_ref_count')
			t_pos_var = x_1.index('t_alt_count')
			n_pos_total = x_1.index('n_depth')
			n_pos_ref = x_1.index('n_ref_count')
			n_pos_var = x_1.index('n_alt_count')
			mutation_status = x_1.index('Mutation_Status')
			sequence_source = x_1.index('Sequence_Source')
			sequencer = x_1.index('Sequencer')
			codon = x_1.index('HGVSc')
			protein = x_1.index('HGVSp_Short')
			domain = x_1.index('DOMAINS')

			header_count = x_1[t_pos_total: n_pos_var+1]
			del x_1[t_pos_total: n_pos_var+1]
			x_1[-1:] = header_count
			x_1.append('Variant_allele_ratio')
			x_1.append('Significant_call(VAF>=0.05)')
			x_1.append('Tier')
			x_1.append('Tumor_Normal_Mode')
			x_1.append('Gene_Codon_ProteinChange')
			x_1.append('CGC_Mutation_Type')
			x_1.append('CGC_Tumor_Types_Somatic')
			x_1.append('CGC_Tumor_Types_Germline')
			x_1.append('CGC_Translocation_Partner')
			x_1.append('CGC_Other_Diseases')
			x_1.append('TARGET_Category')
			x_1.append('TARGET_Rationale')
			x_1.append('TARGET_Alteration_Types')
			x_1.append('TARGET_Therapeutic_Agents')
			x_1.append('CancerGene_pathway')
			x_1.append('CancerGene_types')
			x_1.append('dbsnp_filter')
			x_1.append('PON_filter')
			x_1.append('cosmic_filter')
			x_1.append('Pfam domain')
			x_1.append('Up-/Down-stream sequences')
			cmd = '\t'.join(x_1)
			wfile.write('%s\n' % cmd)
			wfile2.write('%s\n' % cmd)

		else:
			exonic = check_region(x)

			x_1 = x.replace('\n','').split('\t')
			x_2 = x.replace('\n','').split('\t')
			for no, k in enumerate(x_1):
				if k == '':
					x_1[no] = '-'
			gene = x_1[0]
			mutation_type = x_1[9]
			chr = x_1[4]
			if mutation_type == 'DEL':
				pos = str(int(x_1[5])-1)
			else:
				pos = x_1[5]
			ref = x_1[10]
			alt = x_1[12]

			if mutation_type == 'DEL':
				exonic_key = '%s_%s' % (x_1[4], int(x_1[5])-1)
			else: 
				exonic_key = '%s_%s' % (x_1[4], x_1[5])


			(pathway, t_type) = cancer_gene(gene)
			variant1 = '%s' % (exonic_key)
			variant2 = '%s_%s_%s' % (exonic_key, ref, alt)

			seq = ext_seq(chr, pos, ref, alt, mutation_type)
			
			x_1[center] = 'ASAN-CCGD'
			x_1[exon_number] = x_1[exon_number].replace('/','//')
			x_1[exon] = x_1[exon].replace('/','//')
			x_1[intron] = x_1[intron].replace('/','//')
			"""
			x_1[t_pos_total] = str(int(t_count_ref) + int(t_count_var))
			x_1[t_pos_ref] = str(t_count_ref)
			x_1[t_pos_var] = str(t_count_var)
			if mode == 'unpaired':
				x_1[n_pos_total] = '0'
				x_1[n_pos_ref] = '0'
				x_1[n_pos_var] = '0'
			else:
				x_1[n_pos_total] = str(int(n_count_ref) + int(n_count_var))
				x_1[n_pos_ref] = str(n_count_ref)
				x_1[n_pos_var] = str(n_count_var)
			"""
			count_content = x_1[t_pos_total: n_pos_var+1]
			del x_1[t_pos_total: n_pos_var+1]
			x_1[-1:] = count_content

			x_1[mutation_status] = 'Somatic'
			x_1[sequence_source] = 'OncoPanel_V2'
			x_1[sequencer] = 'Illumina MiSEQ'
			
			var_ratio = '%5.2f' % float(float(x_2[t_pos_var]) / float(x_2[t_pos_total]))

			x_1.append(str(var_ratio))
			#if (int(t_count_ref) + int(t_count_var)) >= 30 and float(var_ratio) >= 0.1:
			if float(var_ratio) >= 0.05:
				x_1.append('Yes')
			else:
				x_1.append('No')
			if tier_dic.has_key(gene):
				x_1.append(tier_dic[gene])
			else:
				x_1.append('-')
			x_1.append('unpaired')
			x_1.append('%s_%s_%s' % (gene, x_1[codon], x_1[protein]))	
			
			if cgc_dic.has_key(gene):
				x_1.append(cgc_dic[gene]['CGC_Mutation_Type'])
				x_1.append(cgc_dic[gene]['CGC_Tumor_Types_Somatic'])
				x_1.append(cgc_dic[gene]['CGC_Tumor_Types_Germline'])
				x_1.append(cgc_dic[gene]['CGC_Translocation_Partner'])
				x_1.append(cgc_dic[gene]['CGC_Other_Diseases'])
			else:
				x_1.append('-')
				x_1.append('-')
				x_1.append('-')
				x_1.append('-')
				x_1.append('-')
			if target_dic.has_key(gene):
				x_1.append(target_dic[gene]['category'])
				x_1.append(target_dic[gene]['rationale'])
				x_1.append(target_dic[gene]['alteration_types'])
				x_1.append(target_dic[gene]['therapeutic_agents'])
			else:
				x_1.append('-')
				x_1.append('-')
				x_1.append('-')
				x_1.append('-')
			x_1.append(pathway)
			x_1.append(t_type)

			x_1.append(str(db_filter(dbsnp_dic, variant2)))
			if mutation_type == 'DEL' or mutation_type == 'INS':
				x_1.append(str(db_filter(pon_dic_indel, variant1)))
			else:
				x_1.append(str(db_filter(pon_dic_snv, variant2)))
			x_1.append(str(db_filter(cosmic_dic, variant2)))
			if 'Pfam_domain:' in x_2[domain]:
				x_3 = x_2[domain].split(',')
				for y in x_3:
					if 'Pfam_domain:' in y:
						x_1.append(pfam_dic[y.replace('Pfam_domain:','')])
			else:
				x_1.append('-')
				

			x_1.append(seq)
			cmd = '\t'.join(x_1)
			if exonic == 1:	
				wfile.write('%s\n' % cmd)
			wfile2.write('%s\n' % cmd)
	wfile.close()
	wfile2.close()
if __name__=="__main__":__main__()
