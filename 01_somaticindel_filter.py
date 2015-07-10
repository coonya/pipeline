#! /usr/bin/python
import optparse, os, sys, subprocess

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-i', dest='input',help='input file for filtering')
#parser.add_option('-m', dest='mode',default='paired', help='analysis mode. default: paired')

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

def Parameters(line, sample_col):
	line_tab		= line.strip().split('\t')
	t_els			= line_tab[sample_col].split(':')
	t_els2			= line_tab[8].split(':')
	AD			= t_els2.index("AD") 
	DP			= t_els2.index("DP") 
	MM			= t_els2.index("MM") 
	MQS			= t_els2.index("MQS") 
	NQSBQ			= t_els2.index("NQSBQ") 
	NQSMM			= t_els2.index("NQSMM") 
	REnd			= t_els2.index("REnd") 
	RStart			= t_els2.index("RStart") 
	SC			= t_els2.index("SC")
	if line_tab[7] == 'SOMATIC':
		somatic		= line_tab[7]
	else:
		somatic		= 'GERMLINE'

	cov			= line_tab[6]

	t_genotype		= t_els[0]
	t_obs_count_var		= t_els[AD].split(',')[0]
	t_obs_count_ref		= t_els[AD].split(',')[1]
	t_obs_count_total	= t_els[DP]
	t_av_mm_var		= t_els[MM].split(',')[0]
	t_av_mm_ref		= t_els[MM].split(',')[1]
	t_av_mapq_var		= t_els[MQS].split(',')[0]
	t_av_mapq_ref		= t_els[MQS].split(',')[1]
	t_nqs_av_qual_var	= t_els[NQSBQ].split(',')[0]
	t_nqs_av_qual_ref	= t_els[NQSBQ].split(',')[1]
	t_nqs_mm_rate_var	= t_els[NQSMM].split(',')[0]
	t_nqs_mm_rate_ref	= t_els[NQSMM].split(',')[1]
	t_offset_rend		= t_els[REnd].split(',')[0]
	t_offset_rend_mad	= t_els[REnd].split(',')[1]
	t_offset_rstart		= t_els[RStart].split(',')[0]
	t_offset_rstart_mad	= t_els[RStart].split(',')[1]
	t_strand_count_var1	= t_els[SC].split(',')[0]
	t_strand_count_var2	= t_els[SC].split(',')[1]
	t_strand_count_ref1	= t_els[SC].split(',')[2]
	t_strand_count_ref2	= t_els[SC].split(',')[3]

	return somatic, t_genotype, t_obs_count_var, t_obs_count_ref, t_obs_count_total, t_av_mm_var, t_av_mm_ref, t_av_mapq_var, t_av_mapq_ref, t_nqs_av_qual_var, t_nqs_av_qual_ref, t_nqs_mm_rate_var, t_nqs_mm_rate_ref, t_offset_rend, t_offset_rend_mad, t_offset_rstart, t_offset_rstart_mad, t_strand_count_var1, t_strand_count_var2, t_strand_count_ref1, t_strand_count_ref2, cov



def checkorder(line, mode):
	x_1 = line.strip().split('\t')
	if mode == 'paired':
		if x_1[9][-1] == 'N' and x_1[10][-1] =='T':
			normal_col = 9
			tumor_col = 10
		elif x_1[9][-1] == 'T' and x_1[10][-1] =='N':
			normal_col = 10
			tumor_col = 9
		elif x_1[9] == 'none':
			tumor_col = 10
			normal_col = 'NA'
		elif x_1[10] == 'none':
			tumor_col = 9
			normal = 'NA'

	elif mode == 'unpaired':
		tumor_col = 9
		normal_col = 'NA'

	return normal_col, tumor_col
	
"""
#def pass_filter(mode, somatic, t_obs_count_var,t_av_mm_var, t_av_mm_ref, t_av_mapq_var, t_av_mapq_ref, t_nqs_av_qual_var, t_nqs_av_qual_ref, t_offset_rstart, t_offset_rstart_mad, t_offset_rend, t_offset_rend_mad, cov):
	if int(t_obs_count_var) >= 4\
	or float(t_av_mm_var) <= 4\
	or float(t_av_mm_ref) <= 4\
	or float(t_nqs_mm_rate_var) <= 0.15\
	or float(t_nqs_mm_rate_ref) <= 0.15\
	or float(t_nqs_av_qual_var) >= 20\
	or float(t_nqs_av_qual_ref) >= 20\
	or int(t_offset_rstart) >= 10\
	or int(t_offset_rstart_mad) >=3\
	or int(t_offset_rend) >=10\
	or int(t_offset_rend_mad) >= 3: 

	if int(t_obs_count_var) < 3\
	or float(t_av_mm_var) > 3\
	or float(t_av_mm_ref) > 3\
	or float(t_av_mapq_var) < 28\
	or float(t_av_mapq_ref) < 28\
	or float(t_nqs_av_qual_var) < 19\
	or float(t_nqs_av_qual_ref) < 19\
	or int(t_offset_rstart) < 9\
	or int(t_offset_rstart_mad) < 2\
	or int(t_offset_rend) < 9\
	or :int(t_offset_rend_mad) < 2\
	or 'TCov' in cov\
	or 'NCov' in cov:
		return 0
	else:
		if mode == 'paired':
			if somatic == 'SOMATIC':
				return 1
			else:
				return 0
"""	

def pass_filter(mode, somatic, t_obs_count_var,t_av_mm_var, t_av_mm_ref, t_av_mapq_var, t_av_mapq_ref, t_nqs_mm_rate_var, t_nqs_mm_rate_ref, t_nqs_av_qual_var, t_nqs_av_qual_ref, t_offset_rstart, t_offset_rstart_mad, t_offset_rend, t_offset_rend_mad):


	if int(t_obs_count_var) >= 4\
	and float(t_av_mm_var) <= 4\
	and float(t_av_mm_ref) <= 4\
	and float(t_av_mapq_var) >= 29\
	and float(t_av_mapq_ref) >=29\
	and float(t_nqs_mm_rate_var) <= 0.15\
	and float(t_nqs_mm_rate_ref) <= 0.15\
	and float(t_nqs_av_qual_var) >= 20\
	and float(t_nqs_av_qual_ref) >= 20\
	and int(t_offset_rstart) >= 10\
	and int(t_offset_rstart_mad) >=3\
	and int(t_offset_rend) >=10\
	and int(t_offset_rend_mad) >= 3: 
		if mode == 'paired':
			if somatic == 'SOMATIC':
				return 1
			else:
				return 0
				
		else:
			return 1
	else:
		return 0

def get_pon():
	pon_dic = {}
	pon_file = '/home/D131748/Reference/PON/OPv2_indel_PON_20150121_splited.vcf'
	
	for x in open(pon_file).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]

			if ',' in x_1[4]:
				x_2 = x_1[4].split(',')
				for var in x_2:
					dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)
#					dic_key = '%s_%s' % (chr, pos)
					pon_dic[dic_key] = x
			else:
				var = x_1[4]
				dic_key = '%s_%s_%s_%s' % (chr, pos, ref, var)
#				dic_key = '%s_%s' % (chr, pos)
				pon_dic[dic_key] = x
			if ',' in x_1[3]:
				print 'col 3 has ,.'
	return pon_dic

def get_cosmic():
	cosmic_dic = {}
#	cosmic_file = '/home/D131748/Reference/dbsnp/CosmicAllVariants_v64_02042013_noLimit.resorted2.long_indel_removed.vcf'
	cosmic_file = '/home/D131748/Reference/cosmic/cosmic_v71_indel.b37.vcf'
	
	for x in open(cosmic_file).xreadlines():
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



def __main__():
	check_option(options.input, "You need to use '-i' option with input file\n")

	input_vcf = options.input
	wname = input_vcf.replace('.vcf', '.filtered.vcf') 
        wfile = open(wname,'w') 
	
	pon_dic = get_pon()
#	cosmic_dic = get_cosmic()
	dbsnp_dic = get_dbsnp()
 
        for x in open(input_vcf).xreadlines(): 
                if x[0] != '#': 
			x_1 = x.strip().split('\t')
			chr = x_1[0]
			pos = x_1[1]
			ref = x_1[3]
			var = x_1[4]
			com_key = '%s_%s_%s_%s' % (chr, pos, ref, var)
			com_key2 = '%s_%s' % (chr, pos)




			if mode == 'paired':
				(somatic,\
				t_genotype,		t_obs_count_var,\
				t_obs_count_ref,	t_obs_count_total,\
				t_av_mm_var,		t_av_mm_ref,\
				t_av_mapq_var,		t_av_mapq_ref,\
				t_nqs_av_qual_var,	t_nqs_av_qual_ref,\
				t_nqs_mm_rate_var,	t_nqs_mm_rate_ref,\
				t_offset_rend,		t_offset_rend_mad,\
				t_offset_rstart,	t_offset_rstart_mad,\
				t_strand_count_var1,	t_strand_count_var2,\
				t_strand_count_ref1,	t_strand_count_ref2,\
				cov\
				) = Parameters(x, tumor_col)

#				pass_score = pass_filter('paired', somatic, t_obs_count_var,t_av_mm_var, t_av_mm_ref, t_av_mapq_var, t_av_mapq_ref, t_nqs_av_qual_var, t_nqs_av_qual_ref, t_offset_rstart, t_offset_rstart_mad, t_offset_rend, t_offset_rend_mad, cov)
				pass_score = pass_filter('paired', somatic, t_obs_count_var,t_av_mm_var, t_av_mm_ref, t_av_mapq_var, t_av_mapq_ref, t_nqs_mm_rate_var, t_nqs_mm_rate_ref, t_nqs_av_qual_var, t_nqs_av_qual_ref, t_offset_rstart, t_offset_rstart_mad, t_offset_rend, t_offset_rend_mad)
			elif mode == 'unpaired':
				(somatic,\
				t_genotype,		t_obs_count_var,\
				t_obs_count_ref,	t_obs_count_total,\
				t_av_mm_var,		t_av_mm_ref,\
				t_av_mapq_var,		t_av_mapq_ref,\
				t_nqs_av_qual_var,	t_nqs_av_qual_ref,\
				t_nqs_mm_rate_var,	t_nqs_mm_rate_ref,\
				t_offset_rend,		t_offset_rend_mad,\
				t_offset_rstart,	t_offset_rstart_mad,\
				t_strand_count_var1,	t_strand_count_var2,\
				t_strand_count_ref1,	t_strand_count_ref2,\
				cov\
				) = Parameters(x, tumor_col)
			
#				pass_score = pass_filter('unpaired', somatic, t_obs_count_var,t_av_mm_var, t_av_mm_ref, t_av_mapq_var, t_av_mapq_ref, t_nqs_av_qual_var, t_nqs_av_qual_ref, t_offset_rstart, t_offset_rstart_mad, t_offset_rend, t_offset_rend_mad, cov)
				pass_score = pass_filter('unpaired', somatic, t_obs_count_var,t_av_mm_var, t_av_mm_ref, t_av_mapq_var, t_av_mapq_ref, t_nqs_mm_rate_var, t_nqs_mm_rate_ref, t_nqs_av_qual_var, t_nqs_av_qual_ref, t_offset_rstart, t_offset_rstart_mad, t_offset_rend, t_offset_rend_mad)	

			if pass_score == 1:
#				wfile.write(x)
		
				if pon_dic.has_key(com_key) or dbsnp_dic.has_key(com_key):
#					if cosmic_dic.has_key(com_key):
#						wfile.write(x)
#					else:
#						pass
					pass
				else:
					wfile.write(x)
#			wfile.write(x)
				
		else:
			wfile.write(x)
			if x[0:6] == '#CHROM':
				header_len = len(x.strip().split('\t'))

				#if header_len == 11:
				if header_len >= 11:
					mode = 'paired'
				elif header_len == 10:
					mode = 'unpaired'

				(normal_col, tumor_col) = checkorder(x, mode)

	wfile.close()

if __name__=="__main__":__main__()
