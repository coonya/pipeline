import sys, MySQLdb
from utils import *

def get_sample_list(sample_file):
	dic = {}
	sample_list = []
	normal_list = []
	tumor_list = []
	paired_list = []
	unpaired_list = []

	for x in open(sample_file).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			order = x_1[1]
			sample = x_1[0].split('.')[0]
			t_type = x_1[2].upper()

			sample_list.append(sample)
			
			if t_type == 'NORMAL':
				normal_list.append(sample)
			elif t_type == 'TUMOR':
				tumor_list.append(sample)

			if dic.has_key(str(order)):
				if dic[str(order)].has_key(t_type):
					pass
				else:
					dic[str(order)][t_type] = sample
			else:
				dic[str(order)] = {}
				dic[str(order)][t_type] = sample

	for x in dic:
		if len(dic[x]) == 2:
			paired_list.append(dic[x]['NORMAL'])
			paired_list.append(dic[x]['TUMOR'])
		elif len(dic[x]) == 1 and dic[x].get('TUMOR') != None:
			unpaired_list.append(dic[x]['TUMOR'])
		elif len(dic[x]) == 1 and dic[x].get('NORMAL') != None:
			unpaired_list.append(dic[x]['NORMAL'])

	return dic, sample_list, normal_list, tumor_list, paired_list, unpaired_list


def get_sample_list_agg(project):
	dic = {}
	sample_list = []

	## connet to the mysql server
	db = MySQLdb.connect(host='localhost', user='kdh', passwd='kdh', db='cpcm')
	cursor = db.cursor(MySQLdb.cursors.DictCursor)
	
	cursor.execute("select * from samples_sequenced where project = '%s' and exclusion IS NULL" % (project))
	check_run = cursor.fetchall()
	
	for x in check_run:
		bam1 = x['bam'].replace('.bam','')
		bam2 = bam1.split('_')
		sample = '%s_%s' % (bam2[0], bam2[2])
	
		run_id = x['run_id']
	
		working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % run_id
		working_tmp = '%s/tmp' % working_dir
		recal_bam = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (working_tmp, bam1)
	
		if dic.has_key(sample):
			dic[sample].append(recal_bam)
		else:
			dic[sample] = [recal_bam]
			sample_list.append(sample)
	
	paired_samples = {}
	non_paired_samples = {}
	paired_sample_list = []
	non_paired_sample_list = []

	for y in dic:
		sample_id = y[:-1]
		tumor = '%sT' % sample_id
		normal = '%sN' % sample_id
		
		if dic.has_key(tumor) and dic.has_key(normal):
			paired_samples[tumor] = ['%s' % tumor]
			paired_samples[tumor].append(normal)
			paired_sample_list.append(y)
		else:
			non_paired_samples[y] = y
			non_paired_sample_list.append(y)

	sample_list.sort()
	paired_sample_list.sort()
	non_paired_sample_list.sort()

	return dic, paired_samples, non_paired_samples, sample_list, paired_sample_list, non_paired_sample_list
	

def get_file_name(base_dir, sample_list):
	file_dic = {}

	working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % base_dir
	working_tmp = '%s/tmp' % working_dir
	metrics_dir = '%s/metrics' % working_dir
	mutect_dir = '%s/mutect' % working_dir
	somaticindelocator_dir = '%s/SI' % working_dir
	aggregation_dir = '%s/tmp/aggregation' % working_dir
#	viscapcancer_dir = '%s/VisCapCancer' % working_dir

	for x in sample_list:
		file_dic[x] = {}
		file_dic[x]['bam'] = '%s/%s.bam' % (working_tmp, x)
		file_dic[x]['fastq1'] = '%s/%s_R1.fastq.gz' % (working_tmp, x)
		file_dic[x]['fastq2'] = '%s/%s_R2.fastq.gz' % (working_tmp, x)
		file_dic[x]['sai1'] = '%s/%s_R1.fastq.gz.sai' % (working_tmp, x)
		file_dic[x]['sai2'] = '%s/%s_R2.fastq.gz.sai' % (working_tmp, x)
		file_dic[x]['sam'] = '%s/%s.fastq.gz.initialAlign.sam' % (working_tmp, x)
		file_dic[x]['mergebamalignment'] = '%s/%s.fastq.gz.initialAlign.merged.bam' % (working_tmp, x)
		file_dic[x]['dedup'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.bam' % (working_tmp, x)
		file_dic[x]['dedup_metrics'] = '%s/%s_duplicateMetrics.txt' % (metrics_dir, x)
		file_dic[x]['realign'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.bam' % (working_tmp, x)
		file_dic[x]['recalcsv'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.csv' % (working_tmp, x)
		file_dic[x]['recal'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (working_tmp, x)
		file_dic[x]['map'] = '%s/%s.map' % (aggregation_dir, x)
		file_dic[x]['realign_agg'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.cleaned.bam' % (aggregation_dir, x)
		file_dic[x]['indel_interval'] = '%s/%s_indelCocleanTargets.intervals' % (working_tmp, x)
		file_dic[x]['exon_metrics'] = '%s/%s_exon_targets_hsMetrics.txt' % (metrics_dir, x)
		file_dic[x]['intron_metrics'] = '%s/%s_intron_targets_hsMetrics.txt' % (metrics_dir, x)
		file_dic[x]['full_metrics'] = '%s/%s_full_targets_hsMetrics.txt' % (metrics_dir, x)
		file_dic[x]['exon_coverage'] = '%s/%s_exon_targets_hsIntervals.txt' % (metrics_dir, x)
		file_dic[x]['intron_coverage'] = '%s/%s_intron_targets_hsIntervals.txt' % (metrics_dir, x)
		file_dic[x]['full_coverage'] = '%s/%s_full_targets_hsIntervals.txt' % (metrics_dir, x)
		#file_dic[x]['depthofcoverage'] = '%s/%s.dedup.cleaned.cov' % (viscapcancer_dir, x)
		file_dic[x]['mutect_vcf'] = '%s/%s.mutect.vcf' % (mutect_dir, x)
		file_dic[x]['filtered_mutect_vcf'] = '%s/%s.mutect.filtered.vcf' % (mutect_dir, x)
		file_dic[x]['filtered_mutect_txt'] = '%s/%s.mutect.filtered.txt' % (mutect_dir, x)
		file_dic[x]['mutect_txt'] = '%s/%s.mutect.txt' % (mutect_dir, x)
		file_dic[x]['mutect_wig'] = '%s/%s.mutect.coverage.wig' % (mutect_dir, x)
		file_dic[x]['nonfilter_mutect_vcf'] = '%s/%s.mutect.nonfilter.vcf' % (mutect_dir, x)
		file_dic[x]['nonfilter_mutect_txt'] = '%s/%s.mutect.nonfilter.txt' % (mutect_dir, x)
		file_dic[x]['somaticindelocator'] = '%s/%s.somaticindelocator.vcf' % (somaticindelocator_dir, x)
		file_dic[x]['filtered_somaticindelocator'] = '%s/%s.somaticindelocator.filtered.vcf' % (somaticindelocator_dir, x)
		file_dic[x]['nonfilter_somaticindelocator'] = '%s/%s.somaticindelocator.nonfilter.vcf' % (somaticindelocator_dir, x)
		file_dic[x]['mutect_maf'] = '%s/%s.mutect.filtered.maf' % (mutect_dir, x)
		file_dic[x]['nonfilter_mutect_maf'] = '%s/%s.mutect.nonfilter.maf' % (mutect_dir, x)
		file_dic[x]['somaticindelocator_maf'] = '%s/%s.somaticindelocator.filtered.maf' % (somaticindelocator_dir, x)
		file_dic[x]['nonfilter_somaticindelocator_maf'] = '%s/%s.somaticindelocator.nonfilter.maf' % (somaticindelocator_dir, x)

	return file_dic


def get_file_name_agg(sample_list):
	file_dic = {}

	working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % base_dir
	working_tmp = '%s/tmp' % working_dir
	metrics_dir = '%s/metrics' % working_dir
	mutect_dir = '%s/mutect' % working_dir
	somaticindelocator_dir = '%s/SI' % working_dir
	aggregation_dir = '%s/tmp/aggregation' % working_dir
#	viscapcancer_dir = '%s/VisCapCancer' % working_dir

	for x in sample_list:
		file_dic[x] = {}
		file_dic[x]['bam'] = '%s/%s.bam' % (working_tmp, x)
		file_dic[x]['fastq1'] = '%s/%s_R1.fastq.gz' % (working_tmp, x)
		file_dic[x]['fastq2'] = '%s/%s_R2.fastq.gz' % (working_tmp, x)
		file_dic[x]['sai1'] = '%s/%s_R1.fastq.gz.sai' % (working_tmp, x)
		file_dic[x]['sai2'] = '%s/%s_R2.fastq.gz.sai' % (working_tmp, x)
		file_dic[x]['sam'] = '%s/%s.fastq.gz.initialAlign.sam' % (working_tmp, x)
		file_dic[x]['mergebamalignment'] = '%s/%s.fastq.gz.initialAlign.merged.bam' % (working_tmp, x)
		file_dic[x]['dedup'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.bam' % (working_tmp, x)
		file_dic[x]['dedup_metrics'] = '%s/%s_duplicateMetrics.txt' % (metrics_dir, x)
		file_dic[x]['realign'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.bam' % (working_tmp, x)
		file_dic[x]['recalcsv'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.csv' % (working_tmp, x)
		file_dic[x]['recal'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (working_tmp, x)
		file_dic[x]['map'] = '%s/%s.map' % (aggregation_dir, x)
		file_dic[x]['realign_agg'] = '%s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.cleaned.bam' % (aggregation_dir, x)
		file_dic[x]['indel_interval'] = '%s/%s_indelCocleanTargets.intervals' % (working_tmp, x)
		file_dic[x]['exon_metrics'] = '%s/%s_exon_targets_hsMetrics.txt' % (metrics_dir, x)
		file_dic[x]['intron_metrics'] = '%s/%s_intron_targets_hsMetrics.txt' % (metrics_dir, x)
		file_dic[x]['full_metrics'] = '%s/%s_full_targets_hsMetrics.txt' % (metrics_dir, x)
		file_dic[x]['exon_coverage'] = '%s/%s_exon_targets_hsIntervals.txt' % (metrics_dir, x)
		file_dic[x]['intron_coverage'] = '%s/%s_intron_targets_hsIntervals.txt' % (metrics_dir, x)
		file_dic[x]['full_coverage'] = '%s/%s_full_targets_hsIntervals.txt' % (metrics_dir, x)
		#file_dic[x]['depthofcoverage'] = '%s/%s.dedup.cleaned.cov' % (viscapcancer_dir, x)
		file_dic[x]['mutect_vcf'] = '%s/%s.mutect.vcf' % (mutect_dir, x)
		file_dic[x]['filtered_mutect_vcf'] = '%s/%s.mutect.filtered.vcf' % (mutect_dir, x)
		file_dic[x]['filtered_mutect_txt'] = '%s/%s.mutect.filtered.txt' % (mutect_dir, x)
		file_dic[x]['mutect_txt'] = '%s/%s.mutect.txt' % (mutect_dir, x)
		file_dic[x]['mutect_wig'] = '%s/%s.mutect.coverage.wig' % (mutect_dir, x)
		file_dic[x]['nonfilter_mutect_vcf'] = '%s/%s.mutect.nonfilter.vcf' % (mutect_dir, x)
		file_dic[x]['nonfilter_mutect_txt'] = '%s/%s.mutect.nonfilter.txt' % (mutect_dir, x)
		file_dic[x]['somaticindelocator'] = '%s/%s.somaticindelocator.vcf' % (somaticindelocator_dir, x)
		file_dic[x]['filtered_somaticindelocator'] = '%s/%s.somaticindelocator.filtered.vcf' % (somaticindelocator_dir, x)
		file_dic[x]['nonfilter_somaticindelocator'] = '%s/%s.somaticindelocator.nonfilter.vcf' % (somaticindelocator_dir, x)
		file_dic[x]['mutect_maf'] = '%s/%s.mutect.filtered.maf' % (mutect_dir, x)
		file_dic[x]['nonfilter_mutect_maf'] = '%s/%s.mutect.nonfilter.maf' % (mutect_dir, x)
		file_dic[x]['somaticindelocator_maf'] = '%s/%s.somaticindelocator.filtered.maf' % (somaticindelocator_dir, x)
		file_dic[x]['nonfilter_somaticindelocator_maf'] = '%s/%s.somaticindelocator.nonfilter.maf' % (somaticindelocator_dir, x)

	return file_dic



def print_sample(sample_list, normal_list, tumor_list, paired_list, unpaired_list):
	logger = logging.getLogger('root')

	try:
		print '##########################################################################################'
		logger.info('Total input samples : %s - [%s]' % (len(sample_list), ','.join(sample_list)))
		logger.info('Normal samples      : %s - [%s]' % (len(normal_list), ','.join(normal_list)))
		logger.info('Tumor samples       : %s - [%s]' % (len(tumor_list), ','.join(tumor_list)))
		logger.info('Paired samples      : %s - [%s]' % (len(paired_list)/2, ','.join(paired_list)))
		logger.info('Unpaired samples    : %s - [%s]' % (len(unpaired_list), ','.join(unpaired_list)))

		print_info('Total input samples : %s - [%s]' % (len(sample_list), ','.join(sample_list)))
		print_info('Normal samples      : %s - [%s]' % (len(normal_list), ','.join(normal_list)))
		print_info('Tumor samples       : %s - [%s]' % (len(tumor_list), ','.join(tumor_list)))
		print_info('Paired samples      : %s - [%s]' % (len(paired_list)/2, ','.join(paired_list)))
		print_info('Unpaired samples    : %s - [%s]' % (len(unpaired_list), ','.join(unpaired_list)))
		print '##########################################################################################\n'
	except:
		print sys.exc_info()
		logger.error('WARNING	PIPELINE_OPv2 - error\n')
