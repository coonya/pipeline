#! /usr/bin/python
import optparse, os, sys, subprocess, re, MySQLdb, glob, time, glob

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-1', dest='dir1', help='analyzed dir1 e.g) 131016_M02071_0004_000000000-A50PV')
parser.add_option('-2', dest='dir2', help='analyzed dir2 e.g) 131016_M02071_0004_000000000-A50PV')

parser.add_option('-g', '--buildver',dest='buildver', default='b37', help='reference genome version (b37, hg19). default: b37')
(options, args) = parser.parse_args()
def check_dir(check_dir):
	if os.path.exists(check_dir):
		pass
	else:
		try:
			os.makedirs(check_dir)
		except:
			pass

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def java_run(options):
	java_command = '/usr/bin/java %s -jar' % options
	return java_command

def java7_run(options):
	java_command = '/home/D131748/apps/JAVA/jre1.7.0_60/bin/java %s -jar' % options
	return java_command



def pbs_header(nodes, ppn, depend, no_wjob='', job_id=''):
	if depend == 0:
		cmd = '#PBS -l nodes=%d:ppn=%d\n#PBS -l walltime=10:00:00\n#PBS -j oe\n#PBS -q batch\ncd $PBS_O_WORKDIR\n' % (nodes, ppn)
	elif depend ==1:
		if no_wjob == 1:
			cmd = '#PBS -l nodes=%d:ppn=%d\n#PBS -l walltime=10:00:00\n#PBS -j oe\n#PBS -q batch\n#PBS -W depend=afterok:%s.master\ncd $PBS_O_WORKDIR\n' % (nodes, ppn, job_id)
		elif no_wjob == 2:
			jobs = job_id.split(',')
			jobs_order = ''
			for x in jobs:
				if jobs_order == '':
					jobs_order = '%s.master' % x
				else:
					jobs_order = jobs_order + ':%s.master' % x

			cmd = '#PBS -l nodes=%d:ppn=%d\n#PBS -l walltime=10:00:00\n#PBS -j oe\n#PBS -q batch\n#PBS -W depend=afterok:%s\ncd $PBS_O_WORKDIR\n' % (nodes, ppn, jobs_order)

	return cmd


def file_name(dir1, dir2):
	sample_dic = {}
	working_dir = '/home/D131748/Research/OncoPanel/aggregation/%s_%s' % (options.dir1.replace('/',''), options.dir2)
	metrics_dir = '%s/metrics' % working_dir
	working_tmp = '%s/tmp' % working_dir
	mutect_dir = '%s/mutect' % working_dir
	somaticindelocator_dir = '%s/somaticindelocator' % working_dir
	viscapcancer_dir = '%s/VisCapCancer' % working_dir

	bam2 = {}

	sample_file1 = '/home/D131748/Research/OncoPanel/lane_level/%s/configure/sample.txt' % (dir1)
	sample_file2 = '/home/D131748/Research/OncoPanel/lane_level/%s/configure/sample.txt' % (dir2)
	

	for x in open(sample_file1).xreadlines():
		if x[0] != '#':
			x_1 = x.strip().split('\t')
			sample = x_1[0].split('.')[0].split('_')
			sample2 = '%s_%s' % (sample[0], sample[2][:-1])
			type = sample[2][-1]

			if bam2.has_key(sample2):
				"""
				if type == 'T':
					bam2[sample2]['tumor'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', '')))
				elif type == 'N':
					bam2[sample2]['normal'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', '')))
				"""

				if type == 'T':
					if bam2[sample2].has_key('tumor'):
						bam2[sample2]['tumor'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', '')))
					else:
						bam2[sample2]['tumor'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', ''))]
				elif type == 'N':
					if bam2[sample2].has_key('normal'):
						bam2[sample2]['normal'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', '')))
					else:
						bam2[sample2]['normal'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', ''))]

			

			else:
				bam2[sample2] = {}

				if type == 'T':
					if bam2[sample2].has_key('tumor'):
						bam2[sample2]['tumor'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', '')))
					else:
						bam2[sample2]['tumor'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', ''))]
				elif type == 'N':
					if bam2[sample2].has_key('normal'):
						bam2[sample2]['normal'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', '')))
					else:
						bam2[sample2]['normal'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir1, x_1[0].replace('.bam', ''))]

	for y in open(sample_file2).xreadlines():
		if y[0] != '#':
			y_1 = y.strip().split('\t')
			sample = y_1[0].split('.')[0].split('_')
			sample2 = '%s_%s' % (sample[0], sample[2][:-1])
			type = sample[2][-1]

			if bam2.has_key(sample2):
				"""
				if type == 'T':
					bam2[sample2]['tumor'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
				elif type == 'N':
					bam2[sample2]['normal'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
				"""



				if type == 'T':
					if bam2[sample2].has_key('tumor'):
						bam2[sample2]['tumor'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
					else:
						bam2[sample2]['tumor'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', ''))]
				elif type == 'N':
					if bam2[sample2].has_key('normal'):
						bam2[sample2]['normal'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
					else:
						bam2[sample2]['normal'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', ''))]



			else:
				if bam2.has_key(sample2):

					if type == 'T':
						if bam2[sample2].has_key('tumor'):
							bam2[sample2]['tumor'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
						else:
							bam2[sample2]['tumor'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', ''))]
					elif type == 'N':
						if bam2[sample2].has_key('normal'):
							bam2[sample2]['normal'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
						else:
							bam2[sample2]['normal'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', ''))]
				else:
					bam2[sample2] = {}	
					if type == 'T':
						if bam2[sample2].has_key('tumor'):
							bam2[sample2]['tumor'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
						else:
							bam2[sample2]['tumor'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', ''))]
					elif type == 'N':
						if bam2[sample2].has_key('normal'):
							bam2[sample2]['normal'].append('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', '')))
						else:
							bam2[sample2]['normal'] = ['/home/D131748/Research/OncoPanel/lane_level/%s/tmp/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (dir2, y_1[0].replace('.bam', ''))]
	


	for l in bam2:
		sample_dic[l] = {}
		sample_dic[l]['dedup_N'] = '%s/%sN.aggregation.dedup.bam' % (working_tmp, l)
		sample_dic[l]['dedup_T'] = '%s/%sT.aggregation.dedup.bam' % (working_tmp, l)
		sample_dic[l]['dedup_metrics_N'] = '%s/%sN_duplicateMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['dedup_metrics_T'] = '%s/%sT_duplicateMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['indel_interval'] = '%s/%s_indelCocleanTargets.intervals' % (working_tmp, l)
		sample_dic[l]['map'] = '%s/%s.map' % (working_tmp, l)
		sample_dic[l]['aggregation_N'] = '%s/%sN.aggregation.dedup.cleaned.bam' % (working_tmp, l)
		sample_dic[l]['aggregation_T'] = '%s/%sT.aggregation.dedup.cleaned.bam' % (working_tmp, l)
	
		sample_dic[l]['exon_metrics_N'] = '%s/%sN_exon_targets_hsMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['intron_metrics_N'] = '%s/%sN_intron_targets_hsMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['full_metrics_N'] = '%s/%sN_full_targets_hsMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['exon_coverage_N'] = '%s/%sN_exon_targets_hsIntervals.txt' % (metrics_dir, l)
		sample_dic[l]['intron_coverage_N'] = '%s/%sN_intron_targets_hsIntervals.txt' % (metrics_dir, l)
		sample_dic[l]['full_coverage_N'] = '%s/%sN_full_targets_hsIntervals.txt' % (metrics_dir, l)
		sample_dic[l]['depthofcoverage_N'] = '%s/%sN.dedup.cleaned.cov' % (viscapcancer_dir, l)

		sample_dic[l]['exon_metrics_T'] = '%s/%sT_exon_targets_hsMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['intron_metrics_T'] = '%s/%sT_intron_targets_hsMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['full_metrics_T'] = '%s/%sT_full_targets_hsMetrics.txt' % (metrics_dir, l)
		sample_dic[l]['exon_coverage_T'] = '%s/%sT_exon_targets_hsIntervals.txt' % (metrics_dir, l)
		sample_dic[l]['intron_coverage_T'] = '%s/%sT_intron_targets_hsIntervals.txt' % (metrics_dir, l)
		sample_dic[l]['full_coverage_T'] = '%s/%sT_full_targets_hsIntervals.txt' % (metrics_dir, l)
		sample_dic[l]['depthofcoverage_T'] = '%s/%sT.dedup.cleaned.cov' % (viscapcancer_dir, l)


		sample_dic[l]['mutect_vcf'] = '%s/%s.mutect.vcf' % (mutect_dir, l)
		sample_dic[l]['filtered_mutect_vcf'] = '%s/%s.mutect.filtered.vcf' % (mutect_dir, l)
		sample_dic[l]['mutect_txt'] = '%s/%s.mutect.txt' % (mutect_dir, l)
		sample_dic[l]['mutect_wig'] = '%s/%s.mutect.coverage.wig' % (mutect_dir, l)
		sample_dic[l]['somaticindelocator'] = '%s/%s.somaticindelocator.vcf' % (somaticindelocator_dir, l)
		sample_dic[l]['filtered_somaticindelocator'] = '%s/%s.somaticindelocator.filtered.vcf' % (somaticindelocator_dir, l)
		sample_dic[l]['mutect_maf'] = '%s/%s.mutect.filtered.maf' % (mutect_dir, l)
		sample_dic[l]['somaticindelocator_maf'] = '%s/%s.somaticindelocator.filtered.maf' % (somaticindelocator_dir, l)
		
	return sample_dic, bam2

def get_start_jobid():
	find_jobid = subprocess.Popen("qmgr -c 'p s'|awk '/next_job_number/{print $5}'", shell=True, stdout=subprocess.PIPE)
	jobid = find_jobid.stdout.readline().split('\n')[0]

	return jobid

def check_files(bam_dic, dir1, dir2):
	files1 = glob.glob('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/*.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % dir1)
	files2 = glob.glob('/home/D131748/Research/OncoPanel/lane_level/%s/tmp/*.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % dir2)
	print files1
	print files2


def get_sample_list(bam_dic):
	paired_list = []
	unpaired_list = []

	for x in bam_dic:
		if bam_dic.has_key('normal'):
			paired_list.append(x)
		else:
			unpaired_list.append(x)


	return paired_list, unpaired_list

def check_dup(bam_dic):
	dup_list_tumor = []
	dup_list_normal = []
	dup_list_all = []
	matched_list = []
	tumor_only_list = []
	normal_only_list = []

	print bam_dic

	for x in bam_dic:
		if bam_dic[x].has_key('normal'):
			if len(bam_dic[x]['normal']) >= 2:
				dup_list_normal.append(x)
				if x not in dup_list_all:
					dup_list_all.append(x)
		if bam_dic[x].has_key('tumor'):
			if len(bam_dic[x]['tumor']) >= 2:
				dup_list_tumor.append(x)
				if x not in dup_list_all:
					dup_list_all.append(x)
	for y in dup_list_normal:
		tumor = '%sT' % y[:-1]
		if tumor in dup_list_tumor:
			matched_list.append(y[:-1])
		else:
			normal_only_list.append(y)
	for z in dup_list_tumor:
		normal = '%sN' % z[:-1]
		if normal not in dup_list_normal:
			tumor_only_list.append(z)
	
	return dup_list_normal, dup_list_tumor, dup_list_all, matched_list, tumor_only_list, normal_only_list
		

class picard:
	def __init__(self):
#		self.basecalls_dir = 'BASECALLS_DIR=%s/%sData/Intensities/BaseCalls' % (cwd, options.basecall)
		self.tmp_dir = 'TMP_DIR=/data/scratch'
		self.read_structure = 'READ_STRUCTURE=151T6B151T'
		self.lane = 'LANE=1'
		self.validation_stringency = 'VALIDATION_STRINGENCY=SILENT'
		self.compression_level = 'COMPRESSION_LEVEL=1'
		self.picard_path = '/home/D131748/apps/ETC/Picard/picard-tools-1.92'
		self.ref = '/home/D131748//Reference/Human/b37/human_g1k_v37.fasta'

	"""
	def dedup(self, sample1, sample2, output, metrics_file):
		cmd = []
		
		cmd.append(java_run('-Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G -Xmx4G'))
		cmd.append('%s/MarkDuplicates.jar' % self.picard_path)
		cmd.append('INPUT=%s' % sample1)
		cmd.append('INPUT=%s' % sample2)
		cmd.append('OUTPUT=%s' % output)
		cmd.append('METRICS_FILE=%s' % metrics_file)
		cmd.append(self.tmp_dir)
		cmd.append('MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200')
		cmd.append('CREATE_INDEX=true')
		cmd.append('CREATE_MD5_FILE=false')
		cmd.append(self.compression_level)
		cmd.append(self.validation_stringency)

		exec_command = ' '.join(cmd)

		return exec_command
	"""

	def dedup(self, sample_list, output, metrics_file):
		cmd = []
		
		cmd.append(java_run('-Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G -Xmx4G'))
		cmd.append('%s/MarkDuplicates.jar' % self.picard_path)
		for x in sample_list:
			cmd.append('INPUT=%s' % x)
		cmd.append('OUTPUT=%s' % output)
		cmd.append('METRICS_FILE=%s' % metrics_file)
		cmd.append(self.tmp_dir)
		cmd.append('MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200')
		cmd.append('CREATE_INDEX=true')
		cmd.append('CREATE_MD5_FILE=false')
		cmd.append(self.compression_level)
		cmd.append(self.validation_stringency)

		exec_command = ' '.join(cmd)

		return exec_command


	def CalculateHsMetrics(self, sample, output, coverage, region):
		cmd = []

		bait_interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/Oncopanel_v2_baits_1kGenomeB37.interval_list'

		if region == 'exon':
			interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/OPv2_Exons.interval_list'
		elif region == 'full':
			interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/OPv2_Exon_Intron_GWAS_SNP.merged.interval_list'
		elif region == 'intron':
			interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/OPv2_Introns.interval_list'

		cmd.append(java_run('-Xmx4G -Xms2G'))
		cmd.append('%s/CalculateHsMetrics.jar' % self.picard_path)
		cmd.append('INPUT=%s' % sample)
		cmd.append(self.tmp_dir)
		cmd.append('BAIT_INTERVALS=%s' % bait_interval)
		cmd.append('OUTPUT=%s' % output)
		cmd.append('PER_TARGET_COVERAGE=%s' % coverage)
		cmd.append('TARGET_INTERVALS=%s' % interval)
		cmd.append('REFERENCE_SEQUENCE=%s' % self.ref)
		exec_command = ' '.join(cmd)

		return exec_command


class gatk:
	def __init__(self):
		self.gatk_path = '/home/D131748/apps/ETC/GATK/GenomeAnalysisTK-1.6-5-g557da77'
		
		self.interval = '/home/D131748//Reference/GATK_bundle/2.5/b37/known_indel.intervals'
		self.known1 = '/home/D131748//Reference/GATK_bundle/2.5/b37/Mills_and_1000G_gold_standard.indels.b37.vcf'
		self.known2 = '/home/D131748//Reference/GATK_bundle/2.5/b37/1000G_phase1.indels.b37.vcf'
		self.dbsnp = '/home/D131748//Reference/GATK_bundle/2.5/b37/dbsnp_137.b37.vcf'
		self.ref = '/home/D131748//Reference/Human/b37/human_g1k_v37.fasta'
		self.bait_interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/130823_OPv2_All_True_Targets.no_overlaps.interval_list'

	def RealignerTargetCreator(self, indel_interval, tumor, normal=''):
		cmd = []
		
		cmd.append(java_run('-XX:GCHeapFreeLimit=10 -XX:GCTimeLimit=50 -Djava.io.tmpdir=/data/scratch -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T RealignerTargetCreator' % gatk_path)

		if normal != '':
			cmd.append('-I %s -I %s' % (normal, tumor))

		else:
			cmd.append('-I %s' % tumor)
		

		cmd.append('-o %s' % indel_interval)
		cmd.append('-R %s' % self.ref)
		cmd.append('-nt 1')
		cmd.append('-dcov 250')
		exec_command = ' '.join(cmd)

		return exec_command
	

	def IndelRealigner_agg(self, map, file_dic, indel_interval, tumor, normal=''): 
		cmd = []
		
		wname = map
		wfile = open(wname, 'w')
		wfile2 = open('%s/breakmer_samples.txt' % config_dir, 'a')

## generate map file
		cmd.append(java_run('-Djava.io.tmpdir=/data/scratch -Dsamjdk.use_async_io=true -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T IndelRealigner' % self.gatk_path)

		if normal != '':
			normal_sample = os.path.basename(normal).split('.')[0]
			tumor_sample = os.path.basename(tumor).split('.')[0]
			sample_id = tumor_sample[:-1]

			wfile2.write('%s\t%s\n%s\t%s\n' % (normal_sample, normal, tumor_sample, tumor))

			cmd.append('-I %s -I %s' % (normal, tumor))
			wfile.write('%s\t%s\n' % (os.path.basename(normal), file_dic[sample_id]['aggregation_N']))
			wfile.write('%s\t%s\n' % (os.path.basename(tumor), file_dic[sample_id]['aggregation_T']))
			wfile.close()
			
			cmd.append('--nWayOut %s' % wname)

		else:
			tumor_sample = os.path.basename(tumor).split('.')[0]
			print tumor_sample
			sample_id = tumor_sample[:-1]
			wfile2.write('%s\t%s\n' % (tumor_sample, tumor))

			cmd.append('-I %s' % tumor)
			cmd.append('-o %s' % file_dic[sample_id]['aggregation_T'])
	
		cmd.append('-targetIntervals %s' % indel_interval)
		cmd.append('-R %s' % self.ref)
		cmd.append('-maxInMemory 1000000')
		cmd.append('-model USE_READS')
		exec_command = ' '.join(cmd)

		wfile2.close()
		return exec_command


	def depthofcoverage(self, input, output):
		cmd = []

		cmd.append(java_run('-Xmx4G -Djava.io.tmpdir=/data/scratch'))
		cmd.append('%s/GenomeAnalysisTK.jar -T DepthOfCoverage' % self.gatk_path) 
		cmd.append('--omitDepthOutputAtEachBase')
		cmd.append('--interval_merging OVERLAPPING_ONLY')
		cmd.append('-L %s'% self.bait_interval)
		cmd.append('-I %s' % input)
		cmd.append('-o %s' % output)
		cmd.append('-R %s' % self.ref)
		exec_command = ' '.join(cmd)

		return exec_command


class somatic:
	def __init__(self):
		self.gatk_path = '/home/D131748/apps/ETC/GATK/GenomeAnalysisTK-1.6-5-g557da77'

		self.ref = '/home/D131748//Reference/Human/b37/human_g1k_v37.fasta'
		self.bait_interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/Oncopanel_v2_baits_1kGenomeB37.interval_list'

		self.dbsnp = '/home/D131748//Reference/GATK_bundle/2.5/b37/dbsnp_137.b37.vcf'
		self.dbsnp_common = '/home/D131748//Reference/GATK_bundle/2.5/b37/dbsnp_137_common.b37.vcf'
#		self.cosmic = '/home/D131748//Reference/GATK_bundle/2.5/b37/cosmic_v54_120711_b37.vcf'
		self.cosmic = '/home/D131748/Reference/cosmic/cosmic_v71_snv.b37.vcf'


	def mutect(self, type, file_dic, output_txt, output_vcf, tumor, output_wig='', normal= ''):
		cmd = []

		tumor_sample = os.path.basename(tumor)

#		mutect_path = '/home/D131748/apps/Cancer/Mutect/Mutect-1.1.5'
		mutect_path = '/home/D131748/apps/Cancer/Mutect/Mutect-1.1.7'

		cmd.append(java7_run('-Xmx4G -Djava.io.tmpdir=/data/scratch'))
		cmd.append('%s/mutect-1.1.7.jar --analysis_type MuTect' % mutect_path)
		cmd.append('--reference_sequence %s' % self.ref)
#		cmd.append('--dbsnp %s' % self.dbsnp_common)
#		cmd.append('--cosmic %s' % self.cosmic)
		cmd.append('--intervals %s' % self.bait_interval)

		if normal != '':
			cmd.append('--input_file:normal %s' % normal)
		
		cmd.append('--input_file:tumor %s' % tumor)
		
		if type == 0:
#			if normal != '':
#				PON = '/home/D131748/Reference/PON/OPv2_mutect_PON_20140624.vcf'
#				cmd.append('--normal_panel %s' % PON)
#			else:
#				PON = '/home/D131748/Reference/PON/OPv2_mutect_PON_tumoronly_20140624.vcf'
#				cmd.append('--normal_panel %s' % PON)
			
			cmd.append('--out %s' % output_txt)
			cmd.append('--vcf %s' % output_vcf)
			cmd.append('--coverage_file %s' % output_wig)
		elif type == 1:
			cmd.append('--out %s' % output_txt)
			cmd.append('--vcf %s' % output_vcf)

		cmd.append('--enable_extended_output')
		exec_command = ' '.join(cmd)

		if type == 0:
			if normal != '':
				cmd15_1 = 'grep -v REJECT %s > %s\ngrep -v REJECT %s > %s\n' % (output_vcf, output_vcf.replace('.vcf', '.filtered.vcf'), output_txt, output_txt.replace('.txt', '.filtered.txt'))

			else:
				cmd15_1 = 'grep -v REJECT %s |grep -v \"rs\" > %s\ngrep -v REJECT %s > %s\n' % (output_vcf, output_vcf.replace('.vcf', '.filtered.vcf'), output_txt, output_txt.replace('.txt', '.filtered.txt'))

			exec_command = '%s\n%s' % (exec_command, cmd15_1)


		return exec_command



	def somaticindelocator(self, t_type, output, tumor, normal=''):
		cmd = []

		src_path = '/home/D131748//src'
		
		if normal != '':
			cmd = 'java -Xmx4G -jar %s/GenomeAnalysisTK.jar -T SomaticIndelDetector -R %s -I:normal %s -I:tumor %s -o %s -L %s --window_size 500 -filter T_COV\<2\|\|N_COV\<0\|\|T_INDEL_F\<0.05\|\|T_INDEL_CF\<0.1' % (self.gatk_path, self.ref, normal, tumor, output, self.bait_interval)
			cmd2 = '01_somaticindel_filter.py -i %s' % output
			
			exec_command = '%s\n%s' % (cmd, cmd2)


		else:
			if t_type == 0:
				cmd = 'java -Xmx4G -jar %s/GenomeAnalysisTK.jar -T SomaticIndelDetector -R /%s --unpaired -I:tumor %s -o %s -L %s --window_size 500 -filter COV\<2\|\|INDEL_F\<0.05\|\|INDEL_CF\<0.1' % (self.gatk_path, self.ref, tumor, output, self.bait_interval)
				cmd2 = '01_somaticindel_filter.py -i %s' % (output)
		
				exec_command = '%s\n%s' % (cmd, cmd2)


			elif t_type == 1:
				cmd = 'java -Xmx4G -jar %s/GenomeAnalysisTK.jar -T SomaticIndelDetector -R /%s --unpaired -I:tumor %s -o %s -L %s --window_size 500 -filter COV\<2\|\|INDEL_F\<0.05\|\|INDEL_CF\<0.1' % (self.gatk_path, self.ref, tumor, output.replace('.vcf', '.nonfilter.vcf'), self.bait_interval)
				exec_command = cmd

		return exec_command



class sv:
	def k_mer(self, project):
		cmd = []

		sample_file = '%s/kmer_sample.txt' % sv_out_dir

		cmd = 'sh %s/setup_analysis.sh %s %s %s %s/kmer_region.op2.config' % (k_mer_path, sv_out_dir, options.basecall.replace('/',''), sample_file, k_mer_path)


		return cmd


class viscapcancer:
	def viscapcancer(self, viscapcancer_dir):
		cmd = []
		script = '/home/D131748/apps/CCGD/VisCapCancer/VisCapCancer.R'
		cfg = '/home/D131748/apps/CCGD/VisCapCancer/VisCapCancer.CCGD.OPv2.cfg'

		cmd.append('Rscript %s %s %s %s' % (script, cfg, viscapcancer_dir, viscapcancer_dir))
		exec_command = ' '.join(cmd)


		return exec_command


class annotation:

	def vcf2maf(self, vcf, maf, t_type):
		cmd = []
		new_id = os.path.basename(vcf).split('.')[0]
#		patient = sample.split('_')[2]
#		cpcm = sample.split('_')[0]
#		new_id = '%s_%s' % (cpcm,patient)

		if t_type == 0:
			cmd.append('vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s --normal-id %s' % (vcf, maf, new_id, new_id.replace('T', 'N')))

		elif t_type == 1:		
			cmd.append('vcf2maf.pl --input-vcf %s --output-maf %s --tumor-id %s --normal-id none' % (vcf, maf, new_id))

		exec_command = '\n'.join(cmd)

		return exec_command


	def filter_exonic(self, input, caller):

		#cmd = '01_filter_exonic.py -i %s -c %s' % (input, caller)
		#cmd = '01_add_annotation.py -i %s -c %s' % (input, caller)
		cmd = '01_add_annotation.py -i %s' % (input)

		return cmd


#def __main__():

#################### check directories #######################
check_option(options.dir1, "You need to use '-1' option with basecall dir\n")
check_option(options.dir2, "You need to use '-2' option with basecall dir\n")


##################### directories ###########################	
cwd = os.getcwd()
working_dir = '/home/D131748/Research/OncoPanel/aggregation/%s_%s' % (options.dir1.replace('/',''), options.dir2)
config_dir = '%s/configure' % working_dir
metrics_dir = '%s/metrics' % working_dir
log_dir = '%s/log' % working_dir
working_tmp = '%s/tmp' % working_dir
mutect_dir = '%s/mutect' % working_dir
viscapcancer_dir = '%s/VisCapCancer' % working_dir
somaticindelocator_dir = '%s/somaticindelocator' % working_dir
picard_path = '/home/D131748/apps/ETC/Picard/picard-tools-1.92'
gatk_path = '/home/D131748/apps/ETC/GATK/GenomeAnalysisTK-1.6-5-g557da77'
mutect_path = '/home/D131748/apps/Cancer/Mutect/Mutect-1.1.4'


###################### check directory ######################
check_dir(working_dir)
check_dir(log_dir)
check_dir(metrics_dir)
check_dir(config_dir)
check_dir(working_tmp)
check_dir(mutect_dir)
check_dir(somaticindelocator_dir)
check_dir(viscapcancer_dir)

now = time.strftime('%H_%M_%S')

##################### reference genome version ##################
if options.buildver == 'hg19':
	ref = '~/Reference/Human/Genome/hg19/ucsc.hg19.fasta'
elif options.buildver == 'b37':
	ref = '/home/D131748//Reference/Human/b37/human_g1k_v37.fasta'
	dbsnp = '/home/D131748//Reference/GATK_bundle/2.5/b37/dbsnp_137.b37.vcf'
	dbsnp_common = '/home/D131748//Reference/GATK_bundle/2.5/b37/dbsnp_137_common.b37.vcf'
#	cosmic = '/home/D131748//Reference/GATK_bundle/2.5/b37/cosmic_v54_120711_b37.vcf'
	cosmic = '/home/D131748/Reference/cosmic/cosmic_v71_snv.b37.vcf'
	interval = '/home/D131748//Reference/GATK_bundle/2.5/b37/known_indel.intervals'
	known1 = '/home/D131748//Reference/GATK_bundle/2.5/b37/Mills_and_1000G_gold_standard.indels.b37.vcf'
	known2 = '/home/D131748//Reference/GATK_bundle/2.5/b37/1000G_phase1.indels.b37.vcf'
	bait_interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/Oncopanel_v2_baits_1kGenomeB37.interval_list'
	exon_interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/OPv2_Exons.interval_list'
	exon_intron_interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/OPv2_Exon_Intron_GWAS_SNP.merged.interval_list'
	intron_interval = '/home/D131748//Reference/Target_intervals/OncoPanel_v2/b37/OPv2_Introns.interval_list'
else:
	sys.exit('wrong reference genome version')


os.chdir(log_dir)
log_files = glob.glob('*')
if len(log_files) != 0:
	for log_file in log_files:
		os.remove(log_file)

#################### main #######################################
qsub_id = {}
qsub_dic = {}


duke1 = picard()	
duke2 = gatk()
duke3 = somatic()
duke4 = annotation()
duke5 = viscapcancer()


(file_dic, samples) = file_name(options.dir1, options.dir2)

#check_files(file_dic, options.dir1, options.dir2)

(paired_list, unpaired_list) = get_sample_list(samples)
(dup_list_normal, dup_list_tumor, dup_list_all, matched_list, tumor_only_list, normal_only_list) = check_dup(samples)
jobid = int(get_start_jobid())



### deduplication

if len(dup_list_normal) != 0:
	for z in dup_list_normal:
		sample = z

		wname = '01_dedup_%sN.sh' % sample
		wfile = open(wname, 'w')
		wfile.write('%s\n' % pbs_header(1, 1, 0))
	
		qsub_dic[sample] = {}
		qsub_dic[sample]['dedup_N'] = duke1.dedup(samples[z]['normal'], file_dic[sample]['dedup_N'], file_dic[sample]['dedup_metrics_N'])
		wfile.write(qsub_dic[sample]['dedup_N'])
		wfile.close()
		q_id = 'dedup_%sN' % sample
		qsub_id[q_id] = jobid
		jobid += 1
		subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)

if len(dup_list_tumor) != 0:
	for z in dup_list_tumor:
		sample = z

		wname = '01_dedup_%sT.sh' % sample
		wfile = open(wname, 'w')
		wfile.write('%s\n' % pbs_header(1, 1, 0))
	
		qsub_dic[sample] = {}
		qsub_dic[sample]['dedup_T'] = duke1.dedup(samples[z]['tumor'], file_dic[sample]['dedup_T'], file_dic[sample]['dedup_metrics_T'])
		wfile.write(qsub_dic[sample]['dedup_T'])
		wfile.close()
		q_id = 'dedup_%sT' % sample
		qsub_id[q_id] = jobid
		jobid += 1
		subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)

print dup_list_all
if len(dup_list_all) != 0:
	for k in dup_list_all:

		q_id = 'aggregation_%s' % k
		qsub_id[q_id] = jobid
		wname = '02_aggregation_%s.sh' % k
		wfile = open(wname, 'w')

		normal_sample = file_dic[k]['dedup_N']
		tumor_sample = file_dic[k]['dedup_T']

		indel_interval = file_dic[k]['indel_interval']
		map_file = file_dic[k]['map']

		aggregation_T = file_dic[k]['aggregation_T']
		exon_metrics_T = file_dic[k]['exon_metrics_T']
		exon_coverage_T = file_dic[k]['exon_coverage_T']
		intron_metrics_T = file_dic[k]['intron_metrics_T']
		intron_coverage_T = file_dic[k]['intron_coverage_T']
		full_metrics_T = file_dic[k]['full_metrics_T']
		full_coverage_T = file_dic[k]['full_coverage_T']

		aggregation_N = file_dic[k]['aggregation_T']
		exon_metrics_N = file_dic[k]['exon_metrics_N']
		exon_coverage_N = file_dic[k]['exon_coverage_N']
		intron_metrics_N = file_dic[k]['intron_metrics_N']
		intron_coverage_N = file_dic[k]['intron_coverage_N']
		full_metrics_N = file_dic[k]['full_metrics_N']
		full_coverage_N = file_dic[k]['full_coverage_N']

	
		if k in dup_list_normal:
			if k in normal_only_list:
				qsub_dic[k] = {}
				qsub_dic[k]['aggregation'] = [duke2.RealignerTargetCreator(indel_interval, tumor_sample, normal_sample)]
				qsub_dic[k]['aggregation'].append(duke2.IndelRealigner_agg(map_file, file_dic, indel_interval, tumor_sample, normal_sample))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_N, exon_metrics_N, exon_coverage_N, 'exon'))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_N, intron_metrics_N, intron_coverage_N, 'intron'))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_N, full_metrics_N, full_coverage_N, 'full'))
	#			qsub_dic[k]['aggregation'].append(duke2.depthofcoverage(file_dic[normal_sample]['aggregation'], file_dic[normal_sample]['depthofcoverage']))
	#			qsub_dic[k]['aggregation'].append(duke2.depthofcoverage(file_dic[tumor_sample]['aggregation'], file_dic[tumor_sample]['depthofcoverage']))
	
				cmd = '\n'.join(qsub_dic[k]['aggregation'])
	
				wfile.write('%s\n' % pbs_header(1, 1, 1, 1, '%s' % (qsub_id['dedup_%sN' % k])))
				wfile.write(cmd)
				wfile.close()
				subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)
	
			else:

				qsub_dic[k] = {}
				qsub_dic[k]['aggregation'] = [duke2.RealignerTargetCreator(indel_interval, tumor_sample, normal_sample)]
				qsub_dic[k]['aggregation'].append(duke2.IndelRealigner_agg(map_file, file_dic, indel_interval, tumor_sample, normal_sample))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_T, exon_metrics_T, exon_coverage_T, 'exon'))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_T, intron_metrics_T, intron_coverage_T, 'intron'))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_T, full_metrics_T, full_coverage_T, 'full'))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_N, exon_metrics_N, exon_coverage_N, 'exon'))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_N, intron_metrics_N, intron_coverage_N, 'intron'))
				qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_N, full_metrics_N, full_coverage_N, 'full'))
	#			qsub_dic[k]['aggregation'].append(duke2.depthofcoverage(file_dic[normal_sample]['aggregation'], file_dic[normal_sample]['depthofcoverage']))
	#			qsub_dic[k]['aggregation'].append(duke2.depthofcoverage(file_dic[tumor_sample]['aggregation'], file_dic[tumor_sample]['depthofcoverage']))
	
				cmd = '\n'.join(qsub_dic[k]['aggregation'])
	
				wfile.write('%s\n' % pbs_header(1, 1, 1, 2, '%s,%s' % (qsub_id['dedup_%sT' % k], qsub_id['dedup_%sN' % k])))
				wfile.write(cmd)
				wfile.close()
				subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)
	
		else:
			qsub_dic[k] = {}
			qsub_dic[k]['aggregation'] = [duke2.RealignerTargetCreator(indel_interval, tumor_sample)]
			qsub_dic[k]['aggregation'].append(duke2.IndelRealigner_agg(map_file, file_dic, indel_interval, tumor_sample))
			qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_T, exon_metrics_T, exon_coverage_T, 'exon'))
			qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_T, intron_metrics_T, intron_coverage_T, 'intron'))
			qsub_dic[k]['aggregation'].append(duke1.CalculateHsMetrics(aggregation_T, full_metrics_T, full_coverage_T, 'full'))
#			qsub_dic[k]['aggregation'].append(duke2.depthofcoverage(aggregation, file_dic[k]['depthofcoverage']))

			cmd = '\n'.join(qsub_dic[k]['aggregation'])

			wfile.write('%s\n' % pbs_header(1, 1, 1, 1, '%s' % (qsub_id['dedup_%sT' % k])))
			wfile.write(cmd)
			wfile.close()
			subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)

		jobid += 1

#################### Calling Mutation step ###############################
	qsub_id_calling_mutation = []

	for k in dup_list_all:

		q_id = 'CallingMutation_%s' % k
		qsub_id[q_id] = jobid
		wname = '03_CallingMutation_%s.sh' % k
		wfile = open(wname, 'w')
		qsub_id_calling_mutation.append(str(jobid))


		mutect_txt = file_dic[k]['mutect_txt']
		mutect_vcf = file_dic[k]['mutect_vcf']
		aggregation_T = file_dic[k]['aggregation_T']
		aggregation_N = file_dic[k]['aggregation_N']
		mutect_wig = file_dic[k]['mutect_wig']
		filtered_mutect_vcf = file_dic[k]['filtered_mutect_vcf']
		mutect_maf = file_dic[k]['mutect_maf']

		somaticindelocator = file_dic[k]['somaticindelocator']
		filtered_somaticindelocator = file_dic[k]['filtered_somaticindelocator']
		somaticindelocator_maf = file_dic[k]['somaticindelocator_maf']


		if k in dup_list_normal:
	
			qsub_dic[k] = {}
			qsub_dic[k]['CallingMutation'] = [duke3.mutect(0, file_dic, mutect_txt, mutect_vcf, aggregation_T, mutect_wig, aggregation_N)]
			qsub_dic[k]['CallingMutation'].append(duke3.somaticindelocator(0, somaticindelocator, aggregation_T, aggregation_N))
	
			qsub_dic[k]['CallingMutation'].append(duke4.vcf2maf(filtered_mutect_vcf, mutect_maf, 0))        ## vcf, maf, t_type
			qsub_dic[k]['CallingMutation'].append(duke4.vcf2maf(filtered_somaticindelocator, somaticindelocator_maf, 0))        ## vcf, maf, t_type
			
			qsub_dic[k]['CallingMutation'].append(duke4.filter_exonic(mutect_maf, 'mutect'))  ## vcf, maf, t_type
			qsub_dic[k]['CallingMutation'].append(duke4.filter_exonic(somaticindelocator_maf, 'somaticindelocator'))  ## vcf, maf, t_type
		
		else:
			qsub_dic[k] = {}
			qsub_dic[k]['CallingMutation'] = [duke3.mutect(0, file_dic, mutect_txt, mutect_vcf, aggregation_T, mutect_wig)]
			qsub_dic[k]['CallingMutation'].append(duke3.somaticindelocator(0, somaticindelocator, aggregation_T))
	
			qsub_dic[k]['CallingMutation'].append(duke4.vcf2maf(filtered_mutect_vcf, mutect_maf, 1))        ## vcf, maf, t_type
			qsub_dic[k]['CallingMutation'].append(duke4.vcf2maf(filtered_somaticindelocator, somaticindelocator_maf, 1))        ## vcf, maf, t_type
			
			qsub_dic[k]['CallingMutation'].append(duke4.filter_exonic(mutect_maf, 'mutect'))  ## vcf, maf, t_type
			qsub_dic[k]['CallingMutation'].append(duke4.filter_exonic(somaticindelocator_maf, 'somaticindelocator'))  ## vcf, maf, t_type
		
		jobid += 1
	
		cmd = '\n'.join(qsub_dic[k]['CallingMutation'])
	
		wfile.write('%s\n' % pbs_header(1, 1, 1, 1, '%s' % (qsub_id['aggregation_%s' % k])))
		wfile.write(cmd)
		wfile.close()
		subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)

if len(qsub_id_calling_mutation) > 1:
	hold_id = ','.join(qsub_id_calling_mutation)
else:
	hold_id = qsub_id_calling_mutation[0]

"""
wname = '04_CallingCNVs.sh'
wfile = open(wname,'w')
wfile.write('%s\n' % pbs_header(1, 1, 1, 2, hold_id))
wfile.write(duke5.viscapcancer(viscapcancer_dir))
wfile.close()
subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)
"""
###### R10 Generating report ###############
wname = '05_Generating_report.sh'
wfile = open(wname,'w')

if len(qsub_id_calling_mutation) > 1:
	wfile.write('%s\n' % pbs_header(1, 1, 1, 2, hold_id))
else:
	wfile.write('%s\n' % pbs_header(1, 1, 1, 1, hold_id))

wfile.write('02_generating_excel-2.py -d %s' % working_dir)
wfile.close()
subprocess.call('qsub %s/%s' % (log_dir, wname), shell=True, stdout=subprocess.PIPE)
