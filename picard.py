import time
from utils import *

class picard:
	def __init__(self, params, basecall, file_dic):
		self.basecall = basecall
		self.basecalls_dir = 'BASECALLS_DIR=/home/D131748/Illumina/%s/Data/Intensities/BaseCalls' % basecall
		self.tmp_dir = 'TMP_DIR=/data/scratch'
		self.read_structure = 'READ_STRUCTURE=151T6B151T'
		self.lane = 'LANE=1'
		self.validation_stringency = 'VALIDATION_STRINGENCY=SILENT'
		self.compression_level = 'COMPRESSION_LEVEL=1'

		self.working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % basecall
		self.metrics_dir = '%s/metrics' % self.working_dir
		self.config_dir = '%s/configure' % self.working_dir
		self.barcode_file = '%s/barcode.txt' % self.config_dir
		self.sample_file = '%s/sample.txt' % self.config_dir
		self.param_file = '%s/param.txt' % self.config_dir

		self.picard_path = params['picard_path']
		self.ref = params['ref']
		self.file_dic = file_dic
		
		self.bait_interval = params['target_interval']
		self.exon_interval = params['exon_interval']
		self.intron_interval = params['intron_interval']
		self.full_interval = params['full_interval']

		self.today = time.strftime('%Y/%m/%d')


	def ExtractIlluminaBarcodes(self):
		cmd = []
		
		cmd.append(java_run('-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m'))
		cmd.append('%s/ExtractIlluminaBarcodes.jar' % self.picard_path)
		cmd.append(self.basecalls_dir)
		cmd.append(self.lane)
		cmd.append(self.read_structure)
		cmd.append('BARCODE_FILE=%s' % self.barcode_file)
		cmd.append('METRICS_FILE=%s/barcodeMetrics.txt' % self.metrics_dir)
		cmd.append(self.tmp_dir)
		cmd.append('NUM_PROCESSORS=16')
		cmd.append('MAX_MISMATCHES=1')
		cmd.append('MIN_MISMATCH_DELTA=1')
		cmd.append('MAX_NO_CALLS=2')
		cmd.append('MINIMUM_BASE_QUALITY=0')
		cmd.append('COMPRESS_OUTPUTS=true')
		exec_command = ' '.join(cmd)

		return exec_command


	def IlluminaBasecallsToSam(self):
		cmd = []
		
		run_barcode = self.basecall.split('_')[3].split('/')[0]

		cmd.append(java_run('-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx10G'))
		cmd.append('%s/IlluminaBasecallsToSam.jar' % self.picard_path)
		cmd.append(self.basecalls_dir)
		cmd.append(self.lane)
		cmd.append('SEQUENCING_CENTER=ASAN-CCGD')
		cmd.append('RUN_BARCODE=%s' % run_barcode)
		cmd.append('RUN_START_DATE=%s' % self.today)
		cmd.append('NUM_PROCESSORS=16')
		cmd.append('LIBRARY_PARAMS=%s' % self.param_file)
		cmd.append('READ_GROUP_ID=%s' % run_barcode)
		cmd.append(self.read_structure)
		cmd.append('ADAPTERS_TO_CHECK=INDEXED')
		cmd.append(self.compression_level)
		exec_command = ' '.join(cmd)
		
		return exec_command


	def SamToFastq(self, sample):
		bam = self.file_dic[sample]['bam']
		fastq1 = self.file_dic[sample]['fastq1']
		fastq2 = self.file_dic[sample]['fastq2']

		cmd = []

		cmd.append(java_run('-Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx256m -Xms256m'))
		cmd.append('%s/SamToFastq.jar' % self.picard_path)
		cmd.append('INPUT=%s' % bam)
		cmd.append('FASTQ=%s' % fastq1)
		cmd.append(self.tmp_dir)
		cmd.append('SECOND_END_FASTQ=%s' % fastq2)
		cmd.append('RE_REVERSE=true')
		cmd.append('INCLUDE_NON_PF_READS=true')
		cmd.append('CLIPPING_ATTRIBUTE=XT')
		cmd.append('CLIPPING_ACTION=2')
		exec_command  = ' '.join(cmd)

		return exec_command

	
	def MergeBamAlignment(self, sample):
		bam = self.file_dic[sample]['bam']
		sam = self.file_dic[sample]['sam']
		output = self.file_dic[sample]['mergebamalignment']

		cmd = []

		cmd.append(java_run('-Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms5G  -Xmx5G'))
		cmd.append('%s/MergeBamAlignment.jar' % self.picard_path)
		cmd.append('UNMAPPED_BAM=%s' % bam)
		cmd.append('ALIGNED_BAM=%s' % sam)
		cmd.append('OUTPUT=%s' % output)
		cmd.append('REFERENCE_SEQUENCE=%s' % self.ref)
		cmd.append(self.tmp_dir)
		cmd.append('PAIRED_RUN=true')
		cmd.append('MAX_INSERTIONS_OR_DELETIONS=-1')
		cmd.append(self.compression_level)
		cmd.append('ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=X1 ATTRIBUTES_TO_RETAIN=XM ATTRIBUTES_TO_RETAIN=XO ATTRIBUTES_TO_RETAIN=XG')
		cmd.append('EXPECTED_ORIENTATIONS=FR')
		cmd.append('MAX_RECORDS_IN_RAM=3000000')
		cmd.append(self.validation_stringency)
		cmd.append('CLIP_ADAPTERS=false')
		cmd.append('IS_BISULFITE_SEQUENCE=false')
		cmd.append('ALIGNED_READS_ONLY=false')
		cmd.append('PROGRAM_RECORD_ID=bwa')
		cmd.append('PROGRAM_GROUP_VERSION=0.5.9')
		cmd.append('PROGRAM_GROUP_NAME=bwa')
		cmd.append('PROGRAM_GROUP_COMMAND_LINE=.')
		exec_command = ' '.join(cmd)
		
		return exec_command


	def dedup(self, sample):
		input = self.file_dic[sample]['mergebamalignment']
		output = self.file_dic[sample]['dedup']
		metrics_file = self.file_dic[sample]['dedup_metrics']

		cmd = []

		cmd.append(java_run('-Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2G  -Xmx4G'))
		cmd.append('%s/MarkDuplicates.jar' % self.picard_path)
		cmd.append('INPUT=%s' % input)
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


	def CalculateHsMetrics(self, sample, region):
		input = self.file_dic[sample]['realign_agg']
		cmd = []

		bait_interval = self.bait_interval

		if region == 'exon':
			interval = self.exon_interval
			output = self.file_dic[sample]['exon_metrics']
			coverage = self.file_dic[sample]['exon_coverage']
		elif region == 'full':
			interval = self.full_interval
			output = self.file_dic[sample]['full_metrics']
			coverage = self.file_dic[sample]['full_coverage']
		elif region == 'intron':
			interval = self.intron_interval
			output = self.file_dic[sample]['intron_metrics']
			coverage = self.file_dic[sample]['intron_coverage']

		cmd.append(java_run('-Xmx4G -Xms2G'))
		cmd.append('%s/CalculateHsMetrics.jar' % self.picard_path)
		cmd.append('INPUT=%s' % input)
		cmd.append(self.tmp_dir)
		cmd.append('BAIT_INTERVALS=%s' % bait_interval)
		cmd.append('OUTPUT=%s' % output)
		cmd.append('PER_TARGET_COVERAGE=%s' % coverage)
		cmd.append('TARGET_INTERVALS=%s' % interval)
		cmd.append('REFERENCE_SEQUENCE=%s' % self.ref)
		exec_command = ' '.join(cmd)

		return exec_command

	### Aggregation step

class picard_agg:
	def __init__(self, project, params):
		self.tmp_dir = 'TMP_DIR=/data/scratch'
		self.validation_stringency = 'VALIDATION_STRINGENCY=SILENT'
		self.compression_level = 'COMPRESSION_LEVEL=1'

		self.working_dir = '/home/D131748/Research/OncoPanel/aggregation/Projects/%s' % project
		self.metrics_dir = '%s/metrics' % self.working_dir

		self.picard_path = params['picard_path']
	
		self.ref = params['ref']
		self.bait_interval = params['target_interval']
		self.exon_interval = params['exon_interval']
		self.intron_interval = params['intron_interval']
		self.full_interval = params['full_interval']


	def dedup_agg(self, sample_id, sample_list):
		output = '%s/tmp/%s.agg.dedup.bam' % (self.working_dir, sample_id)
		metrics_file = '%s/%s_duplicateMetrics.txt' % (self.metrics_dir, sample_id)

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

	def CalculateHsMetrics(self, sample, region):
		input = '%s/tmp/%s.agg.dedup.cleaned.bam' % (self.working_dir, sample)

		cmd = []

		bait_interval = self.bait_interval

		if region == 'exon':
			interval = self.exon_interval
			output = '%s/%s_exon_targets_hsMetrics.txt' % (self.metrics_dir, sample)
			coverage = '%s/%s_exon_targets_hsIntervals.txt' % (self.metrics_dir, sample)
		elif region == 'full':
			interval = self.full_interval
			output = '%s/%s_full_targets_hsMetrics.txt' % (self.metrics_dir, sample)
			coverage = '%s/%s_full_targets_hsIntervals.txt' % (self.metrics_dir, sample)
		elif region == 'intron':
			interval = self.intron_interval
			output = '%s/%s_intron_targets_hsMetrics.txt' % (self.metrics_dir, sample)
			coverage = '%s/%s_intron_targets_hsIntervals.txt' % (self.metrics_dir, sample)

		cmd.append(java_run('-Xmx4G -Xms2G'))
		cmd.append('%s/CalculateHsMetrics.jar' % self.picard_path)
		cmd.append('INPUT=%s' % input)
		cmd.append(self.tmp_dir)
		cmd.append('BAIT_INTERVALS=%s' % bait_interval)
		cmd.append('OUTPUT=%s' % output)
		cmd.append('PER_TARGET_COVERAGE=%s' % coverage)
		cmd.append('TARGET_INTERVALS=%s' % interval)
		cmd.append('REFERENCE_SEQUENCE=%s' % self.ref)
		exec_command = ' '.join(cmd)

		return exec_command


