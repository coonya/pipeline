from utils import *

class gatk:

        def __init__(self, params, basecall, file_dic):
                self.basecall = basecall
                self.working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % basecall
                self.config_dir = '%s/configure' % self.working_dir

		self.file_dic = file_dic

		self.gatk_path = params['gatk_path']
	
		self.interval = params['known_interval']
		self.known1 = params['known1']
		self.known2 = params['known2']
		self.dbsnp = params['dbsnp']
		self.ref = params['ref']
		self.bait_interval = params['target_interval']
		self.depth_interval = params['depth_interval']

	def IndelRealigner(self, sample):
		input = self.file_dic[sample]['dedup']
		output = self.file_dic[sample]['realign']

		cmd = []

		cmd.append(java_run('-Xmx4G -Djava.io.tmpdir=/data/scratch -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10'))
		cmd.append('%s/GenomeAnalysisTK.jar -T IndelRealigner' % self.gatk_path)
		cmd.append('-I %s' % input)
		cmd.append('-R %s' % self.ref)
		cmd.append('-targetIntervals %s' % self.interval)
		cmd.append('-o %s' % output)
		cmd.append('-known %s -known %s' % (self.known1, self.known2))
		cmd.append('--consensusDeterminationModel KNOWNS_ONLY')
		cmd.append('-LOD 0.4')
		cmd.append('-compress 1')
		cmd.append('-maxInMemory 1000000')
		exec_command = ' '.join(cmd)
		
		return exec_command


	def CountCovariates(self, sample):
		input = self.file_dic[sample]['realign']
		output = self.file_dic[sample]['recalcsv']

		cmd = []

		cmd.append(java_run('-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T CountCovariates' % self.gatk_path)
		cmd.append('-l INFO')
		cmd.append('-R %s' % self.ref)
		cmd.append('-knownSites %s' % self.dbsnp)
		cmd.append('-I %s' % input)
		cmd.append('-recalFile %s' % output)
		cmd.append('-nt 1')
		cmd.append('-cov ReadGroupCovariate -cov CycleCovariate -cov DinucCovariate -cov QualityScoreCovariate -OQ')
		exec_command = ' '.join(cmd)
		
		return exec_command

	
	def TableRecalibration(self, sample):
		input = self.file_dic[sample]['realign']
		output = self.file_dic[sample]['recal']
		recal_csv = self.file_dic[sample]['recalcsv']

		cmd = []

		cmd.append(java_run('-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T TableRecalibration' % self.gatk_path)
		cmd.append('-l INFO')
		cmd.append('-R %s' % self.ref)
		cmd.append('-I %s' % input)
		cmd.append('--out %s' % output)
		cmd.append('-recalFile %s' % recal_csv)
		cmd.append('-OQ --generate_md5')
		exec_command = ' '.join(cmd)

		return exec_command


	def depthofcoverage(self, sample):
		
		cmd = []

		cmd.append(java_run('-Xmx4G -Djava.io.tmpdir=/data/scratch'))
		cmd.append('%s/GenomeAnalysisTK.jar -T DepthOfCoverage' % gatk_path) 
		cmd.append('--omitDepthOutputAtEachBase')
		cmd.append('--interval_merging OVERLAPPING_ONLY')
		cmd.append('-L %s' % self.depth_interval)
		cmd.append('-I %s' % input)
		cmd.append('-o %s' % output)
		cmd.append('-R %s' % self.ref)
		exec_command = ' '.join(cmd)

		return exec_command

	
	### aggregation steps for lane level

	def RealignerTargetCreator(self, tumor, normal=''): ## output, tumor, normal
		indel_interval = self.file_dic[tumor]['indel_interval']
		tumor =	self.file_dic[tumor]['recal']

		if normal != '':
			normal = self.file_dic[normal]['recal']
		else:
			tumor_sample = os.path.basename(tumor).split('.')[0]
			normal_sample = os.path.basename(tumor).split('.')[0].replace('T', 'N')
			normal = tumor.replace(tumor_sample, normal_sample)

		cmd = []
		
		cmd.append(java_run('-XX:GCHeapFreeLimit=10 -XX:GCTimeLimit=50 -Djava.io.tmpdir=/data/scratch -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T RealignerTargetCreator' % self.gatk_path)

		cmd.append('-I %s -I %s' % (normal, tumor))

		cmd.append('-o %s' % indel_interval)
		cmd.append('-R %s' % self.ref)
		cmd.append('-nt 1')
		cmd.append('-dcov 250')
		exec_command = ' '.join(cmd)

		return exec_command
	

	def IndelRealigner_agg(self, tumor1, normal1=''): 
		map = self.file_dic[tumor1]['map']
		indel_interval = self.file_dic[tumor1]['indel_interval']
		tumor = self.file_dic[tumor1]['recal']
		tumor_agg = self.file_dic[tumor1]['realign_agg']

		cmd = []

		wfile = open(map, 'w')
		wfile2 = open('%s/breakmer_samples.txt' % self.config_dir, 'a')

		## generate map file
		cmd.append(java_run('-Djava.io.tmpdir=/data/scratch -Dsamjdk.use_async_io=true -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T IndelRealigner' % self.gatk_path)

		if normal1 != '':
			normal_sample = os.path.basename(normal1).split('.')[0]
			normal = self.file_dic[normal1]['recal']
			normal_agg = self.file_dic[normal1]['realign_agg']
			tumor_sample = os.path.basename(tumor).split('.')[0]

			wfile2.write('%s\t%s\n%s\t%s\n' % (normal_sample, normal_agg, tumor_sample, tumor_agg))

			cmd.append('-I %s -I %s' % (normal, tumor))

			wfile.write('%s\t%s\n' % (os.path.basename(normal), normal_agg))
			wfile.write('%s\t%s\n' % (os.path.basename(tumor), tumor_agg))

			cmd.append('--nWayOut %s' % map)

		else:
			wfile = open(map, 'w')
			tumor_sample = os.path.basename(tumor).split('.')[0]
			normal_sample = os.path.basename(tumor).split('.')[0].replace('T', 'N')
			normal = tumor.replace(tumor_sample, normal_sample)
			normal_agg = '%s/aggregation/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.cleaned.bam' % (os.path.dirname(tumor), normal_sample)

			wfile2.write('%s\t%s\n%s\t%s\n' % (normal_sample, normal_agg, tumor_sample, tumor_agg))

			cmd.append('-I %s -I %s' % (normal, tumor))

			wfile.write('%s\t%s\n' % (os.path.basename(normal), normal_agg))
			wfile.write('%s\t%s\n' % (os.path.basename(tumor), tumor_agg))

			cmd.append('--nWayOut %s' % map)


		cmd.append('-targetIntervals %s' % indel_interval)
		cmd.append('-R %s' % self.ref)
		cmd.append('-maxInMemory 1000000')
		cmd.append('-model USE_READS')
		exec_command = ' '.join(cmd)
		wfile.close()
		wfile2.close()

		return exec_command

### aggregation steps for final

class gatk_agg:

        def __init__(self, project, params):
                self.working_dir = '/home/D131748/Research/OncoPanel/aggregation/Projects/%s' % project
                self.working_tmp = '%s/tmp' % self.working_dir
                self.config_dir = '%s/configure' % self.working_dir

		self.gatk_path = params['gatk_path']
	
		#self.interval = params['known_interval']
		#self.known1 = params['known1']
		#self.known2 = params['known2']
		#self.dbsnp = params['dbsnp']
		self.ref = params['ref']
		#self.bait_interval = params['target_interval']
		#self.depth_interval = params['depth_interval']


	def RealignerTargetCreator2(self, tumor, normal=''): ## output, tumor, normal
		indel_interval = '%s/%s_indelCocleanTargets.intervals' % (self.working_tmp, tumor)
		if normal == '':
			normal_id = tumor.replace('T', 'N')
			normal = '%s/%s.agg.dedup.bam' % (self.working_tmp, normal_id)
		else:
			normal = '%s/%s.agg.dedup.bam' % (self.working_tmp, normal)

		tumor =	'%s/%s.agg.dedup.bam' % (self.working_tmp, tumor)
		cmd = []

		cmd.append(java_run('-XX:GCHeapFreeLimit=10 -XX:GCTimeLimit=50 -Djava.io.tmpdir=/data/scratch -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T RealignerTargetCreator' % self.gatk_path)

		cmd.append('-I %s -I %s' % (normal, tumor))

		cmd.append('-o %s' % indel_interval)
		cmd.append('-R %s' % self.ref)
		cmd.append('-nt 1')
		cmd.append('-dcov 250')
		exec_command = ' '.join(cmd)

		return exec_command
	

	def IndelRealigner_agg2(self, tumor1, normal1=''): 
		map = '%s/%s.map' % (self.working_tmp, tumor1)
		indel_interval = '%s/%s_indelCocleanTargets.intervals' % (self.working_tmp, tumor1)
		tumor =	'%s/%s.agg.dedup.bam' % (self.working_tmp, tumor1)
		tumor_agg = tumor.replace('.bam', '.cleaned.bam')

		wfile = open(map, 'w')
		wfile2 = open('%s/breakmer_samples.txt' % self.config_dir, 'a')

		wfile.write('%s\t%s\n' % (os.path.basename(tumor), tumor_agg))

		if normal1 != '':
			normal = '%s/%s.agg.dedup.bam' % (self.working_tmp, normal1)
			normal_agg = normal.replace('.bam', '.cleaned.bam')

			wfile.write('%s\t%s\n' % (os.path.basename(normal), normal_agg))
			wfile2.write('%s\t%s\n%s\t%s\n' % (normal1, normal_agg, tumor1, tumor_agg))
		else:
			normal_id = tumor1.replace('T', 'N')
			normal = '%s/%s.agg.dedup.bam' % (self.working_tmp, normal_id)
			normal_agg = normal.replace('.bam', '.cleaned.bam')

			wfile.write('%s\t%s\n' % (os.path.basename(normal), normal_agg))
			wfile2.write('%s\t%s\n' % (tumor1, tumor_agg))

		cmd = []

		cmd.append(java_run('-Djava.io.tmpdir=/data/scratch -Dsamjdk.use_async_io=true -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4G -Xms2G'))
		cmd.append('%s/GenomeAnalysisTK.jar -T IndelRealigner' % self.gatk_path)
		cmd.append('-I %s -I %s' % (normal, tumor))
		cmd.append('--nWayOut %s' % map)
		cmd.append('-targetIntervals %s' % indel_interval)
		cmd.append('-R %s' % self.ref)
		cmd.append('-maxInMemory 1000000')
		cmd.append('-model USE_READS')
		exec_command = ' '.join(cmd)

		return exec_command


