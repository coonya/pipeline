from utils import *

class somatic:
	def __init__(self, params, file_dic):
		self.file_dic = file_dic
		self.gatk_path = params['gatk_path']
		self.mutect_path = params['mutect_path']

		self.ref = params['ref']
		self.bait_interval = params['target_interval']
		self.dbsnp = params['dbsnp']
        	self.dbsnp_common = params['dbsnp_common']
		self.cosmic = params['cosmic']
		self.PON = params['PON']

	def mutect(self, tumor, normal=''):
		input_tumor = self.file_dic[tumor]['realign_agg']
		output_txt = self.file_dic[tumor]['mutect_txt']
		output_vcf = self.file_dic[tumor]['mutect_vcf']
		output_wig = self.file_dic[tumor]['mutect_wig']

		filter_vcf = self.file_dic[tumor]['filtered_mutect_vcf']
		filter_txt = self.file_dic[tumor]['filtered_mutect_txt']
		
		cmd = []

		cmd.append(java7_run('-Xmx4G -Djava.io.tmpdir=/data/scratch'))
		cmd.append('%s/mutect-1.1.7.jar --analysis_type MuTect' % self.mutect_path)
		cmd.append('--reference_sequence %s' % self.ref)
		cmd.append('--dbsnp %s' % self.dbsnp_common)
		cmd.append('--cosmic %s' % self.cosmic)
		cmd.append('--intervals %s' % self.bait_interval)

		if normal != '':
			input_normal = self.file_dic[normal]['realign_agg']

		else:
			tumor_name = os.path.basename(tumor).split('.')[0]
			normal_name = tumor_name.split('.')[0].replace('T', 'N')
			input_normal = input_tumor.replace(tumor_name, normal_name)
		
		cmd.append('--input_file:normal %s' % input_normal)
		cmd.append('--input_file:tumor %s' % input_tumor)
		
		cmd.append('--normal_panel %s' % self.PON)
		cmd.append('--out %s' % output_txt)
		cmd.append('--vcf %s' % output_vcf)
		cmd.append('--coverage_file %s' % output_wig)
		cmd.append('--enable_extended_output')
		exec_command = ' '.join(cmd)

		cmd15_1 = 'grep -v REJECT %s > %s\ngrep -v REJECT %s > %s\n' % (output_vcf, filter_vcf, output_txt, filter_txt)

		exec_command = '%s\n%s' % (exec_command, cmd15_1)


		return exec_command



	def SI(self, tumor, normal=''):
		input_tumor = self.file_dic[tumor]['realign_agg']
		output = self.file_dic[tumor]['somaticindelocator']

		cmd = []

		cmd.append('java -Xmx4G -jar %s/GenomeAnalysisTK.jar' % self.gatk_path)
		cmd.append('-T SomaticIndelDetector')
		cmd.append('-R %s' % self.ref)
		cmd.append('-o %s' % output)
		cmd.append('-L %s' % self.bait_interval)
		cmd.append('--window_size 500')
	
		cmd2 = '/home/D131748/src/01_somaticindel_filter.py -i %s' % output

		if normal != '':
			input_normal = self.file_dic[normal]['realign_agg']

			cmd.append('-filter T_COV\<2\|\|N_COV\<0\|\|T_INDEL_F\<0.05\|\|T_INDEL_CF\<0.1')
			
			cmd1 = ' '.join(cmd)

		else:
			tumor_name = os.path.basename(tumor).split('.')[0]
			normal_name = tumor_name.split('.')[0].replace('T', 'N')
			input_normal = input_tumor.replace(tumor_name, normal_name)
		
			cmd.append('-filter T_COV\<2\|\|T_INDEL_F\<0.05\|\|INDEL_CF\<0.1')

		cmd.append('-I:normal %s' % input_normal)
		cmd.append('-I:tumor %s' % input_tumor)
		cmd1 = ' '.join(cmd)
		
		exec_command = '%s\n%s' % (cmd1, cmd2)

		return exec_command


class somatic_agg:
	def __init__(self, project, params):
		self.gatk_path = params['gatk_path']
		self.mutect_path = params['mutect_path']

		self.ref = params['ref']
		self.bait_interval = params['target_interval']
		#self.dbsnp = params['dbsnp']
        	self.dbsnp_common = params['dbsnp_common']
		self.cosmic = params['cosmic']
		self.PON = params['PON']
		
		self.working_dir = '/home/D131748/Research/OncoPanel/aggregation/Projects/%s' % project
		self.working_tmp = '%s/tmp' % self.working_dir
		self.mutect_dir = '%s/mutect' % self.working_dir
		self.SI_dir = '%s/SI' % self.working_dir

	def mutect(self, tumor, normal):
		input_tumor = '%s/%s.agg.dedup.cleaned.bam' % (self.working_tmp, tumor)
		input_normal = '%s/%s.agg.dedup.cleaned.bam' % (self.working_tmp, normal)

		output_txt = '%s/%s.agg.mutect.txt' % (self.mutect_dir, tumor)
		output_vcf = '%s/%s.agg.mutect.vcf' % (self.mutect_dir, tumor)
		output_wig = '%s/%s.agg.mutect.wig' % (self.mutect_dir, tumor)

		filter_vcf = '%s/%s.agg.mutect.filtered.vcf' % (self.mutect_dir, tumor)
		filter_txt = '%s/%s.agg.mutect.filtered.txt' % (self.mutect_dir, tumor)
		
		filter_vcf2 = '%s/%s.agg.mutect.filtered.dbsnp_cosmic.vcf' % (self.mutect_dir, tumor)

		cmd = []

		cmd.append(java7_run('-Xmx4G -Djava.io.tmpdir=/data/scratch'))
		cmd.append('%s/mutect-1.1.7.jar --analysis_type MuTect' % self.mutect_path)
		cmd.append('--reference_sequence %s' % self.ref)
		cmd.append('--dbsnp %s' % self.dbsnp_common)
		cmd.append('--cosmic %s' % self.cosmic)
		cmd.append('--intervals %s' % self.bait_interval)
		cmd.append('--input_file:normal %s' % input_normal)
		cmd.append('--input_file:tumor %s' % input_tumor)
		
		cmd.append('--normal_panel %s' % self.PON)
		cmd.append('--out %s' % output_txt)
		cmd.append('--vcf %s' % output_vcf)
		cmd.append('--coverage_file %s' % output_wig)
		cmd.append('--enable_extended_output')
		exec_command = ' '.join(cmd)

		cmd2 = ['grep -v REJECT %s > %s\ngrep -v REJECT %s > %s' % (output_vcf, filter_vcf, output_txt, filter_txt)]
		cmd2.append('filter_dbsnp_cosmic.py -i %s' % filter_vcf)

		cmd3 = '\n'.join(cmd2)

		exec_command = '%s\n%s' % (exec_command, cmd3)

		return exec_command



	def SI(self, tumor, normal):
		input_tumor = '%s/%s.agg.dedup.cleaned.bam' % (self.working_tmp, tumor)
		input_normal = '%s/%s.agg.dedup.cleaned.bam' % (self.working_tmp, normal)

		output = '%s/%s.agg.SI.vcf' % (self.SI_dir, tumor)

		filter_vcf = '%s/%s.agg.SI.filtered.vcf' % (self.SI_dir, tumor)

		cmd = []

		cmd.append('java -Xmx4G -jar %s/GenomeAnalysisTK.jar' % self.gatk_path)
		cmd.append('-T SomaticIndelDetector')
		cmd.append('-R %s' % self.ref)
		cmd.append('-I:normal %s' % input_normal)
		cmd.append('-I:tumor %s' % input_tumor)
		cmd.append('-o %s' % output)
		cmd.append('-L %s' % self.bait_interval)
		cmd.append('--window_size 500')
		cmd.append('-filter T_COV\<2\|\|N_COV\<0\|\|T_INDEL_F\<0.05\|\|T_INDEL_CF\<0.1')

		cmd1 = ' '.join(cmd)
		
		cmd2 = ['/home/D131748/src/01_somaticindel_filter.py -i %s' % output]
		cmd2.append('fix_SI_DP.py -i %s' % filter_vcf)

		cmd3 = '\n'.join(cmd2)

		exec_command = '%s\n%s' % (cmd1, cmd3)

		return exec_command


