import os
from utils import *

class annotation:
	def __init__(self, params, file_dic):
		self.file_dic = file_dic

	def vcf2maf(self, sample, type=''):
		cmd = []
		sample_id = sample
		cpcm = sample_id.split('_')[0]
		patient = sample_id.split('_')[2]
		tumor_id = '%s_%s' % (cpcm, patient)
		normal_id = '%s_%s' % (cpcm, patient.replace('T','N'))

		vcf1 = self.file_dic[sample_id]['filtered_mutect_vcf']
		maf1 = self.file_dic[sample_id]['mutect_maf']
		vcf2 = self.file_dic[sample_id]['filtered_somaticindelocator']
		maf2 = self.file_dic[sample_id]['somaticindelocator_maf']


		cmd.append('vcf2maf.pl')
		cmd.append('--input-vcf %s' % vcf1)
		cmd.append('--output-maf %s' % maf1)
	
		cmd.append('--tumor-id %s' % tumor_id)
		if type != '':
			cmd.append('--normal-id %s' % normal_id)
		
		exec_cmd1 = ' '.join(cmd)
	
		cmd = []

		cmd.append('vcf2maf.pl')
		cmd.append('--input-vcf %s' % vcf2)
		cmd.append('--output-maf %s' % maf2)
	
		cmd.append('--tumor-id %s' % tumor_id)
		if type != '':
			cmd.append('--normal-id %s' % normal_id)
		
		exec_cmd2 = ' '.join(cmd)

		exec_cmd3 = '%s\n%s' % (exec_cmd1, exec_cmd2)
		return exec_cmd3
	
	def maf_filter(self, sample):
		cmd = []

		maf_mutect = self.file_dic[sample]['mutect_maf']
		maf_SI = self.file_dic[sample]['somaticindelocator_maf']
		
		cmd.append('01_add_annotation.py -i %s' % maf_mutect)
		cmd.append('01_add_annotation.py -i %s' % maf_SI)
		
		exec_cmd = '\n'.join(cmd)
		return exec_cmd



class annotation_agg:
	def __init__(self, project, params):
		self.working_dir = '/home/D131748/Research/OncoPanel/aggregation/Projects/%s' % project
		self.working_tmp = '%s/tmp' % self.working_dir
		self.mutect_dir = '%s/mutect' % self.working_dir
		self.SI_dir = '%s/SI' % self.working_dir

	def vcf2maf(self, sample):
		cmd = []

		mutect_vcf = '%s/%s.agg.mutect.filtered.dbsnp_cosmic.vcf' % (self.mutect_dir, sample)
		SI_vcf = '%s/%s.agg.SI.filtered.fixDP.vcf' % (self.SI_dir, sample)
	
		cmd.append('vcf2maf.py')
		cmd.append('-i %s' % mutect_vcf)
		
		exec_cmd1 = ' '.join(cmd)
	
		cmd = []

		cmd.append('vcf2maf.py')
		cmd.append('-i %s' % SI_vcf)
		
		exec_cmd2 = ' '.join(cmd)


		exec_cmd3 = '%s\n%s' % (exec_cmd1, exec_cmd2)
		return exec_cmd3
	
	def maf_filter(self, sample):
		mutect_maf = '%s/%s.agg.mutect.filtered.dbsnp_cosmic.maf' % (self.mutect_dir, sample)
		mutect_maf2 = mutect_maf.replace('.maf', '.fixDP.maf')
		SI_maf = '%s/%s.agg.SI.filtered.fixDP.maf' % (self.SI_dir, sample)


		cmd = []
		cmd.append('fix_mutect_maf_DP.py -i %s' % mutect_maf)	
		cmd.append('01_add_annotation.py -i %s' % mutect_maf2)
		cmd.append('01_add_annotation.py -i %s' % SI_maf)
		
		exec_cmd = '\n'.join(cmd)
		return exec_cmd

