import subprocess
from utils import *

def pbs():
	cmd = []
	cmd.append('#!/bin/sh')
	cmd.append('#PBS -l nodes=1:ppn=8')
	cmd.append('#PBS -q batch')
	cmd.append('#PBS -l walltime=10:00:00')
	cmd.append('cd $PBS_O_WORKDIR')
	cmd.append('python /data/ASAN-CCGD/apps/CCGD/Breakmer/BreaKmer/breakmer.py breakmer.cfg')
	
	cmd2 = '\n'.join(cmd)

	return cmd2


class breakmer():
	def __init__(self, basecall, file_dic):
		self.basecall = basecall
		self.basecalls_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % basecall
		self.breakmer = '/data/ASAN-CCGD/apps/CCGD/Breakmer/BreaKmer/breakmer.py'
		self.breakmer_sample = '%s/configure/breakmer_samples.txt' % self.basecalls_dir
		self.config_file = '/data/ASAN-CCGD/apps/CCGD/Breakmer/BreaKmer/breakmer.cfg'

	def prep(self):
		cmd = []

		for x in open(self.breakmer_sample).xreadlines():
			x_1 = x.strip().split('\t')
			
			sample_analysis_dir = '%s/breakmer/%s' % (self.basecalls_dir, x_1[0])
			runlog = '%s/runlog' % sample_analysis_dir
			breakmer_data = '%s/data' % sample_analysis_dir
			sample_bam = '%s/sample.bam' % breakmer_data
			sample_bai = '%s/sample.bam.bai' % breakmer_data
			breakmer_config = '%s/breakmer.cfg' % sample_analysis_dir

			cmd.append('mkdir -p %s' % runlog)
			cmd.append('mkdir -p %s' % breakmer_data)
			cmd.append('ln -s %s %s' % (x_1[1], sample_bam))
			cmd.append('ln -s %s %s' % (x_1[1].replace('.bam','.bai'), sample_bai))
			
			cmd.append('sed "s@<sample_bam_path>@%s@g" %s > %s' % (sample_bam, self.config_file, breakmer_config))
			cmd.append('sed -i "s/<sample_id>/%s/g" %s' % (x_1[0], breakmer_config))
			cmd.append('sed -i "s@<analysis_dir>@%s@g" %s' % (sample_analysis_dir, breakmer_config))

		cmd2 = '\n'.join(cmd)
		return cmd2
		

	def prep2(self):
		cmd = []
		cmd.append('sh /data/ASAN-CCGD/apps/CCGD/Breakmer/BreaKmer/setup.sh')
		cmd.append('%s' % self.basecalls_dir)
		cmd.append('%s/configure/breakmer_samples.txt' % self.basecalls_dir)
		
		cmd2 = ' '.join(cmd)

		return cmd2


	def run_breakmer(self, sample):
		sample_analysis_dir = '%s/breakmer/%s' % (self.basecalls_dir, sample)
		cmd = 'python %s %s/breakmer.cfg' % (self.breakmer, sample_analysis_dir)
		
		return cmd



class breakmer_agg():
	def __init__(self, project, params):
		self.working_dir = '/home/D131748/Research/OncoPanel/aggregation/Projects/%s' % project
		self.working_tmp = '%s/tmp' % self.working_dir
		self.config_dir = '%s/configure' % self.working_dir
		self.breakmer_dir = '%s/breakmer' % self.working_dir

		self.breakmer = '/data/ASAN-CCGD/apps/CCGD/Breakmer/BreaKmer/breakmer.py'
		self.breakmer_sample = '%s/breakmer_samples.txt' % self.config_dir
		self.config_file = '/data/ASAN-CCGD/apps/CCGD/Breakmer/BreaKmer/breakmer.cfg'

	def prep(self):
		cmd = []

		for x in open(self.breakmer_sample).xreadlines():
			x_1 = x.strip().split('\t')
			
			sample_analysis_dir = '%s/%s' % (self.breakmer_dir, x_1[0])
			runlog = '%s/runlog' % sample_analysis_dir
			breakmer_data = '%s/data' % sample_analysis_dir
			sample_bam = '%s/sample.bam' % breakmer_data
			sample_bai = '%s/sample.bam.bai' % breakmer_data
			breakmer_config = '%s/breakmer.cfg' % sample_analysis_dir

			cmd.append('mkdir -p %s' % runlog)
			cmd.append('mkdir -p %s' % breakmer_data)
			cmd.append('ln -s %s %s' % (x_1[1], sample_bam))
			cmd.append('ln -s %s %s' % (x_1[1].replace('.bam','.bai'), sample_bai))
			
			cmd.append('sed "s@<sample_bam_path>@%s@g" %s > %s' % (sample_bam, self.config_file, breakmer_config))
			cmd.append('sed -i "s/<sample_id>/%s/g" %s' % (x_1[0], breakmer_config))
			cmd.append('sed -i "s@<analysis_dir>@%s@g" %s' % (sample_analysis_dir, breakmer_config))

		cmd2 = '\n'.join(cmd)
		return cmd2
		

	def prep2(self):
		cmd = []
		cmd.append('sh /data/ASAN-CCGD/apps/CCGD/Breakmer/BreaKmer/setup.sh')
		cmd.append('%s' % self.working_dir)
		cmd.append('%s/breakmer_samples.txt' % self.config_dir)
		
		cmd2 = ' '.join(cmd)

		return cmd2


	def run_breakmer(self, sample):
		sample_analysis_dir = '%s/%s' % (self.breakmer_dir, sample)
		cmd = 'python %s %s/breakmer.cfg' % (self.breakmer, sample_analysis_dir)
		
		return cmd
