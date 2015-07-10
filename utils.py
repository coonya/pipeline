import sys, os, glob, optparse, shutil, subprocess, logging
from bwa import *

def logo():

	logo = """
###############################################################################
### @@@@@@@ ###################################################################
## @@@    @@ ###     #####     #####     ##############     ###################
## @@ ###    ## @@@@@ ### @@@@@ ### @@@@@ ## @@ # @@ # @@@@@@ #################
## @@ ######## @@ # @@ # @@ # @@ # @@   @@ ## @@ @@ # @@    @@ ################
## @@ ###    # @@ # @@ # @@ # @@ # @@ # @@ ### @@@ # @@ #### @@ ###############
## @@@    @@ # @@ # @@ # @@ # @@ # @@ # @@ ### @@ ### @@    @@@@ ##############
### @@@@@@@ ### @@@@@ ### @@@@@ ## @@ # @@ ## @@ ##### @@@@@ # @@ #############
###############################################################################
######################################### OncoPanelV2 Analysis Pipeline v2.0 ##
###############################################################################
"""
        print logo

def check_working_dir(working_dir):
	if os.path.exists(working_dir):
		print "'%s' directory already exist." % working_dir
		print "Would you like to delete existing '%s' directory?" % working_dir
		answer = raw_input('[y/n/q] :')
		if answer.upper() == 'Y':
			shutil.rmtree(working_dir)
		elif answer.upper() == 'N':
			pass
		elif answer.upper() == 'Q':
			sys.exit('You choose quit!')
		else:
			sys.exit('You typeed wrong botton or bottons. Please run again!!')
	else:
		pass

def check_log(log_dir):
	os.chdir(log_dir)
	log_files = glob.glob('*')
	if len(log_files) != 0:
		for log_file in log_files:
			os.remove(log_file)


def read_config(config_file):
	params = {}
	for x in open(config_file).xreadlines():
		if x[0] != '#' and x != '\n':
			var, val = x.strip().split('=')
			params[var] = val
	return params

def java_run(options):
	java_command = '/usr/bin/java %s -jar' % options
	return java_command


def java7_run(options):
        java_command = '/home/D131748/apps/JAVA/jre1.7.0_60/bin/java %s -jar' % options
        return java_command


def check_dir(check_dir):
	if os.path.exists(check_dir):
		pass
	else:
		try:
			os.makedirs(check_dir)
		except:
			pass

def check_file(file):
	try:
		files = glob.glob(file)
		if file not in files:
			sys.exit('%s is not in configuration directory. Please make sure that %s is located in proper directory\n' % (file, file))
	except:
		sys.exit('%s is not in configuration directory. Please make sure that %s is located in proper directory\n' % (file, file))

def print_info(msg):
	print 'INFO	PIPELINE_OPv2 - %s' % (msg)


def get_start_jobid():
	find_jobid = subprocess.Popen("qmgr -c 'p s'|awk '/next_job_number/{print $5}'", shell=True, stdout=subprocess.PIPE)
	jobid = find_jobid.stdout.readline().split('\n')[0]
	
	return jobid

def pbs_header1(ppn, content):
	cmd = []
	cmd.append('#PBS -l nodes=1:ppn=%d' % ppn)
	cmd.append('#PBS -l walltime=48:00:00')
	cmd.append('#PBS -j oe')
	cmd.append('#PBS -q batch')
	cmd.append('cd $PBS_O_WORKDIR\n')
	cmd.append('%s' % content)

	cmd2 = '\n'.join(cmd)

        return cmd2

def pbs_header2(ppn, multi, job_id, content):
	cmd = []

	if multi == 'no':
		cmd.append('#PBS -l nodes=1:ppn=%d' % ppn)
		cmd.append('#PBS -l walltime=48:00:00')
		cmd.append('#PBS -j oe')
		cmd.append('#PBS -q batch')
		cmd.append('#PBS -W depend=afterok:%s.master' % job_id)
		cmd.append('cd $PBS_O_WORKDIR\n')
		cmd.append('%s' % content)

	elif multi == 'yes':
		jobs = job_id.split(',')
		jobs_order = []

		for x in jobs:
			jobs_order.append('%s.master' % x)
		jobs = ':'.join(jobs_order)
		cmd.append('#PBS -l nodes=1:ppn=%d' % ppn)
		cmd.append('#PBS -l walltime=48:00:00')
		cmd.append('#PBS -j oe')
		cmd.append('#PBS -q batch')
		cmd.append('#PBS -W depend=afterok:%s' % jobs)
		cmd.append('cd $PBS_O_WORKDIR\n')
		cmd.append('%s' % content)
	
	cmd2 = '\n'.join(cmd)

	return cmd2




def setup_logger(log_dir, name) :
  output_path = log_dir

  logger = logging.getLogger(name)
  logger.setLevel(logging.DEBUG)

  # FileHandler
  fh = logging.FileHandler(os.path.join(output_path,'log.txt'),mode='w')
  fh.setLevel(logging.DEBUG)
  # ConsoleHandler
  ch = logging.StreamHandler()
  ch.setLevel(logging.ERROR)

  formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p')
  fh.setFormatter(formatter)
  ch.setFormatter(formatter)

  logger.addHandler(fh)
  logger.addHandler(ch)

def add_logger(file, tool):
	cmd = []
	cmd.append('python /home/D131748/src/pipeline/cmd_run.py -t %s -s 0' % tool)
	cmd.append('sh %s' % file)
	cmd.append('pyton /home/D131748/src/pipeline/cmd_run.py -t %s -s 1' % tool)

	cmd2 = '\n'.join(cmd)
	
	return cmd2

def cmd_run(file):
	cmd = 'sh %s' % file
	
	return cmd


def bwa_aln(sample):
	bwa.aln1(sample)
	bwa.aln2(sample)

def excel_report(working_dir):

	cmd = '02_generating_excel.py -d %s' % working_dir
	
	return cmd


def copy_unmateched_normal(tumor_sample, pwd):
	tumor_path = '%s/tmp' % pwd 
	normal_sample = tumor_sample.split('.')[0].replace('T', 'N')
	
	cmd = []
	cmd.append('ln -s /home/D131748/src/pipeline/13A055_0005_46242189N.aggregation.dedup.bam %s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bam' % (tumor_path, normal_sample))
	cmd.append('ln -s /home/D131748/src/pipeline/13A055_0005_46242189N.aggregation.dedup.bai %s/%s.fastq.gz.initialAlign.merged.dedup.realign.recal.bai' % (tumor_path, normal_sample))
	cmd2 = '\n'.join(cmd)

	return cmd2

def copy_unmateched_normal2(tumor_sample, pwd):
	tumor_path = pwd 
	normal_sample = tumor_sample.split('.')[0].replace('T', 'N')
	
	cmd = []
	cmd.append('ln -s /home/D131748/src/pipeline/13A055_0005_46242189N.aggregation.dedup.bam %s/%s.agg.dedup.bam' % (tumor_path, normal_sample))
	cmd.append('ln -s /home/D131748/src/pipeline/13A055_0005_46242189N.aggregation.dedup.bai %s/%s.agg.dedup.bai' % (tumor_path, normal_sample))
	cmd2 = '\n'.join(cmd)

	return cmd2

