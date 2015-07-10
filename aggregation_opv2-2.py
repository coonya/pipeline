#! /usr/bin/python
import optparse, os, sys, subprocess, MySQLdb, re
from subprocess import *
from utils import *
from sample import *
from picard import *
from gatk import *

## Global parameters
usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v2.0 mede by coonya')
parser.add_option('-p', dest='project',help='project name')
#parser.add_option('-c', dest='caller',help='somatic mutation caller e.g.) mutect, somaticindelocator')

(options, args) = parser.parse_args()

def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def realign(tumor, normal=''):
	cmd = []

	if normal != '':
		cmd.append(gatk.RealignerTargetCreator2(tumor, normal))
		cmd.append(gatk.IndelRealigner_agg2(tumor, normal))

	else:
		cmd.append(copy_unmateched_normal2(tumor, tmp_dir))
		cmd.append(gatk.RealignerTargetCreator2(tumor))
		cmd.append(gatk.IndelRealigner_agg2(tumor))

	exec_command = '\n'.join(cmd)

	return exec_command


## Check option
check_option(options.project, "You need to use -p option.")

### path
working_dir = '/home/D131748/Research/OncoPanel/aggregation/Projects/%s' % options.project
log_dir = '%s/log' % working_dir
metrics_dir = '%s/metrics' % working_dir
config_dir = '%s/configure' % working_dir
mutect_dir = '%s/mutect' % working_dir
SI_dir = '%s/SI' % working_dir
tmp_dir = '%s/tmp' % working_dir

logo()

check_working_dir(working_dir)
check_dir(log_dir)
check_dir(metrics_dir)
check_dir(config_dir)
check_dir(mutect_dir)
check_dir(SI_dir)
check_dir(tmp_dir)

config_file = '/home/D131748/src/pipeline/pipeline_config.cfg'
params = read_config(config_file)


picard = picard_agg(options.project, params)
gatk = gatk_agg(options.project, params)

## the number of cores
cpu_no = {
	'dedup': 2,
	'gatk': 1}

## get samples
(dic, paired_samples, non_paired_samples, sample_list, paired_sample_list, non_paired_sample_list) = get_sample_list_agg(options.project)

print sample_list
print paired_sample_list
print non_paired_sample_list
print paired_samples

## merge and deduplication
os.chdir(log_dir)
jobno = int(get_start_jobid())
wait_id1 = jobno

qsub_list = []
qsub_id = {}

### Merge and deduplication
for x in sample_list:
	wname = '%s/01.Merge_deduplication.%s.%s.sh' % (log_dir, x, wait_id1)
	wfile = open(wname,'w')
	wfile.write('%s\n' % pbs_header1(cpu_no['dedup'], picard.dedup_agg(x, dic[x])))
	wfile.close()

	qsub_key = 'dedup_%s' % x
	qsub_id[qsub_key] = wait_id1
	wait_id1 += 1

	qsub_name = 'qsub %s' %  wname
	#qsub_list.append('qsub %s/01.Merge_deduplication.%s.%s.sh' % (log_dir, x, wait_id1))

	subprocess.call(qsub_name, shell=True, stdout=subprocess.PIPE)

### Realignment around novel indel site
for x in paired_samples:
	tumor_sample = x
	normal_sample = x.replace('T', 'N')

	wait_id2 = []
	for y in paired_samples[x]:
		wait_qsub_key = 'dedup_%s' % y
		wait_id2.append(str(qsub_id[wait_qsub_key]))
	wait_id3 = ','.join(wait_id2)


	wname = '%s/02.agg_realignment.%s.%s.sh' % (log_dir, x, wait_id1)
	wfile = open(wname,'w')


	wfile.write('%s\n' % pbs_header2(cpu_no['gatk'], 'yes', wait_id3,  realign(tumor_sample, normal_sample)))
	wfile.close()

	qsub_key = 'gatk_%s' % x
	qsub_id[qsub_key] = wait_id1
	wait_id1 += 1

	qsub_name = 'qsub %s' % wname
	subprocess.call(qsub_name, shell=True, stdout=subprocess.PIPE)

for x in non_paired_samples:
	tumor_sample = x
	normal_sample = x.replace('T', 'N')

	wait_qsub_key = 'dedup_%s' % y
	wait_id3 = str(qsub_id[wait_qsub_key])


	wname = '%s/02.agg_realignment.%s.%s.sh' % (log_dir, x, wait_id1)
	wfile = open(wname,'w')


	wfile.write('%s\n' % pbs_header2(cpu_no['gatk'], 'yes', wait_id3,  realign(tumor_sample)))
	wfile.close()

	qsub_key = 'gatk_%s' % x
	qsub_id[qsub_key] = wait_id1
	wait_id1 += 1

	qsub_name = 'qsub %s' % wname
	subprocess.call(qsub_name, shell=True, stdout=subprocess.PIPE)


"""
	########### bwa mapping ########################
	wait_id2 = wait_id1 + len(sample_list)
	wfile = open('%s/03.bwa_aln.%s.%s.sh' % (log_dir, x, wait_id2),'w')
	wfile.write('%s\n' % pbs_header2(cpu_no['aln'], 'no', wait_id1, mapping(x)))
	wfile.close()

	qsub_key = 'bwa_aln_%s' % x
	qsub_id[qsub_key] = wait_id2
	qsub_list.append('qsub %s/03.bwa_aln.%s.%s.sh' % (log_dir, x, wait_id2))



	############ post mapping procedure #############
	cmd = []
	cmd.append(bwa.sampe(x))
	cmd.append(picard.MergeBamAlignment(x))
	cmd.append(picard.dedup(x))
	cmd.append(gatk.IndelRealigner(x))
	cmd.append(gatk.CountCovariates(x))
	cmd.append(gatk.TableRecalibration(x))
	cmd2 = '\n'.join(cmd)

	wait_id3 = wait_id2 + len(sample_list)
	wfile = open('%s/04.post_mapping.%s.%s.sh' % (log_dir, x, wait_id3),'w')
	wfile.write('%s\n' % pbs_header2(cpu_no['postmapping'], 'no', wait_id2, cmd2))
	wfile.close()

	qsub_key = 'postmapping_%s' % x
	qsub_id[qsub_key] = wait_id3
	qsub_list.append('qsub %s/04.post_mapping.%s.%s.sh' % (log_dir, x, wait_id3))

qsub_list.sort()
for x in qsub_list:
	subprocess.call(x, shell=True, stdout=subprocess.PIPE)

"""
