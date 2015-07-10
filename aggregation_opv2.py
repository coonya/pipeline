#! /usr/bin/python
import optparse, os, sys, subprocess, MySQLdb, re, time
from subprocess import *
from utils import *
from sample import *
from picard import *
from gatk import *
from somatic import *
from annotation import *
from breakmer import *

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


def callingmutation(tumor, normal):
	cmd = []
	cmd.append(somatic.mutect(tumor, normal))
	cmd.append(somatic.SI(tumor, normal))
	cmd.append(anno.vcf2maf(tumor))
	cmd.append(anno.maf_filter(tumor))

	exec_command = '\n'.join(cmd)

	return exec_command


def calculate_metrics(sample):
	cmd = []

	cmd.append(picard.CalculateHsMetrics(sample, 'exon'))
	cmd.append(picard.CalculateHsMetrics(sample, 'full'))
	cmd.append(picard.CalculateHsMetrics(sample, 'intron'))

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
cbio_dir = '%s/cbioportal/case_lists' % working_dir

logo()

check_working_dir(working_dir)
check_dir(log_dir)
check_dir(metrics_dir)
check_dir(config_dir)
check_dir(mutect_dir)
check_dir(SI_dir)
check_dir(tmp_dir)
check_dir(cbio_dir)

config_file = '/home/D131748/src/pipeline/pipeline_config.cfg'
params = read_config(config_file)


picard = picard_agg(options.project, params)
gatk = gatk_agg(options.project, params)
somatic = somatic_agg(options.project, params)
anno = annotation_agg(options.project, params)
breakmer = breakmer_agg(options.project, params)

## the number of cores
cpu_no = {
	'dedup': 2,
	'gatk': 1,
	'calling_mutation': 1,
	'SV': 8,
	'calculate_metrics': 1,
	'report':1}

time.sleep(1)

## get samples
(dic, paired_samples, non_paired_samples, sample_list1, paired_sample_list, non_paired_sample_list) = get_sample_list_agg(options.project)

## merge and deduplication
os.chdir(log_dir)
jobno = int(get_start_jobid())
wait_id1 = jobno

qsub_list = []
qsub_id = {}

### Merge and deduplication
for x in sample_list1:
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

print "Submitting job for merge and deduplication steps."

sample_list = []
qsub_wait = []

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
	sample_list.append(x)
	qsub_id[qsub_key] = wait_id1
	qsub_wait.append(str(wait_id1))
	wait_id1 += 1

	qsub_name = 'qsub %s' % wname
	subprocess.call(qsub_name, shell=True, stdout=subprocess.PIPE)

for x in non_paired_samples:
	tumor_sample = x
	normal_sample = x.replace('T', 'N')

	wait_qsub_key = 'dedup_%s' % x
	wait_id3 = str(qsub_id[wait_qsub_key])


	wname = '%s/02.agg_realignment.%s.%s.sh' % (log_dir, x, wait_id1)
	wfile = open(wname,'w')


	wfile.write('%s\n' % pbs_header2(cpu_no['gatk'], 'yes', wait_id3,  realign(tumor_sample)))
	wfile.close()

	qsub_key = 'gatk_%s' % x
	sample_list.append(x)
	qsub_id[qsub_key] = wait_id1
	qsub_wait.append(str(wait_id1))
	wait_id1 += 1

	qsub_name = 'qsub %s' % wname
	subprocess.call(qsub_name, shell=True, stdout=subprocess.PIPE)

print "Submitting jobs for realignment steps."
qsub_final = []

### calling mutation
for x in sample_list:
	print x
	tumor_sample = x
	normal_sample = x.replace('T', 'N')

	wait_qsub_key = 'gatk_%s' % x
	wait_id3 = qsub_id[wait_qsub_key]

	wname = '%s/03.calling_mutations.%s.%s.sh' % (log_dir, x, wait_id1)
	wfile = open(wname, 'w')

	wfile.write('%s\n' % pbs_header2(cpu_no['calling_mutation'], 'no', wait_id3, callingmutation(tumor_sample, normal_sample)))
	wfile.close()

	qsub_key = 'callingmutation_%s' % x
	qsub_id[qsub_key] = wait_id1
	qsub_final.append(str(wait_id1))
	wait_id1 += 1
	
	qsub_name = 'qsub %s' % wname
	subprocess.call(qsub_name, shell=True, stdout=subprocess.PIPE)

print "Submitting jobs for calling mutation steps."


### calling sv
qsub_id_all = ','.join(qsub_wait)
	
wname = '%s/04.01.callingSV_prep.%s.sh' % (log_dir, wait_id1)
wfile = open(wname, 'w')
wfile.write('%s\n' % pbs_header2(cpu_no['SV'], 'yes', qsub_id_all, breakmer.prep()))
wfile.close()
subprocess.call('qsub %s' % wname, shell=True, stdout=subprocess.PIPE)

qsub_prep_id = wait_id1
wait_id1 += 1

print "Submitting job for preparing breakmer information."

for z in sample_list1:
	print z
        wname = '%s/04.02.callingSV.%s.%s.sh' % (log_dir, z, wait_id1)
	wfile = open(wname, 'w')
        wfile.write('%s\n' % pbs_header2(cpu_no['SV'], 'no', qsub_prep_id, breakmer.run_breakmer(z)))
        wfile.close()

        subprocess.call('qsub %s' % wname, shell=True, stdout=subprocess.PIPE)

        qsub_key = 'SV_%s' % z
        qsub_id[qsub_key] = wait_id1

        qsub_final.append(str(wait_id1))
	wait_id1 += 1

print "Submitting jobs for calling SV steps."


### calculate metrics


for z in sample_list1:
        wname = '%s/05.calculate_metrics.%s.%s.sh' % (log_dir, z, wait_id1)
	wfile = open(wname, 'w')

        wfile.write('%s\n' % pbs_header2(cpu_no['calculate_metrics'], 'yes', qsub_id_all, calculate_metrics(z)))
        wfile.close()

        subprocess.call('qsub %s' % wname, shell=True, stdout=subprocess.PIPE)

        qsub_key = 'calculate_metrics_%s' % z
        qsub_id[qsub_key] = wait_id1

        qsub_final.append(str(wait_id1))
	wait_id1 += 1

print "Submitting jobs for calculating metrics."

### generating excel report
qsub_id_last = ','.join(qsub_final)
wname = '%s/06.report.%s.sh' % (log_dir, options.project)
wfile = open(wname, 'w')
wfile.write('%s\n' % pbs_header2(cpu_no['report'], 'yes', qsub_id_last, excel_report(working_dir)))
wfile.close()
subprocess.call('qsub %s' % wname, shell=True, stdout=subprocess.PIPE)

print "Complete to submit jobs to AMC cluster."
