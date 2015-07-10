#! /usr/bin/python
import optparse
import sys
import logging
import time
from utils import *
from info_ext import *
from bwa import *
from sample import *
from picard import *
from gatk import *
from somatic import *
from breakmer import *
from annotation import *

usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v3.0 made by coonya')
parser.add_option('-d', dest='basecall',help='basecall DIR')
(options, args) = parser.parse_args()



def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())



def mapping(sample):
	cmd = []

	cmd.append(bwa.aln1(sample))
	cmd.append(bwa.aln2(sample))
	
	exec_command = '\n'.join(cmd)
	
	return exec_command



def aggregation(tumor, normal=''):
	cmd = []

	if normal != '':
		cmd.append(gatk.RealignerTargetCreator(tumor, normal))
		cmd.append(gatk.IndelRealigner_agg(tumor, normal))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'exon'))
		#cmd.append(picard.CalculateHsMetrics(normal, 'exon'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'full'))
		#cmd.append(picard.CalculateHsMetrics(normal, 'full'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'intron'))
		#cmd.append(picard.CalculateHsMetrics(normal, 'intron'))
	else:
		cmd.append(copy_unmateched_normal(tumor, working_dir))
		cmd.append(gatk.RealignerTargetCreator(tumor))
		cmd.append(gatk.IndelRealigner_agg(tumor))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'exon'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'full'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'intron'))
	
	exec_command = '\n'.join(cmd)

	return exec_command



def calculate_metrics(sample):
	cmd = []

	cmd.append(picard.CalculateHsMetrics(sample, 'exon'))
	cmd.append(picard.CalculateHsMetrics(sample, 'full'))
	cmd.append(picard.CalculateHsMetrics(sample, 'intron'))

	exec_command = '\n'.join(cmd)

	return exec_command

	

def callingmutation(tumor, normal=''):
	cmd = []

	if normal != '':
		cmd.append(somatic.mutect(tumor, normal))
		cmd.append(somatic.SI(tumor, normal))
		cmd.append(anno.vcf2maf(tumor))
		cmd.append(anno.maf_filter(tumor))

	else:
		cmd.append(somatic.mutect(tumor))
		cmd.append(somatic.SI(tumor))
		cmd.append(anno.vcf2maf(tumor))
		cmd.append(anno.maf_filter(tumor))

	exec_command = '\n'.join(cmd)

	return exec_command


	
cpu_no = {
	'basecall': 16,
	'SamToFastq': 2,
	'aln': 16,
	'postmapping': 2,
	'agg': 1,
	'callingmutation': 1,
	'SV': 8,
	'calculate_metrics': 1,
	'report': 1}



########## path or file ################
## config file
config_file = '/home/D131748/src/pipeline/pipeline_config.cfg'

## path for pipeline
working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % options.basecall
log_dir = '%s/log' % working_dir
config_dir = '%s/configure' % working_dir
config_dir = '%s/configure' % working_dir
mutect_dir = '%s/mutect' % working_dir
SI_dir = '%s/SI' % working_dir
metrics_dir = '%s/metrics' % working_dir
sample_file = '%s/sample.txt' % config_dir
agg_dir = '%s/tmp/aggregation' % working_dir

logo()



########## checking options or directories or files #####################
check_option(options.basecall, "You need to use '-d' option with basecall DIR\n")
check_working_dir(working_dir)
check_dir(log_dir)
check_dir(mutect_dir)
check_dir(SI_dir)
check_dir(metrics_dir)
check_dir(agg_dir)



################# Generation of logger ################
os.chdir(log_dir)
setup_logger(log_dir, 'root')
logger = logging.getLogger('root')



############# extract parameters ###################
params = read_config(config_file)
logger.info('A extraction of parameters is done.')



############### extract information ##################
info_ext = extract(options.basecall).ext_info()
logger.info('A extraction of information is done.')



################# R1 : get sample information ########
(sample_dic, sample_list, normal_list, tumor_list, paired_list, unpaired_list) = get_sample_list(sample_file)
print sample_dic
file_dic = get_file_name(options.basecall, sample_list)
print_sample(sample_list, normal_list, tumor_list, paired_list, unpaired_list)



############## Generation of picard object
picard = picard(params, options.basecall, file_dic)
bwa = bwa(params, file_dic)
gatk = gatk(params, options.basecall, file_dic)
somatic = somatic(params, file_dic)
breakmer = breakmer(options.basecall, file_dic)
anno = annotation(params, file_dic)



################ R2 : Extract Illumina Barcodes and Illumina Basecalls To Sam ####################
jobno = int(get_start_jobid())
qsub_id = {'basecall': jobno}
wait_id = qsub_id['basecall']
wait_id1 = jobno

cmd = []
cmd.append(picard.ExtractIlluminaBarcodes())
cmd.append(picard.IlluminaBasecallsToSam())

wfile = open('%s/01.BaseCall.%s.sh' % (log_dir, wait_id),'w')
cmd2 = '\n'.join(cmd)
wfile.write('%s\n' % pbs_header1(cpu_no['basecall'], cmd2))
wfile.close()
subprocess.call('qsub %s/01.BaseCall.%s.sh' % (log_dir, wait_id), shell=True, stdout=subprocess.PIPE)

print_info('Submitting job for Illumina Basecalls to SAM is done.')
qsub_list = []
sample_list.sort()

for x in sample_list:

	########### Sam to Fastq #####################
	wait_id1 += 1
	wfile = open('%s/02.SamToFastq.%s.%s.sh' % (log_dir, x, wait_id1),'w')
	wfile.write('%s\n' % pbs_header2(cpu_no['SamToFastq'], 'no', wait_id, picard.SamToFastq(x)))
	wfile.close()

	qsub_key = 'SamToFastq_%s' % x
	qsub_id[qsub_key] = wait_id1
	qsub_list.append('qsub %s/02.SamToFastq.%s.%s.sh' % (log_dir, x, wait_id1))



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



############# aggregation ###########################
print_info('Submitting job for aggregation.')

wait_id4 = wait_id3
qsub_id_callingmutation = []
qsub_list = []
agg_num = 0

for z in range(1, len(sample_dic)+1):
	y = str(z)

	### normal-tumor paired mode
	if len(sample_dic[y]) == 2:
		wait_id4 += 1

		tumor = sample_dic[y]['TUMOR']
		normal = sample_dic[y]['NORMAL']
	
		wfile = open('%s/05.aggregation.%s.%s.sh' % (log_dir, tumor, wait_id4),'w')

		qsub_wait = '%s,%s' % (qsub_id['postmapping_%s' % normal], qsub_id['postmapping_%s' % tumor])
		wfile.write('%s\n' % pbs_header2(cpu_no['agg'], 'yes', qsub_wait, aggregation(tumor, normal)))
	
		wfile.close()
		#q_dic[z] = 'qsub %s/05.aggregation.%s.%s.sh' % (log_dir, tumor, wait_id4)
		qsub_list.append('qsub %s/05.aggregation.%s.%s.sh' % (log_dir, tumor, wait_id4))

		qsub_key = 'agg_%s' % y
		qsub_id[qsub_key] = wait_id4
		agg_num += 1


	### tumor only mode
	elif len(sample_dic[y]) == 1 and sample_dic[y].keys()[0] == 'TUMOR':
		tumor = sample_dic[y]['TUMOR']

		wait_id4 += 1
		wfile = open('%s/05.aggregation.%s.%s.sh' % (log_dir, tumor, wait_id4),'w')
	
		qsub_wait = qsub_id['postmapping_%s' % tumor]
		wfile.write('%s\n' % pbs_header2(cpu_no['agg'], 'no', qsub_wait, aggregation(tumor)))
	
		wfile.close()
		qsub_list.append('qsub %s/05.aggregation.%s.%s.sh' % (log_dir, tumor, wait_id4))
		#q_dic[z] = 'qsub %s/05.aggregation.%s.%s.sh' % (log_dir, tumor, wait_id4)

		qsub_key = 'agg_%s' % y
		qsub_id[qsub_key] = wait_id4
		agg_num += 1

for x in qsub_list:
	subprocess.call(x, shell=True, stdout=subprocess.PIPE)



################ calling mutation #####################
print_info('Submitting job for calling mutation.')
wait_id5 = wait_id4

qsub_list = []
for z in range(1, len(sample_dic)+1):
	y = str(z)


	### normal-tumor paired mode
	if len(sample_dic[y]) == 2:
		tumor = sample_dic[y]['TUMOR']
		normal = sample_dic[y]['NORMAL']

		#wait_id5 = wait_id4 + agg_num
		wait_id5 = wait_id5 + 1

		wfile = open('%s/06.CallingMutation.%s.%s.sh' % (log_dir, tumor, wait_id5),'w')
	
		qkey = 'agg_%s' % y
		wfile.write('%s\n' % pbs_header2(cpu_no['callingmutation'], 'no', qsub_id[qkey], callingmutation(tumor, normal)))
		
		wfile.close()
		qsub_list.append('qsub %s/06.CallingMutation.%s.%s.sh' % (log_dir, tumor, wait_id5))
		#q_dic[z+len(sample_dic)] = 'qsub %s/06.CallingMutation.%s.%s.sh' % (log_dir, tumor, wait_id5)

		qsub_key = 'callingmutation_%s' % y
		qsub_id[qsub_key] = wait_id5
		qsub_id_callingmutation.append(str(wait_id5))


	### tumor only mode
	elif len(sample_dic[y]) == 1 and sample_dic[y].keys()[0] == 'TUMOR':
	
		tumor = sample_dic[y]['TUMOR']

		#wait_id5 = wait_id4 + agg_num
		wait_id5 = wait_id5 + 1

		wfile = open('%s/06.CallingMutation.%s.%s.sh' % (log_dir, tumor, wait_id5),'w')
	
		qkey = 'agg_%s' % y
		wfile.write('%s\n' % pbs_header2(cpu_no['callingmutation'], 'no', qsub_id[qkey], callingmutation(tumor)))

		wfile.close()
		qsub_list.append('qsub %s/06.CallingMutation.%s.%s.sh' % (log_dir, tumor, wait_id5))
		#q_dic[z+len(sample_dic)] = 'qsub %s/06.CallingMutation.%s.%s.sh' % (log_dir, tumor, wait_id5)

		qsub_key = 'callingmutation_%s' % y
		qsub_id[qsub_key] = wait_id5
		qsub_id_callingmutation.append(str(wait_id5))

for x in qsub_list:
	subprocess.call(x, shell=True, stdout=subprocess.PIPE)



###### calling SV #################
print_info('Submitting job for calling SV.')

qsub_id_all = ','.join(qsub_id_callingmutation)
wait_id6 = wait_id5 + 1
wfile = open('%s/07.01.callingSV_prep.%s.sh' % (log_dir, wait_id6),'w')
wfile.write('%s\n' % pbs_header2(cpu_no['SV'], 'yes', qsub_id_all, breakmer.prep()))
wfile.close()
subprocess.call('qsub %s/07.01.callingSV_prep.%s.sh' % (log_dir, wait_id6), shell=True, stdout=subprocess.PIPE)

wait_id7 = wait_id6
qsub_id_sv = []

for z in sample_list:
	wait_id7 += 1
	wfile = open('%s/07.02.callingSV.%s.%s.sh' % (log_dir, z, wait_id7),'w')
	wfile.write('%s\n' % pbs_header2(cpu_no['SV'], 'no', wait_id6, breakmer.run_breakmer(z)))
	wfile.close()
	
	subprocess.call('qsub %s/07.02.callingSV.%s.%s.sh' % (log_dir, z, wait_id7), shell=True, stdout=subprocess.PIPE)

	qsub_key = 'SV_%s' % z
	qsub_id[qsub_key] = wait_id5

	qsub_id_sv.append(str(wait_id7))


###### calculate metrics #################
print_info('Submitting job for calculating metrics.')

qsub_id_all = ','.join(qsub_id_sv)
wait_id8 = wait_id7

qsub_id_cal_metrics = []

for z in sample_list:
	wait_id8 += 1
	wfile = open('%s/08.calculate_metrics.%s.%s.sh' % (log_dir, z, wait_id8),'w')


	qkey = 'SV_%s' % z
	wfile.write('%s\n' % pbs_header2(cpu_no['calculate_metrics'], 'no', qsub_id[qkey], calculate_metrics(z)))
	wfile.close()
	
	subprocess.call('qsub %s/08.calculate_metrics.%s.%s.sh' % (log_dir, z, wait_id8), shell=True, stdout=subprocess.PIPE)

	qsub_key = 'calculate_metrics_%s' % z
	qsub_id[qsub_key] = wait_id8

	qsub_id_cal_metrics.append(str(wait_id8))


###### Generating excel report #################
print_info('Submitting job for generating reeport.')
qsub_id_all = ','.join(qsub_id_cal_metrics)
wfile = open('%s/09.report.%s.sh' % (log_dir, options.basecall.replace('/','')),'w')
wfile.write('%s\n' % pbs_header2(cpu_no['report'], 'yes', qsub_id_all, excel_report(working_dir)))
wfile.close()
subprocess.call('qsub %s/09.report.%s.sh' % (log_dir, options.basecall.replace('/','')), shell=True, stdout=subprocess.PIPE)
print_info('Complete to submit jobs to AMC clusters.')
