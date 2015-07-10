import MySQLdb, re, datetime, os, sys
from utils import *

class extract:
	def __init__(self, basecall):
		self.basecall = '/home/D131748/Illumina/%s' % basecall
		self.working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % basecall
		self.working_tmp = '%s/tmp' % self.working_dir
		self.conf = '%s/configure' % self.working_dir
		print self.working_dir
		print self.working_tmp
		print self.conf

		check_dir(self.working_tmp)
		check_dir(self.conf)


	def ext_info(self):
		## connet to the mysql server
		db = MySQLdb.connect(host='localhost', user='kdh', passwd='kdh', db='cpcm')
		cursor = db.cursor(MySQLdb.cursors.DictCursor)
		
		## read files for extration of information
		try:
			run_info = open('%s/runParameters.xml' % (self.basecall)).read()
		except:
			run_info = open('%s/RunParameters.xml' % (self.basecall)).read()
##	for resequencing
		
		sample_sheet = open('%s/SampleSheet.csv' % (self.basecall)).read()
		
		application = re.search('Workflow,(.*)\r\n', sample_sheet).group(1)
		if ',' in application:
			application = application.replace(',','')

		print "Note - Workflow mode was selected as %s." % application

		if application == 'Resequencing':
			sample_info = open('%s/SampleSheet.csv' % (self.basecall)).read().split('\r\n')[23:-1]
		elif application == 'GenerateFASTQ':
			sample_info = open('%s/SampleSheet.csv' % (self.basecall)).read().split('\r\n')[21:-1]
		else:
			print 'Warning - It has an error in extracting sample sheet and run information.'
			sys.exit()

		sample_info2 = open('%s/SampleSheet.csv' % (self.basecall)).read()
		
		## extract instrument information
		instrument = re.search('<ScannerID>(.*)</ScannerID>',run_info).group(1)
		
		## extract flowcell information
		(flowcell_id, flowcell_PN, flowcell_exp) = re.search('<FlowcellRFIDTag>\s*<SerialNumber>(.*)</SerialNumber>\s*<PartNumber>(\d+)</PartNumber>\s*<ExpirationDate>(.*)</ExpirationDate>',run_info).groups()
		flowcell_exp2 = flowcell_exp.split('T')[0]
	
		## extract reegent information
		(PR2_SN, PR2_PN, PR2_exp) = re.search('<PR2BottleRFIDTag>\s*<SerialNumber>(.*)</SerialNumber>\s*<PartNumber>(\d+)</PartNumber>\s*<ExpirationDate>(.*)</ExpirationDate>',run_info).groups()
		(reagentkit_SN, reagentkit_PN, reagentkit_exp) = re.search('<ReagentKitRFIDTag>\s*<SerialNumber>(.*)</SerialNumber>\s*<PartNumber>(\d+)</PartNumber>\s*<ExpirationDate>(.*)</ExpirationDate>', run_info).groups()
		PR2_exp2 = PR2_exp.split('T')[0]
		reagentkit_exp2 = reagentkit_exp.split('T')[0]
	
		## extract run id information
		run_id = re.search('<RunID>(.*)</RunID>',run_info).group(1)
		run_num = re.search('<RunNumber>(\d+)</RunNumber>',run_info).group(1)
		
		## extract read length
		R1_len = re.search('<RunInfoRead Number="1" NumCycles="(\d+)" IsIndexedRead="N" />',run_info).group(1)
		index_len = re.search('<RunInfoRead Number="2" NumCycles="(\d+)" IsIndexedRead="Y" />',run_info).group(1)
		R2_len = re.search('<RunInfoRead Number="3" NumCycles="(\d+)" IsIndexedRead="N" />',run_info).group(1)
		
		## date run
		run_date = re.search('<RunStartDate>(\d+)</RunStartDate>',run_info).group(1)
		application_ver= re.search('<ApplicationVersion>(.*)</ApplicationVersion>',run_info).group(1)
		FPGA_ver = re.search('<FPGAVersion>(.*)</FPGAVersion>',run_info).group(1)
		MCS_ver = re.search('<MCSVersion>(.*)</MCSVersion>',run_info).group(1)
		RTA_ver = re.search('<RTAVersion>(.*)</RTAVersion>',run_info).group(1)
		experiment_name = re.search('<ExperimentName>(.*)</ExperimentName>',run_info).group(1)
		
		## generate param.txt, barcode.txt and sample.txt files
		dic = {}
		dic_project = {}
		dic_type = {}
		dic_run = {}

		#conf = '%s' % config_dir

		for x in sample_info:
			print x
			x_1 = x.split(',')
			library_name = x_1[0].rstrip()
			index_seq = x_1[5]
			if application == 'Resequencing':
				sample_project = x_1[7]
				sample_type = x_1[8]
			elif application == 'GenerateFASTQ':
				sample_project = x_1[6]
				sample_type = x_1[7]

			if sample_project == '' or sample_type == '':
				print 'Please make sure project number and type of samples in sheet\n'
				sys.exit()

			dic[index_seq] = library_name
			dic_project[index_seq] = sample_project
			dic_type[index_seq] = sample_type

		try:
			investigator = re.search('Investigator Name,(.*)\r\n', sample_info2).group(1)
		except:
			investigator = 'NA'
		date_sequenced = re.search('Date,(.*)\r\n', sample_info2).group(1)
		if ',' in investigator:
			investigator = investigator.replace(',','')
		if ',' in date_sequenced:
			date_sequenced = date_sequenced.replace(',','')

		i = datetime.datetime.now()
		date_analysed = '%s-%s-%s' % (i.year, i.month, i.day)

		wfile = open('%s/barcode.txt' % self.conf,'w')
		wfile2 = open('%s/param.txt' % self.conf,'w')
		wfile3 = open('%s/sample.txt' % self.conf,'w')

		wfile.write('barcode_sequence_1\tBarcode_sequence2\tbarcode_name\tlibrary_name\n')
		wfile2.write('OUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\tBARCODE_1\n%s/undetermined.bam\tundetermined\tundetermined\tN\n' % self.working_tmp)
		wfile3.write('#bam\ttype\tclasification_chr\n')
	
		cursor.execute('select * from barcodes')
		recs = cursor.fetchall()
		adaptor_num = 1
		
		for y in recs:
			if dic.has_key(y['barcode_seq']):
				wfile.write('%s\t\t%s\t%s_%s\n' % (y['barcode_seq'], y['barcode_id'], dic[y['barcode_seq']].rstrip().split('_')[0], dic[y['barcode_seq']].rstrip().split('_')[2]))
				wfile2.write('%s/%s.bam\t%s_%s\t%s_%s\t%s\n' % (self.working_tmp, dic[y['barcode_seq']].rstrip(), dic[y['barcode_seq']].rstrip().split('_')[0], dic[y['barcode_seq']].rstrip().split('_')[2], dic[y['barcode_seq']].rstrip().split('_')[0], dic[y['barcode_seq']].rstrip().split('_')[2], y['barcode_seq']))
				wfile3.write('%s.bam\t%s\t%s\n' % (dic[y['barcode_seq']].rstrip(), dic_project[y['barcode_seq']], dic_type[y['barcode_seq']]))
				cursor.execute("select * from samples_sequenced where run_id = '%s' and bam = '%s.bam'" % (run_id, dic[y['barcode_seq']].rstrip()))
				check_run = cursor.fetchall()
				
				if len(check_run) == 0:
					cmd_sample = "INSERT INTO samples_sequenced(`run_id`, `bam`, `sample_project`, `sample_type`, `date_sequenced`, `date_analysed`, `investigator`) VALUES('%s', '%s.bam', '%s', '%s', '%s', '%s', '%s')" % (run_id, dic[y['barcode_seq']].rstrip(), dic_project[y['barcode_seq']], dic_type[y['barcode_seq']].strip().upper(), date_sequenced, date_analysed, investigator)
					cursor.execute(cmd_sample)
				else:
					pass

			else:
				adaptor_name = 'Unused%s' % adaptor_num
				wfile.write('%s\t\t%s\tUnused\n' % (y['barcode_seq'], y['barcode_id']))
				wfile2.write('%s/%s.bam\t%s\t%s\t%s\n' % (self.working_tmp, adaptor_name.rstrip(), adaptor_name.rstrip(), adaptor_name.rstrip(), y['barcode_seq']))
				adaptor_num += 1
	
		wfile.close()
		wfile2.close()
		wfile3.close()

		cursor.execute("select * from instrument where instrument_name = '%s'" % instrument)
		recs = cursor.fetchall()
		if len(recs) == 0:
			cursor.execute("INSERT INTO instrument (`instrument_name`) VALUES('%s')" % instrument)
			cursor.execute("select * from instrument where instrument_name = '%s'" % instrument)
			recs = cursor.fetchall()
			instrument_id = recs[0]['id']
		else:
			instrument_id = recs[0]['id']
		
		cursor.execute("select * from sequencing_run_information where run_id = '%s'" % run_id)
		run = cursor.fetchall()
		
		if len(run) == 0:
			cmd = "INSERT INTO sequencing_run_information(`run_number`, `run_id`, `R1_length`, `index_length`, `R2_length`, `date_run`, `instrument_id`, `flowcell_id`, `flowcell_PN`, `flowcell_exp`, `PR2_SN`, `PR2_PN`, `PR2_exp`, `reagentkit_SN`, `reagentkit_PN`, `reagentkit_exp`, `application_ver`, `FPGA_ver`, `MCS_ver`, `RTA_ver`, `experiment_name`, `investigator`) VALUES('%d', '%s', '%d', '%d', '%d', '%s', '%d', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')" % (int(run_num), run_id, int(R1_len), int(index_len), int(R2_len), date_sequenced, instrument_id, flowcell_id, flowcell_PN, flowcell_exp2, PR2_SN, PR2_PN, PR2_exp2, reagentkit_SN, reagentkit_PN, reagentkit_exp2, application_ver, FPGA_ver, MCS_ver, RTA_ver, experiment_name, investigator)
			cursor.execute(cmd)
		else:
			pass
	
		db.commit()
		cursor.close()

		print 'Configuration files are generated in %s/configure directory!!!' % self.working_tmp


		barcode_file = '%s/barcode.txt' % self.conf
		sample_file = '%s/sample.txt' % self.conf
		param_file = '%s/param.txt' % self.conf

		check_file(barcode_file)
		check_file(sample_file)
		check_file(param_file)


