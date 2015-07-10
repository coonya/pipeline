
class bwa:
	def __init__(self, params, file_dic):
		self.params = params
		self.bwa_path = self.params['bwa_path']
		self.ref = self.params['ref']

		self.file_dic = file_dic
	

	def aln1(self, sample):
		fastq = self.file_dic[sample]['fastq1']
		sai = self.file_dic[sample]['sai1']

		cmd = []

		cmd.append('%s/bwa aln' % self.bwa_path)
		cmd.append('-t 16')
		cmd.append('-q 5')
		cmd.append('-l 32')
		cmd.append('-k 2')
		cmd.append('-o 1')
		cmd.append('-f %s' % sai)
		cmd.append('%s' % self.ref)
		cmd.append('%s' % fastq)

		cmd_f = ' '.join(cmd)
		
		#print 'CMD : %s' % cmd_f

		return cmd_f

	def aln2(self, sample):
		fastq = self.file_dic[sample]['fastq2']
		sai = self.file_dic[sample]['sai2']

		cmd = []

		cmd.append('%s/bwa aln' % self.bwa_path)
		cmd.append('-t 16')
		cmd.append('-q 5')
		cmd.append('-l 32')
		cmd.append('-k 2')
		cmd.append('-o 1')
		cmd.append('-f %s' % sai)
		cmd.append('%s' % self.ref)
		cmd.append('%s' % fastq)

		cmd_f = ' '.join(cmd)
		
		#print 'CMD : %s' % cmd_f

		return cmd_f


	def sampe(self, sample):
		fastq1 = self.file_dic[sample]['fastq1']
		fastq2 = self.file_dic[sample]['fastq2']
		sai1 = self.file_dic[sample]['sai1']
		sai2 = self.file_dic[sample]['sai2']
		output = self.file_dic[sample]['sam']
		
		cmd = []

		cmd.append('%s/bwa sampe' % self.bwa_path)
		cmd.append('-f %s' % output)
		cmd.append('-P %s' % self.ref)
		cmd.append('%s' % sai1)
		cmd.append('%s' % sai2)
		cmd.append('%s' % fastq1)
		cmd.append('%s' % fastq2)

		cmd_f = ' '.join(cmd)
		
		#print 'CMD: %s' % cmd_f

		return cmd_f


