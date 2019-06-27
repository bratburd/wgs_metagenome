
#run over all samples in a folder of folders
import sys
import os

usage = """<folder of input folders>"""
argv = sys.argv[1:]
if len(argv) == 1:
	input_folder = argv.pop(0)
	#output = argv.pop(0)
else:
	sys.exit(usage)

folders = os.listdir(input_folder)
fastq_dict = {}
for f in folders:
	prefix = f.split('_')[0]
	if prefix not in fastq_dict:
		fastq_dict[prefix] = {}
	if 'R1' in f: #hopefully don't ever have R numbers in sample names
		fastq_dict[prefix]['R1'] = f
	elif 'R2' in f:
		fastq_dict[prefix]['R2'] = f




#adapter trimming anf base pair length trimming with fastp
for p in fastq_dict:
	os.system('fastp -i {0}{1} -I {0}{2} -o fastptrimmed/{3}_R1.fastq.gz -O fastptrimmed/{3}_R2.fastq.gz -j fastp_qc/{3}.json -h fastp_qc/{3}.html'.format(input_folder, fastq_dict[p]['R1'], fastq_dict[p]['R2'],p))
	os.system('mv {0}/{1} fastqgz_raw/'.format(input_folder, p))
	
	#alternatively run with bbduk below
	#os.system('{0} in1={1}{2} in2={1}{3} out1={4}{2} out2={4}{3} ref={5} ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 qtrim=rl trimq=10 maxns=5 entropy=0.2 entropywindow=50 entropyk=5'.format(bbduk_path,output_folder_cat,fastq_dict[p]['R1'],fastq_dict[p]['R2'],output_folder_trimmed,ref))
