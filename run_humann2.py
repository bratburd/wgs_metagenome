#run humann2 over all samples
import sys
import os

usage = """<folder of input folders>"""
argv = sys.argv[1:]
if len(argv) == 1:
	input_folder = argv.pop(0)
else:
	sys.exit(usage)

output_folder_humann = '/home/jenny/Desktop/wls_resistance/metagenome2019/humann/'
files = os.listdir(input_folder)

#specify number of cores to use???

#extract bugs list file and metaphlan log file for later use, then zip output and delete unzipped version
for f in files:
	os.system("humann2 --input {0}{1} --output {2} --metaphlan /home/jenny/tools/metaphlan2/ --diamond /home/jenny/tools/diamond0822/".format(input_folder, f, output_folder_humann))
