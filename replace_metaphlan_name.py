#rename metaphlan output with filename as sample name in header!
import sys
import os
import shutil

usage = """<folder of input folders>
	<output folder>"""
argv = sys.argv[1:]
if len(argv) == 1:
	input_folder = argv.pop(0)
	#output = argv.pop(0)
else:
	sys.exit(usage)

allfiles = os.listdir(input_folder)
#os.system('mkdir metaphlancpmrename')
for f in allfiles:
	newname = f.split('_metaphlan')[0]
	first_row=True
	with open (input_folder + f) as sourcefile:
		lns = sourcefile.readlines()
		lns[0] = lns[0].replace('_Abundance-RELAB','')
		#lns[0] = lns[0].replace('Metaphlan2_Analysis-RELAB',newname)
	with open(input_folder + f,'w') as sourcefile:	
		sourcefile.writelines(lns)






#for ln in sourcefile:
#			if first_row:
#				ln = ln.replace('Metaphlan2_Analysis-CPM',newname)
#				first_row=False
#			sourcefile.write(ln)
				



	#from_file = open(input_folder + f) 
	#line = from_file.readline()
	# make any changes to line here
	#replacementline = line.replace('Metaphlan2_Analysis',newname)
	#to_file = open(input_folder + f,mode="w")
	#to_file.write(replacementline)
	#shutil.copyfileobj(from_file,to_file)
	#from_file.close()
	#to_file.close()
