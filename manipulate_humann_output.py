#need to create tables to import into R from Humann output
#importing relative abundance table combined
#convert into R long file

#rename sample names to something more readable and to match another metadata file
#for each taxonomic level
#remove possible contaminants that are less than x% of community
#recalculate rel. abundance or just add to other?


import os, sys
usage = 'usage error input'

argv = sys.argv[1:]

if len(argv) == 3:
	input_file = argv.pop(0) 
	metadata = argv.pop(0)
	output = argv.pop(0)
else:
	sys.exit(usage)

minimum_cutoff = 0.002 #filter out taxa below this cutoff as contaminants/spurious
#cutoffs From Karstens et al. 2016
mean_cutoff = 0 #0.5 need to include all for PCA/NMDS
samples = []

#this would probably be better with pandas

relab_dict = {}
#each taxonomic level
relab_dict['k'],relab_dict['p'],relab_dict['c'],relab_dict['o'],relab_dict['f'],relab_dict['g'],relab_dict['s'],relab_dict['t'] = \
{},{},{},{},{},{},{},{}
meta_dict = {}

contam_dict ={}
contam_dict['k'],contam_dict['p'],contam_dict['c'],contam_dict['o'],contam_dict['f'],contam_dict['g'],contam_dict['s'],contam_dict['t'] = \
{},{},{},{},{},{},{},{}
with open (input_file,'r') as relab_combined:
	for ln in relab_combined:
		if '#SampleID' in ln:
			samples = ln.split()
			samples.pop(0)
			for s in samples:
				for i in relab_dict:
					if s not in relab_dict[i]:
						relab_dict[i][s] = {}
						contam_dict[i][s] = 0.0
						relab_dict[i][s]['Other'] = float(0.0)
					else:
						print("Sample names duplicated: " + s)
		else:
			abundances = ln.split()				
			taxa = abundances.pop(0)
			taxa_final = taxa.split('|')[-1]
			taxa_level = taxa_final.split('__')[0]
			abundances = [float(i) for i in abundances]
			#calculate mean abundance, so those below mean abundance are added to 'Other'
			mean_abundance = sum(abundances)/len(abundances)
			for i in range(len(abundances)):
				current_percent = abundances[i]
				if current_percent < minimum_cutoff: #eventually will need to update to remove and renormalize
					contam_dict[taxa_level][samples[i]]+= current_percent
					relab_dict[taxa_level][samples[i]][taxa_final] = 0.0

				#if mean_abundance < mean_cutoff: #should be elif if going to use both
				#	relab_dict[taxa_level][samples[i]]['Other'] += current_percent
				else:
					if taxa_final not in relab_dict[taxa_level][samples[i]]:
						relab_dict[taxa_level][samples[i]][taxa_final] = 0.0
					relab_dict[taxa_level][samples[i]][taxa_final] = current_percent

#for d in relab_dict:
#	df = pd.DataFrame.from_dict(d)
#	df.columns = d.keys()
#	rows = 


categories = []
with open (metadata,'r') as meta:
	for ln in meta:
		if 'SampleID' in ln:
			categories = ln.split()
		else:
			meta_info = ln.split()
			sampleID = meta_info.pop(0)
			#meta_dict[sampleID] = []
			meta_dict[sampleID] = meta_info
categories.pop(0)
for x in relab_dict:
	with open(output + '_' + x + '.tsv','w') as out:
		c = '\t'.join(categories)
		out.write('{0}\t{1}\t{2}\t{3}\n'.format('SampleID','Taxa','Percent',c))
		for i in relab_dict[x]:
			for j in relab_dict[x][i]:
				if contam_dict[x][i] != 0.0:
					new_percent = 100*(relab_dict[x][i][j] / (1-contam_dict[x][i]))
				else:
					new_percent = 100*relab_dict[x][i][j]
				m = '\t'.join(meta_dict[i])
				out.write('{0}\t{1}\t{2}\t{3}\n'.format(i,j,new_percent,m))
		
				


			
		
