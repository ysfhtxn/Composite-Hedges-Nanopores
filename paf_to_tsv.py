import os

# input .paf files path
input_paf_paths = [
    r'/home/lijy/workspace/dna_storage/sequencing_qc/fqy_output/',
]


head = ['Query_sequence_name', 'Query_sequence_length', 'Query_start', 'Query_end', 'Relative_strand',
        'Target_sequence_name', 'Target_sequence_length', 'Target_start_on_original_strand',
        'Target_end_on_original_strand', 'Number_of_residue_matches', 'Alignment_block_length', 'Mapping_quality',
        'NM', 'ms', 'AS', 'nn', 'tp', 'cm','s1','s2','de','rl','cg']

for inputpath in input_paf_paths:
	filelist = os.listdir(inputpath)
	file_raw_list = list(filter(lambda filename: filename[-4:] == '.paf', filelist))
	for i in range(len(file_raw_list)):
		curfile = file_raw_list[i]
		print('reading %s' % os.path.join(inputpath, curfile))
		with open(os.path.join(inputpath, curfile), 'r+', encoding='utf-8') as f:
			old = f.read()
			f.seek(0)
			for j in range(len(head)):
				if j != len(head)-1:
					f.write(head[j] + '\t')
				else:
					f.write(head[j] + '\n')
			f.write(old)