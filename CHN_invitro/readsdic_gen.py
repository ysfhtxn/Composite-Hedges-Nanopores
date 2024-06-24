import os
import pandas as pd
import pickle
import tqdm

# PAF file header
columns = ['Query_sequence_name', 'Query_sequence_length', 'Query_start', 'Query_end', 'Relative_strand',
        'Target_sequence_name', 'Target_sequence_length', 'Target_start_on_original_strand',
        'Target_end_on_original_strand', 'Number_of_residue_matches', 'Alignment_block_length', 'Mapping_quality',
        'NM', 'ms', 'AS', 'nn', 'tp', 'cm','s1','s2','de','rl','cg']

# input paf files path
input_paf_paths = r'/data/biolab-nvme-pool1/lijy/mapping_res/va_qscore/paf_smeseqs_230703_2_barcode_custom/smeseqs_230703_2_q16/'

# output dictionary file
output_file = r'/data/biolab-nvme-pool1/lijy/smeseqs/smeseqs_230703_2_q16.P'

reads_hash_dic = {}

# depend on minimap2 results to filter reads
mapping_qulity_threshold = 3

filelist = os.listdir(input_paf_paths)
file_raw_list = list(filter(lambda filename: filename[-len('.paf'):] == '.paf', filelist))[:]


with tqdm.tqdm(total=len(file_raw_list), desc=f"processing...") as pbar:
    for curfile in file_raw_list:
        data = pd.read_csv(os.path.join(input_paf_paths, curfile), sep='\t', on_bad_lines='skip')
        mapping_q = data['Mapping_quality'].to_list()
        query_s = data['Query_start'].to_list()
        query_e = data['Query_end'].to_list()
        query_len = data['Query_sequence_length'].to_list()
        query_name = data['Query_sequence_name'].to_list()
        target_name = data['Target_sequence_name'].to_list()
        relative_strand = data['Relative_strand'].to_list()

        
        for j in range(len(mapping_q)):
            if mapping_q[j] >= mapping_qulity_threshold:
                if query_name[j] not in reads_hash_dic:
                    reads_hash_dic[query_name[j]] = {
                    'Q_s':query_s[j], 'Q_e':query_e[j], 'Q_l':query_len[j],
                    'T_n':target_name[j], 'M_q':mapping_q[j], 'R_s':relative_strand[j]
                    }
                elif query_name[j] in reads_hash_dic and reads_hash_dic[query_name[j]]['M_q'] < mapping_q[j]:
                    reads_hash_dic[query_name[j]]['M_q'] = mapping_q[j]
        pbar.update()


print('writing to file: %s' % output_file)

cnt=0
for _ in reads_hash_dic:
    cnt += 1
print('reads cnt:', cnt)

df = open(output_file, 'wb')
pickle.dump(reads_hash_dic, df)
df.close()