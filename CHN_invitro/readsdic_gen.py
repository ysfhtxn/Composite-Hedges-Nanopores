import os
import pandas as pd
import pickle
import tqdm

# input paf files path
input_paf_paths = [
    r'/home/lijy/workspace/dna_storage/sequencing_qc/fqy_output/'
    ]

# output dictionary file
output_file = r'./fqy_smeseqs_q15_barcode_k12w8m15n2s15_mq5.P'

reads_hash_dic = {}

# depend on minimap2 results to filter reads
mapping_qulity_threshold = 5

for inputpath in input_paf_paths:
    print('processing %s' % inputpath)
    filelist = os.listdir(inputpath)
    file_raw_list = list(
        filter(lambda filename: filename[-4:] == '.paf', filelist))
    file_raw_list = list(
        filter(lambda filename: filename[-4:] == '.paf', filelist))
    with tqdm.trange(len(file_raw_list), desc='reading pafs') as tbar:

        for i in range(len(file_raw_list)):
            curfile = file_raw_list[i]
            data = pd.read_csv(os.path.join(inputpath, curfile), sep='\t', on_bad_lines='skip')
            mapping_q = data['Mapping_quality'].to_list()
            query_s = data['Query_start'].to_list()
            query_e = data['Query_end'].to_list()
            query_len = data['Query_sequence_length'].to_list()
            query_name = data['Query_sequence_name'].to_list()
            target_name = data['Target_sequence_name'].to_list()
            relative_strand = data['Relative_strand'].to_list()

            for j in range(len(mapping_q)):
                if mapping_q[j] >= mapping_qulity_threshold:
                    reads_hash_dic[query_name[j]] = {
                        'Q_s':query_s[j], 'Q_e':query_e[j], 'Q_l':query_len[j],
                        'T_n':target_name[j], 'M_q':mapping_q[j], 'R_s':relative_strand[j]
                    }
            tbar.update()


print('writing to file: %s' % output_file)

cnt=0
for _ in reads_hash_dic:
    cnt += 1
print('reads cnt:', cnt)

df = open(output_file, 'wb')
pickle.dump(reads_hash_dic, df)
df.close()