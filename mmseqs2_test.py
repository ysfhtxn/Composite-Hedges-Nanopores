from Bio import SeqIO
import os

c, min_seq_id, cov_mode = 0.85, 0.9, 0

grouping_res_path = f'/data/biolab-nvme-pool1/zhaoxy/guppy_basecalled/smeseqs_final/merge_smeseqs_time/merge_0_0_120mins.fasta'

grouping_res = {}
# grouping_res_path = f'/data/biolab-nvme-pool1/zhaoxy/guppy_basecalled/smeseqs_final/merge_smeseqs_time/merge_{indx}_0_{time_node}mins.fasta'

if len([item.seq for item in SeqIO.parse(grouping_res_path, 'fasta')]) == 0:
    print('No seq found, skip this index.\n')


mmseqs_command = f'mmseqs easy-cluster {grouping_res_path} clusterRes tmp --min-seq-id {min_seq_id} -c {c} --cov-mode {cov_mode} --rescore-mode 3 --threads 24 > test.log 2>&1'

os.system(command = mmseqs_command)

print('yes')