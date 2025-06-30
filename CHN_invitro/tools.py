import subprocess, os, re, copy, random
from Bio import SeqIO
from scipy.spatial import distance

def read_file(path):
    with open(path) as handle:
        record = SeqIO.parse(handle,"fasta")
        
        all_record = []
        for lin in record:
            all_record.append((str(lin.id), str(lin.seq)))
        # print(len(all_record))
    return all_record

def assembly(fasta_path, result_path):
    p = subprocess.Popen(('muscle -align {} -output {}').format(fasta_path, result_path),
                        shell = True,
                        stdout = subprocess.PIPE,
                        stderr = subprocess.STDOUT)
    for line in p.stdout.readlines():
        print('',end = '')

def batch_assembly(fasta_path, result_path, batch_size=50):
    # batch_size = 50
    batch_output = []
    all_record = read_file(fasta_path)
    random.shuffle(all_record)
    for i in range((len(all_record)//batch_size)+1):
        selected_elements = all_record[batch_size*i : batch_size*(i+1)]
        temp_batch_path = os.path.split(result_path)[0] + r'/temp_batch.fasta'
        temp_batch_output_path = os.path.split(result_path)[0] + r'/temp_batch_output.fasta'
        with open(temp_batch_path, 'w')as f:
            for item in selected_elements:
                f.write(f'>{item[0]}\n')
                f.write(f'{item[1]}\n')
        assembly(temp_batch_path, temp_batch_output_path)
        consensus = get_consensus(temp_batch_output_path)
        batch_output.append(consensus)

    batches_path = os.path.split(result_path)[0] + r'/temp_batches.fasta'
    batches_output_path = os.path.split(result_path)[0] + r'/temp_batches_output.fasta'
    with open(batches_path, 'w')as f:
        for idx, item in enumerate(batch_output):
            f.write(f'>batch_consensus_{idx}\n')
            f.write(f'{item}\n')
    assembly(batches_path, batches_output_path)
    consensus = get_consensus(batches_output_path)
    return consensus

    
def get_consensus(result_path, sigma=4):
    assembly_seqs = [str(i.seq) for i in SeqIO.parse(result_path, "fasta")]
    seq_len = max([len(item) for item in assembly_seqs])
    ratio_list = []
    ratio_div = []
    for idx in range(len(assembly_seqs)):
        assembly_seqs[idx] = assembly_seqs[idx]+'-'*(seq_len-len(assembly_seqs[idx])) if len(assembly_seqs[idx]) < seq_len else assembly_seqs[idx]
    for m in range(min([len(item) for item in assembly_seqs])):
        temp = []
        for n in range(len(assembly_seqs)):
            temp.append(assembly_seqs[n][m])
        lis = str(temp)
        ratio_list.append([lis.count('A'),lis.count('C'),lis.count('G'),lis.count('T')])
        ratio_div.append(lis.count('-'))

    minus_num = seq_len - 243 
    t = copy.deepcopy(ratio_div)
    max_index = []
    for _ in range(minus_num):
        number = max(t)
        index = t.index(number)
        t[index] = 0
        max_index.append(index)

    consensus = []
    for idx , ratio in enumerate(ratio_list):
        if idx in max_index:
            continue
        res_A = distance.jensenshannon([1,0,0,0],ratio)
        res_C = distance.jensenshannon([0,1,0,0],ratio)
        res_G = distance.jensenshannon([0,0,1,0],ratio)
        res_T = distance.jensenshannon([0,0,0,1],ratio)
        res_M = distance.jensenshannon([0.5,0.5,0,0],ratio)
        res_K = distance.jensenshannon([0,0,0.5,0.5],ratio)
        res_R = distance.jensenshannon([0.5,0,0.5,0],ratio)
        res_Y = distance.jensenshannon([0,0.5,0,0.5],ratio)

        res = [res_A,res_C,res_G,res_T,res_M,res_K,res_R,res_Y]

        min_index = res[:sigma].index(min(res[:sigma]))

        index_to_char = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'M', 5: 'K', 6: 'R', 7: 'Y',}

        if min_index in index_to_char:
            consensus.append(index_to_char[min_index])

    consensus_str = ''.join(m for m in consensus)
    return consensus_str

def mmseqs(grouping_res_path, c, min_seq_id, cov_mode, batch_size):
    os.system(command = 'rm -rf ./tmp/*')

    grouping_res = {}
    
    if len([item.seq for item in SeqIO.parse(grouping_res_path, 'fasta')]) == 0:
        print('No seq found, skip this index.\n')
        return 'A'*243
    
    mmseqs_command = f'mmseqs easy-cluster {grouping_res_path} clusterRes tmp --min-seq-id {min_seq_id} -c {c} --cov-mode {cov_mode} --rescore-mode 3 --threads 24 > mmseqs2.log 2>&1'

    os.system(command = mmseqs_command)
    
    for item in SeqIO.parse(grouping_res_path, 'fasta'):
        grouping_res[str(item.id)] = str(item.seq)

    f = open('clusterRes_cluster.tsv', 'r')
    all_lines = [item.split('\t') for item in f.readlines()]
    f.close()

    all_cluster = list(set([item[0] for item in all_lines]))

    all_cluster_dict = {}

    for cluster_id in all_cluster:
        all_cluster_dict[cluster_id] = 0
                
    all_cluster_dict = dict(sorted(all_cluster_dict.items(), key=lambda item: item[1]))
        
    all_eligible_cluster = {}
    for cluster_id in all_cluster:
        cluster_lib = [item[-1].split('\n')[0] for item in all_lines if item[0] == cluster_id]
        all_cluster_dict[cluster_id] = cluster_lib
        if len(cluster_lib) > 10 :
            all_eligible_cluster[cluster_id] = cluster_lib
    
    
    all_eligible_cluster = dict(sorted(all_eligible_cluster.items(), key=lambda item: len(item[1]), reverse=True))

    # delete the former files
    if os.path.exists('./test/'):
        for filename in os.listdir('./test/'):
            file_path = os.path.join('./test/', filename)
            if os.path.isfile(file_path):
                os.remove(file_path)

    all_cluster_num = 10
    
    idx = 0
    for key, value in all_eligible_cluster.items():
        with open(f'./test/test_cluster_{idx}.fasta', 'w')as f:
            for i in value:
                f.write(f'>{i}\n')
                f.write(f'{grouping_res[i]}\n')
        idx += 1
        if idx >= all_cluster_num:
            break
    
    seqs_path = r'./test/'
    
    cluster_num = 0
    for item in os.listdir('./test/'):
        item_path = os.path.join('./test/', item)
        if os.path.isfile(item_path):
            cluster_num += 1

    all_res = []

    correct_count = 0
    for idx in range(cluster_num):
        print(f'Process {idx} working...')
        seq_name = seqs_path + f'test_cluster_{idx}.fasta'
        result_index = str(re.findall(r'\d', seq_name)[0])

        result_dir = r'./test_result/'
        result_path = result_dir + f'align_result_{result_index}.fasta'
        assembly_result = batch_assembly(fasta_path = seq_name, result_path = result_path, batch_size=batch_size)
        print(f'Process {idx} done.')
        all_res.append(assembly_result)

    with open('in.fasta', 'w')as f:
        for i, data in enumerate(all_res):
            f.write(f'>{i}\n')
            f.write(f'{data}\n')

    consensus_str = get_consensus('in.fasta',sigma=8)
    return consensus_str
