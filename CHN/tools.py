import subprocess, os, copy, random
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

def batch_assembly(fasta_path, result_path):
    batch_size = 50
    batch_output = []
    all_record = read_file(fasta_path)
    random.shuffle(all_record)

    selected_elements = all_record
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
    consensus = get_consensus(batches_output_path, sigma=4)
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