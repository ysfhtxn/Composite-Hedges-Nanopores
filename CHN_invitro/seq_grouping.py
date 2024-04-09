import pickle, os, tqdm
from Bio import SeqIO

def get_all_anchor(anchor):
    all_anchor = [anchor]
    # delete one
    subanchors = []
    for i in range(1,len(anchor)):
        subanchor = anchor[:i]+anchor[i+1:]
        subanchors.append(subanchor)
    all_anchor += list(set(subanchors))
    # insert one
    subanchors = []
    for i in range(1,len(anchor)):
        for insert in ['A','C','G','T']:
            subanchor = anchor[:i]+insert+anchor[i:]
            subanchors.append(subanchor)
    all_anchor += list(set(subanchors))
    return all_anchor

def search_anchor(dna, anchor, ref_position):
    # all_anchor = get_all_anchor(anchor)
    all_anchor = [anchor]
    
    true_position = [0,0,0]
    anchor_len = [0,0,0]
    
    ps_offset = 0 
    
    for ps in ref_position:
        ps += ps_offset
        for anchor in all_anchor:
            is_found = False
            search_path = [0] + [num if num == 0 else (-1) ** num * (num // 2) for num in range(2,(10+1)*2)] 
            for walk in search_path:
                if dna[ps+walk:ps+walk+len(anchor)] == anchor:
                    true_position[ref_position.index(ps-ps_offset)] = ps+walk
                    anchor_len[ref_position.index(ps-ps_offset)] = len(anchor)
                    is_found = True
                    ps_offset += walk
                    break
            if is_found:
                break
    
    return sum(1 for num in true_position if num != 0), abs(ps_offset), anchor_len

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = dna_sequence[::-1]
    complement_sequence = ''.join([complement_dict[base] for base in reverse_sequence])
    return complement_sequence

if __name__ == "__main__":
    basecalled_dir = r'your/fastq/dir/of/basecalled/reads/'

    pickle_path = r'your/pickle/path/of/pickle/file'

    bcd_path = r'your/fasta/path/of/barcodes'

    all_bcd_id = [(str(item.id)) for item in SeqIO.parse(bcd_path, 'fasta')]

    with open(pickle_path, 'rb') as f:
        pickle_data = pickle.load(f)

    filelist = os.listdir(basecalled_dir)
    fastq_list = list(
                filter(lambda filename: os.path.split(filename)[-1].split('.')[-1] == 'fastq', filelist))

    ref_num = 122

    st_flag = 'TTTCT'
    ed_flag = 'TAGAGC'

    anchor_5mer = r'your_anchor'

    ref_position = [97, 142, 187] # which is determined by the DNA structure described in the paper
    ref_dna_len = 243 # which is determined by the DNA structure described in the paper

    all_seqs = [[] for _ in range(ref_num)]

    print('Seqs Grouping working...')
    for fastq_index in tqdm.tqdm(range(len(fastq_list))):
        seqs = [(str(lin.id), str(lin.seq)) for lin in SeqIO.parse(basecalled_dir+fastq_list[fastq_index], "fastq")]
        for seq_id, seq_string in seqs:
            if seq_id not in pickle_data:
                continue
            if pickle_data[seq_id]['R_s'] == '+':
                sequence = seq_string
            if pickle_data[seq_id]['R_s'] == '-':
                sequence = reverse_complement(seq_string)

            target_name = pickle_data[seq_id]['T_n']

            st_index = sequence.find(st_flag)
            ed_index = sequence.find(ed_flag)

            if  (st_index != -1) and (ed_index != -1) and len(sequence[st_index:ed_index+len(ed_flag)]) in range(ref_dna_len-20, ref_dna_len+20):
                target_index = all_bcd_id.index(target_name)

                dna = sequence[st_index:ed_index+len(ed_flag)]

                anchor_counts, anchor_offsets, anchor_lens = search_anchor(dna, anchor_5mer, ref_position)

                if anchor_counts > 1:            
                    all_seqs[target_index].append((seq_id, sequence[st_index:ed_index+len(ed_flag)]))

    print('Seqs Grouping Done...')

    print('\n')

    print('Grouping result write working...')

    grouping_res_dir = r'your/fasta/path/of/grouping/res/'

    for idx in tqdm.tqdm(range(ref_num)):
        with open(grouping_res_dir+f'grouping_res_{idx}.fasta', 'w') as f:
            for item in all_seqs[idx]:
                f.write(f'>{item[0]}\n')
                f.write(f'{item[1]}\n')

    print('Grouping result write Done...')