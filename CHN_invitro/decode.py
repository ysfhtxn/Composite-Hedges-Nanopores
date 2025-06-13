from itertools import product
from reedsolo import RSCodec
import time, os
import CHNcodec, tools
from CHNcodec import chn


if __name__ == '__main__':
    resolution = 2
    sigma = 8
    step = 2
    N, K = 40, 36

    rsc = RSCodec(N - K, nsize=N)

    CHN = chn(resolution=resolution, sigma=sigma, step=step)
    CHN.alphabet_show()

    txt_decoded_file_path = r'your/txt/decode/file/path'

    # mmseqs easy-cluster parameters
    cs= [0.85, 0.90, 0.80, 0.95, 0.75]
    min_seq_ids = [0.9, 0.95, 0.85]
    cov_modes = [0,3,1,2]
    mmseqs_paras_combinations = list(product(cs, min_seq_ids, cov_modes))

    decode_fail_num = 0
    fail_index = []

    txt_decode_res = []

    grouping_res_dir = r'your/fasta/dir/of/grouping/res/'
    
    filelist = os.listdir(grouping_res_dir)
    file_list = list(
                filter(lambda filename: os.path.split(filename)[-1].split('.')[-1] == 'fastq', filelist))
    grouping_res_num = len(file_list)

    ref_num = 7 # which can be changed

    print('Decoding...')
    start_time = time.time()
    for grouping_res_indx in range(ref_num):
        grouping_res_path = os.path.join(grouping_res_dir, f'grouping_res_{grouping_res_indx}.fasta')
        for mmseqs_paras_combination in mmseqs_paras_combinations:
            c, min_seq_id, cov_mode = mmseqs_paras_combination
            
            consensus_str = tools.mmseqs(grouping_res_path, c, min_seq_id, cov_mode, batch_size=50)
            if consensus_str == 'A'*243:
                continue

            consensus_payload_with_anchors = consensus_str[57:-11]
            consensus_payload = consensus_payload_with_anchors[:40]+consensus_payload_with_anchors[45:85]+consensus_payload_with_anchors[90:130]+consensus_payload_with_anchors[135:]
            data = CHN.combine_to_letters(consensus_payload)
            res = CHN.decode_c(data)
            # print('decode res nums: %d' % len(res))
            decode_vbits, min_penalty = CHN.hypo_backtrack(res)
            # print(len(decode_vbits),min_penalty)
            # print('same min penalty num: %d' % len(decode_vbits))
            # print('penalty: %f' % min_penalty)

            # randomly choose a min penalty decode_vbits
            if decode_vbits == []:
                rs_decode_res = bytearray(K)
            for idx, item in enumerate(decode_vbits):
                try:
                    msg = CHNcodec.pack_vbits(item, CHN.pattern)
                    rs_decode_res = rsc.decode(bytearray(msg))[0]
                    break
                except:
                    rs_decode_res = bytearray(K)
                    continue
            if bytearray(K) == rs_decode_res:
                print(f'{grouping_res_indx} Failed')
                decode_fail_num += 1
                fail_index.append(grouping_res_indx)
                continue
            txt_decode_res.append(rs_decode_res)

    end_time = time.time()
    print('Decoding finished.')
    print(f'Total time cost: {(end_time-start_time)}')

    with open(txt_decoded_file_path, 'wb') as f:
        for item in txt_decode_res:
            f.write(item)

    print(decode_fail_num)
    print(fail_index)
