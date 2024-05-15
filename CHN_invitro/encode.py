import CHNcodec, tqdm, random
from reedsolo import RSCodec
import numpy as np

def file_split(infile, fragment_size=36):
    test_file = open(infile, 'rb').read()
    splited_file = []
    for i in range(0, len(test_file), fragment_size):
        fragment = test_file[i:i + fragment_size]
        if len(fragment) < fragment_size:
            fragment += b' ' * (fragment_size - len(fragment))
        splited_file.append(fragment)
    return splited_file

if __name__ == '__main__':
    NEED_MAP = False

    infile_path = './examples/sme introduction.txt'
    splited_file = file_split(infile_path)

    basic_rsc = RSCodec(40 - 36, nsize=40)
    chn = CHNcodec.chn(resolution=2, sigma=8, step=2)
    dnas = []
    with tqdm.trange(len(splited_file), desc='Encoding') as pbar:
        for msg in splited_file:
            encoded_str = basic_rsc.encode(msg)
            encode_vbits = CHNcodec.unpack_vbits(encoded_str, chn.pattern)
            encoded_str = chn.encode_c(encode_vbits)
            result_string = chn.to_combine(encode_res=encoded_str)
            mapped_str = chn.mapping(result_string, 8)
            # dnas.append(result_string)
            dnas += mapped_str if NEED_MAP else [result_string]
            pbar.update()

    with open('./examples/encoded.fasta', 'w') as f:
        for idx, dna in enumerate(dnas):
            f.write(f'>index_{idx}\n')
            f.write(f'{dna}\n')