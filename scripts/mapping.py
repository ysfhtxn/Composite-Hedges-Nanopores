import os

# reference .fasta file
fapath = r'/home/lijy/workspace/fasta/barcodeljy_384.fasta'

# input basecalled .fastq path
inputpath = r'/data/biolab-nvme-pool1/fanqy/tmp/bonito-train/fastq/'

# output .paf files path
outputpath = r'./fqy_output/'

filelist = os.listdir(inputpath)
os.makedirs(outputpath, exist_ok=True)

# minimap2 path
minimap2_path = r'/home/lijy/workspace/minimap2-2.24/minimap2'

file_raw_list = list(
    filter(lambda filename: filename[-len('.fastq'):] == '.fastq', filelist))

# custom parameters
k_mer_size = 12 #12
minimizer_window_size = 8 #8
minimal_chaining_score = 15 #15
n_discard_chains_minimizers_le = 2
minimal_peak_DP_alignment_score_to_output  = 15

# config file to save custom parameters info
f = open(os.path.join(outputpath,'config.txt'), 'w')
f.write('k_mer_size: %d\n' % k_mer_size)
f.write('minimizer_window_size: %d\n' % minimizer_window_size)
f.write('minimal_chaining_score: %d\n' % minimal_chaining_score)
f.write('n_discard_chains_minimizers_le: %d\n' % n_discard_chains_minimizers_le)
f.write('minimal_peak_DP_alignment_score_to_output: %d\n' % minimal_peak_DP_alignment_score_to_output)
f.close()

# CPU occupation
threads = 8


for i in range(len(file_raw_list)):
    curfile = file_raw_list[i]
    # default preset params
    # command = '%s -cx map-ont %s %s%s > %s%s.paf' % (minimap2_path, fapath,inputpath,curfile,outputpath,curfile)
    
    # default custom params
    command = f"{minimap2_path} -c -k{k_mer_size} -w{minimizer_window_size} -n{n_discard_chains_minimizers_le} -m{minimal_chaining_score} " \
              f"-s{minimal_peak_DP_alignment_score_to_output} -t{threads} {fapath} {inputpath+curfile} > {outputpath+curfile}.paf"

    print('processing %s' % curfile)
    os.system(command)

	# break for check
    # break


