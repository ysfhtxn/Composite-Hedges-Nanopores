# Composite Hedges Nanopores

A package of coding algorithm for rapid readout of digital information storage in DNA, inspired by [HEDGES](https://github.com/whpress/hedges) and [composite DNA letter](https://github.com/leon-anavy/dna-fountain).

>Composite Hedges Nanopores is published in [Nature Communications](https://www.nature.com/articles/s41467-024-53455-3)
>
>Zhao, X., Li, J., Fan, Q. et al. Composite Hedges Nanopores codec system for rapid and portable DNA data readout with high INDEL-Correction. Nat Commun 15, 9395 (2024).
>
>https://doi.org/10.1038/s41467-024-53455-3
---

## Installation

First, clone the repository to a local directory:

```bash
$ git clone https://github.com/ysfhtxn/Composite-Hedges-Nanopores.git
```

Then install required packages using the file [environment.yml](https://github.com/ysfhtxn/Composite-Hedges-Nanopores/blob/main/environment.yml):

```bash
$ conda env create -f environment.yml && conda activate CHN
```

This will create a Python virtual environment in the `Composite-Hedges-Nanopores` folder. The installation should take less than **10 minutes** on a typical desktop pc, but maybe can take longer if an older pip version is used.

---

## Usage

Composite-Hedges-Nanopores integrated the code of this work into the evaluation platform [Chamaeleo](https://github.com/ntpz870817/Chamaeleo). Based on this evaluation platform, Composite-Hedges-Nanopores uses `./CHN/test_file/sme_introduction.txt` to conduct robustness evaluation.

To evaluate, users can simply run `./robustness_test.py` with command:

```bash
$ python robustness_test.py # The output will be displayed on the terminal
```

or

```bash
$ nohup python -u robustness_test.py > robustness_test.log & # The output will be saved to robustness_test.log
```

The output will be displayed like the following :

```bash
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
**************************************************
Run task (1/1).
**************************************************
Create a transcoding pipeline.
Read binary matrix from file: /data/biolab-nvme-pool1/zhaoxy/github/Composite-Hedges-Nanopores/CHN/test_file/sme introduction.txt
The bit size of the encoded file is None bits and the length of final encoded binary segments is None
Encode bit segments to DNA sequences by coding scheme.
2A0C0G0T 0A2C0G0T 0A0C2G0T 0A0C0G2T 1A1C0G0T 0A0C1G1T 1A0C1G0T 0A1C0G1T 
resolution: 2           sigma: 8
encoding...: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 7/7 [00:00<00:00, 368.16it/s]
mapping...: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 7/7 [00:02<00:00,  2.65it/s]
Decode DNA sequences to bit segments by coding scheme.
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 56/56 [00:22<00:00,  2.54it/s]
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 7/7 [00:00<00:00, 51.60it/s]
decoding...: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 7/7 [00:23<00:00,  3.39s/it]
CHN, None, sme introduction.txt, 0.078, 2.66, 45.913, True, 7, 7, 100.0%
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Evaluation log: 
evaluated coding schemes, evaluated files, evaluated error correction, original segment length, perturbation
[CHN], [sme introduction.txt], [None], 288, {nucleotide insertion: 0.001, nucleotide mutation: 0.03, nucleotide deletion: 0.001, sequence loss: 0, iterations: 1}
task id, coding scheme, error-correction, file, payload length, index length, error-correction length, information density, encoding runtime, decoding runtime, error rate, error indices, error bit segments, transcoding state, success rate
task 0, CHN, None, sme introduction.txt, 288, 0, 0, 0.078, 2.66, 45.913, None, None, None, True, 100.0%
```

The robustness testing should take around **40 seconds** on a typical desktop computer.

## What is needed

### MUSCLE

Please make sure that your environment has `MUSCLE`. This module can generate ensembles of high-accuracy alternative alignments. See [https://www.drive5.com/muscle/](https://www.drive5.com/muscle/) for more details and license restrictions.

After installing `MUSCLE`, make sure that it is placed in your system path.

### MMseqs2

For clustering each replica of specific DNA sequence, `MMseqs2` is needed. See [https://github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2) for more details and license restrictions.

After installing `MMseqs2`, make sure that it is placed in your system path.

---

## Files Tree Diagram

```html
├── assembly                          // The temporary output during alignment and assembly while running robustness_test.py
├── Chamaeleo                         // Chamaeleo - a collection focused on different codec methods for DNA storage
│    ├── ...                          // See more details in "https://github.com/ntpz870817/Chamaeleo"
├── CHN                               // CHN codec file and some modified .py file based on Chamaeleo
│    ├── __init__.py                  // 
│    ├── CHNcodec.py                  // The main body of CHN codec
│    ├── data_handle.py               // Modified from Chamaeleo
│    ├── default.py                   // Modified from Chamaeleo
│    ├── pipelines_mod.py             // Modified from Chamaeleo
│    ├── tools.py                     // Some tools for DNA data recovery 
│    ├── utils.py                     // Some functions for CHN codec
├── CHN_invitro                       // Scripts of in vitro DNA storage data recovery
│    ├── test                         // Files generated by intermediate steps during data processing
│    ├── test_result                  // Files generated by intermediate steps during data processing
│    ├── tmp                          // Files generated by intermediate steps during data processing
│    ├── __init__.py                  // 
│    ├── CHNcodec.py                  // The main body of CHN codec
│    ├── encode.py                    // Encoding example files into DNA strands
│    ├── decode.py                    // Recovering stored data files from raw DNA reads
│    ├── mapping.py                   // Minimap2 scripts --- mapping reads
│    ├── readsdic_gen.py              // Generate read info dictionary based on read ID
│    ├── seq_grouping.py              // Grouping raw DNA reads by barcodes and anchors
│    ├── tools.py                     // Some tools for DNA data recovery 
│    ├── utils.py                     // Some functions for CHN codec
│    ├── examples                     // Examples files
│    │    ├── sme introduction.txt    // An example file
│    │    ├── sme logo.jpg            // An example file
├── .gitmodules                       // Upload submodule file
├── environment.yml                   // Modules required for running test.py
├── README.md                         // Description of this repository
├── robustness_test.py                // Functional test for CHN code based on Chamaeleo
```

---
## Reproducing the ***in vitro*** analysis

To evaluate Composite-Hedges-Nanopores, we encoded the files sme introduction.txt and sme logo.jpg, which can be found in the `./CHN_invitro/examples/` folder. The encoded files are then synthesized, amplified, and sequenced. For encoding, run the command is as follows:

```bash
$ cd CHN_invitro && python encode.py
```

In order to reproduce the decoding of sequencing data, it must first be downloaded from the sequence read archive. The sequencing .fastq data has been deposited in the CNSA (https://db.cngb.org/cnsa/) of the CNGBdb with accession CNP0005551. 

```bash
$ cd CHN_invitro && python mapping.py && python readsdic_gen.py && python seq_grouping.py && python decode.py
```

Remember that the paths in the file above must be adjusted before decoding.

---
## Citation
If you think this repository helps or inspires your study, please consider citing it. 😄
```bibtex
@article{Zhao2024,
  author = {Xuyang Zhao, Junyao Li and Qingyuan Fan and Jing Dai and Yanping Long and Ronghui Liu and Jixian Zhai and Qing Pan and Yi Li},
  title = {Composite Hedges Nanopores codec system for rapid and portable DNA data readout with high INDEL-Correction},
  journal = {Nature Communications},
  volume = {15},
  number = {1},
  pages = {9395},
  year = {2024},
  doi = {10.1038/s41467-024-53455-3},
  url = {https://doi.org/10.1038/s41467-024-53455-3},
  abstract = {Reading digital information from highly dense but lightweight DNA medium nowadays relies on time-consuming next-generation sequencing. Nanopore sequencing holds the promise to overcome the efficiency problem, but high indel error rates lead to the requirement of large amount of high quality data for accurate readout. Here we introduce Composite Hedges Nanopores, capable of handling indel rates up to 15.9\% and substitution rates up to 7.8\%. The overall information density can be doubled from 0.59 to 1.17 by utilizing a degenerated eight-letter alphabet. We demonstrate that sequencing times of 20 and 120 minutes are sufficient for processing representative text and image files, respectively. Moreover, to achieve complete data recovery, it is estimated that text and image data require 4× and 8× physical redundancy of composite strands, respectively. Our codec system excels on both molecular design and equalized dictionary usage, laying a solid foundation approaching to real-time DNA data retrieval and encoding.},
  issn = {2041-1723}
}

