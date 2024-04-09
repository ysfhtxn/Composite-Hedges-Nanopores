# Composite Hedges Nanopores

It is a demo README.md which still needs to be improved.

A package of coding algorithm for rapid readout of digital information storage in DNA, inspired by [HEDGES](https://github.com/whpress/hedges) and [composite DNA letter](https://github.com/leon-anavy/dna-fountain).

---

## Installation

First, clone the repository to a local directory:

```bash
git clone https://github.com/ysfhtxn/Composite-Hedges-Nanopores.git
```

Then install required packages using the file [requirements.yaml](https://github.com/ysfhtxn/Composite-Hedges-Nanopores/blob/main/requirements.yaml):

---

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
├── assembly                          // The temporary output during alignment and assembly
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
│    ├── data_recovery.py             // Recovering stored data files from raw DNA reads
│    ├── mapping.py                   // 
│    ├── paf_to_tsv.py                // 
│    ├── readsdic_gen.py              // 
│    ├── seq_grouping.py              // Grouping raw DNA reads by barcodes and anchors
│    ├── tools.py                     // Some tools for DNA data recovery 
│    ├── utils.py                     // Some functions for CHN codec
├── .gitmodules                       // Upload submodule file
├── environment.yml                   // Modules required for running test.py
├── README.md                         // Description of this repository
├── test.py                           // Functional test for CHN code based on Chamaeleo
```

---

If you think this repository helps or inspires your study, please consider refer it. 😄
