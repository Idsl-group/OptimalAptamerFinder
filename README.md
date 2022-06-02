# OptimalAptamerFinder

## Installation
This program is designed to run on a Linux environment. (RNAFold is complicated to install on Windows and haven't tested on Mac.)
1. Install all the packages in the requirements.txt file (all are available via pip)
2. Install RNAFold as indicated in [this guide](https://algosb2019.sciencesconf.org/data/RNAtutorial.pdf)

## Instructions
 Consider we have sequence rounds 06, 10, 12, 15, 16, 18, 20 and 22 form experiments with thophylline as the target molecule.
 
1. Add data to `data/` folder.
   1. Create the directory `data/theophylline_fastq_r2/` and place the fastq files for each rounds in it.
   2. Create the file `data/theo_primers.txt` and write each primer in a separate line.
   3. Create the file `data/theolist_cut_22.txt` and write the name (without .fastq) of each fastq file in a separate line.
```
OptimalAptamerFinder
└───data
    ├───theophylline_fastq_r2
        ├───theo06.fastq
        ├───theo10.fastq
        ...
    └───theo_primers.txt
            "GGTCTCAGCAGGAGTCCTCCT\nGGTCTCAGCAGGAGTCCTCCT"
    └───theolist_cut_22.txt
            "theo06\ntheo10\ntheo12\ntheo15\ntheo16\ntheo18\ntheo20\ntheo22"
└───results
    ├───experimentID
        ├───recommendations
            ├───rec_theophylline_06.txt
            ├───rec_theophylline_06_hist.png
            ...
        └───alignment_records.json
        └───sequencelib.pickle
└───src
    ├───fast12txt.py
    ├───fastqcleaner.py
    ├───Sequence.py
└───README.md
└───requirements.txt
│
```
2. Convert the fastq .files to .txt files using the following command:
```
python3 src/fastq2txt.py theophylline r2
```
3. Run the program and see the output in the `results/` folder by running the following command:
```
python3 src/Sequence.py theophylline
```
To run the program with several binding_targets, separate them by a hyphen.
```
python3 src/Sequence.py theophylline-uricacid
```
