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

To visualize the data we need to look at the 2 main classes of this project: SequenceLibrary and Aptamer

## SequenceLibrary Class
Contains four attributes:
- k: the length of the k-mer
- with_primers: a boolean indicating whether the sequences in the library keep the primers
- info: a dictionary containing information about the library
- sequences: a dictionary of **Aptamer** objects

   ####  info dictionary
   The info dictionary has a nested dictionary structure with the following keys and values:
    ```
    └───info: {
        ├───binding_target_1:{
            ├───R_1: {
                ├─── seqids: list of ids corresponding to unique sequences in the library
                ├─── motif_counts: {
                    ├─── kmer: {
                        ├─── unweighted: {
                            ├─── kmer1: count1
                            ├─── kmer2: count2
                            ...
                        }
                        └─── weighted: {...}
                        ...
                        }
                    }   
                    ├─── loops: number {...}
                    ├─── forward_stems: {...}
                    └─── backward_stems: {...}
                    }
                }
                └─── motif_hvalues: vestigial attribute (Need to make sure that it can be safely removed)
                }
            ├───R_2: {...}
             ...
            └───R_p: {...}
            }
         ...
        ├───binding_target_molecule_n-1: {...}
        └───binding_target_molecule_n: {...}
        }
    ```
  
    #### sequences dictionary
    The sequences dictionary is much simpler. Its keys correspond to sequence ids from the info dictionary and the value
    of each item is an **Aptamer** object. Ech **Aptamer** object has the following attributes:
    ```
    Aptamer() <- seqid
    ├───sequence = nucleotide sequence of the aptamer with the given seqid
    ├───score = final score given to the aptamer
    ├───rounds = {
        'binding_target_1': {
                             'R_1': counts in R_1, 
                             'R_2': counts in R_2 
                               ...
                             },
        'binding_target_2': {},
        ...
        }
    ├───primers = [forward_primer, backward_primer]  (Empty strings if library is done wihtout primers)
    ├───loops = list of loops found in the sequence
    ├───kmers = list of kmers found in the sequence
    ├───forward_stems = list of forward stems found in the sequence
    ├───backward_stems = list of backward stems found in the sequence
    ├───ct = ct structure of the aptamer
    └───dg = 0
    ```
    


dictonary = {'key1':value1,
             'key2':value2,
             ...
             }

class Human():
    def __init__(self, x, y):
        self.name = x
        self.age = y

    def say_hello(self):
        print('Hello, my name is {} and I am {} years old.'.format(self.name, self.age))

human1 = Human('John', 25)
human1.say_hello()