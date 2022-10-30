# OptimalAptamerFinder

## Installation
This program is designed to run on a Linux environment. (RNAFold is complicated to install on Windows and haven't tested on Mac.)
1. Install all the packages in the requirements.txt file (all are available via pip).

    1.1 Create virtual environment (venv) in the root directory of the project as shown in [this guide](https://docs.python.org/3/library/venv.html)

   1.2 Execute the following commands in the terminal from the root directory of the project:
```
source venv/bin/activate
pip install -r requirements.txt
deactivate
```
2. Install RNAFold as indicated in [this guide](https://algosb2019.sciencesconf.org/data/RNAtutorial.pdf)
3. Check that the command `RNAfold --help` works in the terminal. If not, run the following commands before execution:
```
export PATH=${HOME}/Tutorial/Progs:${PATH}
source ~/.bashrc
```

## Execution instructions
 Consider we have sequence rounds 06, 10, 12, 15, 16, 18, 20 and 22 form experiments with thophylline as the target molecule.
 
1. Add data to `data/` folder.
   1. Create the directory `data/theophylline_fastq_r2/` and place the fastq files for each rounds in it.
   2. Create the file `data/theo_primers.txt` and write each primer in a separate line.
   3. Create the file `data/theolist.txt` and write the name (without .fastq) of each fastq file in a separate line.
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


## Project architecture
```
OptimalAptamerFinder
└───data
    ├───theophylline_fastq_r2
        ├───theo06.fastq
        ├───theo10.fastq
        ...
    └───theo_primers.txt
            "GGTCTCAGCAGGAGTCCTCCT\nGGTCTCAGCAGGAGTCCTCCT"
    └───theolist.txt
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


To visualize the data we need to look at the 2 main classes of this project: SequenceLibrary and Aptamer

## fastq2txt.py
This file is only used at the beginning of our experiments to convert the fastq files to txt files that only contain
the section of the DNA chain we are interested in studying. In our case this corresponded to the nucleotides in between
positions 24 and 54. Users might need to change the rel_range variable to match their own experiments.

## Sequence.py
This file is the main file of the project. It contains the Aptamer and SequenceLibrary classes along with the main fucntions
used to control them.

### Aptamer class
The ```Aptamer``` class is used to represent a single aptamer. Each instance has the following attributes:
```
Aptamer() 
 ─ sequence = nucleotide sequence of the aptamer with the given seqid
 ─ score = final score given to the aptamer
 ─ seqid = seqid of the aptamer
 ─ rounds = {
    'binding_target_1': {
                         'R_1': counts in R_1, 
                         'R_2': counts in R_2 
                           ...
                         },
    'binding_target_2': {},
    ...
    }
 ─ primers = [forward_primer, backward_primer]  (Empty strings if library is done wihtout primers)
 ─ loops = list of loops found in the sequence
 ─ kmers = list of kmers found in the sequence
 ─ forward_stems = list of forward stems found in the sequence
 ─ backward_stems = list of backward stems found in the sequence
 ─ ct = ct structure of the aptamer as given by DNAfold
 ─ dg = dg energy as given by DNAfold
```

The ```Aptamer``` class has the following methods: 
- ```get_sequence()```: return the aptamer sequence.
- ```get_ct()```: returns the ct structure of the aptamer.
- ```get_kmers()```: returns the list of kmers found in the aptamer.
- ```get_hairpins()```: returns the forward stems, backward stems and loops found in the aptamer.
- ```add_round(binding_target, rnd, rnd_counts)```: adds the counts of the aptamer in the given round to the rounds 
attribute.
- ```add_ct(ct_dict)```: sets the ct structure and dg energy of the aptamer from the dictionary containing the output 
given by DNAfold.
- ```add_kmers(k)```: sets the kmers attribute with the kmers of length k found in the aptamer.
- ```add_hairpins()```:  uses the ct structure and sequence of the aptamer to find and set the attributes for the 
forward stems, backward stems and loops found in the aptamer.

### SequenceLibrary class
Each ```SequenceLibrary``` instance has four attributes:
```
SequenceLibrary() 
 ─ k = the length of the k-mer
 ─ with_primers = a boolean indicating whether the sequences in the library keep the primers
 ─ info = a dictionary containing information about the library
 ─ sequences = a dictionary of **Aptamer** objects
```

####  info dictionary
The info dictionary has a nested dictionary structure with the following keys and values:
    
```
info = {
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
            }
        ├───R_2: {...}
         ...
        └───R_p: {...}
        }
     ...
    ├───binding_target_molecule_n-1: {...}
    └───binding_target_molecule_n: {...}
    }
# Structures not described {...} are assumed to be similar to the previous ones on the same hierarchical level
```
  
#### sequences dictionary
The sequences dictionary is much simpler. Its keys correspond to sequence ids from the info dictionary and the value
of each item is an **Aptamer** object.
```
sequences = {
    ├───seqid1: Aptamer()
    ├───seqid2: Aptamer()
    ...
    └───seqidn: Aptamer()
    }
```

The ```SequenceLibrary``` class has the following methods: 
- ```get_count_single_run(binding_target, rnd)```: gets the counts of each unique aptamer that appears in round rnd of 
the given binding target experiment (extracts data directly from data directory)
- ```get_primers(binding_target)```: gets and returnss the primers used in the given binding target experiment 
(extracts data directly from data directory)
- ```get_kmer_counts(binding_target, rnd)```: gets the counts of each kmer in the given round of the given binding 
target experiment (Uses all the seqids stored in the info dictionary and uses them to iterate over the sequences 
dictionary)
- ```get_structures(binding_target, rnd)```: gets the ct structrues and dg energies of each aptamer in the given round
og the given biding target experiment. First it creates a text file with all the unique aptamer sequences of the round.
This file is then used as input for DNAfold. The output of DNAfold is then used to set the ct and dg attributes of
each aptamer in the sequences dictionary.
- ```get_structure_counts(binding_target, rnd)```: gets the counts of each loop, forward stem and backward stem in the
given round of the given binding target experiment. 
- ```add_data(binding_target, rnd)```: manages and adds data to the class instance. It calls the 
```get_count_single_run```, ```get_kmer_counts``` ```get_structures```, and ```get_structure_counts``` functions.
- ```filter_dict(unfiltered_dict, percentile)```: removes all element of a dictionary of counts that have values below 
the percentile threshold. Used to remove aptamers that do not contribute much to the averages and speed up the program.
- ```priority_clustering(binding_target, rnd, alignment_threshold, max_num_alignments)```: clusters aptamers based on 
their aligned priority. It performs the following iterative procedure: 
   1. It selects the most common aptamer in the remaining list of aptamers to be ordered. Creates a new cluster and sets 
     this aptamer as the target for the cluster.
   2. It aligns the rest of the remaining aptameters with respect to the target aptamer using local pairwise Striped 
      Smith-Waterman alignment.
   3. All sequences with a score above the alignment threshold are added to the cluster and removed from the list of 
      remaining aptamers.
   4. If the number of identified clusters is less than max_num_alignments steps 1-3 are repeated.
- ```find_clusters_correlation(cluster_motifs, k)```: returns a matrix with the correlation between the k most common 
motifs in each cluster. This correlation score is calculated based on the number of motifs in common between each 
each pair of clusters.
- ```analyze_priority_clusters(alignment_records, binding_target, rnd)```: uses the ```find_clusters_correlation``` with
kmers, stems and loops to combine those priority clusters that have a high correlation score.
- ```print_cluster(sequences, binding_target, rnd, cnum, total_counts, tcs, max_counts)```: reads, adds and displays 
cluster information and computes final scores for each sequence.

### main() function
Takes a list of binding targets as input. Creates the ```SequenceLibrary``` instance and calls the ```add_data``` 
function. Then it calculates the priority clusters by calling the ```get_priority_clusters``` function. Lastly, it calls 
the ```generate_all_recommendations``` which orders the aptamers based on their final_scores and displays the ones with 
the highest scores.
