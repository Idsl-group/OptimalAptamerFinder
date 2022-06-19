import sys
from matplotlib import pyplot as plt
from scipy.stats import norm
sys.path.append("./")
import os, multiprocessing, pickle, random, subprocess
from difflib import SequenceMatcher
import numpy as np
import pandas as pd
import re
from pyfastaq import sequences
import pickle
import time
import json
from itertools import chain
from skbio import TabularMSA, DNA
from skbio.alignment import local_pairwise_align_ssw
import hashlib
from tqdm import tqdm
from termcolor import colored
#from src.nc2hp import
import datetime


# Initialize time that will be used to label this experiment
TIME = str(datetime.datetime.now())
# Create directory for results
if not os.path.exists('results'):
    os.mkdir('results')
if not os.path.exists('results/'+TIME):
    os.mkdir('results/'+TIME)


# Aptamer class to store information about each aptamer
class Aptamer():
    def __init__(self, seqid, sequence, primers):
        """
        Initializes an aptamer object.
        :param seqid: ID of the sequence as given by the input file
        :param sequence: Nucleotide sequence of the aptamer (30-nucleotide region)
        :param primers: Primers used for the aptamer ['ACGTCGT', 'ACGTCGT'] or ['', '']
        """
        self.sequence = sequence
        self.primers = primers
        self.seqid = seqid
        # Initilaize empty attributes
        self.rounds = {}  # Will store counts in each round of each target molecule {binding_target: {round: count}}
        self.ct = ""  # ct structure of the aptamer (if available) as given by DNAFold
        self.dg = 0  # dg of the aptamer (if available) as given by DNAFold
        self.kmers = None # kmers of the aptamer
        self.loops = None # loops structures found in the aptamer (if any)
        self.forward_stems = None # stems to the left of the loop
        self.backward_stems = None # stems to the right of the loop
        self.score = 0 # store final score given to the aptamer
        #self.alignments = {} # Not used ?
        #self.clusters = {} # Not used ?

    def get_sequence(self):
        # Return the sequence of the aptamer
        return self.sequence

    def get_ct(self):
        # Return the ct structure of the aptamer
        return self.ct

    def get_kmers(self):
        # Return the kmers of the aptamer
        return self.kmers

    def get_hairpins(self):
        # Return the hairpins of the aptamer
        return self.forward_stems, self.backward_stems, self.loops

    #def newAttr(self, attr_name, attr_value):
        # Create new class attribute
    #    setattr(self, attr_name, attr_value)

    def add_round(self, binding_target, rnd, rnd_counts):
        """
        Add a new round to the aptamer.
        :param binding_target: binding target molecule of the aptamer
        :param rnd: round of the aptamer
        :param rnd_counts: number of counts of the aptamer in the round for the given binding target
        :return:
        """
        if binding_target not in self.rounds.keys():
            self.rounds[binding_target] = {rnd: rnd_counts}
        else:
            self.rounds[binding_target][rnd] = rnd_counts

    def add_ct(self, ct_dict):
        """
        Add ct structure and dg to the aptamer form dictionary with results from RNAFold
        :param ct_dict: {'ct': ct_structure, 'dg': dg}
        :return:
        """
        try:
            self.ct = ct_dict[self.sequence]['ct']
            self.dg = ct_dict[self.sequence]['dg']
        except KeyError:
            self.ct = ""
            self.dg = 0
            print("KeyError:", self.sequence)

    def add_kmers(self, k=6):
        '''
        Creates a list of kmers of length k in the aptamer
        :param k: length of the kmers
        '''
        kmers = []
        sequence = self.sequence.replace(self.primers[0], '').replace(self.primers[1], '')
        for i in range(0, len(sequence) - k + 1):
            kmers.append(sequence[i:i + k])
        self.kmers = kmers

    def add_hairpins(self):
        '''
        Creates a list of hairpins in the aptamer
        The hairpin structure is store in the forward_stems, backward_stems and loops attributes of the class.
        '''
        ct = self.get_ct()   #'....(((...))....'
        aptamer = self.get_sequence()
        state = 'off'       # Begin in the off state (no structure is being recorded) one of ['off', 'loop', ')stem', '(stem']
        current_loop = ""
        current_forward_stem = ""
        current_backward_stem = ""
        forward_stems = []
        backward_stems = []
        loops = []
        for i in range(0, len(ct)):
            v = ct[i]       # v is the current character in the ct structure one of '.'   ;   '('   ;   ')'
            nc = aptamer[i]         # nc is the current nucleotide in the aptamer
            if v == '.':
                if state == "off":      # No structure is being recorded
                    state == "off"
                if state == 'loop':     # We are recording a loop
                    current_loop += nc
                    state == "loop"     # Stay in the loop state
                elif state == ')stem':  # We are recording a backward stem
                    if len(current_forward_stem) == len(current_backward_stem): # If the length of the forward and backward stems are the same,
                                                                                #  we are done recording the backward stem and can add them to the lists
                        forward_stems.append(current_forward_stem)
                        backward_stems.append(current_backward_stem)
                        current_forward_stem = ""       # Reset the current forward stem
                        current_backward_stem = ""      # Reset the current backward stem
                    state = "off"        # Closed structure, go back to the off state
                elif state == "(stem":   # We are recording a forward stem
                    current_loop += nc
                    state = "loop"       # Begin recording a loop
            elif v == '(':
                if state == "off":      # No structure is being recorded
                    current_forward_stem += nc
                if state == 'loop':     # We are recording a loop
                    current_forward_stem += nc
                    current_loop = ""    # Reset the current loop
                elif state == ')stem':  # We are recording a backward stem
                    if len(current_forward_stem) == len(current_backward_stem):# If the length of the forward and backward stems are the same,
                                                                               #  we are done recording the backward stem and can add them to the lists
                        forward_stems.append(current_forward_stem)
                        backward_stems.append(current_backward_stem)
                        current_forward_stem = ""       # Reset the current forward stem
                        current_backward_stem = ""      # Reset the current backward stem
                    current_forward_stem += nc
                elif state == "(stem":       # We are recording a forward stem
                    current_forward_stem += nc
                state = '(stem'     # Forward stem state
            elif v == ')':
                if state == "off":      # No structure is being recorded
                    current_backward_stem += nc
                if state == 'loop':     # We are recording a loop
                    loops.append(current_loop)     # Add the current loop to the list of loops
                    current_loop = ""    # Reset the current loop
                    current_backward_stem += nc
                elif state == ')stem':      # We are recording a backward stem
                    current_backward_stem += nc
                elif state == "(stem":       # We are recording a forward stem
                    current_forward_stem += nc
                    loops.append("")              # Add an empty loop to the list of loops
                state = ")stem"     # Backward stem state
        # Add the last structures to the lists
        if len(current_forward_stem) == len(current_backward_stem) and len(current_forward_stem) != 0:
            forward_stems.append(current_forward_stem)
            backward_stems.append(current_backward_stem)
        self.forward_stems = forward_stems
        self.backward_stems = backward_stems
        self.loops = loops


# Sequence Library class to store experimental sequences and results
class SequenceLibrary():
    def __init__(self, k=6, with_primers=False):
        """
        Initialize the sequence library
        :param k: kmer size
        :param with_primers: if True, we include primers in structural analysis of each aptamer
        """
        # Initialize the sequence library
        self.k = k
        self.with_primers = with_primers
        # Empty dictionaries to store the sequences and relevant information
        self.sequences = {}
        self.info = {}

    def get_count_single_run(self, binding_target, rnd):
        """
        Get the count of each unique aptamer in a single run of the given binding target
        :param binding_target:
        :param rnd:
        :return: # Returns a dataframe in which each row corresponds to a unique aptamer sequence
        # The columns of the dataframe correspond to:
        # 'seqs' aptamer sequence
        # 'count' number of occurrences of the aptamer sequence
        # 'freq' normalized frequency of the aptamer sequence
        # 'cum_freq' cumulative frequency of the aptamer sequence once the dataframe is ordered in decreasing frequency order
        """
        # Set some parameters to search the directories for the appropriate file
        path = str(os.getcwd()) + f"/data/{binding_target}_fastq_r2" # Path to folder with sequenced libraries of each round
        # Indices of the relevant 30-mer region of the aptamer
        if "r1" in path:
            rel_range = [24, 54]
        elif "r2" in path:
            rel_range = [18, 48]
        #file_prefix = path.split("/data/")[1][:4]
        # Iterate over all aptamer shown in the sequenced library and add them to the list df
        df = []
        for file in os.listdir(path):
            if file.replace(path.split("/data/")[1][:4], "").replace(".fastq", "") in [rnd]:
                # print("Using ", file)
                seq_reader = sequences.file_reader(path + "/" + file)
                for sequence in seq_reader:
                    seq = sequence.seq[rel_range[0]: rel_range[1]]
                    df.append(seq)
            else:
                pass
        # Convert df list into a pandas dataframe
        df = pd.DataFrame(df, columns=['seq'])
        # Create the counts_df which contains the number oc occurrences of each unique aptamer
        counts_df = pd.DataFrame(df['seq'].value_counts()).rename(columns={'seq': 'counts'})
        # Extract som other values from the dataframe
        counts_df['seqs'] = counts_df.index.to_numpy()
        counts_df['freq'] = counts_df['counts'] / counts_df['counts'].sum()
        counts_df = counts_df.reset_index().sort_values(by=['freq'], ascending=False)
        counts_df.drop(columns=['index'], inplace=True)
        freq = counts_df['freq'].tolist()
        cum_freq = []
        for i in range(0, len(freq)):
            cum_freq.append(sum(freq[:i + 1]))
        counts_df['cum_freq'] = cum_freq
        return counts_df

    def get_primers(self, binding_target):
        """
        Get the primers for the given binding target
        :param binding_target:
        :return: list with forward primer and reverse primer
        """
        primer_file = str(os.getcwd()) + f"/data/{binding_target[:4]}_primers.txt"
        with open(primer_file, "r") as f:
            primers = [line.strip() for line in f.readlines()]
        f.close()
        return primers

    def get_kmer_counts(self, binding_target, rnd):
        """
        Get the counts of each kmer in round rnd for the given binding target
        :param binding_target:
        :param rnd:
        :return:
        """
        unweighted = {} # Counts in how many unique aptamers each kmer appears (adds 1 to the count if the aptamer contains the kmer)
        weighted = {} # Counts how many times in total each kmer appears (adds the counts of the aptamer to the count if the aptamer contains the kmer)
        print("Getting kmer counts")
        # Generate weighted (by aptamer counts) and unweighted kmer counts
        for seqid in tqdm(self.info[binding_target][rnd]['seqids']): # Iterate over the sequence ids corresponding to the given bindig_target and the given round
            # Remember that each unique seqid corresponds to a unique aptamer
            kmers = self.sequences[seqid].get_kmers()  # Get list of kmers in the aptamer
            counts = self.sequences[seqid].rounds[binding_target][rnd] # Get counts of the aptamer in the binding_target in the round
            for kmer in kmers:
                if kmer in unweighted: # Add one to the kmer count
                    unweighted[kmer] += 1
                else: # Initialize the kmer count at 1
                    unweighted[kmer] = 1
                if kmer in weighted: # Add counts of the aptamer to the kmer count
                    weighted[kmer] += counts
                else: # Initialize the kmer counts as the counts of the aptamer
                    weighted[kmer] = counts
        # Add weighted and unweighted counts to the info dictionary
        if 'kmer' not in self.info[binding_target][rnd]['motif_counts'].keys():
            self.info[binding_target][rnd]['motif_counts']['kmer'] = {}
        if rnd not in self.info[binding_target][rnd]['motif_counts']['kmer'].keys():
            self.info[binding_target][rnd]['motif_counts']['kmer'] = {'weighted': self.filter_dict(weighted),
                                                                      'unweighted': self.filter_dict(unweighted)}

    def get_structures(self, binding_target, rnd, structure_mode='from_ct'):
        """
        Get the structures of the aptamers for the given binding target in the given round
        :param binding_target:
        :param rnd:
        :param structure_mode: At some point during development we wanted to have different ways to extract the structure
                                eventually we only calculated the structure from the ct file so this isn't necessary anymore.
                                To be removed.
        :return:
        """
        # Create a file (temp_seqs.txt) with all the unique aptamer sequences (one per line)
        seq_file = str(os.getcwd()) + "/temp_seqs.txt"
        with open(seq_file, 'w') as f:
            for seqid in self.info[binding_target][rnd]['seqids']:
                # f.write(">" + seqid + "\n")
                f.write(self.sequences[seqid].sequence + "\n")
        f.close()
        # Use RNAfold to get the structures and mfes of all the aptamers and write them to (temp_ct.txt)
        ct_file = str(os.getcwd()) + "/temp_ct.txt"
        cmd = f"RNAfold --jobs={multiprocessing.cpu_count()} --infile={seq_file} --outfile={ct_file.split('/')[-1]} --noPS -T {37.0} --noconv"
        subprocess.call([cmd], shell=True)
        """
        The structure of temp_ct.txt is the following:
        AGCTGGAATTTCGTTTACGTTACGGGAGTT
        ........((((((.......))))))... (-5.2)
        and so on for each aptamer
        """
        # Create a dictionary from the information in temp_ct.txt
        ct_dict = {}
        with open(ct_file, 'r') as f:
            for i, line in enumerate(f):
                if i % 2 == 0:
                    # If line is even, the line contains the sequence
                    sequence = line.strip()
                else: # If line is odd, the line contains the ct structure and mfe of the aptamer i the previous line
                    ct = line.strip().split(" ")[0]
                    dg = float(line.strip().split(" (")[-1].replace("(", "").replace(")", ""))
                    ct_dict[sequence] = {'ct': ct, 'dg': dg}
        f.close()
        # Add ct structure, mfe and calculate the harpin motifs of all sequences
        for seqid in self.info[binding_target][rnd]['seqids']:
            self.sequences[seqid].add_ct(ct_dict)
            self.sequences[seqid].add_hairpins()

    def get_structure_counts(self, binding_target, rnd):
        """
        Get the counts of each structure in round rnd for the given binding target
        :param binding_target:
        :param rnd:
        :return:
        """
        # Wont't explain since the idea here is the same as in self.get_kmer_counts() but with different hairpins, stems and loops instead of kmers
        hairpins_weighted = {}
        hairpins_unweighted = {}
        stems_weighted = {}
        stems_unweighted = {}
        loops_weighted = {}
        loops_unweighted = {}
        print("Getting structure counts")
        for seqid in tqdm(self.info[binding_target][rnd]['seqids']):
            counts = self.sequences[seqid].rounds[binding_target][rnd]
            hairpins = self.sequences[seqid].hairpins
            stems = self.sequences[seqid].stems
            loops = self.sequences[seqid].loops
            for hairpin in hairpins:
                if hairpin in hairpins_unweighted:
                    hairpins_unweighted[hairpin] += 1
                else:
                    hairpins_unweighted[hairpin] = 1
                if hairpin in hairpins_weighted:
                    hairpins_weighted[hairpin] += counts
                else:
                    hairpins_weighted[hairpin] = counts
            for stem in stems:
                if stem in stems_unweighted:
                    stems_unweighted[stem] += 1
                else:
                    stems_unweighted[stem] = 1
                if stem in stems_weighted:
                    stems_weighted[stem] += counts
                else:
                    stems_weighted[stem] = counts
            for loop in loops:
                if loop in loops_unweighted:
                    loops_unweighted[loop] += 1
                else:
                    loops_unweighted[loop] = 1
                if loop in loops_weighted:
                    loops_weighted[loop] += counts
                else:
                    loops_weighted[loop] = counts
        if 'hairpin' not in self.info[binding_target][rnd]['motif_counts'].keys():
            self.info[binding_target][rnd]['motif_counts']['hairpin'] = {
                'weighted': self.filter_dict(hairpins_weighted), 'unweighted': self.filter_dict(hairpins_unweighted)}
        if 'stem' not in self.info[binding_target][rnd]['motif_counts'].keys():
            self.info[binding_target][rnd]['motif_counts']['stem'] = {'weighted': self.filter_dict(stems_weighted),
                                                                      'unweighted': self.filter_dict(stems_unweighted)}
        if 'loop' not in self.info[binding_target][rnd]['motif_counts'].keys():
            self.info[binding_target][rnd]['motif_counts']['loop'] = {'weighted': self.filter_dict(loops_weighted),
                                                                      'unweighted': self.filter_dict(loops_unweighted)}

    def add_data(self, binding_target, rnd, structure_mode='from_nc', top_k=None):
        """
        Add the data to the info dictionary
        :param binding_target:
        :param rnd:
        :param structure_mode:
        :param top_k: used if we only want to use the top k most common sequences of the library. Never really used it
        :return:
        """
        # Ge primers
        if self.with_primers:
            primers = self.get_primers(binding_target)
        else:
            primers = ["", ""]
        # Get sequence counts
        if top_k is None:
            df = self.get_count_single_run(binding_target, rnd)
        else:
            df = self.get_count_single_run(binding_target, rnd).head(top_k)
        # Get squences and counts as lists. Also initialize an empty list to store the seqids
        sequences = df['seqs'].tolist()
        counts = df['counts'].tolist()
        seqids = []
        # self.info = {} is a dictionary that contains several other dictionaries.
        # On the first level there are the different binding_targets of the study (usually we only compute the library for one binding_target at a
        # time but this allows you to store information of multiple binding_targets at the same time)
        if binding_target not in self.info.keys():
            self.info[binding_target] = {rnd: {'seqids': {}, # Will store seqids and aptamers for given rnd
                                               'motif_counts': {}, # Will store different motif counts for given rnd
                                               }}
        if rnd not in self.info[binding_target].keys():
            self.info[binding_target][rnd] = {'seqids': {},
                                              'motif_counts': {},
                                              }
        # Add the data to the info dictionary
        print(f"Adding data for {binding_target} from round {rnd} to library")
        for i in tqdm(range(0, len(sequences))):
            sequence = primers[0] + sequences[i] + primers[1] # Add primers to the sequence
            seqid = hashlib.sha1(str(sequence).encode('utf8')).hexdigest()[:10] # Generate seqid with a hash
            seqids.append(seqid) # Add seqid to the list of seqids
            if seqid in self.sequences.keys(): # If seqid is already in the library (this can happen if the aptamer appeared
                                                # in an earlier round or previously studied binding_target), add the round counts
                                                # to the Aptamer object
                self.sequences[seqid].add_round(binding_target, rnd, counts[i])
            else:  # If seqid is not in the library, create a new Aptamer object
                aptamer = Aptamer(seqid, sequence, primers)
                aptamer.add_round(binding_target, rnd, counts[i]) # Add the round counts to the Aptamer object
                aptamer.add_kmers(self.k) # Add kmers to the aptamer
                self.sequences[aptamer.seqid] = aptamer # Add aptamer and seqid to the self.sequences dictionary
        self.info[binding_target][rnd]['seqids'] = seqids # Add all seqids (of all aptamers that appear for the given binding_target and rnd to the dictioanry
        # Initialize empty motif counts dicitonary to later hold counts of kmers and
        if 'motif_counts' not in self.info[binding_target][rnd].keys():
            self.info[binding_target][rnd]['motif_counts'] = {}
        # motif_hvalues was initially going to hold the aplificantion values of different motifs but we ended up not using it
        # Ignore it
        if 'motif_hvalues' not in self.info[binding_target][rnd].keys():
            self.info[binding_target][rnd]['motif_hvalues'] = {}
        # Get kmer counts and structures of all aptamers for the given binding_target in the given round
        self.get_kmer_counts(binding_target, rnd)
        self.get_structures(binding_target, rnd, structure_mode=structure_mode)

    def get_motif_hvalues(self, binding_target):
        """
        Calculate amplification values for all motifs for the given binding_target
        This function is never used but I kept it here because it could have been useful
        I encourage you to signore it
        """
        rnds = list(self.info[binding_target].keys())
        for motiftype in self.info[binding_target][rnds[0]]['motif_counts'].keys():
            for weighted in [True, False]:
                if weighted:
                    counts = {rnd: self.info[binding_target][rnd]['motif_counts'][motiftype]['weighted'] for rnd in
                              rnds}
                else:
                    counts = {rnd: self.info[binding_target][rnd]['motif_counts'][motiftype]['unweighted'] for rnd in
                              rnds}
                rnds = list(counts.keys())
                motifs = list(set(chain(*[counts[rnd].keys() for rnd in rnds])))
                motif_values = {}
                for i in range(1, len(rnds)):
                    values = {}
                    prior = list({k: v for k, v in
                                  sorted(counts[rnds[i - 1]].items(), key=lambda item: item[1], reverse=True)}.keys())
                    posterior = list({k: v for k, v in
                                      sorted(counts[rnds[i]].items(), key=lambda item: item[1], reverse=True)}.keys())
                    # Sort prior
                    for motif in tqdm(motifs):
                        try:
                            prior_index = prior.index(motif)
                        except ValueError:
                            prior_index = len(prior)
                        try:
                            posterior_index = posterior.index(motif)
                        except Exception:
                            posterior_index = len(posterior)
                        values[motif] = posterior_index - prior_index
                    motif_values[rnds[i]] = values
                for rnd in rnds[1:]:
                    if motiftype not in self.info[binding_target][rnd]['motif_hvalues'].keys():
                        self.info[binding_target][rnd]['motif_hvalues'][motiftype] = {}
                    if weighted:
                        self.info[binding_target][rnd]['motif_hvalues'][motiftype]['weighted'] = motif_values[rnd]
                    else:
                        self.info[binding_target][rnd]['motif_hvalues'][motiftype]['unweighted'] = motif_values[rnd]

    """
    def calculate_sequence_scores(self, binding_target, rnd, scoretype='counts', motiftype='kmer', weighted=True, scaletype='sum'):
        data = self.info[binding_target][rnd][scoretype][motiftype]['weighted' if weighted else 'unweighted']
        idx = list(data.keys())
        values = list(data.values())
        col = [rnd]
        df = pd.DataFrame(values, index=idx, columns=col)
        if scaletype == 'sum':
            df = df / df.sum(axis=0)
        elif scaletype == 'max':
            df = df / abs(df).max(axis=0)
        elif scaletype == 'standardize':
            df = (df - df.mean(axis=1)) / df.std(axis=1)
        for seqid in self.info[binding_target][rnd]['seqids']:
            self.sequences[seqid].add_score(df, binding_target, rnd, scoretype=scoretype, motiftype=motiftype,
                                            weighted=weighted)
    """

    def filter_dict(self, unfiltered_dict, percentile=0.5):
        """
        Removes all elements from dict that don't have a value above a threshold
        """
        thres = list(set(unfiltered_dict.values()))
        thres = thres[int(len(thres) * (1 - percentile))]
        filtered_dict = {k: v for k, v in unfiltered_dict.items() if v >= thres}
        filtered_dict = {k: v for k, v in sorted(filtered_dict.items(), key=lambda item: item[1], reverse=True)}
        return filtered_dict

    def priority_clustering(self, binding_target, rnd, alignment_threshold=0.8, max_num_alignments=20, primers=["", ""]):
        """
        Cluster sequences based on the aligned priority
        """
        print("Using alignment_threshold", alignment_threshold)
        seqids = [seqid for seqid in self.info[binding_target][rnd]['seqids']]
        counts = [self.sequences[seqid].rounds[binding_target][rnd] for seqid in seqids]
        # Sort seqids in decreasing order of counts
        seqids = [seqid for seqid, count in sorted(zip(seqids, counts), key=lambda item: item[1], reverse=True)]
        total_counts = sum([self.sequences[seqid].rounds[binding_target][rnd] for seqid in seqids])
        alignment_record = {}
        covered = 0
        x, y = [0], [0]
        while (len(seqids) > 1) and len(alignment_record) < max_num_alignments:
            beginning = time.time()
            seqid1 = seqids[0]
            alignment_record[seqid1] = []
            indices = []
            for i, seqid2 in enumerate(seqids):
                seq1 = self.sequences[seqid1].sequence
                seq2 = self.sequences[seqid2].sequence
                if primers[0] not in seq1 or primers[1] not in seq1:
                    seq1 = primers[0] + seq1 + primers[1]
                if primers[0] not in seq2 or primers[1] not in seq2:
                    seq2 = primers[0] + seq2 + primers[1]
                dna1 = DNA(seq1)
                dna2 = DNA(seq2)
                alignment, score, start_end_positions = local_pairwise_align_ssw(dna1, dna2)
                score = score / (2 * len(seq1))
                if score > alignment_threshold:
                    indices.append(i)
                    alignment_record[seqid1].append(seqid2)
                current = time.time()
                eta = (current - beginning) / (i + 1) * (len(seqids) - i)
                print(
                    f"{i}/{len(seqids)} ETA: {round(eta, 1)}s w. avg. time {round((current - beginning) / (i + 1), 4)}s",
                    end="\r")
            print()
            covered += len(indices)
            print(
                f"Covered {round(covered / len(self.info[binding_target][rnd]['seqids']) * 100, 2)}% of all aptamers with {len(alignment_record)} alignments")
            new_seqids = [seqids[i] for i in range(0, len(seqids)) if i not in indices]
            seqids = new_seqids
            y.append(covered / len(self.info[binding_target][rnd]['seqids']))
            x.append(len(alignment_record))
        seqid1 = seqids[0]
        alignment_record[seqid1] = []
        indices = []
        for i, seqid2 in enumerate(seqids):
            indices.append(i)
            alignment_record[seqid1].append(seqid2)
        return alignment_record

    def analyze_priority_clusters(self, alignment_records, binding_target, rnd):
        print("Analyzing priority clusters for", binding_target, rnd)
        try:
            clusters = alignment_records[binding_target][rnd]
        except Exception:
            print(alignment_records.keys())
            print(alignment_records[binding_target].keys())
            clusters = alignment_records[binding_target][rnd]
        num_apts = sum([len(clusters[cluster]) for cluster in clusters])
        total_counts = sum(
            [sum([self.sequences[seqid].rounds[binding_target][rnd] for seqid in clusters[cluster]]) for cluster in
             clusters])
        max_counts = max(
            [max([self.sequences[seqid].rounds[binding_target][rnd] for seqid in clusters[cluster]]) for cluster in
             clusters])
        done = False
        clusters_kmers = []
        #clusters_stems_motifs = []
        clusters_forward_stems = []
        clusters_backward_stems = []
        clusters_loops = []
        for i, cluster_id in enumerate(list(clusters.keys())[:-1]):
            cluster_id = list(clusters.keys())[i]
            print("SUMMARY OF CLUSTER " + str(i))
            cluster_counts = sum([self.sequences[seqid].rounds[binding_target][rnd] for seqid in clusters[cluster_id]])
            uss = round(len(clusters[cluster_id]) / num_apts * 100, 2)
            tcs = round(cluster_counts / total_counts * 100, 2)
            if tcs > 30:  ### param ###
                done = True
            print(f"Unique sequence share: {uss}%   "
                  f"Total sequence share: {tcs}%   "
                  f"({len(clusters[cluster_id])})"
                  )
            sequences = [self.sequences[clusters[cluster_id][s]] for s in range(0, len(clusters[cluster_id]))]
            cluster_kmers, cluster_forward_stems, cluster_backward_stems, cluster_loops = self.print_cluster(sequences, binding_target, rnd, i, total_counts, tcs, max_counts)
            clusters_kmers.append(cluster_kmers)
            clusters_forward_stems.append(cluster_forward_stems)
            clusters_backward_stems.append(cluster_backward_stems)
            clusters_loops.append(cluster_loops)
            print("\n")
        A = self.find_clusters_correlation(clusters_kmers, k=5)
        B = self.find_clusters_correlation(clusters_forward_stems, k=5)
        C = self.find_clusters_correlation(clusters_backward_stems, k=5)
        D = self.find_clusters_correlation(clusters_loops, k=5)
        clusters_keys = list(clusters.keys())[:-1]
        new_clusters = {}
        added = []
        for i in range(0, len(clusters_keys)):
            for j in range(i + 1, len(clusters_keys)):
                ck1 = clusters_keys[i]
                ck2 = clusters_keys[j]
                if ck1 in added or ck2 in added:
                    continue
                if ck1 not in new_clusters.keys():
                    new_clusters[ck1] = clusters[ck1]
                    added.append(ck1)
                a = A[i][j]
                b = B[i][j]
                c = C[i][j]
                d = D[i][j]
                if a <= 2 or b <= 2 or c < 2 or d < 2:
                    new_clusters[ck1].extend(clusters[ck2])
                    added.append(ck2)
        new_clusters[list(clusters.keys())[-1]] = clusters[list(clusters.keys())[-1]]
        clusters = new_clusters
        print("\n\n\n")
        print("Updated clusters:")
        for i, cluster_id in enumerate(list(clusters.keys())[:-1]):
            cluster_id = list(clusters.keys())[i]
            print("SUMMARY OF CLUSTER " + str(i))
            cluster_counts = sum([self.sequences[seqid].rounds[binding_target][rnd] for seqid in clusters[cluster_id]])
            uss = round(len(clusters[cluster_id]) / num_apts * 100, 2)
            tcs = round(cluster_counts / total_counts * 100, 2)
            if tcs > 30:  ### param ###
                done = True
            print(f"Unique sequence share: {uss}%   "
                  f"Total sequence share: {tcs}%   "
                  f"({len(clusters[cluster_id])})"
                  )
            sequences = [self.sequences[clusters[cluster_id][s]] for s in range(0, len(clusters[cluster_id]))]
            cluster_kmers, cluster_forward_stems, cluster_backward_stems, cluster_loops  = self.print_cluster(sequences, binding_target, rnd, i, total_counts, tcs, max_counts)
            print("\n")

    def print_cluster(self, sequences, binding_target, rnd, cnum, total_counts, tcs, max_counts):
        primer_file = str(os.getcwd()) + f"/data/{binding_target[:4]}_primers.txt"
        with open(primer_file, 'r') as f:
            primers = [line.strip() for line in f]
        f.close()
        targ_seq = primers[0] + sequences[0].sequence + primers[1]
        alignments = []
        cts = []
        # Get alignments and cluster kmer counts
        cluster_kmer_counts = {}
        cluster_forward_stem_counts = {}
        cluster_backward_stem_counts = {}
        cluster_loop_counts = {}
        for s, seq in enumerate(sequences):
            seq = primers[0] + sequences[s].sequence + primers[1]
            dna1 = DNA(targ_seq)
            dna2 = DNA(seq)
            alignment, score, start_end_positions = local_pairwise_align_ssw(dna1, dna2)
            al = ''.join([str(x)[2] for x in alignment.loc[1].values])
            al = al[len(primers[0]): -len(primers[1])]
            kmers = sequences[s].get_kmers()
            ct = sequences[s].get_ct()
            idxs = [m.start() for m in re.finditer("-", al)]
            for m in idxs:
                ct = [c for c in ct]
                ct.insert(m, "-")
                ct = ''.join(ct)
            if len(al) != len(ct):
                dif = len(al) - len(ct)
                al = al[dif:]
            for kmer in kmers:
                if kmer not in cluster_kmer_counts.keys():
                    cluster_kmer_counts[kmer] = 0
                cluster_kmer_counts[kmer] += 1
            alignments.append(al)
            cts.append(ct)
            forward_stems, backward_stems, loops = sequences[s].get_hairpins()
            for stem in forward_stems:
                if stem not in cluster_forward_stem_counts.keys():
                    cluster_forward_stem_counts[stem] = 0
                cluster_forward_stem_counts[stem] += 1
            for stem in backward_stems:
                if stem not in cluster_backward_stem_counts.keys():
                    cluster_backward_stem_counts[stem] = 0
                cluster_backward_stem_counts[stem] += 1
            for loop in loops:
                if loop not in cluster_loop_counts.keys():
                    cluster_loop_counts[loop] = 0
                cluster_loop_counts[loop] += 1
        # Sort cluster_kmer_count sin decreasing order
        cluster_kmer_counts = {k: v for k, v in sorted(cluster_kmer_counts.items(), key=lambda item: item[1], reverse=True)}
        cluster_forward_stem_counts = {k: v for k, v in sorted(cluster_forward_stem_counts.items(), key=lambda item: item[1], reverse=True)}
        cluster_backward_stem_counts = {k: v for k, v in sorted(cluster_backward_stem_counts.items(), key=lambda item: item[1], reverse=True)}
        cluster_loop_counts = {k: v for k, v in sorted(cluster_loop_counts.items(), key=lambda item: item[1], reverse=True)}

        for s, al in enumerate(alignments):
            # Get variables
            al = alignments[s]
            ct = cts[s]
            # Get kmer scores
            kmers = sequences[s].get_kmers()
            kmer_dict = self.info[binding_target][rnd]['motif_counts']['kmer']['unweighted']
            # Kmer score
            ks = 100 * sum([kmer_dict[kmer] if kmer in kmer_dict.keys() else 0 for kmer in kmers]) / (
                        len(kmers) * max(kmer_dict.values()))
            # In cluster kmer score
            cks = 100 * sum([cluster_kmer_counts[kmer] for kmer in kmers]) / (
                            len(kmers) * max(cluster_kmer_counts.values()))
            # Secondary structure stability score
            dg = abs(sequences[s].dg)
            # Sequence count score
            sc = 100*sequences[s].rounds[binding_target][rnd]/max_counts
            #Strcutural scores
            forward_stems, backward_stems, loops = sequences[s].get_hairpins()
            # Forward stem score
            fss = 100 * sum([cluster_forward_stem_counts[motif] for motif in forward_stems]) / max((
                            len(forward_stems) * max(cluster_forward_stem_counts.values()),1))
            # Backward stem score
            bss = 100 * sum([cluster_backward_stem_counts[motif] for motif in backward_stems]) / max((
                            len(backward_stems) * max(cluster_backward_stem_counts.values()),1))
            # Loop score
            ls = 100 * sum([cluster_loop_counts[motif] for motif in loops]) / max((
                len(loops) * max(cluster_loop_counts.values()), 1))
            variable_region = ""
            # Generate colored output string
            for i in range(0, len(al)):
                if ct[i] == "(" or ct[i] == ")":
                    if al[i] == 'A':
                        variable_region += colored(al[i], 'red')
                    elif al[i] == 'T':
                        variable_region += colored(al[i], 'red')
                    elif al[i] == 'G':
                        variable_region += colored(al[i], 'red')
                    elif al[i] == 'C':
                        variable_region += colored(al[i], 'red' )
                else:
                    variable_region += al[i]
            out = colored(primers[0], 'cyan') + variable_region + colored(primers[1], 'cyan')
            popularity_scores = [tcs, sc]
            stability_scores = [dg]
            motif_scores = [ks, cks, (fss+bss)/2, ls]
            final_score = sum(popularity_scores)/len(popularity_scores) + sum(stability_scores)/len(stability_scores) + sum(motif_scores)/len(motif_scores)
            self.sequences[sequences[s].seqid].score = final_score
            # Print top 10 sequences in cluster
            if s < 10:
                print(f'seq{str(s):2}| {sequences[s].seqid} | {out}   |||   {tcs:.2f} // {sc:.2f} // {dg:.2f} // {ks:.2f}  // {cks:.2f} // {fss:.2f}  // {bss:.2f} // // {ls:.2f}  --> {final_score:.2f}')
        avgcs = sum([x.score for x in sequences])/len(sequences)
        print(f"Avg. cluster score: {avgcs:.2f}")
        print("...")
        return list(cluster_kmer_counts.keys()), list(cluster_forward_stem_counts.keys()), list(cluster_backward_stem_counts.keys()), list(cluster_loop_counts.keys())

    def find_clusters_correlation(self, clusters_motifs, k =10):
        clusters_motifs = [cluster_motifs[:k] for cluster_motifs in clusters_motifs]
        M = np.zeros((len(clusters_motifs), len(clusters_motifs)))
        for i in range(0, len(clusters_motifs)):
            for j in range(0, len(clusters_motifs)):
                kmers1 = [kmer for kmer in clusters_motifs[i]]
                kmers2 = [kmer for kmer in clusters_motifs[j]]
                score = 0
                for kmer in kmers1:
                    if kmer in kmers2:
                        score -= 1
                        kmers2.remove(kmer)
                    else:
                        score += 1
                M[i, j] = score
        return M


def get_priority_clusters(alignment_threshold=0.6):
    sequencelib_path = f'results/{TIME}/sequencelib.pickle'
    with open(sequencelib_path, 'rb') as f:
        sequencelib = pickle.load(f)
    f.close()
    configs = []
    for binding_target in sequencelib.info.keys():
        for rnd in sequencelib.info[binding_target].keys():
            configs.append({'binding_target': binding_target, 'rnd': rnd})
    alignment_records = {}
    begin = time.time()
    for i, config in enumerate(configs):
        print(f"Config {i} of {len(configs)} --> {config}")
        binding_target = config['binding_target']
        rnd = config['rnd']
        primer_file = str(os.getcwd()) + f"/data/{binding_target[:4]}_primers.txt"
        with open(primer_file, 'r') as f:
            primers = [line.strip() for line in f]
        f.close()
        alignment_record = sequencelib.priority_clustering(binding_target, rnd, max_num_alignments=20, primers=primers, alignment_threshold=alignment_threshold)
        current = time.time()
        time_per_config = (current - begin) / (i + 1)
        time_left = (len(configs) - i - 1) * time_per_config
        d = round(time_left // (24 * 3600))
        time_left -= d * 24 * 3600
        h = round(time_left // 3600)
        time_left -= h * 3600
        m = round(time_left // 60)
        time_left -= m * 60
        s = round(time_left)
        print(f"Clustering ETA: {d}d {h}h {m}m {s}s -> {time_per_config}s per config")
        if binding_target not in alignment_records:
            alignment_records[binding_target] = {}
        alignment_records[binding_target][rnd] = alignment_record
        # Save alignment_records as json
        with open(f'results/{TIME}/alignment_records.json', 'w') as f:
            json.dump(alignment_records, f)
        f.close()
        # Load alignment_records as json
    with open(f'results/{TIME}/alignment_records.json', 'r') as f:
        alignment_records = json.load(f)
    f.close()


def view_clusters(binding_target, rnd):
    with open(f'results/{TIME}/alignment_records.json', 'r') as f:
        alignment_records = json.load(f)
    f.close()
    with open(f'results/{TIME}/sequencelib.pickle', 'rb') as f:
        sequencelib = pickle.load(f)
    f.close()
    binding_target = "caffeine"
    rnd = "12"
    sequencelib.analyze_priority_clusters(alignment_records, binding_target, rnd)
    print("Recommendad sequences:")
    seqids = [seqid for seqid in sequencelib.info[binding_target][rnd]['seqids']]
    for seqid in seqids:
        if isinstance(sequencelib.sequences[seqid].score, dict):
            sequencelib.sequences[seqid].score = 0
    # Sort seqids based on self.sequences[seqid].final_scor in decreasing order
    seqids = sorted(seqids, key=lambda seqid: sequencelib.sequences[seqid].score, reverse=True)
    for i in range(0, 10):
        seqid = seqids[i]
        sequence = sequencelib.sequences[seqid].sequence
        final_score = sequencelib.sequences[seqid].score
        ct = sequencelib.sequences[seqid].ct
        rounds = sequencelib.sequences[seqid].rounds[binding_target]
        print(f"seq{i} # {seqids[i]} | {sequence} | {ct} | {final_score}  --> {rounds}")


def generate_all_recommendations(binding_targets):
    with open(f'results/{TIME}/alignment_records.json', 'r') as f:
        alignment_records = json.load(f)
    f.close()
    with open(f'results/{TIME}/sequencelib.pickle', 'rb') as f:
        sequencelib = pickle.load(f)
    f.close()
    for binding_target in binding_targets: #['ampicillin', 'chloramphenicol', 'liposome' , 'caffeine', 'theophylline', 'lactate', 'uricacid', 'oxytetracycline', 'uricacid']:
        if binding_target == 'caffeine':
            rnds = [12, 15, 20]
        elif binding_target == 'lactate':
            rnds = [14, 18, 19, 20]
        elif binding_target == 'oxytetracycline':
            rnds = [17, 21]
        elif binding_target == 'theophylline':
            rnds = ['06', 10, 12, 15, 16, 18, 20, 22]
        elif binding_target == 'uricacid':
            rnds = [16, 18]
        elif binding_target == 'ampicillin':
            rnds = [1]
        elif binding_target == 'chloramphenicol':
            rnds = [1]
        elif binding_target == 'liposome':
            rnds = [8, 11]
        for rnd in rnds:
            rnd = str(rnd)
            print(rnd, binding_target)
            sequencelib.analyze_priority_clusters(alignment_records, binding_target, rnd)
            # Get counts df
            rnd_ = str(rnds[-1])
            seqids_ = sequencelib.info[binding_target][rnd_]['seqids']
            data_dict = {}
            for seqid in seqids_:
                apt = sequencelib.sequences[seqid]
                seq = apt.sequence
                data_dict[seq] = {rnd_: -1}
                data_dict[seq][rnd_] = apt.rounds[binding_target][rnd_]
            df = pd.DataFrame(data_dict).T
            df.index = seqids_
            # Get recommendations
            print("Recommended sequences:")
            seqids = [seqid for seqid in sequencelib.info[binding_target][rnd]['seqids']]
            for seqid in seqids:
                if isinstance(sequencelib.sequences[seqid].score, dict):
                    sequencelib.sequences[seqid].score = 0
            # Sort seqids based on self.sequences[seqid].final_scor in decreasing order
            seqids = sorted(seqids, key=lambda seqid: sequencelib.sequences[seqid].score, reverse=True)
            # Crete directory recommendations if it doesn't exist
            if not os.path.exists(str(os.getcwd()) + f'/results/{TIME}/recommendations'):
                os.makedirs(str(os.getcwd()) + f'/results/{TIME}/recommendations')
            with open(str(os.getcwd()) + f"/results/{TIME}/recommendations/rec_{binding_target}_{rnd}.txt", 'w') as f:
                for i in range(0, 20):
                    seqid = seqids[i]
                    try:
                        idx = str(seqids_.index(seqid))
                    except Exception:
                        idx = str(-1)
                    sequence = sequencelib.sequences[seqid].sequence
                    final_score = sequencelib.sequences[seqid].score
                    ct = sequencelib.sequences[seqid].ct
                    rounds = sequencelib.sequences[seqid].rounds[binding_target]
                    print(f"seq{i} ({idx:>4s})| {seqids[i]} | {sequence} | {ct} | {final_score}  --> {rounds}")
                    f.write(f"seq{i} ({idx:>4s})| {seqids[i]} | {sequence} | {ct} | {final_score}  --> {rounds} \n")
            f.close()

            # Visualize round score distribution as a histogram
            scores = np.array([sequencelib.sequences[seqid].score for seqid in seqids])[:1000]
            mu, std = norm.fit(scores)
            xmin, xmax = 0, max(scores)
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu, std)
            fig, ax = plt.subplots()
            ax.hist(scores, bins=30, density=True)
            ax.plot(x, p, 'k', linewidth=3, label='score distribution')
            ymax = p.max()
            plt.vlines(x=mu, color='k', linewidth=1, ymin=0, ymax=ymax, label=f'mean score = {mu:.2f}')
            plt.vlines(x=mu + std, color='k', linestyle='dashed', ymin=0, ymax=p.max()/np.exp(1/2), linewidth=1)
            plt.vlines(x=mu - std, color='k', linestyle='dashed', ymin=0, ymax=p.max()/np.exp(1/2), linewidth=1)
            ax.set_title(f"Score distribution in {binding_target}_{rnd}")
            ax.set_xlim(xmin, xmax)
            ax.legend()
            plt.savefig(str(os.getcwd()) + f"/results/{TIME}/recommendations/rec_{binding_target}_{rnd}_hist.png")
            plt.close()


def main(binding_targets=['theophylline']):
    sequencelib = SequenceLibrary(with_primers=False)
    for binding_target in binding_targets:
        with open(f'data/{binding_target[:4]}list.txt') as f:
            rnds = [file.strip().replace(binding_target[:4], '') for file in f.readlines()]
        for rnd in rnds:
            rnd = str(rnd)
            sequencelib.add_data(binding_target, rnd, structure_mode='from_ct', top_k=None)
            cmd = "rm temp_ct.txt"
            subprocess.call([cmd], shell=True)
            cmd = "rm temp_seqs.txt"
            subprocess.call([cmd], shell=True)
    with open(f'results/{TIME}/sequencelib.pickle', 'wb') as f:
        pickle.dump(sequencelib, f)
    f.close()
    get_priority_clusters()
    generate_all_recommendations(binding_targets)
    with open(f'results/{TIME}/sequencelib.pickle', 'wb') as f:
        pickle.dump(sequencelib, f)
    f.close()


if __name__ == '__main__':
    binding_targets = sys.argv[1].split("-")
    main(binding_targets=binding_targets)
