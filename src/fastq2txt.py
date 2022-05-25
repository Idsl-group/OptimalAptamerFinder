import os, sys
from pyfastaq import sequences
import re
import random
def random_nucleotide(match):
    return random.sample(["A", "C", "G", "T"], k=1)[0]

def fastaq2txt(molecule='theophylline', indir=None, r='r1'):
    if r == 'r1':
        rel_range = [24, 54]
    elif r == 'r2':
        rel_range = [18, 48]


    if indir is None:
        outdir = str(os.getcwd()) + "/data/{}_txt_{}/".format(molecule, r)
    else:
        outdir = indir.replace("fastq", "txt")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for infile in os.listdir(indir):
        if infile.endswith(".fastq"):
            seqs = []
            seq_reader = sequences.file_reader(indir + "/" + infile)
            for sequence in seq_reader:
                seq = sequence.seq[rel_range[0]: rel_range[1]]
                if seq.count("N") < 3:
                    seq = re.sub("N", random_nucleotide, seq)
                    seqs.append(seq)
            outfile = infile.replace(".fastq", ".txt")
            with open(outdir + "/" + outfile, 'w') as fout:
                for seq in seqs:
                    fout.write(seq + "\n")

if __name__ == "__main__":
    molecule, r = sys.argv[1], sys.argv[2]
    indir = str(os.getcwd()) + "/data/{}_fastq_{}".format(molecule, r)
    fastaq2txt(molecule, indir=indir, r=r)
