import os, sys
sys.path.append("/home/javier/PycharmProjects/OptimalAptamerFinder")
from pyfastaq import sequences
import re
import random

def random_nucleotide(match):
    return random.sample(["A", "C", "G", "T"], k=1)[0]

def cleanfastq(indir, outdir=None):
    rel_range = [24, 54]

    if outdir is None:
        outdir = indir[:-4]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    i = 0
    for infile in os.listdir(indir):
        newfile = outdir + "/" + str(infile)
        if infile.endswith(".fastq"):
            with open(indir + "/" + infile, "r") as inf:
                lines = inf.readlines()
                with open(newfile, "w") as outf:
                    for i in range(0, len(lines), 4):
                        data = lines[i:i+4]
                        if data[1][rel_range[0]:rel_range[1]].count("N") > 2:
                            continue
                        else:
                            data[1] = re.sub("N", random_nucleotide, data[1])
                            outf.write("".join(data))
                outf.close()
            inf.close()


if __name__ == "__main__":
    indir = str(os.getcwd()) + "/data/theophylline_fastq_r2_bad"
    cleanfastq(indir)