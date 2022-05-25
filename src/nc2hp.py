import numpy as np
import pandas as pd
from itertools import groupby
from termcolor import colored


possible_pairings = [('A', 'C'), ('A', 'T'), ('A', 'G'), ('C', 'A'), ('C', 'T'), ('C', 'G'),
                     ('T', 'A'), ('T', 'C'), ('T', 'G'), ('G', 'A'), ('G', 'C'), ('G', 'T')]
binding_rules = {('A', 'C'):0, ('A', 'T'):1, ('A', 'G'):0, ('C', 'A'):0, ('C', 'T'):0, ('C', 'G'):1,
                 ('T', 'A'):1, ('T', 'C'):0, ('T', 'G'):1, ('G', 'A'):0, ('G', 'C'):1, ('G', 'T'):1,
                 ('A', 'A'):0, ('C', 'C'):0, ('T', 'T'):0, ('G', 'G'):0}


def find_binding_sites(seq):
    rev = seq[::-1]
    bindings = np.zeros((len(rev), len(seq)), dtype=int)
    for i in range(0, len(rev)):
        for j in range(0, len(seq)-i):
            if i == 0 or j == 0:
                bindings[i, j] = binding_rules[(rev[i], seq[j])]
            else:
                m = binding_rules[(rev[i], seq[j])]
                if m == 1:
                    bindings[i, j] = bindings[i - 1, j - 1] + m
                else:
                    bindings[i, j] = 0

    bindings_df = pd.DataFrame(bindings, index=list(rev), columns=list(seq))
    return bindings, bindings_df


class ExactStem():
    def __init__(self, seq, cs):
        rev = seq[::-1]
        self.fw_start = cs[0][1]
        self.rv_start = cs[0][0]
        self.fw_end = cs[-1][1]
        self.rv_end = cs[-1][0]
        self.fw_nc = seq[self.fw_start:self.fw_end+1]
        self.rv_nc = rev[self.rv_start:self.rv_end+1]
        assert len(self.fw_nc) == len(self.rv_nc) == len(cs)
        self.size = len(self.fw_nc)
        self.counts = {nc:self.fw_nc.count(nc) for nc in ['A', 'C', 'G', 'T']}
        self.consecutive_counts = {nc: sum(1 for _ in group) for nc, group in groupby(self.fw_nc)}
        self.gt_pairing = self.gt_pairing_percentage()

    def gt_pairing_percentage(self):
        c = 0
        for i in range(0, len(self.fw_nc)):
            if (self.fw_nc[i] == 'G' and self.rv_nc[i] == 'T') or (self.fw_nc[i] == 'T' and self.rv_nc[i] == 'G'):
                c += 1
        return c / self.size

    def gen_idxs(self):
        idxs =  [(self.rv_start+i, self.fw_start+i) for i in range(0, len(self.fw_nc))]
        assert len(idxs) == len(self.fw_nc)
        return idxs

    def check_criteria(self):
        if self.size < 3:
            return False
        return True

    def get_diagonal_index(self):
        crv, cfw = self.rv_start, self.fw_start
        while cfw != 0 and crv != 0:
            cfw -= 1
            crv -= 1
        return (crv, cfw)

    def get_diagonal_end_index(self, seq):
        l = len(seq)
        lrv, lfw = self.rv_end, self.fw_end
        d = 3 if (lfw + lrv) % 2 == 0 else 4
        while lfw + lrv < l - d:
            lfw += 1
            lrv += 1
        return (lrv, lfw)


def find_exact_stems(seq, bindings):
    rev = seq[::-1]
    # Find exact_stems
    exact_stems = []
    for i in range(0, len(rev)):
        for j in range(0, len(seq)):
            cs = []
            while i < len(rev) and j < len(seq):
                state = "off" if bindings[i, j] == 0 else "on"
                if state == "on":
                    cs.append((i, j))
                elif state == "off":
                    if len(cs) > 0:
                        stem = ExactStem(seq, cs)
                        if stem.check_criteria():
                            exact_stems.append(stem)
                        cs = []
                i += 1
                j += 1
    xsM = np.flip(np.eye(len(bindings), dtype=int), axis=0)*-1
    for stem in exact_stems:
        for idxs in stem.gen_idxs():
            xsM[idxs] = 1
    exact_stem_df = pd.DataFrame(data=xsM, index=list(rev), columns=list(seq))
    return exact_stems, xsM, exact_stem_df


class NonExactStem(ExactStem):
    def __init__(self, seq, xstem1, xstem2):
        rev = seq[::-1]
        if xstem1.fw_start < xstem2.fw_start:
            self.fw_start = xstem1.fw_start
            self.rv_start = xstem1.rv_start
            self.fw_end = xstem2.fw_end
            self.rv_end = xstem2.rv_end
            self.xstem1 = xstem1
            self.xstem2 = xstem2
        else:
            self.fw_start = xstem2.fw_start
            self.rv_start = xstem2.rv_start
            self.fw_end = xstem1.fw_end
            self.rv_end = xstem1.rv_end
            self.xstem1 = xstem2
            self.xstem2 = xstem1
        self.fw_nc = seq[self.fw_start:self.fw_end + 1]
        self.rv_nc = rev[self.rv_start:self.rv_end + 1]
        assert len(self.fw_nc) == len(self.rv_nc)
        self.size = len(self.fw_nc)
        self.counts = {nc: self.fw_nc.count(nc) for nc in ['A', 'C', 'G', 'T']}
        self.consecutive_counts = {nc: sum(1 for _ in group) for nc, group in groupby(self.fw_nc)}
        self.gt_pairing = (self.xstem1.gt_pairing*self.xstem1.size + self.xstem2.gt_pairing*self.xstem2.size) / self.size
        self.ex_non_ex_ratio = (self.xstem1.size + self.xstem2.size) / self.size

    def check_criteria(self):
        if self.size < 3:
            return False
        return True


def find_non_exact_stems(seq, exact_stems):
    # Find non-exact stems
    non_exact_stems = []
    for i in range(0, len(exact_stems)):
        for j in range(i+1, len(exact_stems)):
            xstem1 = exact_stems[i]
            xstem2 = exact_stems[j]
            if xstem1.get_diagonal_index() == xstem2.get_diagonal_index():
                non_exact_stem = NonExactStem(seq, xstem1, xstem2)
                if non_exact_stem.check_criteria():
                    non_exact_stems.append(non_exact_stem)
    return non_exact_stems


class Hairpin(ExactStem):
    def __init__(self, seq, non_exact_stem):
        rev = seq[::-1]
        self.stem = non_exact_stem
        self.loop_fw_start = self.stem.fw_end + 1
        self.loop_rv_start = self.stem.rv_end + 1
        self.loop_fw_end = len(seq) - self.stem.rv_end - 2
        self.loop_rv_end = len(rev) - self.stem.fw_end - 2
        self.fw_nc = seq[self.loop_fw_start:self.loop_fw_end + 1]
        self.rv_nc = rev[self.loop_rv_start:self.loop_rv_end + 1]
        try:
            assert self.fw_nc == self.rv_nc[::-1], f"{self.fw_nc} != {self.rv_nc[::-1]}"
        except Exception:
            print(seq)
            print(" "*self.stem.fw_start + self.stem.fw_nc + "^")
            print(" "*self.loop_fw_start + self.fw_nc + "^")
            print(rev)
            print(" "*self.stem.rv_start + self.stem.rv_nc + "^")
            print(" "*self.loop_rv_start + self.rv_nc + "^")
            assert self.fw_nc == self.rv_nc[::-1], f"{self.fw_nc} != {self.rv_nc[::-1]}"
        self.size = len(self.fw_nc)
        self.counts = {nc: self.fw_nc.count(nc) for nc in ['A', 'C', 'G', 'T']}
        self.consecutive_counts = {nc: sum(1 for _ in group) for nc, group in groupby(self.fw_nc)}
        self.loop_stem_ratio = self.size / (self.stem.size * 2 + self.size)
        self.stem_n_loop = (self.stem.fw_nc + self.stem.rv_nc[::-1], self.fw_nc)

    def view_hairpin(self, seq):
        #print(seq)
        if isinstance(self.stem, NonExactStem):
            out = seq[:self.stem.fw_start]
            out += colored(self.stem.xstem1.fw_nc, "red")
            out += seq[self.stem.xstem1.fw_end + 1:self.stem.xstem2.fw_start]
            out += colored(self.stem.xstem2.fw_nc, "magenta")
            out += colored(self.fw_nc, None, "on_cyan")
            out += colored(self.stem.xstem2.rv_nc[::-1], "magenta")
            out += seq[len(seq) - self.stem.xstem2.rv_start: len(seq) - self.stem.xstem1.rv_end - 1]
            out += colored(self.stem.xstem1.rv_nc[::-1], "red")
            out += seq[len(seq) - self.stem.rv_start:]
        else:
            out = seq[:self.stem.fw_start]
            out += colored(self.stem.fw_nc, "red")
            out += colored(self.fw_nc, None, "on_cyan")
            out += colored(self.stem.rv_nc[::-1], "red")
            out += seq[len(seq) - self.stem.rv_start:]
        print(out)

    def check_criteria(self):
        if self.size > self.stem.size*4:
            return False
        elif self.size > 20:
            return False
        elif self.loop_fw_end + 1 - self.loop_fw_start <=3:
            return False
        return True




def find_hairpins(seq, stems):
    hairpins = []
    for stem in stems:
        hairpin = Hairpin(seq, stem)
        if hairpin.check_criteria():
            hairpins.append(hairpin)
    return hairpins



if __name__ == '__main__':
    seq = "GGAGGCTCTCGGGACGACGATTGTGGTCTATTCATAGGCGTCCGCTGATTACGATTGCAGCATCGGGACG"
    ct = "(((.(((....((((.((....)).)))).......))).)))((((..........))))........."
    bindings, bindings_df = find_binding_sites(seq)
    exact_stems, xsM, exact_stem_df = find_exact_stems(seq, bindings)
    non_exact_stems = find_non_exact_stems(seq, exact_stems)
    hairpins_w_non_exact_stems = find_hairpins(seq, non_exact_stems)
    print(seq)
    for i, hairpin in enumerate(hairpins_w_non_exact_stems):
        print(f"Hairpin w. nxs {i}")
        hairpin.view_hairpin(seq)
    hairpins_w_exact_stems = find_hairpins(seq, exact_stems)
    for i, hairpin in enumerate(hairpins_w_exact_stems):
        print(f"Hairpin w. xs {i}")
        hairpin.view_hairpin(seq)