#!/usr/bin/env python3

import pandas as pd
import gzip
import re
from tqdm import tqdm

class DeMultiSeq:

    def __init__(self, r1, r2, barcode_dict):
        self.r1 = r1
        self.r2 = r2
        self.barcode_dict = barcode_dict

        self.f1 = gzip.open(self.r1, "rt")
        self.f2 = gzip.open(self.r2, "rt")

        # {cell_barcode : multiseq_barcode : umi_set}
        self.data = {}

    def update_data(self, cb, umi, msb):
        if not self.valid_data(cb, umi, msb):
            return

        if cb not in self.data:
            self.data[cb] = {}

        if msb not in self.data[cb]:
            self.data[cb][msb] = set()

        if umi not in self.data[cb][msb]:
            self.data[cb][msb].add(umi)

    def valid_data(self, cb, umi, msb):

        if "N" in cb:
            return False
        if "N" in umi:
            return False
        if "N" in msb:
            return False
        if msb not in self.barcode_dict:
            return False

        return True

    def step(self):
        cell_barcode, umi = self.parse_r1()
        multiseq_barcode = self.parse_r2()

        self.update_data(cell_barcode, umi, multiseq_barcode)
        return True

    def parse_r1(self):
        header, seq, p, qual = [next(self.f1) for _ in range(4)]
        cell_barcode = seq[:16]
        umi = seq[16:26]

        return (cell_barcode, umi)

    def parse_r2(self):
        header, seq, p, qual = [next(self.f2) for _ in range(4)]
        multiseq_barcode = seq[:8]

        return multiseq_barcode

    def output(self):
        frame = []
        for cb in self.data:
            for msb in self.data[cb]:
                frame.append([
                    cb, msb, len(self.data[cb][msb])
                ])
        frame = pd.DataFrame(frame, columns = ["cell_barcode", "multiseq_barcode", "n_umi"])
        frame.to_csv("../output/output.tab", sep="\t", index=False)

    def parse(self):

        n_iter = 0
        while True:
            try:
                self.step()
                n_iter += 1
                if n_iter % 10000 == 0:
                    print(n_iter)
                    # break

            except StopIteration:
                break

        self.output()

def main():
    r1 = "../data/HCC_AA_MS_S6_L001_R1_001.fastq.gz"
    r2 = "../data/HCC_AA_MS_S6_L001_R2_001.fastq.gz"
    barcode_dict = {
        "CACTGTAG" : "Parental",
        "GCCAGTTA" : "LM2-B",
        "TGCCGTGG" : "LM2-C"
    }

    dms = DeMultiSeq(r1, r2, barcode_dict)
    dms.parse()

    pass

if __name__ == '__main__':
    main()

# CACTGTAG Parental
# GCCAGTTA LM2-B
# TGCCGTGG LM2-C
