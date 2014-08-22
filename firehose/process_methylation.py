__author__ = 'Eric T Dawson'
import argparse
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infile", default=None, type=str, required=True,
                        help="An Illumina methylation27 or 450 microarray file for preprocessing.")
    return parser.parse_args()


def process_samples(samples):
    d = OrderedDict()
    d["Gene"] = 0
    count = 1
    for s in samples:
        d[s.strip()] = count
        count += 1
    return "\t".join(d.keys())


def process_line(line):
    tokens = [x.strip() for x in line.split("\t")[1:]]
    count = 1
    vals = []
    gene = None
    if "TCGA" in tokens[0]:
        return process_samples(tokens)
    if "NA" in tokens or "Beta_value" in tokens:
        return None
    for t in tokens:
        if count % 4 == 3:
            pass
        elif count % 4 == 0:
            pass
        elif count % 4 == 2:
            gene = t
            if ";" in gene:
                return None
        else:
            vals.append(t)
        count += 1

    vals.insert(0, gene)
    return "\t".join(vals)


def process_methylation_file(infile, outfile):
    gene_to_ind = {}
    ind = 0
    for line in infile:
        proc = process_line(line)
        if proc is not None:
            gene = proc[0] if "TCGA" not in proc[2] else None
            out_line = "\t".join(proc)
            outfile.write(out_line + "\n")
            if gene is not None:
                index_list = gene_to_ind.get(gene, [])
                index_list.append(ind)

if __name__ == "__main__":
    args = parse_args()
    f = open(args.infile, "r")
    outfile = open("test.methylation.out.txt", "w")
    process_methylation_file(f, outfile)


