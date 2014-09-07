__author__ = 'Eric T Dawson'
import argparse
from collections import OrderedDict

## Pro tip: Move your infile to an ssd or RAMdisk for maximum performance!


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infile", default=None, type=str, required=True,
                        help="An Illumina methylation27 or 450 microarray file for preprocessing.")
    return parser.parse_args()


def process_samples(samples):
    """Takes the header line from the methylation file and generates an ordered set of sample names."""
    d = OrderedDict()
    d["Gene"] = 0
    count = 1
    for s in samples:
        d[s.strip()] = count
        count += 1
    return "\t".join(d.keys())


def process_line(line):
    """"For a given line in the file, this function:
    1. Ignores probes matching to multiple genes or unknown genes
    2. Handles sample processing if the line is the header line.
    3. Removes genomic position/chromosome data
    and returns a tab-separated string of beta-values with the gene name at the first position"""
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


def check_vals(line, proc):
    """Checks a new line of beta values against the existing line in the dictionary. Any values that
    are lesser in the old line than the new one are replaced in-place."""
    assert len(line) == len(proc)
    for i in xrange(0, len(line)):
        if proc[i] > line[i]:
            line[i] = proc[i]


def process_methylation_file(infile, outfile=open("test.methylation.tmp.txt", "w")):
    """Takes in an infile (Illumina 27/450) and performs a series of steps to shrink the file
    and keep only maximum values before writing to outfile. The function intentionally
    avoids loading the infile into memory, as an Illumina450 file can greatly exceed the amount
    of RAM on a typical desktop."""
    gene_to_line = {}
    samples = None
    ind = 0
    for line in infile:
        proc = process_line(line)
        if proc is not None:
            if "TCGA" in proc[2]:
                samples = proc
            else:
                gene = proc[0]
                if gene not in gene_to_line:
                    gene_to_line[gene] = proc
                else:
                    check_vals(gene_to_line[gene], proc)

    for g in gene_to_line:
        outfile.write("".join([gene_to_line[g], "\n"]))
            # gene = proc[0] if "TCGA" not in proc[2] else None
            # out_line = "\t".join(proc)
            # outfile.write(out_line + "\n")
            # if gene is not None:
            #     index_list = gene_to_ind.get(gene, [])
            #     index_list.append(ind)

if __name__ == "__main__":
    args = parse_args()
    f = open(args.infile, "r")
    outfile = open("test.methylation.out.txt", "w")
    process_methylation_file(f, outfile)


