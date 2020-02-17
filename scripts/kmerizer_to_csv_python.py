#!/usr/bin/env python3

import sys
import subprocess
import argparse
import pandas as pd
import numpy as np
parser = argparse.ArgumentParser(description='Get Kmers and Format into Dictionary')
parser.add_argument('-f','--file', type=str,
                    help='Fasta File')
parser.add_argument('-k', "--size", type=int, default=8,
                    help='Kmer Size')
parser.add_argument('-o', "--output", default='kmer_output_table',
                    help="Output CSV File Name")

args = parser.parse_args()

x = pd.read_csv("data/Columns_to_Keep.txt", header=None)
columns_to_keep = list(x[0])




k = args.size
fastaFile = args.file
kmerCmd = 'scripts/kmer-counter-master/kmer-counter --fasta --k=%d %s' % (k, fastaFile)

try:
    output = subprocess.check_output(kmerCmd, shell=True)
    output = output.decode()
    result = {}
    for line in output.splitlines():
        header, counts = line.strip().split('\t')
        header = header[1:]
        kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in counts.split(' ')])
        result[header] = kmers

    diction = pd.DataFrame.from_dict(result, orient="index")
    diction.replace(np.nan, 0, inplace=True)
    mytrans = str.maketrans('ATCG', 'TAGC')
    ## Merge Compliments
    diction.columns = np.sort([diction.columns, [x.translate(mytrans) for x in diction.columns]], axis=0)[0, :]
    diction = diction.groupby(level=0, axis=1).sum()
    diction = diction.reindex(columns_to_keep, axis=1)
    ## Normalize by Length
    diction = diction.div(diction.sum(axis=1), axis=0)
    diction = diction.round(7)




    name = args.output + ".csv"
    diction.to_csv(name)
    
except subprocess.CalledProcessError as error:
    sys.stderr.write("%s\n" % (str(error)))
