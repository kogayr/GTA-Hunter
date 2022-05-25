"""
The script checks cluster of genes based on the the supplied feature table.
Created: 11/06/2018
Modified: 06/02/2019
"""

import argparse
def parser_args():
    
    parser = argparse.ArgumentParser(description="""Clustering of genes via DBSCAN to get RcGTA-like clusters""")
    parser.add_argument("--data", type=str, required=True,
                        help="The text files with the chromosome ID first and genes to cluster.")
    parser.add_argument("--feature", type=str, required=True,
                        help="Supplied feature table from the NCBI.")
    parser.add_argument("-s", "--Cluster", type=int, 
                        dest="size", required=False, default=6,
                        help="""The minimum cluster size to consider as RcGTA-like cluster
                        (default=6).""")
    parser.add_argument("-e", "--Max", type=int,
                        dest="dist", required=False, default=8000,
                        help="""The maximum spatial distance to cluster genes
                        (default=8000).""")
    return parser

import itertools
import operator

def main():
    parser = parser_args()
    args = parser.parse_args()
    f = open(args.data, "r")
    lines = f.readlines()
    f.close()
    a = open(args.feature, "r")
    ids = a.readlines()
    a.close()
    min_size = args.size
    eps =args.dist
    for line in lines:
        line = line.rstrip()
        line = line.split()
        s = 1
        d = len(line)
        dist1 = []
        dist2 = []
        prot = {}
        while s < d:
            for each in ids:
                if line[0] in each and line[s] in each:
                    each = each.rstrip()
                    each = each.split("\t")
                    for w in range(len(each)):
                        if each[w].startswith('NC_') or each[w].startswith('NZ_'):
                            d1 = int(each[w+1])
                            d2 = int(each[w+2])
                            prot[line[s]] = d1
                            break
                        else:
                            continue 
                    dist1.append(d1)
                    dist2.append(d2)
                else:
                    continue
            s = s+1
            continue
        if len(dist1) == 0:
            print("Please supply proper feature table for", line[0])
        else:
            dist1.sort()
            dist2.sort()
            n = 0
            m = 1
            x = len(dist1) - 1
            dist_delta = []
            while n < x:
                b = dist1[m] - dist2[n]
                dist_delta.append(b)
                n = n +1
                m = m +1
            cons = []
            for i in range(len(dist_delta)):
                if dist_delta[i] < eps:
                    cons.append(i)
                    continue
                else:
                    continue
            all_cluster = []
            for a, b in itertools.groupby(enumerate(cons), lambda i_x: i_x[0]-i_x[1]):
                all_cluster.append((list(map(operator.itemgetter(1), b))))
            check = 0
            for i in range(len(all_cluster)):
                cluster_size = len(all_cluster[i]) + 1
                if cluster_size >= min_size:
                    print(line[0] + ' has RcGTA-like cluster; The cluster size is', cluster_size, 'genes')
                    list_of_prots = []
                    for j in range(len(all_cluster[i])):
                        ind = dist1[all_cluster[i][j]]
                        list_of_prots.append(list(prot.keys())[list(prot.values()).index(ind)])
                    last_ind = max(all_cluster[i]) + 1
                    last_prot = dist1[last_ind]
                    list_of_prots.append(list(prot.keys())[list(prot.values()).index(last_prot)])
                    print(*list_of_prots, sep='\t')
                    check = check + 1
            if check == 0:
                print(line[0] + ' has NOT RcGTA-like cluster')

if __name__ == '__main__':
    main()