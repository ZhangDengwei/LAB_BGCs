# -*- coding: utf-8 -*-
# @Author: zhangdw
# @Date:   2021-12-26 16:46:53
# @Last Modified by:   zhangdw
# @Last Modified time: 2022-02-23 22:24:43


# Import modules
import argparse
import os
import sys


# main body
def hmm_stats(inFolder, outPrefix):
    hmmFeature = os.path.join(inFolder, "bgc_features.csv")
    hmmReference = os.path.join(inFolder, "hmm.csv")
    bgc_des = os.path.join(inFolder, "bgc.csv")

    bgc_all = dict()
    with open(bgc_des, 'r') as fin:
        fin.readline()
        for line in fin:
            index, bgc_id, dataset_id, name = line.rstrip("\n").split(",")[:4]
            bgc_all[bgc_id] = name

    if (not os.path.exists(hmmFeature)) or (not os.path.exists(hmmReference)):
        raise FileError("Did not find file 'bgc_features.csv' or 'hmm.csv', please check!!!")
        sys.exit(0)
    # prase hmm found
    hmm_id_name = dict()
    with open(hmmReference, 'r') as fin:
        fin.readline()
        for line in fin:
            index, hmm_id, pfam_accession, pfam_name, db_id, model_length = line.rstrip("\n").split(",")
            hmm_id_name[hmm_id] = pfam_name

    # generate pfam distribution matrix
    dict_bgc_hmm = dict()
    with open(hmmFeature, "r") as fin:
        fin.readline()
        for line in fin:
            index, bgc_id, hmm_id, value = line.rstrip("\n").split(",")
            if bgc_id not in dict_bgc_hmm:
                dict_bgc_hmm[bgc_id] = []
            dict_bgc_hmm[bgc_id].append(hmm_id)

    print(f"{len(bgc_all.keys())-len(dict_bgc_hmm.keys())} BGCs have not found features, excluded from analysis!!!")

    out_matrix = outPrefix + "__hmm_presence_matrix.tsv"
    with open(out_matrix, 'w') as fo:
        hmm_id_list = list(hmm_id_name.keys())
        header = ["bgc_id"] + [hmm_id_name[x] for x in hmm_id_list]
        print(*header, sep="\t", file=fo)
        for k, v in dict_bgc_hmm.items():
            cont = [k]
            for x in hmm_id_list:
                presence = 1 if x in v else 0
                cont.append(presence)
            print(*cont, sep="\t", file=fo)


def main():
    parse = argparse.ArgumentParser(description="statistics for hmm presence/absence distribution calculated by bigslice")

    parse.add_argument("--inFolder", help="input folder containg prased results of bigslice", required=True)
    parse.add_argument("--outPrefix", help="output prefix", required=True)

    args = parse.parse_args()

    hmm_stats(args.inFolder, args.outPrefix)    


if __name__ == "__main__": 
    main()

