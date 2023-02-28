# -*- coding: utf-8 -*-
# @Author: zhangdw
# @Date:   2022-02-21 19:41:28
# @Last Modified by:   zhangdw
# @Last Modified time: 2022-02-21 20:36:00


# Import modules
import argparse
import os
import re


# main body
def get(infolder, out):
    dict_bgc = dict()
    for dirpath, dirnames, filenames in os.walk(infolder):
        if len(dirpath.split("/")) == 7:
            if dirpath.split("/")[-1] not in dict_bgc:
                dict_bgc[dirpath.split("/")[-1]] = []
            for file in filenames:
                if re.search("region.+gbk", file):
                        dict_bgc[dirpath.split("/")[-1]].append(file)

    with open(out, "w") as fo:
        print("Genome\tGBC_Number", file=fo)
        for k, v in dict_bgc.items():
            out_c = [k, len(v)]
            print(*out_c, sep="\t", file=fo)


def main():
    parse = argparse.ArgumentParser(description="get bgc number from antismash output")

    parse.add_argument("--infolder", help="input folder", required=True)
    parse.add_argument("--out", help="output file", required=True)

    args = parse.parse_args()

    get(args.infolder, args.out)


if __name__ == "__main__": 
    main()

