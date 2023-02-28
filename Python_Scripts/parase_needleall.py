# -*- coding: utf-8 -*-
# @Author: zhangdw
# @Date:   2022-06-13 10:53:02
# @Last Modified by:   zhangdw
# @Last Modified time: 2022-06-13 11:16:15


# Import modules
import argparse
import re


# main body
def parase(inPut, outPut):
    pair_dict = dict()
    last_seq_1 = ""
    last_seq_2 = ""
    with open(inPut, "r") as fin, open(outPut, "w") as fo:
        for line in fin:
            if re.search("Aligned_sequences", line):
                last_seq_1 = ""
                last_seq_2 = ""
            elif re.search("1:", line):
                last_seq_1 = re.search("1:\s+(.+)", line.rstrip("\n")).group(1)
                if last_seq_1 not in pair_dict:
                    pair_dict[last_seq_1] = dict()
            elif re.search("2:", line):
                last_seq_2 = re.search("2:\s+(.+)", line.rstrip("\n")).group(1)
                if last_seq_2 not in pair_dict[last_seq_1]:
                    pair_dict[last_seq_1][last_seq_2] = []
            elif re.search("Length", line):
                length = re.search("Length:\s+(.+)", line.rstrip("\n")).group(1)
                pair_dict[last_seq_1][last_seq_2].append(length)
            elif re.search("Identity:", line):
                identity = re.search(".+\((.+)%\)", line.rstrip("\n")).group(1)
                pair_dict[last_seq_1][last_seq_2].append(identity)
            elif re.search("Similarity:", line):
                similarity = re.search(".+\((.+)%\)", line.rstrip("\n")).group(1)
                pair_dict[last_seq_1][last_seq_2].append(similarity)
            elif re.search("Gaps:", line):
                gap = re.search(".+\((.+)%\)", line.rstrip("\n")).group(1)
                pair_dict[last_seq_1][last_seq_2].append(gap)

        header = ["Seq_1", "Seq_2", "Length", "Identity(%)", "Similarity(%)", "Gap"]
        print(*header, file=fo, sep="\t", flush=True)
        for k,v in pair_dict.items():
            seq_1 = k
            for x,y in v.items():
                seq_2 = x
                length, identity, similarity, gap = y
                out_con = [seq_1, seq_2, length, identity, similarity, gap]
                print(*out_con, sep="\t", file=fo, flush=True)


def main():
    parse = argparse.ArgumentParser(description="parase the output of needleall")

    parse.add_argument("--inPut", help="the result of needleall", required=True)
    parse.add_argument("--outPut", help="the output file", required=True)

    args = parse.parse_args()

    parase(args.inPut, args.outPut)


if __name__ == "__main__": 
    main()

