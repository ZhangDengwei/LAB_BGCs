# -*- coding: utf-8 -*-
# @Author: zhangdw
# @Date:   2022-05-13 11:24:20
# @Last Modified by:   zhangdw
# @Last Modified time: 2022-05-19 14:50:28


# Import modules
import argparse
import os
import re
from Bio import SeqIO


# main body
# scan hmmsearch output files
def scanOut(inPath, suffix):
    for root, dirs, files in os.walk(inPath):
        for file in files:
            if file.endswith(suffix):
                file_path = os.path.abspath(os.path.join(root, file))
                yield file_path


# read FASTA file
def scanFasta(inPath, suffix):
    dict_fa = dict()
    for root, dirs, files in os.walk(inPath):
        for file in files:
            if file.endswith(suffix):
                bgc_id = re.search("bgc\d+", file)[0]
                dict_fa[bgc_id] = dict()
                file_path = os.path.abspath(os.path.join(root, file))
                for seq_record in SeqIO.parse(file_path, "fasta"):
                    dict_fa[bgc_id][seq_record.id] = [seq_record.seq.rstrip("*")]
                with open(file_path, "r") as fin:
                    for line in fin:
                        if line.startswith(">"):
                            items = line.split()
                            seq_id = items[0].strip().lstrip(">")
                            startp = items[2].strip()
                            endp = items[4].strip()
                            strand = items[6].strip()
                            dict_fa[bgc_id][seq_id].extend([startp, endp, strand])
    return dict_fa


# parase hmmsearch output files
def paraseHMM(inFile):
    dict_domain = dict()
    with open(inFile, "r") as fin:
        for line in fin:
            if not line.startswith("#"):
                items = line.split()
                cds_id = items[0].strip()
                domain_id = items[2].strip()
                evalue = items[4].strip()
                hit_score = items[5].strip()
                if cds_id not in dict_domain:
                    dict_domain[cds_id] = [domain_id, evalue, hit_score]
                else:
                    if hit_score > dict_domain[cds_id][2]:
                        dict_domain[cds_id] = [domain_id, evalue, hit_score]
                    else:
                        pass
    return dict_domain


def fetch(inPath_hmmsearch, suffix_hmmsearch, inPath_Fasta, suffix_Fasta, out1, out2):
    # out1: summary of hmmsearch results
    # out2: FASTA files containing putative precursors
    with open(out1, "w") as fo1, open(out2, "w") as fo2:
        header = ["bgcID", "ORF", "Start", "End", "Strand", "Domain_Hit", "Hit_Score", "E_value", "Sequence"]
        print(*header, sep="\t", file=fo1, flush=True)
        dict_fa = scanFasta(inPath_Fasta, suffix_Fasta)
        print(dict_fa, flush=True)
        for hmmOut in scanOut(inPath_hmmsearch, suffix_hmmsearch):
            #baseName = os.path.basename(hmmOut)
            bgcID = re.search("bgc\d+", hmmOut)[0]
            dict_domain = paraseHMM(hmmOut)
            if not dict_domain:
                outCon = [bgcID,".",".",".",".",".",".",".","."]
                print(*outCon, sep="\t", file=fo1, flush=True)
            else:
                for k, v in dict_domain.items():
                    cds_seq = dict_fa[bgcID][k][0]
                    cds_startp = dict_fa[bgcID][k][1]
                    cds_endp = dict_fa[bgcID][k][2]
                    cds_strand = dict_fa[bgcID][k][3]
                    outCon = [bgcID, k, cds_startp, cds_endp, cds_strand, v[0], v[2], v[1], cds_seq]
                    print(*outCon, sep="\t", file=fo1, flush=True)

                    fasta_header = ">" + bgcID + "__" + k + "__" + cds_startp + ":" + cds_endp
                    print(fasta_header, file=fo2, flush=True)
                    print(cds_seq, file=fo2, flush=True)
    

def main():
    parse = argparse.ArgumentParser(description="fetch out putative precursors from hmmsearch outputs")

    parse.add_argument("--inPath_hmmsearch", help="input path of hmmsearch", required=True)
    parse.add_argument("--suffix_hmmsearch", help="suffix of hmmsearch results, default: .tblout", default=".tblout", required=False)
    parse.add_argument("--inPath_Fasta", help="input path of FASTA", required=True)
    parse.add_argument("--suffix_Fasta", help="suffix of FASTA files, default: .fas", default=".fas", required=False)
    parse.add_argument("--out1", help="summary file of hmmsearch result", required=True)
    parse.add_argument("--out2", help="FASTA of precursor sequence", required=True)

    args = parse.parse_args()

    fetch(args.inPath_hmmsearch, args.suffix_hmmsearch, args.inPath_Fasta, args.suffix_Fasta, args.out1, args.out2)

if __name__ == "__main__": 
    main()

