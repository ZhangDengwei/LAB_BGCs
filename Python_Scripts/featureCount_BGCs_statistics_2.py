# -*- coding: utf-8 -*-
# @Author: zhangdw
# @Date:   2021-11-26 14:14:51
# @Last Modified by:   zhangdw
# @Last Modified time: 2022-03-19 17:03:01


# Import modules
import argparse
import pandas as pd
import numpy
import re


# main body

def readGFF(gff):
    # read GFF file to identify core genes and biosynthetic genes
    bio_genes = {}
    core_genes = {}

    with open(gff, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                items = line.rstrip('\n').split('\t')
                if items[2] == "CDS":
                    bgc, attribute = items[0], items[8]
                    gene_id = re.search("ID=(.+?);", attribute).group(1)
                    gene_kind = re.search("gene_kind=(.+?);", attribute).group(1)
                    if gene_kind == "biosynthetic-additional":
                        if bgc not in bio_genes:
                            bio_genes[bgc] = []
                        bio_genes[bgc].append(gene_id)
                    elif gene_kind == "biosynthetic":
                        if bgc not in bio_genes:
                            bio_genes[bgc] = []
                        bio_genes[bgc].append(gene_id)
                        if bgc not in core_genes:
                            core_genes[bgc] = []  
                        core_genes[bgc].append(gene_id)
                    else:
                        pass
    return bio_genes, core_genes


def featureCount(featureCount, readT, coverageT, prefix, gff):
    '''
    featureCount: expression matrix from featureCount
    readT: threashold for read number to determine whether a CDS is present, default: >10
    coverageT: threashold for percentage of expressed CDS to determine whether a BGC is present, default: > 0.5
    '''
    d1 = pd.read_table(featureCount, sep='\t', comment="#")

    # transfer read counts to TPM
    print("Calculating TPM...")
    for i in d1.columns[6:]:
        d1[i] = d1[i]*1000/d1['Length']
        sumi = d1[i].sum() + 1 # Sum of some columns might be Zero 
        d1[i] = d1[i].apply(lambda x:x*1000000/sumi)

    bgc_expression = {}
    unique_BGC = []

    # prase GFF file
    bio_genes, core_genes = readGFF(gff)

    with open(featureCount, 'r') as fin:
        fisrt_line = fin.readline()
        second_line = fin.readline()
        # frist six lines are metadat [Geneid, Chr, Start, End, Strand, Length]
        sampleList = second_line.rstrip('\n').split('\t')[6:]
        # creat dictionary for each sample to store information
        for x in sampleList: bgc_expression[x] = {}

        row_index = 0
        for line in fin:
            Geneid, Chr, Start, End, Strand, Length = line.rstrip('\n').split('\t')[:6]
            counts = line.rstrip('\n').split('\t')[6:]
            if Chr not in unique_BGC: unique_BGC.append(Chr)

            try:
                bio_genes_Chr = bio_genes[Chr]
            except KeyError:
                print(f"Warning: Did not find biosynthetic genes for {Chr} !!!")
                bio_genes_Chr = []

            try:
                core_genes_Chr = core_genes[Chr]
            except KeyError:
                print(f"Warning: Did not find core biosynthetic genes for {Chr} !!!")
                core_genes_Chr = []

            for i in range(len(sampleList)):
                #print(int(counts[i]))
                col_index = i + 6
                tpm = float(d1.iloc[row_index, col_index])
                sample = sampleList[i]

                if Chr not in bgc_expression[sample]:
                    # creat table
                    bgc_expression[sample][Chr] = {"CDS_list":[],
                                                    "CDS_present":[],
                                                    "CDS_counts":[],
                                                    "CDS_TPM":[],
                                                    "biosynthesis_genes_present":[],
                                                    "biosynthesis_genes_counts":[],
                                                    "biosynthesis_genes_TPM":[],
                                                    "core_genes_present":[],
                                                    "core_genes_counts":[],
                                                    "core_genes_TPM":[]}
                bgc_expression[sample][Chr]["CDS_list"].append(Geneid)
                if int(counts[i]) >= readT: 
                    bgc_expression[sample][Chr]["CDS_present"].append(1) # 1 for present
                else:
                    bgc_expression[sample][Chr]["CDS_present"].append(0) # 0 for absent
                bgc_expression[sample][Chr]["CDS_counts"].append(int(counts[i]))
                bgc_expression[sample][Chr]["CDS_TPM"].append(tpm)

                # store biosynthetic genes
                if Geneid in bio_genes_Chr:
                    bgc_expression[sample][Chr]["biosynthesis_genes_counts"].append(int(counts[i]))
                    bgc_expression[sample][Chr]["biosynthesis_genes_TPM"].append(tpm)
                    if int(counts[i]) >= readT:
                        bgc_expression[sample][Chr]["biosynthesis_genes_present"].append(1)
                    else:
                        bgc_expression[sample][Chr]["biosynthesis_genes_present"].append(0)
                # store core genes
                if Geneid in core_genes_Chr:
                    bgc_expression[sample][Chr]["core_genes_counts"].append(int(counts[i]))
                    bgc_expression[sample][Chr]["core_genes_TPM"].append(tpm)
                    if int(counts[i]) >= readT:
                        bgc_expression[sample][Chr]["core_genes_present"].append(1)
                    else:
                        bgc_expression[sample][Chr]["core_genes_present"].append(0)

            row_index += 1

    # write out matrix
    print("Output tables...")
    out_bgc_present_1 = open(prefix+"__BGCs_present.tsv", 'w')
    out_bgc_present_2 = open(prefix+"__BGCs_presence_biosynthetic_genes.tsv", 'w')
    out_bgc_count = open(prefix+"__BGCs_counts_considering_all_biosynthetic_genes.tsv", 'w')
    out_bgc_TPM = open(prefix+"__BGCs_TPM_considering_all_biosynthetic_genes.tsv", 'w')
    out_core_counts = open(prefix+"__BGCs_counts_considering_core_genes.tsv", 'w')
    out_core_TPM = open(prefix+"__BGCs_TPM_considering_core_genes.tsv", 'w')

    note = "# 1 refer to BGC present in sample under the criteria:\n" + \
            f"# 1) CDS is present when mapped reads are > {readT};\n" + \
            f"# 2) number of biosynthetic genes found is more than {coverageT} of total biosynthetic genes in BGC\n" + \
            f"# 3) at least one core biosynthetic gene must be found\n" + \
            "# 0 refer to BGC not found in sample\n"
            
    out_bgc_present_1.write(note)
    out_bgc_present_2.write("# [A/B] A refers to biosynthetic_genes present, B refers to total biosynthetic_genes number\n")

    header = ["Samples"] + unique_BGC
    out_bgc_present_1.write("\t".join(header)+"\n")
    out_bgc_present_2.write("\t".join(header)+"\n")
    out_bgc_count.write("\t".join(header)+"\n")
    out_bgc_TPM.write("\t".join(header)+"\n")
    out_core_counts.write("\t".join(header)+"\n")
    out_core_TPM.write("\t".join(header)+"\n")

    for k, v in bgc_expression.items():
        # k: sample
        out_bgc_present_cont_1 = k
        out_bgc_present_cont_2 = k
        out_bgc_count_cont = k
        out_bgc_TPM_cont = k
        out_core_counts_cont = k
        out_core_TPM_cont = k

        for bgc in unique_BGC:
            CDS_list = v[bgc]["CDS_list"]
            biosynthesis_genes_present = v[bgc]["biosynthesis_genes_present"]
            biosynthesis_genes_counts = v[bgc]["biosynthesis_genes_counts"]
            biosynthesis_genes_TPM = v[bgc]["biosynthesis_genes_TPM"]
            core_genes_present = v[bgc]["core_genes_present"]
            core_genes_counts = v[bgc]["core_genes_counts"]
            core_genes_TPM = v[bgc]["core_genes_TPM"]

            # Two criteria to determine whether a BGC is present in a sample
            # 1) expressed biosynthesis genes must exceed coverage threashold, say 50%
            # 2) core gene must be found
            coverage = numpy.mean(biosynthesis_genes_present)
            if coverage >= float(coverageT) and sum(core_genes_present) > 0:
                bgc_present = 1
            else:
                bgc_present = 0
            bgc_present_detail = str(sum(biosynthesis_genes_present)) + "/" + str(len(biosynthesis_genes_present))
            bgc_average_count = numpy.mean(biosynthesis_genes_counts)
            bgc_average_TPM = numpy.mean(biosynthesis_genes_TPM)
            try:
                core_genes_counts = numpy.mean(core_genes_counts)
                core_genes_TPM = numpy.mean(core_genes_TPM)
            except:
                print(f"Warning: Did not find core genes of {bgc} !!!")
                core_genes_counts = "NA"
                core_genes_TPM = "NA"
            out_bgc_present_cont_1 += "\t" + str(bgc_present)
            out_bgc_present_cont_2 += "\t" + str(bgc_present_detail)
            out_bgc_count_cont += "\t" + str(bgc_average_count)
            out_bgc_TPM_cont += "\t" + str(bgc_average_TPM)
            out_core_counts_cont += "\t" + str(core_genes_counts)
            out_core_TPM_cont += "\t" + str(core_genes_TPM)

        out_bgc_present_1.write(out_bgc_present_cont_1+'\n')
        out_bgc_present_2.write(out_bgc_present_cont_2+'\n')
        out_bgc_count.write(out_bgc_count_cont+'\n')
        out_bgc_TPM.write(out_bgc_TPM_cont+'\n')
        out_core_counts.write(out_core_counts_cont+'\n')
        out_core_TPM.write(out_core_TPM_cont+'\n')

    out_bgc_present_1.close()
    out_bgc_present_2.close()
    out_bgc_count.close()
    out_bgc_TPM.close()
    out_core_counts.close()
    out_core_TPM.close()


def main():
    parse = argparse.ArgumentParser(description="Process BGC expression table generated by featureCount, six files would be wrote out")

    parse.add_argument("--gff", help="gff file which is the input for featureCount", required=True)
    parse.add_argument("--featureCount", help="featureCount summary table", required=True)
    parse.add_argument("--readT", help="threashold for read number to determine whether a CDS is present, default: >10",
                        default=10, type=int, required=False)
    parse.add_argument("--coverageT", help="threashold for percentage of expressed CDS to determine whether a BGC is present, default: > 0.5",
                        default=0.5, required=False)
    parse.add_argument("--prefix", help="prefix of output files", required=True)

    args = parse.parse_args()

    # processing

    print("Starting...")
    featureCount(args.featureCount, args.readT, args.coverageT, args.prefix, args.gff)

    print("Finished successfully!")

if __name__ == "__main__": 
    main()

