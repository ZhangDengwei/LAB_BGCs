# -*- coding: utf-8 -*-
# @Author: zhangdw
# @Date:   2022-05-21 17:07:37
# @Last Modified by:   zhangdw
# @Last Modified time: 2022-05-25 14:29:05


# Import modules
import argparse
import re


# main body
def combine(summary_HMM, summary_BAGEL, prefix, Ecutoff):
    out_summary_combine = prefix + "_summary_combined.tsv"
    out_fasta_classI = prefix + "_combined_putative_precursors_classI.fa"
    out_fasta_classII = prefix + "_combined_putative_precursors_classII.fa"
    out_fasta_classIII = prefix + "_combined_putative_precursors_classIII.fa"

    class_I_domain = ["TIGR03651", "BacteriocIIc_cy", "Bacteriocin_IId", "Subtilosin_A", "CclA_1", "AS-BacteriocIIc_cy", 
                      "AS-BacteriocIIc_cy.aligned_c0","AS-BacteriocIIc_cy.aligned_c1","AS-BacteriocIIc_cy.aligned_c2",
                      "AS-Bacteriocin_IId","AS-Bacteriocin_IId.aligned_c0","AS-Bacteriocin_IId.aligned_c1","AS-Bacteriocin_IId.aligned_c2",
                      "AS-Bacteriocin_IId.aligned_c3","AS-Bacteriocin_IId.aligned_c4","AS-Bacteriocin_IId.aligned_c5","AS-Bacteriocin_IId.aligned_c6",
                      "AS-Bacteriocin_IId.aligned_c7","AS-Bacteriocin_IId.aligned_c8","AS-Bacteriocin_IId.aligned_c9",
                      "AS-TIGR03651","AS-TIGR03651.aligned_c0","AS-TIGR03651.aligned_c1","AS-TIGR03651.aligned_c2","AS-TIGR03651.aligned_c3",
                      "AS-TIGR03651.aligned_c4","AS-TIGR03651.aligned_c5","AS-TIGR03651.aligned_c6","AS-TIGR03651.aligned_c7",
                      "AS-TIGR03651.aligned_c8","AS-TIGR03651.aligned_c9","TIGR03601","TIGR03795","AS-TIGR03601",
                      "AS-TIGR03601.aligned_c0","AS-TIGR03601.aligned_c1","AS-TIGR03601.aligned_c2","AS-TIGR03601.aligned_c3","AS-TIGR03601.aligned_c4",
                      "AS-TIGR03601.aligned_c5","AS-TIGR03601.aligned_c6","AS-TIGR03601.aligned_c7","AS-TIGR03601.aligned_c8","AS-TIGR03601.aligned_c9"]
    class_II_domain = ["Antimicrobial17","Bacteriocin_IIi","Lactococcin_972","Bacteriocin_IIc","Bacteriocin_II","Antimicrobial14",
                        "LcnG-beta","Lactococcin","MccV",
                        "AS-Antimicrobial17","AS-Antimicrobial17.aligned_c0","AS-Antimicrobial17.aligned_c1","AS-Antimicrobial17.aligned_c2",
                        "AS-Bacteriocin_IIi","AS-Bacteriocin_IIi.aligned_c0","AS-Bacteriocin_IIi.aligned_c1","AS-Bacteriocin_IIi.aligned_c2",
                        "AS-Lactococcin_972","AS-Lactococcin_972.aligned_c0","AS-Lactococcin_972.aligned_c1","AS-Lactococcin_972.aligned_c2","AS-Lactococcin_972.aligned_c3",
                        "AS-Lactococcin_972.aligned_c4","AS-Lactococcin_972.aligned_c5","AS-Lactococcin_972.aligned_c6","AS-Lactococcin_972.aligned_c7",
                        "AS-Lactococcin_972.aligned_c8","AS-Lactococcin_972.aligned_c9",
                        "AS-Bacteriocin_IIc","AS-Bacteriocin_IIc.aligned_c0","AS-Bacteriocin_IIc.aligned_c1","AS-Bacteriocin_IIc.aligned_c2","AS-Bacteriocin_IIc.aligned_c3",
                        "AS-Bacteriocin_IIc.aligned_c4","AS-Bacteriocin_IIc.aligned_c5","AS-Bacteriocin_IIc.aligned_c6",
                        "AS-Bacteriocin_IIc.aligned_c7","AS-Bacteriocin_IIc.aligned_c8","AS-Bacteriocin_IIc.aligned_c9",
                        "AS-Bacteriocin_II","AS-Bacteriocin_II.aligned_c0","AS-Bacteriocin_II.aligned_c1","AS-Bacteriocin_II.aligned_c2",
                        "AS-Bacteriocin_II.aligned_c3","AS-Bacteriocin_II.aligned_c4","AS-Bacteriocin_II.aligned_c5","AS-Bacteriocin_II.aligned_c6",
                        "AS-Bacteriocin_II.aligned_c7","AS-Bacteriocin_II.aligned_c8","AS-Bacteriocin_II.aligned_c9",
                        "AS-LcnG-beta","AS-LcnG-beta.aligned_c0","AS-LcnG-beta.aligned_c1","AS-LcnG-beta.aligned_c2","AS-LcnG-beta.aligned_c3",
                        "AS-Lactococcin","AS-Lactococcin.aligned_c0","AS-Lactococcin.aligned_c1","AS-Lactococcin.aligned_c2","AS-Lactococcin.aligned_c3",
                        "AS-Lactococcin.aligned_c4","AS-Lactococcin.aligned_c5","AS-Lactococcin.aligned_c6","AS-Lactococcin.aligned_c7",
                        "AS-Lactococcin.aligned_c8","AS-Lactococcin.aligned_c9"]
    class_III_domain = ["Helveticin_J","Cloacin","Linocin_M18","AS-Cloacin","AS-Cloacin.aligned_c0","AS-Cloacin.aligned_c1","AS-Cloacin.aligned_c2","AS-Cloacin.aligned_c3",
                        "AS-Cloacin.aligned_c4","AS-Cloacin.aligned_c5","AS-Cloacin.aligned_c6","AS-Cloacin.aligned_c7","AS-Cloacin.aligned_c8","AS-Cloacin.aligned_c9",
                        "AS-Linocin_M18","AS-Linocin_M18.aligned_c0","AS-Linocin_M18.aligned_c1","AS-Linocin_M18.aligned_c2","AS-Linocin_M18.aligned_c3",
                        "AS-Linocin_M18.aligned_c4","AS-Linocin_M18.aligned_c5","AS-Linocin_M18.aligned_c6","AS-Linocin_M18.aligned_c7","AS-Linocin_M18.aligned_c8",
                        "AS-Linocin_M18.aligned_c9"]
    with open(summary_HMM, "r") as fin1, open(summary_BAGEL, "r") as fin2, open(out_summary_combine, "w") as fo0, open(out_fasta_classII, "w") as fo1, open(out_fasta_classIII, "w") as fo2, open(out_fasta_classI, "w") as fo3:
        dict_hmm = dict()
        dict_bagel = dict()
        for line in fin1:
            if not line.startswith("bgcID"):
                bgcID, ORF, startp, endp, strand, hitDomain, hitScore, evalue, sequence = line.strip("\n").split("\t")
                if bgcID not in dict_hmm:
                    dict_hmm[bgcID] = dict()
                if strand == "1":
                    strand_re = "+"
                elif strand == "-1":
                    strand_re = "-"
                else:
                    strand_re = strand
                unique_pos_id = startp+":"+endp+":"+strand_re
                if (hitDomain != ".") and (float(evalue) < float(Ecutoff)):
                    dict_hmm[bgcID][unique_pos_id] = [hitDomain, hitScore, evalue, sequence]
        for line in fin2:
            if not line.startswith("bgc_file"):
                bgc_file, AOI_start, AOI_end, AOI_class, AOI_geneTable, precursor_class, precursor_start, precursor_end, precursor_strand, precursor_seq = line.strip("\n").split("\t")
                if bgc_file not in dict_bagel:
                    dict_bagel[bgc_file] = dict()
                uniq_id_bagel = precursor_start+":"+precursor_end+":"+precursor_strand
                dict_bagel[bgc_file][uniq_id_bagel] = [precursor_class, precursor_seq]
        all_bgc = list(dict_hmm.keys())

        header = ["bgcID","Locus","Approach","hitDomain","hitScore","Evalue","Class_BAGEL4","Precursor_Class","Sequence"]
        print(*header, sep="\t", file=fo0, flush=True)
        for bgc in all_bgc:
            locus_hmm = list(dict_hmm[bgc].keys())
            locus_bagel = list(dict_bagel[bgc].keys()) if bgc in dict_bagel else []
            locus_common = list(set(locus_hmm).intersection(set(locus_bagel)))
            locus_uniq_hmm = list(set(locus_hmm)^set(locus_common))
            locus_uniq_bagel = list(set(locus_bagel)^set(locus_common))
            if locus_common:
                approach = "both"
                for x in locus_common:
                    hitDomain, hitScore, Evalue, sequence_hmm = dict_hmm[bgc][x]
                    precursor_Class, precursor_seq = dict_bagel[bgc][x]
                    class_bagel4 = precursor_Class.split(";")[0][-1]
                    final_class = ""

                    if (hitDomain in class_I_domain) and (class_bagel4 == "1"):
                        final_class = "Class_I"
                    elif (hitDomain in class_II_domain) and (class_bagel4 == "2"):
                        final_class = "Class_II"
                    elif (hitDomain in class_III_domain) and (class_bagel4 == "3"):
                        final_class = "Class_III"
                    else:
                        final_class = "inconsistent"
                    if final_class:
                        out_con = [bgc, x, approach, hitDomain, hitScore, Evalue, precursor_Class, final_class, sequence_hmm]
                        print(*out_con, sep="\t", file=fo0, flush=True)
                        header_fa = ">"+bgc+"__"+x
                        if final_class == "Class_II":
                            print(header_fa, file=fo1, flush=True)
                            print(sequence_hmm, file=fo1, flush=True)
                        elif final_class == "Class_III":
                            print(header_fa, file=fo2, flush=True)
                            print(sequence_hmm, file=fo2, flush=True)
                        elif final_class == "Class_I":
                            print(header_fa, file=fo3, flush=True)
                            print(sequence_hmm, file=fo3, flush=True)

            if locus_uniq_hmm:
                approach = "hmmsearch"
                for x in locus_uniq_hmm:
                    if x != ".:.:.":
                        hitDomain, hitScore, Evalue, sequence_hmm = dict_hmm[bgc][x]
                        final_class = ""
                        if hitDomain in class_I_domain:
                            final_class = "Class_I"
                        elif hitDomain in class_II_domain:
                            final_class = "Class_II"
                        elif hitDomain in class_III_domain:
                            final_class = "Class_III"
                        if final_class:
                            out_con = [bgc, x, approach, hitDomain, hitScore, Evalue, ".", final_class, sequence_hmm]
                            print(*out_con, sep="\t", file=fo0, flush=True)
                            header_fa = ">"+bgc+"__"+x
                            if final_class == "Class_II":
                                print(header_fa, file=fo1, flush=True)
                                print(sequence_hmm, file=fo1, flush=True)
                            elif final_class == "Class_III":
                                print(header_fa, file=fo2, flush=True)
                                print(sequence_hmm, file=fo2, flush=True)
                            elif final_class == "Class_I":
                                print(header_fa, file=fo3, flush=True)
                                print(sequence_hmm, file=fo3, flush=True)
            if locus_uniq_bagel:
                approach = "bagel4"
                for x in locus_uniq_bagel:
                    precursor_Class, precursor_seq = dict_bagel[bgc][x]
                    class_bagel4 = precursor_Class.split(";")[0][-1]
                    final_class = ""
                    if class_bagel4 == "1":
                        final_class = "Class_I"
                    elif class_bagel4 == "2":
                        final_class = "Class_II"
                    elif class_bagel4 == "3":
                        final_class = "Class_III"
                    if final_class:
                        out_con = [bgc, x, approach, ".", ".", ".", precursor_Class, final_class, precursor_seq]
                        print(*out_con, sep="\t", file=fo0, flush=True)
                        header_fa = ">"+bgc+"__"+x
                        if final_class == "Class_II":
                            print(header_fa, file=fo1, flush=True)
                            print(precursor_seq, file=fo1, flush=True)
                        elif final_class == "Class_III":
                            print(header_fa, file=fo2, flush=True)
                            print(precursor_seq, file=fo2, flush=True)
                        elif final_class == "Class_I":
                            print(header_fa, file=fo3, flush=True)
                            print(sequence_hmm, file=fo3, flush=True)


def main():
    parse = argparse.ArgumentParser(description="combine results of hmmsearch and bagel4")

    parse.add_argument("--summary_HMM", help="summary table of hmmsearch result", required=True)
    parse.add_argument("--summary_BAGEL", help="summary table of BAGEL4 result", required=True)
    parse.add_argument("--prefix", help="prefix of output files", required=True)
    parse.add_argument("--Ecutoff", help="cutoff of evalue for hmmsearch", default = 0.01, required=False)

    args = parse.parse_args()

    combine(args.summary_HMM, args.summary_BAGEL, args.prefix, args.Ecutoff)


if __name__ == "__main__": 
    main()

