# -*- coding: utf-8 -*-
# @Author: zhangdw
# @Date:   2021-12-02 19:17:51
# @Last Modified by:   zhangdw
# @Last Modified time: 2022-07-13 11:40:06


# Import modules
import argparse
import re
import os
import pandas as pd
from Bio import SeqIO



# main body
def bgcCount(html):
    try:
        tables = pd.read_html(html, match="Region")
    except ValueError:
        bgc_count = 0
        tables = None
    else:
        bgc_count = tables[-1].shape[0]
    return(bgc_count, tables)


def bigscape_Class(bgcType):
    dict_class = dict()
    PKS_I = ['t1pks', 'T1PKS']
    PKS_other = ["transatpks", "t2pks", "t3pks", "otherks", "hglks", "transAT-PKS", "transAT-PKS-like", "T2PKS", "T3PKS", "PKS-like", "hglE-KS"]
    NRPS = ["nrps", "NRPS", "NRPS-like", "thioamide-NRP", "NAPAA"]
    RiPPs = ["lantipeptide", "thiopeptide", "bacteriocin", "linaridin", "cyanobactin", "glycocin", "LAP", "lassopeptide", "sactipeptide", "bottromycin",
     "head_to_tail", "microcin", "microviridin", "proteusin", "guanidinotides", "RiPP-like", "lanthipeptide", "lipolanthine", "RaS-RiPP", "fungal-RiPP",
      "thioamitides", "lanthipeptide-class-i", "lanthipeptide-class-ii", "lanthipeptide-class-iii", "lanthipeptide-class-iv", "lanthipeptide-class-v",
       "ranthipeptide", "redox-cofactor", "RRE-containing", "epipeptide", "cyclic-lactone-autoinducer", "spliceotide"]
    Saccharides = ["amglyccycl", "oligosaccharide", "cf_saccharide", "saccharide"]
    Terpene = ["terpene"]
    Others = ["acyl_amino_acids", "arylpolyene", "aminocoumarin", "ectoine", "butyrolactone", "nucleoside", "melanin", "phosphoglycolipid", "phenazine",
     "phosphonate", "other", "cf_putative", "resorcinol", "indole", "ladderane", "PUFA", "furan", "hserlactone", "fused", "cf_fatty_acid", "siderophore",
      "blactam", "fatty_acid", "PpyS-KS", "CDPS", "betalactone", "PBDE", "tropodithietic-acid", "NAGGN", "halogenated", "pyrrolidine", "prodigiosin"]
    for i in PKS_I:
        dict_class[i] = "PKS_I"
    for i in PKS_other:
        dict_class[i] = "PKS_other"
    for i in NRPS:
        dict_class[i] = "NRPS"
    for i in RiPPs:
        dict_class[i] = "RiPPs"
    for i in Saccharides:
        dict_class[i] = "Saccharides"
    for i in Terpene:
        dict_class[i] = "Terpene"
    for i in Others:
        dict_class[i] = "Others"

    re_class = []
    for x in bgcType.split(","):
        re_class.append(dict_class[x])
    re_class_set = set(re_class)
    if len(re_class_set) == 1:
        final_class = list(re_class_set)[0]
    else:
        if re_class_set == {"PKS_I", "NRPS"}:
            final_class = "PKS/NRPS Hybrids"
        else:
            final_class = "Others"
    return final_class

def parsegbkcluster(gbk):
    GCs = []
    length = 0
    CDS_number = 0
    protocluster_number = 0
    contig_edge = ""

    gbkcontents = SeqIO.parse(gbk, "genbank")
    for record in gbkcontents:
        contig_id = record.id
        comment = dict(record.annotations["structured_comment"])
        region_start = dict(comment["antiSMASH-Data"])["Orig. start"]
        region_end = dict(comment["antiSMASH-Data"])["Orig. end"]
        for feature in record.features:
            # Parsing the regionname
            if feature.type == "region":
                if "product" in feature.qualifiers:
                    for product in feature.qualifiers["product"]:
                        GCs.append(product)
                contig_edge = feature.qualifiers['contig_edge'][0]
                coor_start = int(str(feature.location.start).strip("<"))
                coor_end = int(str(feature.location.end).strip(">"))
                length = coor_end
                protocluster_number = str(len(GCs))

            if feature.type == "CDS":
                CDS_number += 1
    GCs.sort()
    return contig_id, region_start, region_end, ",".join(GCs), str(length), protocluster_number, str(CDS_number), contig_edge


def process(inPath, pattern, suffix, outPrefix):
    dict_bgc_count = dict()

    out_bgcCount = outPrefix + "__BGC_Count.tsv"
    out_bgcClass = outPrefix + "__BGC_Class.tsv"
    out_bgcStatistic = outPrefix + "__BGC_Region.tsv"
    dt_overall = None

    with open(out_bgcCount, 'w') as fo1, open(out_bgcClass, 'w') as fo2:
        print(*["Genome", "BGC_Count"], file=fo1, sep="\t", flush=True)
        print(*["Genome", "Contig", "Region_start", "Region_end", "GBK","Type","Class_bigscape_rule","Cluster_number","Contig_edge","Length","CDS_number","GBK_Path"], file=fo2, sep="\t", flush=True)

        for root, dirs, files in os.walk(inPath):
            for file in files:
                if re.search(pattern, file) and file.endswith(suffix):
                    gbk = file
                    gbk_path = os.path.abspath(os.path.join(root, file))
                    genome = root.split("/")[-1]
                    contig_id, region_start, region_end, product, length, cluster_number, CDS_number, contig_edge = parsegbkcluster(gbk_path)
                    bgc_class = bigscape_Class(product)
                    print(*[genome,contig_id,region_start,region_end,gbk,product,bgc_class,cluster_number,contig_edge,length,CDS_number,gbk_path], file=fo2, sep="\t", flush=True)
                elif file == "index.html":
                    bgc_number, dt_regions = bgcCount(os.path.abspath(os.path.join(root, file)))
                    genome_id = root.split("/")[-1]
                    print(*[genome_id, bgc_number], file=fo1, sep="\t", flush=True)

                    if bgc_number != 0:
                        dt_region = dt_regions[-1]
                        dt_region["Genomes"] = genome_id
                        if dt_overall is None:
                            dt_overall = dt_region
                        else:
                            dt_overall = pd.concat([dt_overall, dt_region])

        dt_overall.to_csv(out_bgcStatistic, sep="\t")



def main():
    parse = argparse.ArgumentParser(description="statistics for BGCs")

    parse.add_argument("--inPath", help="input path", required=True)
    parse.add_argument("--pattern", help="keyword of file name, default: region", default="region", required=False)
    parse.add_argument("--suffix", help="suffix of file name, default: .gbk", default=".gbk", required=False)
    parse.add_argument("--outPrefix", help="prefix of output file", required=True)

    args = parse.parse_args()

    process(args.inPath, args.pattern, args.suffix, args.outPrefix)


if __name__ == "__main__": 
    main()

