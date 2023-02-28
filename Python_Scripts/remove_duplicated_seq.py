#!/usr/bin/env python3

"""
Modified from BiG-MAP
"""

# Import statements
import os
import subprocess
import sys
import argparse
from Bio import SeqIO
import json
import re
from pathlib import Path
import fastaparser
import numpy as numpy
from Bio.SeqFeature import FeatureLocation, ExactPosition, BeforePosition, AfterPosition
import pickle
from subprocess import Popen, PIPE
import shutil
from glob import glob
import ntpath
import itertools
import warnings


def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
                                     usage='''
______________________________________________________________________
Note: Only remove redundant FASTA files containing one sequence or
      more than one sequence for genome duplication
______________________________________________________________________
Generic command: python3 remove_duplicated_seq.py [Options]*
-D [input dir(s)] -O [output dir]
Create a redundancy filtered fasta reference file from multiple
anti/gutSMASH outputs. Use BiG-MAP_process conda environment.
Obligatory arguments:
    -D   Specify the path to the directory containing the FASTA files
    -O   Put path to the folder where the MASH filtered gene cluster
         files should be located here. The folder should be an
         existing folder. Default = current folder (.)
Options:
    -tg  Fraction between 0 and 1; the similarity threshold (Mash distance,
         Mash D roughly equal to 1-ANI) that determines when the sequences
         can be considered similar. Default = 0.01.
    -s   Specify the sketch size created by Mash. It is recommended to read
         the Mash instructions when changing this parameter. Default = 5000
    -k   Specify the k-mer size used by Mash. It is recommended to read the
         Mash instructions when changing this parameter. Default = 16
    -ns  Number of FASTA file used for mash distance calculation at a time, 
         otherwise the input FASTA file will be calculate seperately and then
         be combined, default = 20,000
    -p   threads used
    -g   Whether deduplicate genomes [default: "1", True], use "0" to disable
______________________________________________________________________
''')

    parser.add_argument("-D", "--indir", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-tg", "--threshold_GC", help=argparse.SUPPRESS,
                        required=False, default=0.01, type=float)
    parser.add_argument("-k", "--kmer", help=argparse.SUPPRESS,
                        type=int, default=16, required=False)
    parser.add_argument("-s", "--sketch", help=argparse.SUPPRESS,
                        type=int, default=5000, required=False)
    parser.add_argument("-p", "--threads", help=argparse.SUPPRESS,
                        type=int, required=False, default=20)
    parser.add_argument("-ns", "--numberSplit",
                        type=int, default=20000, help=argparse.SUPPRESS, required=False)
    parser.add_argument("-g", "--genome",
                        type=int, default=1, choices=[0,1], help=argparse.SUPPRESS, required=False)
    return (parser.parse_args())


######################################################################
# similarity between gene clusters using MASH
######################################################################
def make_sketch(indir, outdir, kmer, sketch, numberSplit, threads):
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    Parameters
    ----------
    outdir
        string, the path to output directory
    option
        string, either 'GC' for the gene clusters or 'HG' for the housekeeping genes
    returns
    ----------
    """
    outlogfile = os.path.join(outdir, 'log.file')
    with open(outlogfile, "wb") as log_file:
        try:
            # compute how many BGCs
            all_file = os.listdir(indir)
            # use the most frequent suffix
            suffix_list = [Path(x).suffix for x in all_file]
            most_frequent_suffix = max(set(suffix_list), key = suffix_list.count)
            all_file_GC_seq = [os.path.join(indir, i) for i in all_file]
            
            # set the maximal number of BGCs to 20,000, otherwise split those into smaller parts for mash
            if len(all_file_GC_seq) < numberSplit:
                outfile = os.path.join(outdir, 'mash_sketch')
                inp = os.path.join(indir, "*"+most_frequent_suffix)
                cmd_mash = f"mash sketch -o {outfile} -k {kmer} -p {threads} -s {sketch} -a {inp}"
                p = Popen(cmd_mash, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                log_file.write(stderr)
            else:
                # split all input files in pieces
                all_file_GC_seq_partial = [all_file_GC_seq[i:i+numberSplit] for i in range(0, len(all_file_GC_seq), numberSplit)]
                for index, values in list(enumerate(all_file_GC_seq_partial)):
                    out_file = 'mash_sketch_' + str(index)
                    out_folder = os.path.join(outdir, "mash_sketch_folder_" + str(index))
                    outfile = os.path.join(outdir, out_file)
                    # create folder to run mash sketch more conveniently
                    os.mkdir(out_folder)
                    # move fasta to created folder
                    for value in values:
                        shutil.move(value, out_folder)
                    inp = os.path.join(out_folder, "*"+most_frequent_suffix)
                    cmd_mash = f"mash sketch -o {outfile} -k {kmer} -p {threads} -s {sketch} -a {inp}"
                    p = Popen(cmd_mash, shell=True, stdout=PIPE, stderr=PIPE)
                    stdout, stderr = p.communicate()
                    log_file.write(stderr)

        except(subprocess.CalledProcessError):
            # Raise error here for error table
            pass


def calculate_distance(outdir, threads, output_file="mash_output_GC.tab"):
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    Parameters
    ----------
    outdir
        string, the path to output directory
    returns
    ----------
    """

    # read all mash_sketch file
    all_file = os.listdir(outdir)
    all_mash_file = [os.path.join(outdir, i) for i in all_file if i.endswith(".msh")]
    if len(all_mash_file) == 1:
        infile = all_mash_file[0]
        outfile = os.path.join(outdir,output_file)
        try:
            cmd_mash = f"mash dist -p {threads} {infile} {infile} > {outfile}"
            res_download = subprocess.check_output(cmd_mash, shell=True)

        except(subprocess.CalledProcessError):
            # Raise error here for error table
            pass
    else:
        output_file_temp = os.path.join(outdir, "mash_output_GC_tmp")
        output_file_final = os.path.join(outdir, output_file)
        output_list = []
        input_mash_path = all_mash_file
        # combine two differen mash files
        pairwise_mash = list(itertools.permutations(input_mash_path,2))
        # add self-comparison mash files
        pairwise_mash_all = pairwise_mash + [(i, i) for i in input_mash_path]        

        for index, mash in list(enumerate(pairwise_mash_all)):
            out_file_name = "mash_output_GC_" + str(index)
            outfile = os.path.join(outdir,out_file_name)
            output_list.append(outfile)
            try:
                cmd_mash = f"mash dist -p {threads} {mash[0]} {mash[1]} > {outfile}"
                res_download = subprocess.check_output(cmd_mash, shell=True)
            except(subprocess.CalledProcessError):
                # Raise error here for error table
                pass


        mash_output_list = os.path.join(outdir, "mash_output_GC_*")
        #out_mash = f"cat {mash_output_list} > {output_file_temp}"
        fo_final = open(output_file_final, 'w')
        out_mash = f"cat {mash_output_list}"
        #os.system(out_mash)
        #subprocess.run(out_mash, shell=True)
        p = subprocess.Popen(out_mash, shell=True, stdout=fo_final, stderr=PIPE)
        p.wait()
        fo_final.close()
        #cmd_sort = f"sort -k2 {output_file_temp} -S 80% --parallel=40 -o {output_file_final}"
        #os.system(cmd_sort)
        #p2 = subprocess.Popen(cmd_sort, shell=True, stdout=PIPE, stderr=PIPE)
        #p2.wait()
        cmd_remove = f'rm {os.path.join(outdir,"mash_output_GC_*")}'
        subprocess.run(cmd_remove, shell=True)

    return ()


######################################################################
# calculate N50
######################################################################
def calculate_N50(inFASTA):
    """Calculate N50 for a FASTA file
     Returns:
        float: N50 value.
    """
    fasta = open(inFASTA, 'r')
    reader = fastaparser.Reader(fasta, parse_method='quick')
    dict_reader = dict(reader)
    # 1. list of contig sizes
    list_of_lengths = [len(x) for x in dict_reader.values()]
    contig_sizes = numpy.array(list_of_lengths)
    # 2. sort integer array by size
    contig_sizes.sort()
    # 3. reverse the order of the array
    contig_sizes = contig_sizes[::-1]
    # 4. get the total sum (e.g., assembly size)
    total_len = contig_sizes.sum()
    # 5. get half of the total size
    half_total_len = total_len / 2
    # 6. make an array of zeros
    contig_sum_lens = numpy.zeros(contig_sizes.size, dtype=int)
    # 7. interate over the items in the empty array and
    #    fill each index with the sum of all up to that point. 
    for i in range(contig_sizes.size):
        contig_sum_lens[i] = contig_sizes[i:].sum()
    contig_sum_lens
    # 8. calculate which values are > half the total length
    which_contigs_longer_than_half = contig_sum_lens > half_total_len
    # 9. keep the values that answered True above
    contigs_longer_than_half = contig_sizes[which_contigs_longer_than_half]
    # 10. keep the lowest value from above, which is the N50
    n50 = contigs_longer_than_half.min()
    
    return n50


######################################################################
# Extracting clusters using cut-off & calculate medoid or cut-off & N50 (for genome)
######################################################################
def compute_matrix_sum(geneList, pairwiseDistance):
    dic_distance_min = ""
    representative_seq = ""
    for x in geneList:
        geneList_cp = geneList.copy()
        geneList_cp.remove(x)
        distance_all_x = 0

        for y in geneList_cp:
            gene_cis = x + y
            gene_trans = y + x
            try:
                distance_all_x += pairwiseDistance[gene_cis]
            except:
                distance_all_x += pairwiseDistance[gene_trans]
        if not dic_distance_min:
            dic_distance_min = distance_all_x
            representative_seq = x
        else:
            if dic_distance_min <= distance_all_x:
                continue
            else:
                dic_distance_min = distance_all_x
                representative_seq = x

    return representative_seq


def representative_N50(genomeList):
    # input a list of genomes, return the representative one with the highest N50
    n50_list = [calculate_N50(x) for x in genomeList]
    min_N50_index = numpy.argmax(numpy.array(n50_list))
    representative_genome = genomeList[min_N50_index]
    return representative_genome


def calculate_medoid(outdir, cut_off, med={}, input_file="mash_output_GC.tab", genomePattern=1):
    """
    calculates the GCFs based on similarity threshold
    parameters and calculates the medoid of that GCF
    ----------
    outdir
        string, the path to output directory
    cut_off
        float, between 0 and 1
    returns
    ----------
    dict_medoids = {fasta file of medoid: similar fasta files}
    """
    # Parse the input into a dictionary of gene families
    family_by_gene = {}  # store the pairwise gene cluster info, eg. [geneA:cluterA, geneB:clusterB, geneC:clusterC]
    family_by_gene_filtered = {} 
    family_members = {} # store the gene cluster, eg. [clusterA:[geneA,geneB,geneC...]]
    family_distance_matrices = {} 
    dict_medoids = med
    infile = os.path.join(outdir, input_file)


    cluster_genos = dict() # {cluster1:[geneA,geneB,geneC],cluster2:[geneD,geneE,geneF]...}
    cluster_genos_final = dict()
    all_genome_list = []
    cluster_list = dict() # {geneA:geneB, geneC:geneB}, corresponding cluster
    pairwise_distance = dict()

    with open(infile, 'r') as input_f:
        for line in input_f:
            # Skip lines starting with '#'
            if line.startswith('#'):
                continue
            gene1, gene2, distance, number2, overlap = line.strip().split('\t')
            # no1, no2 = overlap.split("/")
            # overlap = float(no1)/float(no2)

            # store distance of a pair of genes
            combine_genes_cis = gene1 + gene2
            combine_genes_trans = gene2 + gene1

            if (combine_genes_cis in pairwise_distance) or (combine_genes_trans in pairwise_distance):
                pass
            else:
                pairwise_distance[combine_genes_cis] = float(distance) # only store cis

            if float(distance) <= cut_off:
                # Case1: gene1 and gene2 are not assigned into cluster, select gene1 as cluster name
                if (gene1 not in cluster_list) and (gene2 not in cluster_list):
                    cluster_name = gene1
                    cluster_list[gene1] = cluster_name
                    cluster_list[gene2] = cluster_name
                    cluster_genos[cluster_name] = [gene1, gene2]
                # Case2: gene1 has been assigned into a cluster but gene2 not, add gene2 to gene1 cluster
                elif  (gene1 in cluster_list) and (gene2 not in cluster_list):
                    cluster_name = cluster_list[gene1]
                    cluster_list[gene2] = cluster_name
                    cluster_genos[cluster_name].append(gene2)

                # Case3: gene2 has been assigned into a cluster but gene1 not, add gene1 to gene2 cluster
                elif  (gene1 not in cluster_list) and (gene2 in cluster_list):
                    cluster_name = cluster_list[gene2]
                    cluster_list[gene1] = cluster_name
                    cluster_genos[cluster_name].append(gene1)

                # Case4: both gene1 and gene2 have been assigned to specific cluster
                elif (gene1 in cluster_list) and (gene2 in cluster_list):
                    cluster_name_gene1 = cluster_list[gene1]
                    cluster_name_gene2 = cluster_list[gene2]
                    # Case4.1: gene1 and gene2 belong to same cluster
                    if cluster_name_gene1 == cluster_name_gene2:
                        pass
                    # Case4.2: gene1 and gene2 belong to different cluster, merge two cluster
                    # only keep the gene1 cluster, remove gene2 cluster, change all items belonging to gene2 cluster name to gene1 cluster
                    else:
                        cluster_genos[cluster_name_gene1].extend(cluster_genos[cluster_name_gene2])
                        for gene in cluster_genos[cluster_name_gene2]:
                            cluster_list[gene] = cluster_name_gene1 # include gene2
                        del cluster_genos[cluster_name_gene2]

    if genomePattern != 1:
        # fish out the representative sequence using the minimal sum distance
        for k,v in cluster_genos.items():
            representative_seq = compute_matrix_sum(v, pairwise_distance)
            cluster_genos_final[representative_seq] = list(set(v))
    else:
        # fish out the representative sequence using the highest N50
        for k,v in cluster_genos.items():
            representative_seq = representative_N50(v)
            cluster_genos_final[representative_seq] = list(set(v))

    return cluster_genos_final



def writeGCFfasta(sim_dict, outdir, outfile):
    """Writes the GCFs reference fasta file and adds n_repr
    parameters
    ----------
    sim_dict
        dict, similarity dictionary for GCs
    outdir
        string, the path to output directory
    outfile
        string, name of the outfile
    returns
    ----------
    outfile = name of the GCF fasta file in outdir
    """
    infiles = sim_dict.keys()
    outfile = os.path.join(outdir, outfile)
    with open(outfile, "w") as fout:
        for fkey in infiles:
            n_repr = len(sim_dict[fkey])
            fkey_contents = fkey.split("/")
            fname = "/".join(fkey_contents)
            with open(fname, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        line = line.strip()
                        if n_repr < 1:
                            fout.write(f"{line}--NR=1\n")
                        else:
                            fout.write(f"{line}--NR={n_repr}\n")
                    else:
                        fout.write(line)
    return (outfile)


def makefastaheadersim(sim_dict):
    """converts the file sim_dict to a fastaheader a similarity dictionary
    parameters
    ----------
    sim_dict
        dict, similarity dictionary for GCs
    returns
    ----------
    ret = {fastaheader: [fastaheaders that are similar]}
    """
    ret = {}
    infiles = sim_dict.keys()
    for fname in infiles:
        sim_fnames = sim_dict[fname]
        n_repr = len(sim_fnames)
        with open(fname, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    fastaheader = f"{line[1:]}--NR={1 if n_repr < 1 else n_repr}"  # stripping '>'
                    ret[fastaheader] = []
        for sf in sim_fnames:
            fkey_contents = sf.split("/")
            with open(sf, "r") as s:
                for line in s:
                    line = line.strip()
                    if line.startswith(">"):
                        sim_fastaheader = line[1:]  # stripping '>'
                        ret[fastaheader].append(sim_fastaheader)
    return (ret)


def writejson(dictionary, outdir, outfile_name):
    """writes results in a dict to json format
    parameters
    ----------
    dictionary
        dict, dicionary containing some results (here GCFs results)
    outdir
        string, the path to output directory
    outfile_name
        string, name for the outfile
    returns
    ----------
    outfile
        name and path of the output file
    """

    outfile = os.path.join(outdir, outfile_name + ".json")
    with open(outfile, "w") as w:
        w.write(json.dumps(dictionary, indent=4))
    return(outfile)



######################################################################
# Writing and purging files
######################################################################
def purge(d, pattern):
    """removes files matching a pattern
    parameters
    ----------
    d
        string, directory path
    pattern
        string, regex
    returns
    ----------
    """
    for f in os.listdir(d):
        if re.search(pattern, f):
            os.remove(os.path.join(d, f))

def purge_folder(d, pattern):
    # remove folder
    for dirpath, dirnames, filenames in os.walk(d):
        if not dirnames and re.search(pattern, dirpath):
            shutil.rmtree(dirpath)


######################################################################
# MAIN
######################################################################
def main():
    warnings.simplefilter('ignore', ResourceWarning)
    args = get_arguments()

    try:
        os.mkdir(args.outdir)
    except:
        pass

    ################################
    # Mash: similarity
    ################################
    print("01. Mash Sketch", flush=True)

    make_sketch(args.indir, args.outdir + os.sep, args.kmer, args.sketch, args.numberSplit, args.threads)

    #checks the output of the mash sketch
    reruns = 0
    total_reruns = 25
    for i in range(total_reruns):
        logfile_name = os.path.join(args.outdir, 'log.file')
        with open (logfile_name, "r") as log_file:
            if "ERROR" in log_file.read():
                make_sketch(args.indir, args.outdir + os.sep, args.kmer, args.sketch, args.numberSplit, args.threads)
                reruns += 1
                print("Encountered error in making sketch file. Retry attempt:", reruns)
                if reruns == total_reruns:
                    sys.exit("Maximum number of reruns is reached. Try decreasing the number of cores (-p flag).")
            else:
                break

    print("02. Mash Distance estimination", flush=True)

    calculate_distance(args.outdir + os.sep, args.threads)

    print("03. Calculate medoid", flush=True)

    GCFs = calculate_medoid(args.outdir + os.sep, args.threshold_GC, genomePattern=args.genome)

    q = os.path.join(args.outdir, "nonredundant_sequence")
    try:
        os.mkdir(q)
    except:
        pass
    de_file = os.path.join(args.outdir, "nonredundant_sequence.list")
    
    with open(de_file, "w") as fin:
        print(*["Representative_sequence", "Sequences"], sep="\t", file=fin)
        for k,v in GCFs.items():  
            print(*[k,v], sep="\t", file=fin)
            
            ln_raw = os.path.abspath(k)
            ln_file = os.path.join(q, Path(k).name)

            cmd_ln = f'ln -s {ln_raw} {ln_file}'
            subprocess.run(cmd_ln, shell=True)
            #os.system(f'ln -s {k} os.path.join({q}, Path({k}).name)')

    print("04. Write out clusters", flush=True)

    fastadict = makefastaheadersim(GCFs)


    # Writing results to outdir
    writejson(fastadict, args.outdir, "BiG-MAP.GCs")
    #fasta_file = writeGCFfasta(GCFs, args.outdir, "BiG-MAP.GCF.fna")


    #############################
    # Cleaning output dir
    #############################
    purge(args.outdir, ".fasta")
    purge(args.outdir, ".txt")
    purge(args.outdir, ".faa")
    #purge_folder(args.outdir, "mash_sketch_")
    #purge(args.outdir, "log.file")


if __name__ == "__main__":
    main()
