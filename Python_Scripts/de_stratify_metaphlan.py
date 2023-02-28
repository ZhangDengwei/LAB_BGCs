#!/usr/bin/env python3

import os
import sys
from sys import argv

script, file1 = argv


def classify():
    with open(file1, "r") as f1, open('merged_kingdom.txt', "w") as kf, open('merged_phylum.txt', "w") as pf, \
        open('merged_class.txt', "w") as cf, open('merged_order.txt', "w") as of, open('merged_family.txt', "w") as ff, \
            open('merged_genus.txt', "w") as gf, open('merged_species.txt', "w") as sf:

        kindom_d, phylum_d, class_d, order_d, family_d, genus_d, species_d = {}, {}, {}, {}, {}, {}, {}
        f1.readline()
        sample_list = f1.readline()
        
        kf.write(sample_list)
        pf.write(sample_list)
        cf.write(sample_list)
        of.write(sample_list)
        ff.write(sample_list)
        gf.write(sample_list)
        sf.write(sample_list)


        for line in f1:
            if line.split('\t')[0].split("|")[-1][0] == 'k':
                kf.write(line)
            elif line.split('\t')[0].split("|")[-1][0] == 'p':
                '''
                if line.split('\t')[0].split('_')[-1] == 'unclassified':
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1], ' '] + line.rstrip('\n').split()[1:]) + '\n'
                else:
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1]] + line.rstrip('\n').split()[1:]) + '\n'
                    '''
                line_cont = line.split("\t")[0].split("|")[-1] + "\t" + "\t".join(line.rstrip('\n').split("\t")[1:]) + '\n'
                pf.write(line_cont)
            elif line.split('\t')[0].split("|")[-1][0] == 'c':
                '''
                if line.split('\t')[0].split('_')[-1] == 'unclassified':
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1], ' '] + line.rstrip('\n').split()[1:]) + '\n'
                else:
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1]] + line.rstrip('\n').split()[1:]) + '\n'
                    '''
                line_cont = line.split("\t")[0].split("|")[-1] + "\t" + "\t".join(line.rstrip('\n').split("\t")[1:]) + '\n'
                cf.write(line_cont)
            elif line.split('\t')[0].split("|")[-1][0] == 'o':
                '''
                if line.split('\t')[0].split('_')[-1] == 'unclassified':
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1], ' '] + line.rstrip('\n').split()[1:]) + '\n'
                else:
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1]] + line.rstrip('\n').split()[1:]) + '\n'
                    '''
                line_cont = line.split("\t")[0].split("|")[-1] + "\t" + "\t".join(line.rstrip('\n').split("\t")[1:]) + '\n'
                of.write(line_cont)
            elif line.split('\t')[0].split("|")[-1][0] == 'f':
                '''
                if line.split('\t')[0].split('_')[-1] == 'unclassified':
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1], ' '] + line.rstrip('\n').split()[1:]) + '\n'
                else:
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1]] + line.rstrip('\n').split()[1:]) + '\n'
                '''
                line_cont = line.split("\t")[0].split("|")[-1] + "\t" + "\t".join(line.rstrip('\n').split("\t")[1:]) + '\n'
                ff.write(line_cont)
            elif line.split('\t')[0].split("|")[-1][0] == 'g':
                '''
                if line.split('\t')[0].split('_')[-1] == 'unclassified':
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1], ' '] + line.rstrip('\n').split()[1:]) + '\n'
                else:
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1]] + line.rstrip('\n').split()[1:]) + '\n'
                '''
                line_cont = line.split("\t")[0].split("|")[-1] + "\t" + "\t".join(line.rstrip('\n').split("\t")[1:]) + '\n'
                gf.write(line_cont)
            elif line.split('\t')[0].split("|")[-1][0] == 's':
                '''
                if line.split('\t')[0].split('_')[-1] == 'unclassified':
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1], ' '] + line.rstrip('\n').split()[1:]) + '\n'
                else:
                    line_cont = "\t".join([line.split("\t")[0].split("|")[-1]] + line.rstrip('\n').split()[1:]) + '\n'
                '''
                line_cont = line.split("\t")[0].split("|")[-1] + "\t" + "\t".join(line.rstrip('\n').split("\t")[1:]) + '\n'
                sf.write(line_cont)

classify()


