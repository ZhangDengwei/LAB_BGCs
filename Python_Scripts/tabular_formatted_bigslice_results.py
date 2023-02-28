# -*- coding: utf-8 -*-
# @Author: cq
# @Date:   2021-11-23 14:00:12
# @Last Modified by:   cq
# @Last Modified time: 2021-11-23 16:11:06


# Import statements
import os
import sys
import argparse
import re
import pandas as pd



def findResults(inputPath):
	# inputPath: results generated by script 'convert_bigslice_data_2_table.py'
	files = os.listdir(inputPath)
	files_csv = [x for x in files if x.endswith(".csv")]

	print("find {0} csv files, as follows: ".format(len(files_csv)))
	for x in files_csv: print(x)

	# find several core files
	core_files = ["bgc.csv",
	              "gcf_membership.csv"]
	for x in core_files:
		if x not in files_csv:
			print("Did not find file {0}, please check...".format(x))
			sys.exit() 

def funcMember(value):
	if value <= 300:
		return "core"
	else:
		return "putative"


def gcfNumber(dataFrame,colName1, colName2):
	'''
	dataFrame: Pandas dataframe
	colName1: col name of gcf membership
	colName2: col name of membership_value [<=300: core, >300: putative]
	'''
	dict_number = {}

	gcf_cluster = set(dataFrame[colName1])

	for gcf in gcf_cluster:
		member_number = dataFrame[dataFrame[colName1]==gcf].shape[0]
		core_number = dataFrame[(dataFrame[colName1]==gcf) & (dataFrame[colName2]<=300)].shape[0]
		putative_number = dataFrame[(dataFrame[colName1]==gcf) & (dataFrame[colName2]>300)].shape[0]

		dict_number[gcf] = [member_number, core_number, putative_number]

	return dict_number


def tabular(inputPath, outFile):
	bgcFile = os.path.join(inputPath, 'bgc.csv')
	gcfFile = os.path.join(inputPath, 'gcf_membership.csv')

	df_bgc = pd.read_csv(bgcFile, header=0, index_col=0, sep=',')
	df_gcf = pd.read_csv(gcfFile, header=0, index_col=0, sep=',')

	df_merge = pd.merge(df_bgc, df_gcf, left_on="id", right_on="bgc_id")

	dict_number = gcfNumber(df_merge, 'gcf_id', 'membership_value')

	# add members type
	df_merge['membership_types'] = df_merge.apply(lambda x: funcMember(x.membership_value), axis=1)
	df_merge['membership_number'] = df_merge.apply(lambda x: dict_number[x.gcf_id][0], axis=1)
	df_merge['membership_core_number'] = df_merge.apply(lambda x: dict_number[x.gcf_id][1], axis=1)
	df_merge['membership_putative_number'] = df_merge.apply(lambda x: dict_number[x.gcf_id][2], axis=1)

	# output
	df_merge.to_csv(outFile, sep="\t", index=False)


def main():
	parse = argparse.ArgumentParser(description="merge 'bgc.csv' and 'gcf_membership.csv'")

	parse.add_argument("--inputPath", help="Folder containing bigslice results with csv format", required=True)
	parse.add_argument("--outFile", help="output file", required=True)

	args = parse.parse_args()

	# processing
	findResults(args.inputPath)

	tabular(args.inputPath, args.outFile)

if __name__ == "__main__": 
	main()