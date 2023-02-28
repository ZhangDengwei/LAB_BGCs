from multiprocessing.pool import ThreadPool
from scipy.spatial import distance
from sklearn.cluster import AgglomerativeClustering
import numpy as np
from pandas import DataFrame, Series
import argparse
import collections
from multiprocessing import Pool
import time
import multiprocessing


def cldist(parameters):
	item = parameters[0]
	dict_domain = parameters[1]
	dict_dist = dict()
	dict_dist[item] = dict()
	for y in dict_domain.keys():
		dist = distance.cosine(dict_domain[item], dict_domain[y])
		dict_dist[item][y] = dist
	return dict_dist


def finalrun(infile, prefix, max_bgcid):
	global dict_domain, dict_dist
	dict_domain = dict()
	dict_dist = dict()
	dict_dist_out = dict()

	with open(infile, "r") as fin:
		fin.readline()
		for line in fin:
			items = line.rstrip("\n").split("\t")
			if (max_bgcid == 0) or (int(items[0]) <= max_bgcid):
				dict_domain[items[0]] = [int(i) for i in items[1:]]
				
	bgcs = dict_domain.keys()
	process = locals()

	p = Pool()
	z = []
	for x in bgcs:
		z.append([x, dict_domain])
	results = p.map(cldist, z)
	p.close()
	p.join()

	for x in results:
		dict_dist.update(x)
	'''
	# sort 
	for k, v in dict_dist.items():
		v_tmp = collections.OrderedDict(sorted(v.items()))
		dict_dist[k] = v_tmp
	dict_dist_sort = collections.OrderedDict(sorted(dict_dist.items()))
	'''

	data = DataFrame(dict_dist)
	out_dist = prefix + "_pairwise.distance.tsv"
	data.to_csv(out_dist, sep="\t", index=True)

	cluster_2 = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage="average", distance_threshold=0.2)
	cluster_3 = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage="average", distance_threshold=0.3)
	cluster_4 = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage="average", distance_threshold=0.4)
	cluster_5 = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage="average", distance_threshold=0.5)
	cluster_6 = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage="average", distance_threshold=0.6)
	cluster_7 = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage="average", distance_threshold=0.7)
	cluster_8 = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage="average", distance_threshold=0.8)
	
	cluster_2.fit_predict(data)
	cluster_3.fit_predict(data)
	cluster_4.fit_predict(data)
	cluster_5.fit_predict(data)
	cluster_6.fit_predict(data)
	cluster_7.fit_predict(data)
	cluster_8.fit_predict(data)

	cluster_frame = DataFrame()
	cluster_frame["BGC"] = Series(data.index)
	cluster_frame["cluster_0.2"] = Series(cluster_2.labels_)
	cluster_frame["cluster_0.3"] = Series(cluster_3.labels_)
	cluster_frame["cluster_0.4"] = Series(cluster_4.labels_)
	cluster_frame["cluster_0.5"] = Series(cluster_5.labels_)
	cluster_frame["cluster_0.6"] = Series(cluster_6.labels_)
	cluster_frame["cluster_0.7"] = Series(cluster_7.labels_)
	cluster_frame["cluster_0.8"] = Series(cluster_8.labels_)
	
	out_cluster = prefix + "clustering.tsv"
	cluster_frame.to_csv(out_cluster, sep="\t", index=False)


def main():
	parse = argparse.ArgumentParser(description="clustering BGCs into family based on cosine distance")

	parse.add_argument("--input", help="BGC features extracted from bigslice", required=True)
	parse.add_argument("--prefix", help="prefix of output file", required=True)
	parse.add_argument("--maxBGCid", help="the maximum BGC id assigned by bigsclice, those BGCs (<= maxBGCid) would be \
											included for clustering, [default=0, including all BGCs]", required=False, default=0, type=int)

	args = parse.parse_args()

	start = time.time()
	finalrun(args.input, args.prefix, args.maxBGCid)
	end = time.time()
	print(end - start)


if __name__ == "__main__": 
    main()

