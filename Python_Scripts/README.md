- remove_duplicated_seq.py

> Deduplicating the bacterial genomes or biosynthetic gene clusters (BGC)

- convert_bigslice_data_2_table.py

> Extracting all tables from sqlite3 database that was generated by BiG-SLiCE. This will creat a `result` folder under the input BiG-SLiCE result folder.

- tabular_formatted_bigslice_results.py

> Combining the results from BiG-SLiCE, fishing out BGCs and corresponding families (GCF) defined by BiG-SCLiCE.  

- bigslice_hmm_statistics.py

> Tabulating the presence/absence matrix of BGC features (biosynthetic-Pfam and sub-Pfam domains) identified by BiG-SCLiCE.

- bgcs_statistics.py

> Summarizing the BGCs detected by antiSMASH 6

- bgc_number.py

> Summarizing the BGC number, for each genome, detected by antiSMASH 6

- clustering_bgc_from_bigsclice.pool_3.py

> Clustering BGCs into families and clans, according to all-to-all cosine distances between BGCs

- featureCount_BGCs_statistics_2.py

> Summarizing the BGC profiles detected in metagenomes/metatransciptome.
> 
> This will generate six outputs:
> 
> - XXX___BGCs_present.tsv, matrix showing presence/absence of BGCs
> 
> - XXX___BGCs_presence_biosynthetic_genes.tsv, matrix showing the presence of the biosynthetic genes of BGCs
> 
> - XXX___BGCs_counts_considering_all_biosynthetic_genes.tsv, matrix showing the read counts of BGCs when condsidering all biosynthetic genes for each BGC
> 
> - XXX___BGCs_TPM_considering_all_biosynthetic_genes.tsv, matrix showing BGC abundance (TPM) when condsidering all biosynthetic genes for each BGC
> 
> - XXX___BGCs_counts_considering_core_genes.tsv, matrix showing read counts of BGCs when only considering core genes
> 
> - XXX___BGCs_TPM_considering_core_genes.tsv, matrix showing BGC abundance (TPM) when only condsidering core genes for each BGC

- de_stratify_metaphlan.py

> Stratifying the metaphlan-formatted abundance table according to different taxonomic level

- hmmsearch.precursor.py

> Fishing out the putative precursors of class II bacteriocins identified by `hmmsearch`

- combine_two_results.py

> Combining the outputs of class II bacteriocin identification using `hmmsearch` and `BAGEL4`

- remodel_cdhit.py

> Summarizing the output of `cd-hit`, generating two outputs
> 
> - XXX__overall.tsv, one cluster per row
> 
> - XXX__one_vs_one.tsv, one query per row

- parase_needleall.py

> Summarizing the output of `needleall`
