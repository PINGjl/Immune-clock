##source activate cpdb

import numpy as np
import pandas as pd
# Update on Aug.15, 2024
import sys
import random
import cellphonedb

random.seed(10)
samples = open("/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/07_celldb/sample.txt","r")
lines = samples.readlines()
samples.close()

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
cpdb_file = '/data/pingjiale/02_Software/01_anaconda/cellphonedb.zip'

for i in range(len(lines)):
	sample = lines[i].strip("\n")
	meta_file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/07_celldb/"+sample+"/cellphonedb_meta.txt"
	counts_file = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/07_celldb/"+sample+"/cellphonedb_count.txt"
	output = "/data/pingjiale/05_Result/01_qz_pbmc/02_scRNA/07_celldb/"+sample+"/"
	
	cpdb_results = cpdb_statistical_analysis_method.call(
		cpdb_file_path = cpdb_file,
		meta_file_path = meta_file,
		counts_file_path = counts_file,
		counts_data = 'gene_name',
		score_interactions = True,
		threshold = 0.1,
		output_path = output)
