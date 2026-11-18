#!/usr/bin/env python

import sys
import random

out_table = sys.argv[1]
out_variant_truth = sys.argv[2]
out_cell_truth = sys.argv[3]

num_variants = 1000
num_cells = 1000
average_ASE_per_cell = 10
Xi = 0.0

variant_truth = []
cell_truth = []
for i in range(0, num_cells):
	cell_truth.append(random.uniform(0, 1) < 0.5)
for i in range(0, num_variants):
	variant_truth.append(random.uniform(0, 1) < 0.5)

with open(out_table, "w") as f:
	for cell in range(0, num_cells):
		ASE = random.randint(0, 2*average_ASE_per_cell)
		for i in range(0, ASE):
			variant = random.randint(0, num_variants-1)
			allele = "REF"
			if cell_truth[cell] == variant_truth[variant]:
				if random.uniform(0, 1) < Xi:
					allele = "ALT"
				else:
					allele = "REF"
			else:
				if random.uniform(0, 1) < Xi:
					allele = "REF"
				else:
					allele = "ALT"
			print("cell" + str(cell) + "\t" + "X:" + str(variant) + ":A:T" + "\t" + allele + "\t" + "1", file=f)

with open(out_variant_truth, "w") as f:
	for variant in range(0, num_variants):
		print("X:" + str(variant) + ":A:T" + "\t" + ("mat" if variant_truth[variant] else "pat"), file=f)

with open(out_cell_truth, "w") as f:
	for cell in range(0, num_cells):
		print("cell" + str(cell) + "\t" + ("mat" if cell_truth[cell] else "pat"), file=f)
