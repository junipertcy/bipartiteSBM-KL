#! bin/sh

# test 1
./biSBM southernWomen.edgelist southernWomen.types southernWomen 2 3 1 3

echo "\n\n"

# test 2
./biSBM test.edgelist test.types test 2 3 1 3

echo "\n\n"

# Should produce output:
# ***** Communities by type/part *****
# 0,1,
# 2,3,4,
# Using biSBM with degree correction.
# edges: 	southernWomen.edgelist
# types: 	southernWomen.types
# K: 	2,3
# KL steps	3
# output: 	southernWomen
# TYPE 0 : 18
# TYPE 1 : 14
# NODES: 32
# Read in twice edges = 178
# >1,0.000000,-367.657617
# >2,0.000000,-367.657617
# >3,0.000000,-367.657617
# Final Score: -367.658
# 
# 
# 
# ***** Communities by type/part *****
# 0,1,
# 2,3,4,
# Using biSBM with degree correction.
# edges: 	test.edgelist
# types: 	test.types
# K: 	2,3
# KL steps	3
# output: 	test
# TYPE 0 : 700
# TYPE 1 : 300
# NODES: 1000
# Read in twice edges = 15392
# >1,0.000000,-65525.804037
# >2,1.000000,-65525.804037
# >3,1.000000,-65525.804037
# Final Score: -65525.8
