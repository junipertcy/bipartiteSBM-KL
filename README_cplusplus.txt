README for biSBM v1.2

Daniel Larremore
Harvard School of Public Health
July 29, 2014
http://danlarremore.com/bipartiteSBM
daniel.larremore@gmail.com

biSBM - a method for community detection in bipartite networks, based on the publication: 
	"Efficiently inferring community structure in bipartite networks"
	Daniel B. Larremore, Aaron Clauset, and Abigail Z. Jacobs 
	Physical Review E 90(1), 012805 (2014).
	http://danlarremore.com/pdf/2014_LCJ_EfficientlyInferringCommunityStructureInBipartiteNetworks_PRE.pdf

/***** MATLAB vs R vs C++ *****/

This file explains how to compile the pure C++ version of the biSBM code. If you prefer using a wrapper in MATLAB or R, see README_MATLAB.txt or README_R.txt. 

/***** Usage *****/

You'll need the files:
	biSBM.cpp
	biSBM.h
	-----
	example_cplusplus.sh
	southernWomen.edgelist
	southernWomen.types
	test.edgelist
	test.types
TO COMPILE:
	g++ -O3 -Wall -g -pedantic -o biSBM biSBM.cpp

TO TEST:
	sh example_cplusplus.sh

EXAMPLE USAGE:
	open example_cplusplus.sh

OPTIONS:
	./biSBM
	with no options will print usage.
	./biSBM <edgeList> <vertexTypes> <outputFOLDER> <Ka>  <Kb> <isDegreeCorrect> <KL steps> <optionalTrueComms> <initialize@truecomms?>

	REQUIRED:

	<edgeList> - file containing edges in format "i \t j \n" for an edge between vertices i and j.
	<vertexTypes> - a file with a 1 or 2 on each line to indicate the type of the corresponding vertex.
	<outputFOLDER> - name of the folder to be created for the output files
	<Ka>, <Kb> - The number of communities for type a (1) and type b (2) vertices.
	<isDegreeCorrect> - boolean degree correction? 0 = no, 1 = yes.
	<KL steps> - The number of random initializations
	
	OPTIONAL:
	
	<optionalTrueComms> - name of file containing true communities.
	<initialize@truecomms?> - boolean. 0 = no, 1 = yes.