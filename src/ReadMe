This program generates a genome graph
1. To compile the program:
	i. Go to the src directory
	ii.make
2. Change the library path to the 'src' directory
	export LD_LIBRARY_PATH= [...]/src

3. Input Fasta file should be changed to the foramt is acceptable by the program run:
	../src/graph fasta_prepare prepared_fasta input1_fasta input2_fasta [...] inputn_fasta
4. Use the prepared_fasta to generate a maf file contains pairwise alignments
	Any aligner can be used 
5. After having the input_maf file, to built the graph run:
	../src/graph build_graph prepared_fasta input_maf output_maf output_encode(noencode) output_fasta ourput_dot output_gfa(nogfa) output_txt 0(1) number_of_threads(optional)
	
	0: generate graphs with longer nodes
	1: generate graphs with shorter nodes
	noencode: if don't need the arithmetically encoded output file
	nogfa: if dont' need the gfa output file


6. To decode the output_encode 
	../src/graph dy_decoding output_encode 0(1) output_decode
	
	0 or 1 determine the center type which was used in buildign the graph.	
