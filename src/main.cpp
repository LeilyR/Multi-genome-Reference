#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <ostream>
#include <vector>
#include "pw_alignment.hpp"
#include "data.hpp"
#include "model.hpp"
#include "encoder.hpp"
#include "test.hpp"
#include "graph.hpp"
#define NO_MAKEFILE
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"
#undef NO_MAKEFILE



#define VERSION "0.0.1" 
#define NUMVERISON 0
#define MIN_READ_VERSION 0

#define TIMING 1


#define TEST 1


#if TEST

int main(int argc, char * argv[]) {
	std::cout << " hello " << std::endl;
	Graph g = Graph();


	//  runNeedleman("CCCGGGGGGTGCA","ATAGTTGCA",2);
	//load data
	all_data data = all_data();

	//Init new Graph //TODO in main or parseData ?
	Graph newGraph = Graph();
	g.parseData(data,newGraph);

	//std::cout << newGraph;
	//std::cout <<" New Graph " <<newGraph << std::endl;




/*
	pw_alignment p(std::string("ATT----TTCTT"), string("AGTGATAT----"), 12, 15, 23, 26,1,1);
	pw_alignment s(std::string("ATT----TTCTT"), string("ACTGATG---AC"),13, 18, 24, 29,2,1);

	for(size_t i=0; i<p.alignment_length(); ++i) {
		char s1;
		char s2;
		p.alignment_col(i, s1, s2);
		std::cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << std::endl;
	  
	}

	pw_alignment p1;
	pw_alignment p2;
	p.split(false,4, p1, p2);
		for(size_t i=0; i<p2.alignment_length(); ++i) {
		char s1;
		char s2;
		p2.alignment_col(i, s1, s2);
		std::cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << std::endl;
	 
	}
	overlap o;
	o.add_alignment(s);
	typedef multimap <size_t, pw_alignment> mmap;
	mmap alignment_map;
	alignment_map.insert(make_pair(12,p));
	alignment_map.insert(make_pair(23,p));
	alignment_map  alignments_on_reference1;
	mmap::const_iterator alignit = alignments_on_reference1.begin();
	for(; alignit != alignments_on_reference1.end(); ++alignit) {
		char s1;
		char s2;
		for(size_t i=0;i<p.alignment_length();++i){
		s.alignment_col(i, s1, s2);
		std::cout << "pos " << i<< " s1 " << s1 << " s2 " << s2 << std::endl;
}

}
*/
	return 0;

}



#endif
void usage() {
	cerr << "Graph program" << std::endl;
	cerr << "Usage:" << std::endl;
}

int do_fasta_prepare(int argc, char * argv[]) {

	if(argc < 4) {
		usage();
		cerr << "Program fasta_prepare" << std::endl; 
		cerr << "Combine fasta files from different species/accessions into a common file" << std::endl;
		cerr << "Parameters: " << std::endl;
		cerr << "output file -- fasta file" << std::endl;
		cerr << "<list of fasta files> -- input file, one file per species/accession, no colon in file name"<< std::endl;
		return 1;
	}


	std::ofstream outf(argv[2]);
	if(outf) {
		std::vector<std::ifstream*> inputs;
		std::vector<std::string> names;
		for(int i=3; i<argc; ++i) {
			std::string istr(argv[i]);
			std::vector<std::string> slashparts;
			strsep(istr, "/", slashparts);
			std::string afterslash = slashparts.at(slashparts.size()-1);
			std::vector<std::string> periodparts;
			strsep(afterslash, ".", periodparts);
			std::string name = periodparts.at(0);
			if(periodparts.size()>1) {
				std::stringstream sstr;
				for(size_t j=0; j<periodparts.size()-1; j++) {
					if(j!=0) sstr << ".";
					sstr << periodparts.at(j);
				}
				name = sstr.str();
			}
			for(size_t j=0; j<name.length(); ++j) {
				if(':' == name.at(j)) {
					cerr << "Error: " << argv[i] << " file name contains a colon" << std::endl;
					exit(1);
				}
			}
			std::ifstream * iin = new ifstream(argv[i]);
			if(!(*iin)) {
				cerr << "Error: cannot read: " << argv[i] << std::endl;
				exit(1);
			}
			inputs.push_back(iin);
			names.push_back(name);
		}

		for(size_t i=0; i<inputs.size(); ++i) {


			std::string str;
			while(getline((*(inputs.at(i))), str)) {
				if(str.length()>0) { // ignore empty lines
				if(str.at(0)=='>') {
					std::string nname = str.substr(1);
					std::stringstream nhead;
					nhead << ">" << names.at(i) << ":" << nname;
					outf << nhead.str() << std::endl;

						
						
/*
	for(size_t col = 0; col < alignment_length(); col++) {
		char c1;
		char c2;
		alignment_col(col, c1, c2);
		std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
	if(col>50) break;	
	}
*/	

				} else {
					outf << str << std::endl;
				}
				}
			}
			inputs.at(i)->close();
			delete inputs.at(i);
		}
		outf.close();
	} else {
		cerr << "Error: cannot write: " << argv[2] << std::endl;
		exit(1);
	}

	outf.close();

	return 0;
}

/*
    src -- The name of one of the source sequences for the alignment. For sequences that are resident in a browser assembly, the form 'database.chromosome' allows automatic creation of links to other assemblies. Non-browser sequences are typically reference by the species name alone.
    start -- The start of the aligning region in the source sequence. This is a zero-based number. If the strand field is '-' then this is the start relative to the reverse-complemented source sequence.
    size -- The size of the aligning region in the source sequence. This number is equal to the number of non-dash characters in the alignment text field below.
    strand -- Either '+' or '-'. If '-', then the alignment is to the reverse-complemented source.
    srcSize -- The size of the entire source sequence, not just the parts involved in the alignment.
    text -- The nucleotides (or amino acids) in the alignment and any insertions (dashes) as well.
*/

void write_maf_record(ostream & out, const std::string & src, size_t start, size_t size, char strand, size_t srcSize, const string & alignment_part) {
	out << "s " << src << " " << start << " " << size << " " << strand << " " << srcSize << " " << alignment_part << std::endl;
}


/*
	write half a pairwise alignment to a maf file 
	(this write the s-line only)

*/
void write_maf_record(ostream & out, const all_data & data, const pw_alignment & al, size_t reference) {
	assert(reference==0 || reference==1);
	// get the requested half alignment
	size_t print_seq;
	size_t print_start;
	size_t print_end;
	std::string print_al;
	// TODO print size not end
	if(reference==0) {
		print_seq = al.getreference1();
		print_start = al.getbegin1();
		print_end = al.getend1();
		print_al = al.get_al_ref1();
	} else {
		print_seq = al.getreference2();
		print_start = al.getbegin2();
		print_end = al.getend2();
		print_al = al.get_al_ref2();
	}
	// translate reference sequence index number to corresponding long name (accession:contig)
	size_t acc = data.accNumber(print_seq);
	std::string accname = data.get_acc(acc);
	std::string seqname = data.get_seq_name(print_seq);
	std::stringstream write_longname;
	write_longname << accname << ':' << seqname;
	char strand = '+';
	size_t size = data.get_seq_size(print_seq);
	size_t al_on_ref_size = print_end - print_start + 1;
	if(print_end < print_start) {
		al_on_ref_size = print_start - print_end + 1;
		strand = '-';
		print_start = size - print_start -1;


		
	} 
	write_maf_record(out, write_longname.str(), print_start, al_on_ref_size, strand, size, print_al);
	

}

/*
	compute star phylogeny multiple sequence alignment
	modifies input alignments by inserting gaps
*/   

void msa_star_alignment(const std::string & center, std::vector<pw_alignment> & alignments) {
	std::vector<std::string> cparts;
	strsep(center, ":", cparts);
	size_t center_ref = atoi(cparts.at(0).c_str());
	size_t center_left = atoi(cparts.at(1).c_str());
	size_t center_length = 0; // length of the cluster center on its reference
	std::vector<size_t> gaps_before;

	for(size_t i=0; i<alignments.size(); ++i) {
		size_t ref1 = alignments.at(i).getreference1();

		size_t ref2 = alignments.at(i).getreference2();
		size_t left1, right1, left2, right2;
		alignments.at(i).get_lr1(left1, right1);
		alignments.at(i).get_lr2(left2, right2);

	//	std::cout << " al " << i << std::endl;
	//	alignments.at(i).print();
	//	std::cout << std::endl;
		// find cluster center and make sure that each alignment goes to cluster center
		// and all cluster centers are identical on the cluster center reference
		// then find the minimal number of gaps before each cluster center reference position
		if(ref1 == center_ref && left1 == center_left) {
			if(center_length==0) {
				center_length = right1 - left1 + 1;
				gaps_before = std::vector<size_t>(center_length + 1, 0);
			} else {
				std::cout << " cl " << center_length << " l1 " << left1 << " r1 " << right1 << std::endl;
				assert(center_length == right1 - left1 +1);
			}
			size_t center_ref_pos=0; // relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0;
			// now go over alignment
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c1 == '-') {
					gapcounter++;
				} else {
					if(gaps_before.at(center_ref_pos) < gapcounter) {
						gaps_before.at(center_ref_pos) = gapcounter;
					}	
					gapcounter = 0;
					center_ref_pos ++;
				}		
			}
			// gaps at end of alignment
			if(gaps_before.at(center_ref_pos) < gapcounter) {
				gaps_before.at(center_ref_pos) = gapcounter;
			}
		
		} else {
			assert(ref2 == center_ref && left2 == center_left);
			if(center_length==0) {
				center_length = right2 - left2 + 1;
				gaps_before = std::vector<size_t>(center_length + 1, 0);
			} else {
				std::cout << " cl " << center_length << " l2 " << left2 << " r2 " << right2 << std::endl;

				assert(center_length == right2 - left2 +1);
			}
			size_t center_ref_pos=0; // relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0;
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c2 == '-') {
					gapcounter++;
				} else {
					if(gaps_before.at(center_ref_pos) < gapcounter) {
						gaps_before.at(center_ref_pos) = gapcounter;
					}	
					gapcounter = 0;
					center_ref_pos ++;
				}
			}
			// gaps at end of alignment
			if(gaps_before.at(center_ref_pos) < gapcounter) {
				gaps_before.at(center_ref_pos) = gapcounter;
			}
		} // if ref	
	} // for i 


	// now we have inserted a lot of gaps in to gaps_before
	// all that gaps are inserted into all alignments to create a multiple alignment where all sequences have the same length
	for(size_t i=0; i<alignments.size(); ++i) {
		size_t ref1 = alignments.at(i).getreference1();
		size_t ref2 = alignments.at(i).getreference2();
		size_t left1, right1, left2, right2;
		alignments.at(i).get_lr1(left1, right1);
		alignments.at(i).get_lr2(left2, right2);

		// for writing new al std::strings
		std::stringstream centerref_al;
		std::stringstream otherref_al;

		if(ref1 == center_ref && left1 == center_left) {
					
			size_t center_ref_pos=0; // we are before that position relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0; // how many gaps did we already see on cluster center ref
			// now go over alignment
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c1 == '-') {
					gapcounter++;
					centerref_al << c1;
					otherref_al << c2;
				} else {
					// add gaps before next center ref symbol
					for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
						centerref_al << '-';
						otherref_al << '-';
					}
					centerref_al << c1;
					otherref_al << c2;
					gapcounter = 0;
					center_ref_pos ++;
				}
			}
			for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
				centerref_al << '-';
				otherref_al << '-';
			}
			pw_alignment newal(centerref_al.str(), otherref_al.str(), 
					alignments.at(i).getbegin1(), alignments.at(i).getbegin2(),
					alignments.at(i).getend1(), alignments.at(i).getend2(),
					alignments.at(i).getreference1(), alignments.at(i).getreference2());
			alignments.at(i) = newal;

		} else {
			assert(ref2 == center_ref && left2 == center_left);
		
			size_t center_ref_pos=0; // relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0;
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c2 == '-') {
					gapcounter++;
					centerref_al << c2;
					otherref_al << c1;
				} else {
					// add gaps before next center ref symbol
					for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
						centerref_al << '-';
						otherref_al << '-';
					}
					centerref_al << c1;
					otherref_al << c2;
					gapcounter = 0;
					center_ref_pos ++;
				}
			
			}
			for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
				centerref_al << '-';
				otherref_al << '-';
			}
			// swap alignment, cluster center is now on reference 1
			pw_alignment newal(centerref_al.str(), otherref_al.str(), 
					alignments.at(i).getbegin2(), alignments.at(i).getbegin1(),
					alignments.at(i).getend2(), alignments.at(i).getend1(),
					alignments.at(i).getreference2(), alignments.at(i).getreference1());
			alignments.at(i) = newal;
		} // else: center on ref 2 of al i
	} // for alignments

// test code
	size_t length = 0;
	for(size_t i=0; i<alignments.size(); ++i) {
		if(length == 0) {
			length = alignments.at(i).alignment_length();
		}
		assert(length == alignments.at(i).alignment_length());
	}

}


	/*
		write the graph as multiple alignment data structure	
	*/

void write_graph_maf(const std::string & graphout, const std::map<string, std::vector<pw_alignment> > & cluster_result_al, const all_data & data) {
	
	std::ofstream gout(graphout.c_str());
	gout << "##maf version=1 scoring=probability" << std::endl;
	gout << "# graph version " << VERSION << std::endl;
	// TODO add more metadata
	gout << "# Coordinates are 0-based.  For - strand matches, coordinates" << std::endl;
	gout << "# in the reverse complement of the 2nd sequence are used." << std::endl;
	gout << "#" << std::endl;
	gout << "# name start alnSize strand seqSize alignment" << std::endl;
	gout << "#" << std::endl;
	size_t cluster_number = 0; 
	for(std::map<std::string, std::vector<pw_alignment> >::const_iterator it = cluster_result_al.begin(); it!=cluster_result_al.end(); ++it) {
		std::string center = it->first;
		std::vector<pw_alignment> als = it->second;
	//	std::cout << "Center: "<< center << " als " << als.size() << std::endl;
		msa_star_alignment(center, als);
		if(als.size()>0) {
			gout << "# cluster " << cluster_number << std::endl;
			gout << "a score=1" << std::endl; // TODO

		

			write_maf_record(gout, data, als.at(0), 0);
			for(size_t i=0; i<als.size(); ++i) {
				write_maf_record(gout, data, als.at(i), 1);
			}	
			cluster_number++;
		}
	}

	gout.close();





}


int do_mc_model(int argc, char * argv[]) {//mco bardar
	typedef mc_model use_model;
//	typedef clustering<use_model> use_clustering;
	typedef initial_alignment_set<use_model> use_ias;
	typedef affpro_clusters<use_model> use_affpro;
	if(argc < 5) {
		usage();
		cerr << "Program: model" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* output maf file for the graph" << std::endl;
		cerr << "* number of threads to use (optional, default 10)" << std::endl;
	}

// TODO read maf or sam

	clock_t read_data_time = clock();
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string graphout(argv[4]);
	size_t num_threads = 1;
	if(argc == 6) {
		num_threads = atoi(argv[5]);
	}
	std::ofstream outs("encode",std::ofstream::binary);
// Read all data
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	overlap ol(data);
	wrapper wrap;
	test_encoder test;

	read_data_time = clock() - read_data_time;
//	encoding_functor functor1(data);


// Find connected components of alignments with some overlap
	clock_t initial_cc_time = clock();
	compute_cc cccs(data);
	std::vector<std::set< const pw_alignment *, compare_pw_alignment> > ccs;
	cccs.compute(ccs); //fill in ccss, order by size(Notice that they are just connected to each other if they have overlap!)
	initial_cc_time = clock() - initial_cc_time;
/*	std::cout << " Found " << ccs.size() << " connected components" << std::endl;
	for(size_t i=0; i<ccs.size(); i++) {
		std::cout << "Connected component "<< i << " contains " << ccs.at(i).size() << " alignments" << std::endl;
	}*/
// Train the model on all data
	use_model m(data);
	m.train(outs);
//	double sum = 0;
//	double c1;
//	double c2;
//	double m1;
//	double m2;
//	size_t len=0;
/*	for(size_t i =0; i <  data.numSequences(); i++){
		const dnastd::string & sequence = data.getSequence(i);
		size_t length = sequence.length();
		len +=length;
	}*/
/*	std::cout<<"length: " << len << std::endl;
	for(size_t i = 0 ; i< data.numAlignments() ; i++){
		const pw_alignment & p = data.getAlignment(i);
		size_t ref1 = p.getreference1();
		size_t ref2 = p.getreference2();
		const dnastd::string & sequence1 = data.getSequence(ref1);
		const dnastd::string & sequence2 = data.getSequence(ref2);
		size_t length1 = sequence1.length();
		size_t length2 = sequence2.length();
		m.cost_function(p, c1,c2,m1,m2,outs);
		sum +=c1+c2;
		len +=length1 + length2;
	}
	std::cout<< "sum is "<<sum <<std::endl;
	std::cout<< "len times two: "<< len*2 <<std::endl;*/
//	counting_functor functor(data);
//	encoder en(data,m);
//	use_clustering clust(ol,data,m);
	encoder en(data,m,wrap);
	std::vector<overlap> cc_overlap(ccs.size(), overlap(data));// An overlap class object is needed for each connected component, thus we cannot use ol
	// base cost to use an alignment (information need for its adress)
	double cluster_base_cost = log2(data.numAlignments());
//	std::cout << " base cost " << cluster_base_cost << std::endl;
// Select an initial alignment std::set for each connected component (in parallel)
	// TODO put the next two into a single data structure
	std::map<std::string, std::vector<string> > global_results;//for each center returns all its cluster members TODO needed?
	std::map<std::string, std::vector<pw_alignment> > alignments_in_a_cluster;//string ---> center of a cluster, vector ---> alignments with that center

	size_t num_clusters = 0;
	size_t num_cluster_members_al = 0;
	size_t num_cluster_seq = 0;
	size_t num_cluster_inputs_al = 0;
#if TIMING
	clock_t ias_time = 0;
	clock_t test_function_time = 0;
	clock_t second_cc_time = 0;
	clock_t ap_time = 0;

#endif


#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for(size_t i=0; i<ccs.size(); ++i) {
#pragma omp critical(print)
{
//		std::cout << " on initial CC " << i << " size " << ccs.at(i).size() << std::endl;
}
		clock_t ias_time_local = clock();
		std::set< const pw_alignment *, compare_pw_alignment> & cc = ccs.at(i);
		use_ias ias(data,cc, m, cluster_base_cost,outs);
		ias.compute(cc_overlap.at(i),outs); //Marion : I have an error that I dont understand
		ias_time_local = clock() - ias_time_local;

		clock_t test_time_local = clock();
		cc_overlap.at(i).test_partial_overlap(); // TODO remove slow test function
		test_time_local = clock() - test_time_local;

	//	std::cout << " number of alignments " << cc_overlap.at(i).size() << std::endl;
		// TODO this can be done a lot faster because there is no partial overlap here
		clock_t second_cc_time_local = clock();
		std::vector< std::set<const pw_alignment *, compare_pw_alignment> > cc_cluster_in; //has shorter length than original sequence pieces
		compute_cc occ(cc_overlap.at(i), data.numSequences());
		occ.compute(cc_cluster_in);
		second_cc_time_local = clock() - second_cc_time_local;
#if TIMING
#pragma omp critical(time)
{
		ias_time += ias_time_local;
		test_function_time += test_time_local;
		second_cc_time += second_cc_time_local;
}
#endif

	//	std::cout<< "cc_cluster_in size: "<< cc_cluster_in.size()<<std::endl;
	//	std::cout << "cc " << i << ": from " << cc.size() << " original als with total gain " << ias.get_max_gain() << " we made " << cc_overlap.at(i).size() << " pieces with total gain " << ias.get_result_gain() <<  std::endl; 
	//	std::cout << " components for clustering: " << std::endl;
		std::vector<std::string> all_centers;
		for(size_t j=0; j<cc_cluster_in.size(); ++j) {

		//	std::cout << " run affpro on " << cc_cluster_in.at(j).size()<<std::endl;
		//	data.numAcc();

			/*
			std::cout << " in al std::set size " << cc_cluster_in.at(j).size() << std::endl;
			size_t nums=0;
			for(std::set<const pw_alignment *, compare_pw_alignment>::iterator ait=cc_cluster_in.at(j).begin(); ait!=cc_cluster_in.at(j).end(); ++ait) {
				std::cout << " nums " << nums << std::endl;
				nums++;
				(*ait)->print();
				std::cout << std::endl;
			}
			*/

			clock_t ap_time_local = clock();
			use_affpro uaf(cc_cluster_in.at(j), m, cluster_base_cost,outs);
			std::map<std::string, std::vector<string> >cluster_result;
			uaf.run(cluster_result);
			ap_time_local = clock() - ap_time_local;
#if TIMING
#pragma omp critical(aptime) 
{
			ap_time += ap_time_local;
}

#endif

			std::map<std::string, std::vector<pw_alignment> > local_al_in_a_cluster;


			size_t counter =0;
			for(std::map<std::string,std::vector<string> >::iterator it=cluster_result.begin();it !=cluster_result.end();it++){
				counter++;
				std::string center = it->first;
			/*	
				std::vector<std::string> clmembers = it->second;
				for(size_t k=0; k<clmembers.size(); ++k) {
					std::cout << "clm " << clmembers.at(k) << " in " << center << std::endl;
				}
			*/

				std::map<std::string, std::vector<pw_alignment> >::iterator it1 = local_al_in_a_cluster.find(center);
				if(it1 == local_al_in_a_cluster.end()){
					local_al_in_a_cluster.insert(make_pair(center, std::vector<pw_alignment>()));
					it1 = local_al_in_a_cluster.find(center);
				}
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int ref = atoi(center_parts.at(0).c_str());
				unsigned int left = atoi(center_parts.at(1).c_str());

				// look at all input alignments and determine which cluster it belongs to
				for(std::set<const pw_alignment*, compare_pw_alignment>::iterator it2 = cc_cluster_in.at(j).begin(); it2!=cc_cluster_in.at(j).end(); ++it2){
					const pw_alignment * al = *it2;
			//		std::cout<< "alignment in cc cluster in at " << j << std::endl;
			//		al->print();
			//		double g1;
			//		double g2;
			//		m.gain_function(*al,g1,g2,outs);
			//		std::cout<< "g1: "<<g1 << " g2: "<< g2 <<std::endl;
					size_t ref1 = al->getreference1();
					size_t ref2 = al->getreference2();
					size_t left1;
					size_t left2;
					size_t right1;
					size_t right2;
					al->get_lr1(left1,right1);
					al->get_lr2(left2,right2);
					if(ref1 == ref && left1 == left){
						for(size_t k = 0; k < it->second.size();k++){
							std::vector<std::string> member_parts;
							strsep(it->second.at(k), ":" , member_parts);
							unsigned int mem_ref = atoi(member_parts.at(0).c_str());
							unsigned int mem_left = atoi(member_parts.at(1).c_str());
							if(ref2 == mem_ref && left2 == mem_left){
								it1->second.push_back(*al);
							}else continue;
						}
					}else if(ref2 == ref && left2 == left){
							for(size_t k = 0; k < it->second.size();k++){
								std::vector<std::string> member_parts;
								strsep(it->second.at(k), ":" , member_parts);
								unsigned int mem_ref = atoi(member_parts.at(0).c_str());
								unsigned int mem_left = atoi(member_parts.at(1).c_str());
								if(ref1 == mem_ref && left1 == mem_left){
									it1->second.push_back(*al);
								}else continue;
							}
					}else continue;
				} // for cluster_in std::set
			}
#pragma omp critical 
{
			num_clusters += cluster_result.size();
			num_cluster_inputs_al += cc_cluster_in.at(j).size();
			
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it = local_al_in_a_cluster.begin(); it!=local_al_in_a_cluster.end(); ++it) {
				alignments_in_a_cluster.insert(*it);
				num_cluster_members_al += 1 + it->second.size();
				
			
			}
}

		/*	std::set<std::string> intermediate_center;
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
				if(it->second.size()==0){
					intermediate_center.insert(it->first);
				}
			}
			for(std::set<std::string>::iterator it = intermediate_center.begin();it !=intermediate_center.end();it++){
				std::string cent = *it;
				std::map<std::string, std::vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(cent);
				alignments_in_a_cluster.erase(it1);				
			}*/	
			
		//	std::cout<<"counter: "<<counter<<std::endl;
			
#pragma omp critical
{
			for(std::map<std::string,std::vector<string> >::iterator it=cluster_result.begin();it != cluster_result.end();it++){
				std::string center = it->first;
				num_cluster_seq+= 1 + it->second.size();
				std::map<std::string,std::vector<string> >::iterator it1 = global_results.find(center);
				if(it1 == global_results.end()){
					global_results.insert(make_pair(center, std::vector<std::string>()));
					it1 = global_results.find(center);
				}
				for(size_t i =0; i < it->second.size();i++){
					it1->second.push_back(it->second.at(i));
				}
			}
}
		} // for cluster_in std::set  
		std::cout << std::endl;



/*	for(std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.begin(); it != alignments_in_a_cluster.end();it++){
		for(size_t h=0; h< it->second.size();h++){
			it->second.at(h).print();
		}
	}*/	
	
				
	/*	for(std::set< const pw_alignment *, compare_pw_alignment>::iterator it=cc.begin();it!=cc.end();it++ ){
			const pw_alignment *al = *it;
			size_t index = al ->getreference1();

		}*/
//#pragma omp critical(print) 
//{
	//	std::cout << " initial CC " << i << " done " << std::endl << flush;
//}
	} // for connected components



	clock_t graph_ma_time = clock();
	// TODO better separation of the different applications of our program: create/read model, compress/decompress sequences, create graph
	// write graph result in maf format
	write_graph_maf(graphout, alignments_in_a_cluster, data);//includes clusters with no associated member.
	graph_ma_time = clock() - graph_ma_time;

	clock_t arithmetic_encoding_time = clock();
	std::set<std::string> intermediate_center;
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
		if(it->second.size()==0){
			intermediate_center.insert(it->first);
		}
	}
	for(std::set<std::string>::iterator it = intermediate_center.begin();it !=intermediate_center.end();it++){
		std::string cent = *it;
		std::map<std::string, std::vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(cent);
		alignments_in_a_cluster.erase(it1);				
	}
/*	std::cout<< "al cent size: "<<  alignments_in_a_cluster.size() << std::endl;
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
		std::cout << " al cent: "<< it->first <<std::endl;
	}*/
	//Defining weights of global clustering results! 
	std::map<std::string, unsigned int>weight;
	size_t max_bit = 8;	
	size_t max_members = 0;//Returns the largest cluster
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){	
        	if(it->second.size()> max_members){		
			max_members = it->second.size();
		}
		
	}
	
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){
		std::string cluster = it ->first;
//		std::cout<< "global result size: "<< it->second.size() << std::endl;
		std::map<std::string, std::vector<pw_alignment> >::iterator cent = alignments_in_a_cluster.find(cluster);
//		if(it->second.size() !=1){
		if(cent != alignments_in_a_cluster.end()){//Removing the centers with no associated member
			std::map<std::string, unsigned int >::iterator it1 = weight.find(cluster);
			if(it1 == weight.end()){
				weight.insert(make_pair(cluster,0));
				it1 = weight.find(cluster);
			}
			it1->second = (unsigned int)(it->second.size()*((m.get_powerOfTwo().at(max_bit)-1)/(double)max_members));
			if(it1 ->second == 0){
				it1->second = 1;
			}
		}else continue;
	}
	size_t no_al = 0;
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin();it !=global_results.end();it++){
			no_al += it->second.size();
	}
//	std::cout<<"no of clustered alignments: "<< no_al <<std::endl;	
	std::map<std::string, std::vector<string> >membersOfCluster;//first string represents center and vector of strings are associated members
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
		std::map<std::string, std::vector<string> >::iterator it1 = membersOfCluster.find(it->first);
		for(size_t j =0; j < it->second.size();j++){
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			pw_alignment p = it->second.at(j);
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			if(it1 == membersOfCluster.end()){
				membersOfCluster.insert(make_pair(it->first,std::vector<std::string>()));
				it1 = membersOfCluster.find(it->first);
			}
			std::stringstream sample1;
			std::stringstream sample2;
			sample1 << ref1 << ":" << left_1;
			sample2 << ref2 << ":" << left_2;
			it1->second.push_back(sample1.str());
			it1->second.push_back(sample2.str());
		}
	}
	std::map<std::string, string>member_of_cluster;// first string is a associated one and the second one is its center
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
		for(size_t j =0; j < it->second.size();j++){
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			pw_alignment p = it->second.at(j);
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			std::stringstream sample1;
			std::stringstream sample2;
			sample1 << ref1 << ":" << left_1;
			sample2 << ref2 << ":" << left_2;
			if(sample1.str() != it->first){
				member_of_cluster.insert(make_pair(sample1.str(),it->first));
			}else{
				member_of_cluster.insert(make_pair(sample2.str(), it->first));
			}
		}
		member_of_cluster.insert(make_pair(it->first, it->first));		
	}
	std::cout<< "size of memeber_of_cluster is : "<< member_of_cluster.size()<<std::endl;


//	c.calculate_similarity();
//	c.update_values();
//	c.update_clusters();

//	initial_alignment_std::set<model> ias(data, m);
//	ias.compute(o);
//	std::cout << "There are " << o.size() << " alignment parts without partial overlap" << std::endl;
//	std::ifstream in("encode",std::ifstream::binary);
//	m.std::set_patterns(in);
//	for(std::map<std::string, std::vector<unsigned int> >::const_iterator it = m.get_high(0).begin(); it != m.get_high(0).end(); it++){
//		for(size_t k =0; k <5; k++){
//							std::cout<< it->second.at(k)<< " ; ";
//						}
//						std::cout<< " "<< std::endl;

//	}
/*	for(size_t i = 0 ; i < data.numAlignments(); i++){
		const pw_alignment & p = data.getAlignment(i);
		size_t left1;
		size_t left2;
		size_t right1;
		size_t right2;
		p.get_lr1(left1,right1);
		p.get_lr2(left2,right2);
		if(p.getreference1() == 1 && p.getreference2() == 25 && left2 == 70424){
			std::cout<< "left is 70424 and ref 2 is 25"<<std::endl;
			std::cout<< "left 1 is " << left1 << std::endl;
		//	p.print();
			double g1;
			double g2;
			double c1;
			double c2;
			double m1;
			double m2;
			m.cost_function(p,c1,c2,m1,m2,outs);
			std::cout << "m1: "<< m1 << " m2: "<< m2<< " c1: "<< c1<< " c2: "<< c2 << std::endl;
			m.gain_function(p,g1,g2,outs);
			std::cout<< "g1: "<<g1 << " g2: "<< g2 <<std::endl;
		}
		if (p.getreference1() == 25&&p.getreference2()==1&&left1 == 70424){
			std::cout<< "left is 70424 and ref 1 is 25"<<std::endl;
			std::cout<< "left 2 is " << left2 << std::endl;
		//	p.print();
			double g1;
			double g2;
			double c1;
			double c2;
			double m1;
			double m2;
			m.cost_function(p,c1,c2,m1,m2,outs);
			std::cout << "m1: "<< m1 << " m2: "<< m2<< " c1: "<< c1<< " c2: "<< c2 << std::endl;
			m.gain_function(p,g1,g2,outs);
			std::cout<< "g1: "<<g1 << " g2: "<< g2 <<std::endl;
		}
	}

	*/
//Data compression:
	std::cout<< "weight size: "<< weight.size()<<std::endl;
//	en.arithmetic_encoding_alignment(weight,member_of_cluster,alignments_in_a_cluster,outs);
//	en.write_to_stream(alignments_in_a_cluster,outs);
	std::ofstream al_encode("align_encode",std::ofstream::binary);
	dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
//	en.arithmetic_encoding_seq(outs);
//	en.calculate_high_in_partition(weight,alignments_in_a_cluster);
//	en.arithmetic_encoding_centers(alignments_in_a_cluster,outs);
//	en.arithmetic_encoding_alignment(weight,member_of_cluster,alignments_in_a_cluster,outs,*enc);

//	en.test_al_encoding(weight,member_of_cluster,alignments_in_a_cluster,outs,*enc);
	delete enc;
	outs.close();
//	test.encode();



	std::ifstream in("encode",std::ifstream::binary);
//	en.read_from_stream(in);
	dlib::entropy_decoder_kernel_1  dec;
//	en.arithmetic_decoding_alignment(in,dec);
//	test.compare();
//	en.arithmetic_decoding_centers(in);


//	en.test_al_decoding(in,dec);
//	test.compare();
	arithmetic_encoding_time = clock() - arithmetic_encoding_time;

	std::cout << "Clustering summary: " << std::endl;
	std::cout << "Input: " << num_cluster_inputs_al << " pw alignments on " << num_cluster_seq <<" sequence pieces " <<std::endl;
	std::cout << "Output: " << num_clusters << " clusters containing " << num_cluster_members_al << " alignments " << std::endl;
#if TIMING
	std::cout << "Time overview: " << std::endl;
	std::cout << "Read data " << (double)read_data_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute initial connected components " << (double)initial_cc_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute initial alignments std::set + remove partial overlap " << (double)ias_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Test functions " << (double)test_function_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute second connected components " << (double)second_cc_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Affinity propagation clustering " << (double)ap_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write MSA-maf graph " << (double)graph_ma_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write compressed file " << (double)arithmetic_encoding_time/CLOCKS_PER_SEC << std::endl;
#endif



	return 0;
}
int do_model_seq(int argc, char * argv[]){//mc o bardashtam
	typedef model use_model;
	if(argc < 5) {
		usage();
		cerr << "Program: model" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* output maf file for the graph" << std::endl;
		cerr << "* number of threads to use (optional, default 10)" << std::endl;
	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string graphout(argv[4]);
	size_t num_threads = 1;
	if(argc == 6) {
		num_threads = atoi(argv[5]);
	}

//Reading data:
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	overlap ol(data);
//	wrapper wrapp;
//	counting_functor functor(data);
//	encoding_functor functor1(data);

//std::ofstream outs("encode",std::ofstream::binary);
//Train all the sequences:
	use_model m(data);
//	encoder en(data,m,wrapp);
	m.train();//outs o ham bardashtam
	double sum = 0;
	double g1;
	double g2;
	for(size_t i = 0 ; i< data.numAlignments() ; i++){
		const pw_alignment & p = data.getAlignment(i);
		m.gain_function(p, g1,g2);
		sum +=g1+g2;
	}
	std::cout<< "sum is "<<sum <<std::endl;
// Find connected components of alignments with some overlap
//	compute_cc cccs(data);
//	std::vector<std::set< const pw_alignment *, compare_pw_alignment> > ccs;
//	cccs.compute(ccs);
//test
//	en.arithmetic_encoding_seq(outs);
//	en.arithmetic_decoding_seq();
//	en.arithmetic_encoding_alignment();
	
	return 0;
}


#if !TEST
	
int main(int argc, char * argv[]) {

	if(argc < 2) {
		usage();
		exit(1);
	}
	std::string program(argv[1]);

	if(0==program.compare("fasta_prepare")) {
		return do_fasta_prepare(argc, argv);
	}
	else if(0==program.compare("model"))
		return do_mc_model(argc, argv);
	
	 
/*	else if(0==program.compare("model")){
		return do_model_seq(argc,argv);
	}*/
	else {
		usage();
	}

	
	return 1;
}



#endif
#include "model.cpp"























