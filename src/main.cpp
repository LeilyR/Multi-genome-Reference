#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include<ostream>
#include<vector>
//#include "/ebio/abt6/lrabbani/Downloads/dlib/dlib/entropy_encoder.h"
#include "pw_alignment.hpp"
#include "data.hpp"
#include "model.hpp"



using namespace std;
//using namespace dlib;


#define TEST 0


#if TEST

int main(int argc, char * argv[]) {
	cout << " hello " << endl;

	pw_alignment p(string("ATT----TTCTT"), string("AGTGATAT----"), 12, 15, 23, 26,1,1);
	pw_alignment s(string("ATT----TTCTT"), string("ACTGATG---AC"),13, 18, 24, 29,2,1);

	for(size_t i=0; i<p.alignment_length(); ++i) {
		char s1;
		char s2;
		p.alignment_col(i, s1, s2);
		cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << endl;
	  
	}

	pw_alignment p1;
	pw_alignment p2;
	p.split(false,4, p1, p2);
		for(size_t i=0; i<p2.alignment_length(); ++i) {
		char s1;
		char s2;
		p2.alignment_col(i, s1, s2);
		cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << endl;
	 
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
		cout << "pos " << i<< " s1 " << s1 << " s2 " << s2 << endl;
}

}
	return 0;

}



#endif
void usage() {
	cerr << "Graph program" << endl;
	cerr << "Usage:" << endl;
}

int do_fasta_prepare(int argc, char * argv[]) {

	if(argc < 4) {
		usage();
		cerr << "Program fasta_prepare" << endl; 
		cerr << "Combine fasta files from different species/accessions into a common file" << endl;
		cerr << "Parameters: " << endl;
		cerr << "output file -- fasta file" << endl;
		cerr << "<list of fasta files> -- input file, one file per species/accession, no colon in file name"<< endl;
		return 1;
	}


	std::ofstream outf(argv[2]);
	if(outf) {
		vector<ifstream*> inputs;
		vector<string> names;
		for(int i=3; i<argc; ++i) {
			string istr(argv[i]);
			vector<string> slashparts;
			strsep(istr, "/", slashparts);
			string afterslash = slashparts.at(slashparts.size()-1);
			vector<string> periodparts;
			strsep(afterslash, ".", periodparts);
			string name = periodparts.at(0);
			if(periodparts.size()>1) {
				stringstream sstr;
				for(size_t j=0; j<periodparts.size()-1; j++) {
					if(j!=0) sstr << ".";
					sstr << periodparts.at(j);
				}
				name = sstr.str();
			}
			for(size_t j=0; j<name.length(); ++j) {
				if(':' == name.at(j)) {
					cerr << "Error: " << argv[i] << " file name contains a colon" << endl;
					exit(1);
				}
			}
			ifstream * iin = new ifstream(argv[i]);
			if(!(*iin)) {
				cerr << "Error: cannot read: " << argv[i] << endl;
				exit(1);
			}
			inputs.push_back(iin);
			names.push_back(name);
		}

		for(size_t i=0; i<inputs.size(); ++i) {


			string str;
			while(getline((*(inputs.at(i))), str)) {
				if(str.length()>0) { // ignore empty lines
				if(str.at(0)=='>') {
					string nname = str.substr(1);
					stringstream nhead;
					nhead << ">" << names.at(i) << ":" << nname;
					outf << nhead.str() << endl;

						
						
/*
	for(size_t col = 0; col < alignment_length(); col++) {
		char c1;
		char c2;
		alignment_col(col, c1, c2);
		cout <<col <<"\t"<< c1<<"\t"<<c2<<endl;
	if(col>50) break;	
	}
*/	

				} else {
					outf << str << endl;
				}
				}
			}
			inputs.at(i)->close();
			delete inputs.at(i);
		}
		outf.close();
	} else {
		cerr << "Error: cannot write: " << argv[2] << endl;
		exit(1);
	}

	outf.close();

	return 0;
}
 


int do_mc_model(int argc, char * argv[]) {
	typedef model use_model;
	typedef clustering<use_model> use_clustering;
	typedef initial_alignment_set<use_model> use_ias;
	typedef affpro_clusters<use_model> use_affpro;
	if(argc < 4) {
		usage();
		cerr << "Program: model" << endl;
		cerr << "Parameters:" << endl;
		cerr << "* fasta file from fasta_prepare" << endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << endl;
		cerr << "* number of threads to use (optional, default 10)" << endl;
	}

	string fastafile(argv[2]);
	string maffile(argv[3]);
	size_t num_threads = 1;
	if(argc == 5) {
		num_threads = atoi(argv[4]);
	}
	
// Read all data
	all_data data(fastafile, maffile);
	overlap ol(data);
	counting_functor functor(data);


// Find connected components of alignments with some overlap
	compute_cc cccs(data);
	vector<set< const pw_alignment *, compare_pw_alignment> > ccs;
	cccs.compute(ccs);
	cout << " Found " << ccs.size() << " connected components" << endl;
	for(size_t i=0; i<ccs.size(); i++) {
		cout << "Connected component "<< i << " contains " << ccs.at(i).size() << " alignments" << endl;
	}
// Train the model on all data
	use_model m(data);
	m.train();
//	mc_model m(data);
//	entropy_encoder_kernel_1 k();
//	m.markov_chain();
//	m.markov_chain_alignment();
	use_clustering clust(ol,data,m);

	vector<overlap> cc_overlap(ccs.size(), overlap(data));
	// base cost to use an alignment (information need for its adress)
	double cluster_base_cost = log2(data.numAlignments());
	cout << " base cost " << cluster_base_cost << endl;
// Select an initial alignment set for each connected component (in parallel)
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for(size_t i=0; i<ccs.size(); ++i) {
		set< const pw_alignment *, compare_pw_alignment> & cc = ccs.at(i);
		use_ias ias(data, cc, m, cluster_base_cost);
		ias.compute(cc_overlap.at(i));

		// TODO this can be done a lot faster because there is no partial overlap here
		vector< set<const pw_alignment *, compare_pw_alignment> > cc_cluster_in;
		compute_cc occ(cc_overlap.at(i), data.numSequences());
		occ.compute(cc_cluster_in);
#pragma omp critical
{
		cout << "cc " << i << ": from " << cc.size() << " original als with total gain " << ias.get_max_gain() << " we made " << cc_overlap.at(i).size() << " pieces with total gain " << ias.get_result_gain() <<  endl; 
		cout << " components for clustering: " << endl;
		for(size_t j=0; j<cc_cluster_in.size(); ++j) {
			cout << " run affpro on " << cc_cluster_in.at(j).size();

			use_affpro uaf(cc_cluster_in.at(j), m, cluster_base_cost);
			uaf.run();
		}
		cout << endl;


}

		
	
				
	//	for(set< const pw_alignment *, compare_pw_alignment>::iterator it=cc.begin();it!=cc.end();it++ ){
	//		const pw_alignment *al = *it;
	//		al->getreference1();
	//		k.getstream();

	//	}
	}

//	c.calculate_similarity();
//	c.update_values();
//	c.update_clusters();

//	initial_alignment_set<model> ias(data, m);
//	ias.compute(o);
//	cout << "There are " << o.size() << " alignment parts without partial overlap" << endl;

	return 0;
}


#if !TEST
	
int main(int argc, char * argv[]) {

	if(argc < 2) {
		usage();
		exit(1);
	}
	string program(argv[1]);

	if(0==program.compare("fasta_prepare")) {
		return do_fasta_prepare(argc, argv);
	} else if(0==program.compare("model")) {
		return do_mc_model(argc, argv);
	}

	else {
		usage();
	}


	return 1;
}



#include "model.cpp"
#endif























