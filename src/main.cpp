#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include<ostream>
#include<vector>
#include "pw_alignment.hpp"
#include "data.hpp"
#include "model.hpp"
#include "encoder.hpp"
#include "test.hpp"
#define NO_MAKEFILE
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"




using namespace std;
// using namespace dlib;

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
	typedef mc_model use_model;
//	typedef clustering<use_model> use_clustering;
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
ofstream outs("encode",std::ofstream::binary);
// Read all data
	all_data data(fastafile, maffile);
	overlap ol(data);
	wrapper wrap;
//	encoding_functor functor1(data);


// Find connected components of alignments with some overlap
	compute_cc cccs(data);
	vector<set< const pw_alignment *, compare_pw_alignment> > ccs;
	cccs.compute(ccs); //fill in ccss, order by size(Notice that they are just connected to each other if they have overlap!)
/*	cout << " Found " << ccs.size() << " connected components" << endl;
	for(size_t i=0; i<ccs.size(); i++) {
		cout << "Connected component "<< i << " contains " << ccs.at(i).size() << " alignments" << endl;
	}*/
// Train the model on all data
	use_model m(data);
	m.train(outs);
//	counting_functor functor(data);
//	encoder en(data,m);
//	use_clustering clust(ol,data,m);
	encoder en(data,m,wrap);
	test_encoder test;


	vector<overlap> cc_overlap(ccs.size(), overlap(data));// vase in ke vase har connected component i yek overlap lazem darim ke dar natije bishtar az yeki overlap mikhaim vase hamin ol o estefade nakaradim
	// base cost to use an alignment (information need for its adress)
	double cluster_base_cost = log2(data.numAlignments());
//	cout << " base cost " << cluster_base_cost << endl;
// Select an initial alignment set for each connected component (in parallel)
	map<string, vector<string> > global_results;//for each center returns all its cluster members
	map<string, vector<pw_alignment> > alignments_in_a_cluster;//string ---> center of a cluster, vector ---> alignments with that center
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for(size_t i=0; i<ccs.size(); ++i) {
		set< const pw_alignment *, compare_pw_alignment> & cc = ccs.at(i);
		use_ias ias(data,cc, m, cluster_base_cost,outs);
		ias.compute(cc_overlap.at(i),outs);
	//	cout << " number of alignments " << cc_overlap.at(i).size() << endl;
		// TODO this can be done a lot faster because there is no partial overlap here
		vector< set<const pw_alignment *, compare_pw_alignment> > cc_cluster_in; //has shorter length than original sequence pieces
		compute_cc occ(cc_overlap.at(i), data.numSequences());
		occ.compute(cc_cluster_in);
	//	cout<< "cc_cluster_in size: "<< cc_cluster_in.size()<<endl;
	//	cout << "cc " << i << ": from " << cc.size() << " original als with total gain " << ias.get_max_gain() << " we made " << cc_overlap.at(i).size() << " pieces with total gain " << ias.get_result_gain() <<  endl; 
	//	cout << " components for clustering: " << endl;
		vector<string> all_centers;
		for(size_t j=0; j<cc_cluster_in.size(); ++j) {
		//	cout << " run affpro on " << cc_cluster_in.at(j).size()<<endl;
			data.numAcc();
			use_affpro uaf(cc_cluster_in.at(j), m, cluster_base_cost,outs);
			map<string, vector<string> >cluster_result;
			uaf.run(cluster_result);
			size_t counter =0;
			for(map<string,vector<string> >::iterator it=cluster_result.begin();it !=cluster_result.end();it++){
				counter++;
				string center = it->first;
				map<string, vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(center);
				if(it1 == alignments_in_a_cluster.end()){
					alignments_in_a_cluster.insert(make_pair(center, vector<pw_alignment>()));
					it1 = alignments_in_a_cluster.find(center);
				}
				vector<string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int ref = atoi(center_parts.at(0).c_str());
				unsigned int left = atoi(center_parts.at(1).c_str());
				for(set<const pw_alignment*, compare_pw_alignment>::iterator it2 = cc_cluster_in.at(j).begin(); it2!=cc_cluster_in.at(j).end(); ++it2){
					const pw_alignment * al = *it2;
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
							vector<string> member_parts;
							strsep(it->second.at(k), ":" , member_parts);
							unsigned int mem_ref = atoi(member_parts.at(0).c_str());
							unsigned int mem_left = atoi(member_parts.at(1).c_str());
							if(ref2 == mem_ref && left2 == mem_left){
								it1->second.push_back(*al);
							}else continue;
						}
					}else if(ref2 == ref && left2 == left){
							for(size_t k = 0; k < it->second.size();k++){
								vector<string> member_parts;
								strsep(it->second.at(k), ":" , member_parts);
								unsigned int mem_ref = atoi(member_parts.at(0).c_str());
								unsigned int mem_left = atoi(member_parts.at(1).c_str());
								if(ref1 == mem_ref && left1 == mem_left){
									it1->second.push_back(*al);
								}else continue;
							}					
					}else continue;
				}
			}
		/*	set<string> intermediate_center;
			for(map<string, vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
				if(it->second.size()==0){
					intermediate_center.insert(it->first);
				}
			}
			for(set<string>::iterator it = intermediate_center.begin();it !=intermediate_center.end();it++){
				string cent = *it;
				map<string, vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(cent);
				alignments_in_a_cluster.erase(it1);				
			}*/	
			
		//	cout<<"counter: "<<counter<<endl;
			
#pragma omp critical
{
			for(map<string,vector<string> >::iterator it=cluster_result.begin();it != cluster_result.end();it++){
				string center = it->first;
				map<string,vector<string> >::iterator it1 = global_results.find(center);
				if(it1 == global_results.end()){
					global_results.insert(make_pair(center, vector<string>()));
					it1 = global_results.find(center);
				}
				for(size_t i =0; i < it->second.size();i++){
					it1->second.push_back(it->second.at(i));
				}
			}
}
		}
		cout << endl;



/*	for(map<string, vector<pw_alignment> >::iterator it=alignments_in_a_cluster.begin(); it != alignments_in_a_cluster.end();it++){
		for(size_t h=0; h< it->second.size();h++){
			it->second.at(h).print();
		}
	}*/	
	
				
	/*	for(set< const pw_alignment *, compare_pw_alignment>::iterator it=cc.begin();it!=cc.end();it++ ){
			const pw_alignment *al = *it;
			size_t index = al ->getreference1();

		}*/
	}
	set<string> intermediate_center;
	for(map<string, vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
		if(it->second.size()==0){
			intermediate_center.insert(it->first);
		}
	}
	for(set<string>::iterator it = intermediate_center.begin();it !=intermediate_center.end();it++){
		string cent = *it;
		map<string, vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(cent);
		alignments_in_a_cluster.erase(it1);				
	}
/*	cout<< "al cent size: "<<  alignments_in_a_cluster.size() << endl;
	for(map<string, vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
		cout << " al cent: "<< it->first <<endl;
	}*/
	//Defining weights of global clustering results! 
	map<string, unsigned int>weight;
	size_t max_bit = 8;	
	size_t members = 0;//Returns the largest cluster
	for(map<string, vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){	
        	if(it->second.size()> members){		
			members = it->second.size();
		}
		
	}
	
	for(map<string, vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){
		string cluster = it ->first;
	//	cout<< "global result size: "<< it->second.size() << endl;
		map<string, vector<pw_alignment> >::iterator cent = alignments_in_a_cluster.find(cluster);
//		if(it->second.size() !=1){
		if(cent != alignments_in_a_cluster.end()){//Removing the centers with no associated member
			map<string, unsigned int >::iterator it1 = weight.find(cluster);
			if(it1 == weight.end()){
				weight.insert(make_pair(cluster,0));
				it1 = weight.find(cluster);
			}
			it1->second = (unsigned int)(it->second.size()*((m.get_powerOfTwo().at(max_bit)-1)/(double)members));
			if(it1 ->second == 0){
				it1->second = 1;
			}
		}else continue;
	}

/*	for(map<string,unsigned int>::iterator it=weight.begin();it !=weight.end();it++){
		cout<<"weight cent: "<< it ->first <<endl;
	}*/
	map<string, vector<string> >membersOfCluster;//first string represents center and vector of strings are associated members(It can be removed and replaced by global result!)
/*	for(map<string, vector<string> >::iterator it= global_results.begin(); it !=global_results.end(); it++){
		for(size_t i =0; i < it->second.size();i++){
			membersOfCluster.insert(make_pair(it->second.at(i),it->first));
		}
	}
	set<string> intermediate;
	for(map<string, vector<string> >::iterator it = global_results.begin(); it!=global_results.end();it++){
		if(it->second.size()==1){
			intermediate.insert(it->first);
		}
	}
	for(set<string>::iterator it = intermediate.begin();it !=intermediate.end();it++){
		string cent = *it;
		map<string, string>::iterator it1 = membersOfCluster.find(cent);
		membersOfCluster.erase(it1);				
	}*/ 
//Replaced by:
	for(map<string, vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
		map<string, vector<string> >::iterator it1 = membersOfCluster.find(it->first);
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
				membersOfCluster.insert(make_pair(it->first,vector<string>()));
				it1 = membersOfCluster.find(it->first);
			}
			stringstream sample1;
			stringstream sample2;
			sample1 << ref1 << ":" << left_1;
			sample2 << ref2 << ":" << left_2;
			it1->second.push_back(sample1.str());
			it1->second.push_back(sample2.str());
		}
	}
//	c.calculate_similarity();
//	c.update_values();
//	c.update_clusters();

//	initial_alignment_set<model> ias(data, m);
//	ias.compute(o);
//	cout << "There are " << o.size() << " alignment parts without partial overlap" << endl;
//	ifstream in("encode",std::ifstream::binary);
//	m.set_patterns(in);
//	for(map<string, vector<unsigned int> >::const_iterator it = m.get_high(0).begin(); it != m.get_high(0).end(); it++){
//		for(size_t k =0; k <5; k++){
//							cout<< it->second.at(k)<< " ; ";
//						}
//						cout<< " "<< endl;

//	}
//Data compression:
	cout<< "weight size: "<< weight.size()<<endl;
//	en.arithmetic_encoding_seq(outs);
	en.arithmetic_encoding_alignment(weight,membersOfCluster,alignments_in_a_cluster,outs);
//	test.encode();
	outs.close();
//	test.decode();
	ifstream in("encode",std::ifstream::binary);
	en.arithmetic_decoding_alignment(in);
//	en.arithmetic_decoding_seq();



	return 0;
}
int do_mc_model_seq(int argc, char * argv[]){
	string fastafile(argv[2]);
	string maffile(argv[3]);
//	size_t num_threads = 1;
	typedef mc_model use_model;
//	typedef initial_alignment_set<use_model> use_ias;
//	typedef affpro_clusters<use_model> use_affpro;

//Reading data:(Use an empty maf file!)
	all_data data(fastafile, maffile);
	wrapper wrapp;
//	counting_functor functor(data);
//	encoding_functor functor1(data);

ofstream outs("encode",std::ofstream::binary);
//Train all the sequences:
	use_model m(data);
	encoder en(data,m,wrapp);
	overlap ol(data);
	m.train(outs);
// Find connected components of alignments with some overlap
	compute_cc cccs(data);
	vector<set< const pw_alignment *, compare_pw_alignment> > ccs;
	cccs.compute(ccs);
//test
	en.arithmetic_encoding_seq(outs);
	en.arithmetic_decoding_seq();
//	en.arithmetic_encoding_alignment();
	
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
	}
	else if(0==program.compare("model"))
		return do_mc_model(argc, argv);
	
	 
/*	else{
		return do_mc_model_seq(argc,argv);
	}*/
	else {
		usage();
	}

	
	return 1;
}



#include "model.cpp"
#endif























