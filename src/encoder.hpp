#ifndef ENCODER_HPP
#define ENCODER_HPP

#include "data.hpp"
#include "model.hpp"
#include <iostream>
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"

#include <map>
#include <vector>

using namespace std;

class encoder{
	public:
	encoder( all_data& , mc_model&);
	~encoder();
	void arithmetic_encoding_seq();
	void arithmetic_decoding_seq_test(int t);
	void calculating_clusters_high(map<string, unsigned int> & weight);
	void add_high_to_stream();
	void arithmetic_encoding_alignment(map<string, string> & cluster_members, map<string, vector<pw_alignment> > & alignmentOfCluster);
	void arithmetic_decoding_alignment();
	void encoding_seq_test();
	void arithmetic_decoding_seq();
	void set_center_high(ifstream& in);

	private:
	all_data & data;
	mc_model & model;
	map<string,vector<double> > all_the_patterns;
	map<string, unsigned int> total;
	map< size_t, vector<unsigned int> > lower_bound;
	vector<map< size_t, vector<unsigned int> > >upper_bound;
	map<string, vector<unsigned int> >cluster_high;
	




};


#endif

