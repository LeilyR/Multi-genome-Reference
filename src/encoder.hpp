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
	encoder( all_data& , mc_model&, wrapper &);
	~encoder();
	void arithmetic_encoding_seq(ofstream&);
	void arithmetic_decoding_seq_test(int t);
	void calculating_clusters_high(map<string, unsigned int> & weight);
	void arithmetic_encoding_alignment(map<string, unsigned int> & weight, map<string, string > & cluster_members, map<string, vector<pw_alignment> > & alignmentOfCluster,ofstream&);
	void arithmetic_decoding_alignment(ifstream&);
	void arithmetic_enc_centers(map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream &);
	void arithmetic_dec_centers(ifstream&, dlib::entropy_decoder_kernel_1&);
	void encoding_seq_test();
	void arithmetic_decoding_seq();
	void setOfAlignments(map<string,vector<pw_alignment> > & );
	void add_acc_to_stream(ofstream &);
	void set_acc_from_stream(ifstream &);
	void set_pattern_from_stream(ifstream & );
	void arithmetic_enc_centId(map<string, vector<pw_alignment> > & alignmentOfCluster, ofstream &);
	void arithmetic_dec_centId(ifstream &, dlib::entropy_decoder_kernel_1&);
	const map<string, vector<unsigned int> > & get_center_high() const;
	void partition_centers(map<string, vector<pw_alignment> > & , map<string, unsigned int>&);
	void calculate_high_in_partition(map<string, unsigned int> &, map<string,vector<pw_alignment> > &);
	void partition_high();
	void add_center_to_stream(ofstream & );
	void set_center_from_stream (ifstream &);
	void arithmetic_encoding_centId(map<string, vector<pw_alignment> > &, ofstream &);
	void arithmetic_decoding_centId(ifstream&, dlib::entropy_decoder_kernel_1&);
	void arithmetic_encoding_centers(map<string, vector<pw_alignment> > &,	dlib::entropy_encoder_kernel_1 &,ofstream &);
	void arithmetic_decoding_centers(ifstream &,dlib::entropy_decoder_kernel_1 & );
	void add_partition_high_to_stream(ofstream & );
	void set_partition_high_from_stream(ifstream&);
	private:
	all_data & data;
	mc_model & model;
	wrapper& wrappers;
	map<string,vector<double> > all_the_patterns;
	map<string, unsigned int> total;
	map< size_t, vector<unsigned int> > lower_bound;
	vector<map< size_t, vector<unsigned int> > >upper_bound;
	map<string, vector<unsigned int> >cluster_high;
	vector<multimap<size_t , pw_alignment *> > AlignmentsFromClustering; //vector---> goes over all the references, size_t --> left of an alignment on the sequence
	map<size_t, vector<string> > acc_of_center;
	vector< map<string, vector<unsigned int> > >cluster_high_partition;//vector ---> partitions , map<> --> center of each partition with their high and low value.
//	map<size_t, vector<string> > first_pattern_after_al; // size_t -----> sequence id, string -----> is the last pattern is covered by an alignment during the encoding.
	//plan is saving the right on the stream and fill in this map in decoding. we need those right values for the first base after al.
	map<string,string> decoded_centers; // string --->dna sequence, string ---> ref:left(ino estefade nemikonam)
//	vector<string> centerId;//keep the centers in the order that has been decoded
	map<size_t , vector<string> > partition; // size_t ---> partition number, vector---> centers of a partition
	vector<unsigned int>partitionHigh; //high value of each partition.
	vector<map <string, string> >decoded_center_in_partition;//string ---> ref:left, string --->dna sequence



};



#endif

