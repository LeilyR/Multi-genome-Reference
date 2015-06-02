#ifndef ENCODER_HPP
#define ENCODER_HPP

#include "data.hpp"
#include "model.hpp"
#include <iostream>
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"

#include <map>
#include <vector>


class encoder{
	public:
	encoder( all_data& , mc_model&, wrapper &);
	~encoder();
	void arithmetic_encoding_seq(std::ofstream&);
	void arithmetic_decoding_seq_test(int t);
	void calculating_clusters_high(std::map<std::string, unsigned int> & weight);
	void arithmetic_encoding_alignment(std::map<std::string, unsigned int> & weight, std::map<std::string, std::string> & cluster_members, std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster,std::ofstream&,dlib::entropy_encoder_kernel_1 &);
	void arithmetic_decoding_alignment(std::ifstream&,dlib::entropy_decoder_kernel_1&);
	void arithmetic_enc_centers(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream &);
	void arithmetic_dec_centers(std::ifstream&, dlib::entropy_decoder_kernel_1&);
	void encoding_seq_test();
	void arithmetic_decoding_seq();
	void setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & );
	void add_acc_to_stream(std::ofstream &);
	void set_acc_from_stream(std::ifstream &);
	void set_pattern_from_stream(std::ifstream & );
	void arithmetic_enc_centId(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster, std::ofstream &);
	void arithmetic_dec_centId(std::ifstream &, dlib::entropy_decoder_kernel_1&);
	const std::map<std::string, std::vector<unsigned int> > & get_center_high() const;
	void partition_centers(std::map<std::string, std::vector<pw_alignment> > & , std::map<std::string, unsigned int>&);
	void calculate_high_in_partition(std::map<std::string, unsigned int> &, std::map<std::string,std::vector<pw_alignment> > &);
	void partition_high();
	void add_center_to_stream(std::ofstream & );
	void set_center_from_stream (std::ifstream &);
	void arithmetic_encoding_centId(std::map<std::string, std::vector<pw_alignment> > &, std::ofstream &);
	void arithmetic_decoding_centId(std::ifstream&, dlib::entropy_decoder_kernel_1&);
	void arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > &,std::ofstream &,dlib::entropy_encoder_kernel_1 & );
	void arithmetic_decoding_centers(std::ifstream &,dlib::entropy_decoder_kernel_1 &);
	void add_partition_high_to_stream(std::ofstream & );
	void set_partition_high_from_stream(std::ifstream&);
	void test_al_encoding(std::map<std::string, unsigned int> & , std::map<std::string, std::string > & , std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &, dlib::entropy_encoder_kernel_1 &);
	void test_al_decoding(std::ifstream &, dlib::entropy_decoder_kernel_1 &);
	void test_al_decoding( std::string center,size_t part);
	void write_to_stream(std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &);
	void read_from_stream(std::ifstream&);
	std::string associatedMember(std::string & , size_t&, size_t&);
	private:
	all_data & data;
	mc_model & model;
	wrapper& wrappers;
	std::map<std::string,std::vector<double> > all_the_patterns;
	std::map<std::string, unsigned int> total;
	std::map< size_t, std::vector<unsigned int> > lower_bound;
	std::vector<std::map< size_t, vector<unsigned int> > >upper_bound;
	std::map<std::string, std::vector<unsigned int> >cluster_high;
	std::vector<std::multimap<size_t , pw_alignment *> > AlignmentsFromClustering; //vector---> goes over all the references, size_t --> left of an alignment on the sequence
	std::map<size_t, std::vector<std::string> > acc_of_center;
	std::vector< std::map<std::string, vector<unsigned int> > >cluster_high_partition;//vector ---> partitions , map<> --> center of each partition with their high and low value.
	std::map<std::string,std::string> decoded_centers; // string --->dna sequence, string ---> ref:left(ino estefade nemikonam)
//	std::vector<std::string> centerId;//keep the centers in the order that has been decoded
	std::map<size_t , std::vector<std::string> > partition; // size_t ---> partition number, vector---> centers of a partition
	std::vector<unsigned int>partitionHigh; //high value of each partition.
	std::vector<std::map <std::string, std::string> >decoded_center_in_partition;//string ---> ref:left, string --->dna sequence



};



#endif

