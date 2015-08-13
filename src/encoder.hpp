#ifndef ENCODER_HPP
#define ENCODER_HPP

#include "data.hpp"
#include "dynamic_mc.hpp"
#include "pw_alignment.hpp"
#include <iostream>
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"

#include <map>
#include <vector>
#include <cassert>


template<typename T>
class encoder{
	public:
	encoder(all_data& , T &, wrapper &);
	~encoder();
	void arithmetic_encoding_seq(std::ofstream&);
	void calculating_clusters_high(std::map<std::string, unsigned int> & weight);
	void encoding_seq_test();
	void setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & );
	void add_acc_to_stream(std::ofstream &);
	void set_acc_from_stream(std::ifstream &);
	void set_pattern_from_stream(std::ifstream & );
	const std::map<std::string, std::vector<unsigned int> > & get_center_high() const;
	void partition_centers(std::map<std::string, vector<pw_alignment> > & , std::map<std::string, unsigned int>& );
	void calculate_high_in_partition(std::map<std::string, unsigned int> & , std::map<std::string,std::vector<pw_alignment> > & );
	void partition_high();
	void add_center_to_stream(std::ofstream &);
	void set_center_from_stream (std::ifstream & );
	void arithmetic_encoding_centId(std::map<std::string, std::vector<pw_alignment> > &, std::ofstream &);
	void arithmetic_decoding_centId(std::ifstream&, dlib::entropy_decoder_kernel_1&);
	void arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster,std::ofstream & ,dlib::entropy_encoder_kernel_1 & );
	void arithmetic_decoding_centers(std::ifstream &,dlib::entropy_decoder_kernel_1 & );
	void add_partition_high_to_stream(std::ofstream & );
	void set_partition_high_from_stream(std::ifstream& );
	void al_encoding(std::map<std::string, unsigned int> & , std::map<std::string, std::string > & , std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &, dlib::entropy_encoder_kernel_1 &);
	void al_decoding(std::ifstream & , dlib::entropy_decoder_kernel_1 &, std::ofstream &);
	void write_to_stream(std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &);
	void read_from_stream(std::ifstream&);
	const std::multimap<size_t, pw_alignment*> & get_alignment(std::map<std::string,std::vector<pw_alignment> > & , size_t);
	void size_calculator( size_t &, std::map<std::string, std::vector<pw_alignment> > &, std::map<std::string, std::string>&);
	private:
	all_data & data;
	T & model;
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

template<typename T>
class decoder{
	public:
	decoder(all_data& , T &, decoding_wrapper &);
	~decoder();
	std::string associatedMember(std::string &  , size_t& , size_t& );
	void set_acc_from_stream(std::ifstream &);
	void set_pattern_from_stream(std::ifstream & );
	void set_center_from_stream (std::ifstream & );
	void arithmetic_decoding_seq();
	void arithmetic_decoding_seq_test(int t);
	void arithmetic_decoding_centId(std::ifstream&, dlib::entropy_decoder_kernel_1&);
	void arithmetic_decoding_centers(std::ifstream &,dlib::entropy_decoder_kernel_1 & );
	void set_partition_high_from_stream(std::ifstream& );
	void al_decoding(std::ifstream & , dlib::entropy_decoder_kernel_1 &, std::ofstream &);
	void read_from_stream(std::ifstream&);
	private:
	all_data & data;
	T & model;
	decoding_wrapper& wrappers;
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

