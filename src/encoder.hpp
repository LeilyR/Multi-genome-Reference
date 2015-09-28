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
	void partition_centers(std::map<std::string, vector<pw_alignment> > & , std::map<std::string, unsigned int>& );
	void partitioning_long_centers(std::vector<std::vector<std::string> > &, std::map<std::vector<std::string>, unsigned int > &);
	void calculate_high_in_partition(std::map<std::string, unsigned int> & , std::map<std::string,std::vector<pw_alignment> > & );
	void partition_high();
	void calculate_long_center_high(std::vector<std::vector<std::string> > &, std::map<std::vector<std::string>, unsigned int> &);
	void setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & );
	void add_partition_high_to_stream(std::ofstream & );
	void add_center_to_stream(std::ofstream &);
	void add_long_center_to_the_stream(std::ofstream &);//add both long centers info and their partition info to the stream
	void arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > & alignmentOfCluster,std::ofstream & ,dlib::entropy_encoder_kernel_1 & );
	void encode_flags_of_reverse_parts(pw_alignment & , unsigned int &,unsigned int&, dlib::entropy_encoder_kernel_1 &);
	void al_encoding(std::map<std::string, unsigned int> & , std::map<std::string, std::string > & , std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &, dlib::entropy_encoder_kernel_1 &);
	void al_encode_with_long_center(std::map<std::string, unsigned int> &, std::map<std::vector<std::string> , unsigned int> & ,std::map<std::string, std::vector<pw_alignment> > &, std::vector<std::map<size_t, std::vector<std::string> > > &, std::map<std::string, std::string > & , std::vector<std::vector<std::string> > &, std::ofstream & ,dlib::entropy_encoder_kernel_1 & );//TODO
	void encoding_seq_test();
	void write_to_stream(std::map<std::string, std::vector<pw_alignment> > & ,std::ofstream &);
	void size_calculator( size_t &, std::map<std::string, std::vector<pw_alignment> > &, std::map<std::string, std::string>&);
	private:
	all_data & data;
	T & model;
	wrapper& wrappers;
	std::map<std::string,std::vector<double> > all_the_patterns;
	std::map<std::string, unsigned int> total;
	std::map< size_t, std::vector<unsigned int> > lower_bound;
	std::vector<std::map< size_t, vector<unsigned int> > >upper_bound;
	std::vector<std::multimap<size_t , pw_alignment *> > AlignmentsFromClustering; //vector---> goes over all the references, size_t --> left of an alignment on the sequence
	std::map<size_t, std::vector<std::string> > acc_of_center;
	std::vector< std::map<std::string, vector<unsigned int> > >cluster_high_partition;//vector ---> partitions , map<> --> center of each partition with their high and low value.
	std::map<std::string,std::string> decoded_centers; // string --->dna sequence, string ---> ref:left(ino estefade nemikonam)
//	std::vector<std::string> centerId;//keep the centers in the order that has been decoded
	std::map<size_t , std::vector<std::string> > partition; // size_t ---> partition number, vector---> centers of a partition
	std::vector<unsigned int>partitionHigh; //high value of each partition.
	std::vector<std::map <std::string, std::string> >decoded_center_in_partition;//string ---> ref:left, string --->dna sequence
	
	//long centers: 
	std::map<size_t , std::vector<std::vector<std::string> > > long_center_partition; // size_t is the number of partion and the vector includes all the long centers in that partition.
	std::vector<std::map<std::vector<std::string> , vector<unsigned int> > > long_center_high; //vector ---> partitions , map<> --> long centers of each partition with their high and low value.
	std::vector<unsigned int>HighOfPartition; //high value of each partition for long centers.




};

template<typename T>
class decoder{
	public:
	decoder(all_data& , T &, decoding_wrapper &);
	~decoder();
	std::string associatedMember(std::string &  , size_t& , size_t& );
	void set_partition_high_from_stream(std::ifstream& );
	void set_center_from_stream (std::ifstream & );
	void set_acc_from_stream(std::ifstream &);
	void set_long_centers_from_stream(std::ifstream &);
	void arithmetic_decoding_seq();
	void arithmetic_decoding_seq_test(int t);
	void arithmetic_decoding_centers(std::ifstream &,dlib::entropy_decoder_kernel_1 & );
	void decode_flag_of_reverse_parts(std::ifstream & , dlib::entropy_decoder_kernel_1 & , std::ofstream &, bool& ,bool& ,bool&,unsigned int & );
	void al_decoding(std::ifstream & , dlib::entropy_decoder_kernel_1 &, std::ofstream &);
	void al_decode_with_long_center(std::ifstream & , dlib::entropy_decoder_kernel_1 &, std::ofstream &);
	void insert_character(size_t &, size_t &);
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

	//long center:
	std::vector<std::map<std::vector<std::string>, std::vector<unsigned int> > > long_centers_high_value; //vector ---> partitions, map->first ---> long_center, map ->second ---> its high and low
	std::vector<unsigned int> long_center_partitions_high_value; // Each row of the vector represents the high value of its correspondence partition 
};



#endif

