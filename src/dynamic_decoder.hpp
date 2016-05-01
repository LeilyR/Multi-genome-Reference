#ifndef DYNAMIC_DECODER_HPP
#define DYNAMIC_DECODER_HPP

//#include "data.hpp"
#include "dynamic_mc.hpp"
#include <iostream>
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"

#include <map>
#include <vector>
#include <cassert>

template<typename T>
class dynamic_decoder{
	public:
	dynamic_decoder(T &, decoding_wrapper &);
	~dynamic_decoder();
	void set_center_from_stream(std::ifstream & in);
	void set_partition_high_from_stream(std::ifstream & in);
	void set_center_flags_from_the_stream(std::ifstream & in);
	std::string first_sequence_pattern()const;
	std::string first_alignment_pattern()const;
	void arithmetic_decoding_centers(std::ifstream & in, dlib::entropy_decoder_kernel_1 & dec);
	void decoding_directions(std::ifstream &, dlib::entropy_decoder_kernel_1 &, bool&, bool&, bool&, bool&, uint32_t & );
	std::string associatedMember(std::string & center, size_t & modificationPattern, size_t & position);
	void al_decoding(std::ifstream & , dlib::entropy_decoder_kernel_1 &, std::ofstream &);


//Long centers
	void set_long_centers_from_stream(std::ifstream & );
	void set_long_center_partition_high_from_stream(std::ifstream& );
	void decoding_long_centers(std::ifstream &, dlib::entropy_decoder_kernel_1 & );
	void al_decode_with_long_center(std::ifstream & , dlib::entropy_decoder_kernel_1 &, std::ofstream &);


	private: 
	T & model;
	decoding_wrapper& wrappers;
	std::vector< std::vector<std::pair<size_t, vector<unsigned int> > > >cluster_high_partition;// size_t is center id , vector<unsigned int> is its boundries
	std::map<size_t , std::vector<size_t> > acc_of_center;
	std::vector<unsigned int>partitionHigh; //high value of each partition
	std::vector<uint32_t> center_flags;
	std::vector<std::map <size_t, std::string> >decoded_center_in_partition;

//Long centers
	std::map<size_t, std::vector<int> > acc_of_long_center;
	std::vector<std::map <int, std::string> >decoded_long_center_in_partition;//string ---> ref:left, string --->dna sequence
	std::vector<std::vector<std::pair<std::vector<int>, std::vector<unsigned int> > > >long_centers_high_value; //vector ---> partitions, map->first ---> long_center, map ->second ---> its high and low
	std::vector<unsigned int> long_center_partitions_high_value; // Each row of the vector represents the high value of its correspondence partition 
	std::vector<unsigned int> flag_value;




};

#endif
