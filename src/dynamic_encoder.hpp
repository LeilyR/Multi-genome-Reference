#ifndef DYNAMIC_ENCODER_HPP
#define DYNAMIC_ENCODER_HPP

#include "data.hpp"
#include "pw_alignment.hpp"
#include "dynamic_mc.hpp"
#include <iostream>
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"

#include <map>
#include <vector>
#include <cassert>


template<typename T>
class dynamic_encoder{
	public:
	dynamic_encoder(all_data& , T &, wrapper &);
	~dynamic_encoder();

//Original centers which are obtained from clustering:
	void partition_centers(std::map<std::string,std::vector<pw_alignment> > & , std::map<std::string, unsigned int> & );
	void calculate_high_in_partition(std::map<std::string,std::vector<pw_alignment> > & , std::map<std::string, unsigned int> & );
	void partition_high();
	void add_center_to_the_stream(std::ofstream &);
	void add_partition_high_to_the_stream(std::ofstream &);
	void setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & );
	void calculate_end_of_center_flags(std::map<std::string, std::vector<pw_alignment> > &, std::vector<uint32_t> & widths,std::ofstream & outs);
	void add_center_flags_to_the_stream(std::vector<uint32_t> & bits , unsigned char & model_index, std::ofstream &);
	void arithmetic_encoding_centers(std::map<std::string, std::vector<pw_alignment> > & , std::ofstream & ,dlib::entropy_encoder_kernel_1 & );
	void encoding_directions(unsigned int & cent_ref, unsigned int & cent_left, const pw_alignment & p, std::ofstream &, dlib::entropy_encoder_kernel_1 & );
	void al_encode(std::map<std::string, unsigned int> &, std::map<std::string, std::string > &, std::map<std::string, std::vector<pw_alignment> > &, std::ofstream &, dlib::entropy_encoder_kernel_1 & );

//Long centers which are obtained from collapsing the original one:
	void partitioning_long_centers(std::map<std::vector<std::string>, std::vector<pw_alignment> > & , std::map<std::string ,std::vector<pw_alignment> > & , std::map<std::vector<std::string>, unsigned int> & );
	void calculate_long_center_high(std::map<std::vector<std::string>, std::vector<pw_alignment> > &,  std::map<std::string ,std::vector<pw_alignment> > & , std::map<std::vector<std::string>, unsigned int> &);
	void add_long_center_to_the_stream(std::ofstream &);
	void find_short_center(const pw_alignment* , std::map<std::string, std::string > & , std::vector<std::string> & , size_t & , size_t & );
	void encoding_long_center(std::map<std::string,std::vector<pw_alignment> > &, std::ofstream &, dlib::entropy_encoder_kernel_1 &);
	void make_fully_reverse_center(std::vector<std::string> &, std::vector<std::string> & );
	void al_encode_with_long_center(std::vector<std::multimap<size_t , std::string> >&,std::map<std::vector<std::string> , unsigned int> & ,std::map<std::string, std::vector<pw_alignment> > &, std::vector<std::map<size_t, std::vector<std::string> > > &, std::map<std::string, std::string > &, std::ofstream & ,dlib::entropy_encoder_kernel_1 &,std::map<std::vector<std::string>, std::vector<pw_alignment> >& );


	private:
	all_data & data;
	T & model;
	wrapper& wrappers;

	std::map<size_t , std::vector<std::string> > partition; // size_t ---> partition number, vector---> centers of a partition
	std::vector<unsigned int>partitionHigh; //high value of each partition.
	std::vector< std::map<std::string, vector<unsigned int> > >cluster_high_partition;//vector ---> partitions , string --> centers of each partition with their high and low value.
	std::vector<std::multimap<size_t , pw_alignment*> > AlignmentsFromClustering;

//Long centers
	std::vector<std::vector<std::pair<std::vector<std::string> , vector<unsigned int> > > > long_center_high; //vector ---> partitions , vector<> --> long centers of each partition with their high and low value.
	std::map<size_t , std::vector<std::vector<std::string> > > long_center_partition; // size_t is the number of partion and the vector includes all the long centers nad original centers in that partition.
	std::vector<unsigned int>HighOfPartition; //high value of each partition for long centers.






};
#endif

