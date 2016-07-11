
#ifndef MERGING_CENTERS_HPP
#define MERGING_CENTERS_HPP

#include "pw_alignment.hpp"
#include"data.hpp"
#include"model.hpp"
#include "suffix_tree.hpp"

#include <map>
#include <vector>
#include <cassert>
#include <omp.h>


class merging_centers{
	public:
		merging_centers(all_data & d, finding_centers & cent):data(d),centers(cent), all_current_centers(data.numSequences()){
		
		}
		~merging_centers(){

		}
		void updating_centers(std::vector<size_t> & , size_t &);
		void merg_gain_value( suffixTree & );
		void find_highest_gain(std::vector<size_t> & , int &, std::vector<size_t> & );
		void adding_new_centers(std::vector<std::vector<std::string> >&, std::vector<std::map<size_t , std::vector<std::string> > >&);// Saves new centers in the 'long_centers' vector in main. It gives a unique id to each of them(Their row in the outer vector). vector<string> includes are the centers in the new one. they are saved in the form of ref:left.
		void index_centers(std::map<std::string , std::vector<pw_alignment> > & );
		void add_long_centers_to_map(std::vector<std::vector<std::string> > & , std::map<vector<std::string>, std::vector<pw_alignment> > &); //Initially fills in the new_centers map
		void create_alignment(std::vector<std::vector<std::string> > &, std::map<vector<std::string>, vector<pw_alignment> >&, std::map<std::string , std::vector<pw_alignment> > &,  std::vector<std::map<size_t , std::vector<std::string> > > &,  std::vector<std::map<size_t , std::string> > & ); //Creates all the als using long centers 
		void remove_fully_reverse_refs(vector<vector<std::string> > & ,std::map<std::vector<std::string>, std::vector<pw_alignment> > &);
		void find_new_centers(size_t &, std::vector<std::string> & , size_t &, std::vector<std::map<size_t, std::vector<std::string> > >& );//Finds long centers on a certain sequence	
		void find_long_center_length(std::vector<std::string> & , std::map<std::string,std::vector<pw_alignment> > & , size_t & ,size_t & ,size_t &,size_t &, std::vector<std::map<size_t , std::string> > &);

	private:
		all_data & data;
		finding_centers & centers;
		std::map<std::string, int> center_index;
		std::map<std::vector<size_t> , size_t> merged_centers;//merged centers, vector<size_t> ----> original indices of megerd centers , size_t ----> its new index
		std::map<std::vector<size_t>, int >gains; // vector<size_t> ----> series of centers, int ----> gain of it. It includes the last combination of centers
		std::vector<std::vector<std::vector<size_t> > >all_current_centers;
};

#endif
