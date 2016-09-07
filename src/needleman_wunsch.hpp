#ifndef NEEDLEMAN_WUNSCH_HPP
#define NEEDLEMAN_WUNSCH_HPP

#include "pw_alignment.hpp"
#include "data.hpp"
#include "dynamic_mc.hpp"
#include "iostream"
#include <algorithm>
#include <map>
#include <vector>
#include <cassert>


template<typename T>
class needleman{
	public:
	//	needleman(all_data & , T & );
		needleman(all_data & , T  &, std::string &  , std::string &);
		~needleman();
		void compute_matrix(size_t & read_id, size_t & ref_id);
		std::string get_first_context();
		void count_first_row_delete(size_t &, char &, std::string & );
		std::string get_previous_contexts();
		void update_keep(size_t & , std::string & , std::string & , char &, std::string &);
		void count_keep(size_t & ,std::string & );
		void update_delete(size_t & , std::string & , std::string &  , char & , std::string &);
		void count_delete(size_t & ,std::string & );
		void add_minimums(size_t &  , size_t &  , std::vector<double> &  , std::vector<double> &  , std::vector<std::string> & , size_t &  , size_t & , bool & );
		void add_path(size_t &, size_t &, size_t &);
		void delete_cost(const std::map< std::string, std::vector<double>  >& , std::string & , std::string & , double &);
		void keep_cost(const std::map< std::string, std::vector<double>  >& , std::string &,  std::string & , double &);
		void find_the_best_path(size_t &, std::string &, std::string & );
		void print_score_matrix();
		void print_context_matrix();
		void print_path();

	private:
		T & model;
		const all_data & data;
		std::vector<std::vector<double> > mod_matrix;//TODO can be removed
		std::vector<std::vector<double> > score_matrix;
		std::vector<std::vector<std::string> > context_matrix;
		std::vector<std::vector<size_t> > keep_number;
		std::vector<std::vector<size_t> > delete_number;
		std::vector<std::vector<size_t> > path;//For each cell shows how we reached it. 1->match/mismatch(diagonal) 2->deletion(horizontal) 3->insertion(vertical), only the first cell gets zero.
		std::string read;
		std::string ref;
		

};

#endif


