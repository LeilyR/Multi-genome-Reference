#ifndef DYNAMIC_MC_HPP
#define DYNAMIC_MC_HPP

#include "data.hpp"


#define BITS_PER_SEQUENCE_HIGH_VALUE 12
#define MAX_SEQUENCE_MARKOV_CHAIN_LEVEL 8


class dynamic_mc_model{//TODO level 0;

	public:
		dynamic_mc_model(all_data & );
		~dynamic_mc_model();
		void train_sequence_model();	
//		void markov_chain(size_t &, size_t);
		void calculate_sequence_high(size_t & accession, const std::map<std::string, std::vector<size_t> > & counts);
		void make_all_patterns(const size_t &);
		void write_parameters(std::ofstream & );
		void set_patterns(std::ifstream &);
		void recursive_al_model();
		void make_all_alignment_pattern(const size_t & level);
		void calculate_alignment_high(size_t & , size_t &);
		void write_alignment_parameters(std::ofstream &);
		void set_alignment_pattern(std::ifstream &);
		void count_context_OneToTwo(const pw_alignment &  ,size_t , abstract_context_functor & )const;
		void count_context_TwoToOne(const pw_alignment &  ,size_t , abstract_context_functor & )const;
		void modification(char , int & , int & , int &, int & )const;
		char modification_character(int , int , int , int )const;
		void train(std::ofstream &);
		void cost_function(const pw_alignment& , double & , double & , double & , double & )const ;
		void gain_function(const pw_alignment& , double & , double & )const ;
		size_t get_acc_level(size_t&)const;//use it in mc model class as the level of model
		void set_acc_level(size_t &); //sets accessions level after reading them from the file
		size_t get_al_level(size_t& , size_t&)const;//use it in mc model class as the level of model
		void set_al_level(size_t& , size_t & , size_t & );
		const std::map<std::string, std::vector<unsigned int> > & get_high(size_t acc)const;
		std::string get_firstPattern(size_t &)const;
		std::string get_firstAlignmentPattern(size_t& , size_t& )const;
		void get_encoded_member(pw_alignment & , size_t center_ref, size_t center_left, encoding_functor & ,std::ofstream&)const;
		std::vector<unsigned int> get_high_at_position(size_t seq_index, size_t position) const;
		const std::map<std::string, std::vector<unsigned int> >& get_highValue(size_t acc1, size_t acc2)const;


	private:

		all_data & data;	
		std::vector<size_t>accLevel; // accession i uses a markov chain of that level, 0 means simple model without context
		std::vector<std::vector<double> > create_cost;//accession(cost of each base)

		std::vector<std::map<std::string, std::vector<double> > >sequence_successive_bases;//accession(string ---> context , vector<double> ---> probability of the next base)
		std::vector<std::vector<std::map <std::string, std::vector<double> > > >mod_cost; //accession1(accession2(context, next context))
		std::vector<size_t> powersOfTwo;
		std::set<std::string> all_alignment_patterns; // all possible alignment patterns
		std::map<std::string,std::vector<double> > all_the_patterns;
		std::vector<std::map<std::string, std::vector<unsigned int> > > high;//sequences patterns
		std::vector<std::vector<std::map<std::string , std::vector<unsigned int> > > >highValue;//alignments patterns
		std::vector<vector<size_t> >al_level;



		double total_sequence_cost(const std::map<std::string, std::vector<size_t> > & all_context_counts) const;
		void context_count(std::map<std::string, std::vector<size_t> > & countmap, const std::string & context, size_t base, size_t num) const;





};

#endif
