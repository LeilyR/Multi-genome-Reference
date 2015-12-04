#ifndef DYNAMIC_MC_HPP
#define DYNAMIC_MC_HPP

#include "data.hpp"


#define MIN_ENCODING_WIDTH_BITS 1
#define MAX_ENCODING_WIDTH_BITS 20
#define NUM_ENCODING_TYPES 40 // (max - min + 1)*2, 40 types encoded 0 to 39, 40 is used as flag for no model
#define MAX_SEQUENCE_MARKOV_CHAIN_LEVEL 5
#define MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL 3 // context is three modification instructions + current base on known strand
#define TOTAL_FOR_ENCODING 4294967295 // 2^32-1
#define ENCODING_MIN_WIDTH 3 // min difference between low and high value for arithmetic encoding
#define NUM_DELETE_DYN 5 // number of delete encodings, 2^0, ... , 2^4
#define NUM_KEEP_DYN 10 // number of keep encodings, 2^0, ..., 2^9
#define NUM_MODIFICATIONS 25 // NUM_DELETE_DYN + NUM_KEEP_DYN + 10 ( for change to ACTGN, insert ACTGN)




template<typename reader> 
class alignment_contexts {


	public:
	alignment_contexts(reader & r, size_t level);
	~alignment_contexts();

	void read_alignment_1_2(const pw_alignment & p);
	void read_alignment_2_1(const pw_alignment & p);

	private:
	reader & rd;
	size_t level;



	void next_modification_1_2(const pw_alignment & al, size_t at, char & modification, char & on_char, size_t & num_cols);
	void next_modification_2_1(const pw_alignment & al, size_t at, char & modification, char & on_char, size_t & num_cols);

};



class reader_counter {
public:
	inline void see_context(const std::string & context, size_t modification);

	std::map<std::string, std::vector<size_t> > countsmap;
};


class reader_costs {
public:
//	inline void see_context(const std::string & context, size_t modification);

	double cost_sum;

};


typedef alignment_contexts<reader_counter> context_counter;
typedef alignment_contexts<reader_costs> context_cost;

class dynamic_mc_model {

	public:
		dynamic_mc_model(all_data &, size_t num_threads = 1);
		~dynamic_mc_model();
		void train_sequence_model();
		void train_alignment_model();	
		void compute_sequence_model();
		void compute_alignment_model();

		void write_parameters(std::ofstream & ) const;
		void set_patterns(std::ifstream &);
	

	
		double get_alignment_base_cost(const size_t & acc1, const size_t & acc2) const;


		static char modification_character(int modify_base, int num_delete, int insert_base, int num_keep);
		static std::string tranlsate_modification_character(char enc);
		static size_t modification_length(char mod); 
		static void modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep);

		static void bits_write(const std::vector<bool> & bits, std::ofstream & out);
		static void bits_read(std::ifstream & in, std::vector<bool> & bits, size_t num_bits_to_read);

		static void ints_to_bits(const std::vector<uint32_t> & ints, std::vector<bool> & bits, size_t bits_per_int);
		static void bits_to_ints(const std::vector<bool> & bits, std::vector<uint32_t> & ints, size_t num_ints, size_t bits_per_int);


// TODO check all those functions
//		void markov_chain(size_t &, size_t);
		void calculate_sequence_high(size_t & accession, const std::map<std::string, std::vector<size_t> > & counts);
		void make_all_patterns(const size_t &);
		void recursive_al_model();
		void make_all_alignment_pattern(const size_t & level);
		void calculate_alignment_high(size_t & , size_t &);
		void write_alignment_parameters(std::ofstream &);
		void set_alignment_pattern(std::ifstream &);
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
		size_t num_threads;


// model parameters which are stored in encoded file:
		std::vector<std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > sequence_models; // acc -> context -> (model_type, bits)
		// low/high of data and flags for sequence data. In each accession we have base encodings from 0 to all_bases_total
		// alignment begin encodings from all_bases_total to alignment_flag
		// sequence end encodings from alignment_flag to TOTAL_FOR_ENCODING 	
		std::vector<uint32_t> acc_sequence_all_bases_total;
		std::vector<uint32_t> acc_sequence_alignment_flag;

		std::vector<std::vector<uint32_t> > alignments_all_modification_total;

		std::vector<std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > > alignment_models; // acc_from -> acc_to -> context -> (model_type, bits) 	
// model parameters which are computed from stored parameters:
		std::vector<std::map< std::string, std::vector<uint32_t> > > sequence_model_widths; // acc -> pattern (length is MAX_SEQUENCE_MARKOV_CHAIN_LEVEL  ) -> encoding high value per base 

		std::vector<std::vector<std::map<std::string, std::vector<uint32_t> > > > alignment_model_widths; //  acc-from -> acc-to -> pattern (length is MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL +1 ) -> encoding high value per base 


		std::vector<double> alignment_begin_cost; 
		std::vector<std::vector<double> > alignment_end_cost;
		double alignment_adress_cost; // estimated cluster center adressing cost
		std::vector<std::vector< double > > alignment_base_cost; // summary base cost




		// TODO delete	
		std::vector<size_t> sequence_level; // accession i uses a markov chain of that level, 0 means simple model without context
		std::vector< std::map<std::string, std::vector<uint32_t> > > sequence_high_values; // acc -> pattern (length in sequence_level) -> encoding high value per base 		
		std::vector<uint32_t> alignment_begin_high; // acc -> high value for al begin in that accession
		std::vector<uint32_t> sequence_end_high; // acc -> sequence end flag high value, this should be equal to total. Using this flag twice (encoding empty sequence) means end of accession
		 
		/*
			This is the function to convert statistical models (probability distribution function on a finite set)
			to encoding widths (distance between low and high value for arithmetic encoding

			We try to approximate the float distribution using integer ranges
			distri.at(i) =approx= encoding_widths.at(i)/total
			each width is >= min and <=max
			The sum of all widths is total

		*/
// good function
		static void distribution_to_sqroot_bits(const std::vector< double > & distri, size_t bits, std::vector<uint32_t> & sqbits); 
		static void distribution_to_bits(const std::vector<double> & distri, size_t bits, std::vector<uint32_t> & vbits);
		static void sqroot_bits_to_width_encoding(const std::vector<uint32_t> & sqbits, uint32_t target_total, std::vector<uint32_t> & widths);
		static void bits_to_width_encoding(const std::vector<uint32_t> & bits, uint32_t target_total, std::vector<uint32_t> & widths);

		static unsigned char get_model_index(bool use_sqroot, size_t num_bits);
		static void get_model_type(unsigned char model_index, bool & use_sqroot, size_t & num_bits);

		static void counts_to_distribution(const std::vector<size_t> & counts, std::vector<double> & distri);
		static double counts_to_bits_eval(const std::vector<size_t> & counts, uint32_t target_total, unsigned char & model_index,  std::vector<uint32_t> & bits);
		static double count_cost_eval(const std::vector<size_t> & counts, const std::vector<uint32_t> & enc_widths);

		static double context_selector(const std::map< std::string, std::vector<size_t> > & context_counts, const std::vector<unsigned char> & alphabet, uint32_t target_total, 
			std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & cpref_bits, size_t max_length );
		static void context_subtree_cost(const std::string & cur_context, const std::map< std::string, std::vector<size_t> > & context_counts, const std::vector<unsigned char> & alphabet, size_t max_length, 
			std::map< std::string, std::vector<size_t> > & sum_counts);

		static void model_bits_to_full_high_values(const std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & bitmodel, uint32_t target_total, size_t target_pattern_length, 
			const std::vector<unsigned char> & alphabet, std::map< std::string, std::vector<uint32_t> > & hv_results);

		static void alignment_counts_add(std::map<std::string, std::vector<size_t> > & countsmap, const std::map<std::string, std::vector<size_t> > & newcounts);

// TODO Static
		void read_paramters(std::ifstream & in, 	std::vector<std::string> & acc_names, 
			std::vector<std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > & sequence_models,
			std::vector<uint32_t> & acc_sequence_all_bases_total,
			std::vector<uint32_t> & acc_sequence_alignment_flag,
			std::vector<std::vector<uint32_t> > & alignments_all_modification_total,
			std::vector<std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > > & alignment_models) const;

		static void check_read_error(std::ifstream & in);

		static void get_sequence_alphabet(std::vector<unsigned char> & alphabet);
		static void get_alignment_alphabet(std::vector<unsigned char> & alphabet);


		void check_parameters_in_file(const std::string & file) const;
		static void mapmodel_write(std::ofstream & outs, const std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & mapmodel, const std::vector<unsigned char> & alphabet, size_t max_level);
		static void mapmodel_read(std::ifstream & in,  std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & mapmodel, const std::vector<unsigned char> & alphabet, size_t max_level);

/*
		std::vector<std::map<std::string, std::vector<double> > >sequence_successive_bases;//accession(string ---> context , vector<double> ---> probability of the next base)
		std::vector<std::vector<std::map <std::string, std::vector<double> > > >mod_cost; //accession1(accession2(context, next context))
		std::vector<size_t> powersOfTwo;
		std::set<std::string> all_alignment_patterns; // all possible alignment patterns
		std::map<std::string,std::vector<double> > all_the_patterns;
		std::vector<std::map<std::string, std::vector<unsigned int> > > high;//sequences patterns
		std::vector<std::vector<std::map<std::string , std::vector<unsigned int> > > >highValue;//alignments patterns
		std::vector<vector<size_t> >al_level;

*/

		double total_sequence_cost(const std::map<std::string, std::vector<size_t> > & all_context_counts) const;
		void context_count(std::map<std::string, std::vector<size_t> > & countmap, const std::string & context, size_t base, size_t num) const;
		
		static void binary_write(std::ofstream & out, size_t v);
		static void binary_write(std::ofstream & out, uint32_t v);
		static void binary_write(std::ofstream & out, const std::string & str);
		static void binary_read(std::ifstream & in, size_t & v);
		static void binary_read(std::ifstream & in, uint32_t & v);
		static void binary_read(std::ifstream & in, std::string & str);
		static void binary_read(std::ifstream & in, unsigned char & c);
		

};

#endif
