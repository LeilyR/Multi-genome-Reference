#ifndef DYNAMIC_MC_HPP
#define DYNAMIC_MC_HPP

#include "data.hpp"
#include "pw_alignment.hpp"

#define MIN_ENCODING_WIDTH_BITS 1
#define MAX_ENCODING_WIDTH_BITS 20
#define NUM_ENCODING_TYPES 40 // (max - min + 1)*2, 40 types encoded 0 to 39, 40 is used as flag for no model
#define MAX_SEQUENCE_MARKOV_CHAIN_LEVEL 5
#define MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL 2 // context is two modification instructions + current base on known strand
#define TOTAL_FOR_ENCODING 16384 //2^14 
#define ENCODING_MIN_WIDTH 3 // min difference between low and high value for arithmetic encoding
#define NUM_DELETE_DYN 6 // number of delete encodings, 2^0, ... , 2^4, 2^5
#define NUM_KEEP_DYN 12 // number of keep encodings, 2^0, ..., 2^9, 2^10,2^11
#define NUM_MODIFICATIONS 28 // NUM_DELETE_DYN + NUM_KEEP_DYN + 10 ( for change to ACTGN, insert ACTGN)

// TODO distance between members of a long center: base_cost / (modification cost per base in al - creation cost per base   )

template<typename reader>
class sequence_contexts {
	public:
	sequence_contexts(reader & r, size_t level);
	~sequence_contexts();

	void read_sequence(const std::string & d);
	void read_sequence(const dnastring & d, const size_t & from, const size_t & to);

	private:
	reader & rd;
	const size_t level;
	std::string astring;

};


template<typename reader> 
class alignment_contexts {


	public:
	alignment_contexts(reader & r, size_t level);
	~alignment_contexts();

	void read_alignment_1_2(const pw_alignment & p);
	void read_alignment_2_1(const pw_alignment & p);
	void read_alignment_1_2(const pw_alignment & p, wrapper & wrap, std::ofstream &);
	void read_alignment_2_1(const pw_alignment & p, wrapper & wrap, std::ofstream &);

	void read_alignment_1_2(const pw_alignment & p, const size_t & from, const size_t & to);
	void read_alignment_2_1(const pw_alignment & p, const size_t & from, const size_t & to);

	private:
	reader & rd;
	const size_t level;



	void next_modification_1_2(const pw_alignment & al, size_t at, char & modification, char & on_char, size_t & num_cols);
	void next_modification_2_1(const pw_alignment & al, size_t at, char & modification, char & on_char, size_t & num_cols);

};



class a_reader_counter {
public:
	inline void see_context(const std::string & context, size_t modification);

	std::map<std::string, std::vector<size_t> > countsmap;
};


class a_reader_costs {
public:
	a_reader_costs(const std::map<std::string, std::vector<double> > & mcost);
	~a_reader_costs();
	inline void see_context(const std::string & context, size_t modification);

	double cost_sum;
	const std::map<std::string, std::vector<double> > & model_costs;

};

class s_reader_counter {
public:
	s_reader_counter();
	~s_reader_counter();

	inline void see_context(const std::string & context, size_t modification);

	std::map<std::string, std::vector<size_t> > countsmap;
	
private:
	const size_t numchars;

};

class s_reader_costs {
public:
	s_reader_costs(const std::map<std::string, std::vector<double> > & mcost);
	~s_reader_costs();
	inline void see_context(const std::string & context, size_t modification);

	double cost_sum;
	const std::map<std::string, std::vector<double> > & model_costs;

};

class a_reader_encode{
public:
	a_reader_encode(const std::map<std::string, std::vector<uint32_t> > & mhigh);
	~a_reader_encode();
	inline void see_context(const std::string & context , size_t modification);

	const std::map<std::string, std::vector<uint32_t> > & model_high;
	std::vector<std::vector<unsigned int> > low_high_values;

};

typedef alignment_contexts<a_reader_encode> acontext_encode;
typedef alignment_contexts<a_reader_counter> acontext_counter;
typedef alignment_contexts<a_reader_costs> acontext_cost;
typedef sequence_contexts<s_reader_counter> scontext_counter;
typedef sequence_contexts<s_reader_costs> scontext_cost;
// TODO average costs for modify, insert, delete
class dynamic_mc_model {
	public:
		dynamic_mc_model(all_data &, wrapper &, size_t num_threads = 1);
		~dynamic_mc_model();

		void train(std::ofstream &);



		void write_parameters(std::ofstream & ) const;
		void set_patterns(std::ifstream &);
	

		// XXX attention the base cost is now included in modify cost / gain cost	
		double get_alignment_base_cost(const size_t & acc1, const size_t & acc2) const;


		static char modification_character(int modify_base, int num_delete, int insert_base, int num_keep);
		static std::string translate_modification_character(char enc);
		static size_t modification_length(char mod); 
		static void modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep);


		void cost_function(const pw_alignment& , double & , double & , double & , double & )const;
		void cost_function(const pw_alignment& p , double & m1 , double & m2, size_t & acc1 , size_t acc2)const;

		inline	void gain_function(const pw_alignment& p, double & g1, double & g2)const{
			double c1;
			double c2;
			double m1;
			double m2;
			cost_function(p, c1, c2, m1,m2);
			g1 = c2 - m1;
			g2 = c1 - m2;
		}
		double get_the_gain(const pw_alignment&, std::string &)const;
		const std::vector<uint32_t> & get_seq_high_at_position(size_t & ref, size_t & position)const;
		const std::vector<uint32_t> & get_center_high_at_position(size_t & ref, size_t & left , size_t & position)const;
		void calculate_center_flags(vector<size_t> & counts , uint32_t & target_total , unsigned char & model_index, std::vector<uint32_t> & bits, std::vector<uint32_t> & widths);
		void calculate_center_high_values(uint32_t & target_total);
		void get_seq_flag(size_t & acc, std::vector<uint32_t> & al_begin , std::vector<uint32_t> & seq_acc_end)const;
		void get_end_al_flag(size_t& acc1, size_t & acc2, std::vector<uint32_t> & al_end);
		void get_center_flag(std::vector<uint32_t> & bits, unsigned char model_index, uint32_t & target_total, std::vector<uint32_t> & widths)const;
		const std::map< std::string, std::vector<uint32_t>  > & get_sequence_model_highs(size_t & acc)const;
		const std::map< std::string, std::vector<uint32_t>  > & get_center_model_highs(size_t & acc)const;
		const std::map< std::string, std::vector<uint32_t>  > & get_alignment_model_highs(size_t & acc1, size_t & acc2)const;
		void arith_encode_al(const pw_alignment & p, unsigned int & cent_ref, unsigned int & cent_left, std::vector<std::vector< unsigned int> > & low_high, wrapper & wrap, std::ofstream & encout);
		void arith_encode_long_al(const pw_alignment & p, size_t & acc1, size_t & acc2, unsigned int & cent_ref, unsigned int & cent_left, std::vector<std::vector< unsigned int> > & low_high);
		void write_al_high_onstream(std::ofstream & al_high_out);
		size_t get_acc()const;
		const std::map<std::string, std::vector<double>  > get_al_cost(size_t & acc1, size_t & acc2)const; 
		size_t estimate_gap_in_long_centers()const;



	


	private:

		all_data & data;
		wrapper& wrappers;
		size_t num_threads;


// model parameters which are stored in encoded file:
		std::vector<std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > sequence_models; // acc -> context (of all lengths) -> (model_type, bits)
		// low/high of data and flags for sequence data. In each accession we have base encodings from 0 to all_bases_total
		// alignment begin encodings from all_bases_total to alignment_flag
		// sequence end encodings from alignment_flag to TOTAL_FOR_ENCODING 	
		std::vector<uint32_t> acc_sequence_all_bases_total;
		std::vector<uint32_t> acc_sequence_alignment_flag;

		std::vector<std::vector<uint32_t> > alignments_all_modification_total;// End of alignment flag ---> acc1(acc2)

		std::vector<std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > > alignment_models; // acc_from -> acc_to -> context (variable length) -> (model_type, bits) 	
// model parameters which are computed from stored parameters:
		std::vector<std::map< std::string, std::vector<uint32_t>  > > sequence_model_highs; // acc -> pattern (length is MAX_SEQUENCE_MARKOV_CHAIN_LEVEL  ) -> encoding high value per base 
		std::vector<std::map< std::string, std::vector<uint32_t>  > > center_model_highs; 
		std::vector<std::map< std::string, std::vector<double>  > > sequence_model_costs;

		std::vector<std::vector<std::map<std::string, std::vector<uint32_t>  > > > alignment_model_highs; //  acc-from -> acc-to -> pattern (length is MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL +1 ) -> encoding high value per base 
		std::vector<std::vector<std::map<std::string, std::vector<double>  > > > alignment_model_costs;

		
		std::vector<double> alignment_begin_cost; 
		std::vector<std::vector<double> > alignment_end_cost;
		double alignment_adress_cost; // estimated cluster center adressing cost
		std::vector<std::vector< double > > alignment_base_cost; // summary base cost


// TODO to compute global base/subst cost we need count data and model. as we dont need these for decoding they are not stored in the compressed file. They are needed for estimating gap_in_long_centers
		double substitution_cost(const std::map<std::string, std::pair<unsigned char, std::vector<uint32_t> > > & mod_model, std::map<std::string, std::vector<size_t> > & counts, const size_t & target_total, size_t & numsubst) const;
		double global_base_cost;
		double global_substitution_cost;

		 
		void train_sequence_model();
		void train_alignment_model();	
		void compute_sequence_model();
		void compute_alignment_model();
		void compute_sequence_model_decoding();
		void compute_alignment_model_decoding();


		/*
			This is the function to convert statistical models (probability distribution function on a finite set)
			to encoding widths (distance between low and high value for arithmetic encoding

			We try to approximate the float distribution using integer ranges
			distri.at(i) =approx= encoding_widths.at(i)/total
			each width is >= min and <=max
			The sum of all widths is total

		*/
		static double counts_to_bits_eval(const std::vector<size_t> & counts, uint32_t target_total, unsigned char & model_index, std::vector<uint32_t> & bits);
		static void distribution_to_sqroot_bits(const std::vector< double > & distri, size_t bits, std::vector<uint32_t> & sqbits); 
		static void distribution_to_bits(const std::vector<double> & distri, size_t bits, std::vector<uint32_t> & vbits);
		static void sqroot_bits_to_width_encoding(const std::vector<uint32_t> & sqbits, uint32_t target_total, std::vector<uint32_t> & widths);
		static void bits_to_width_encoding(const std::vector<uint32_t> & bits, uint32_t target_total, std::vector<uint32_t> & widths);
		static double distri_to_bits(const std::vector<double> & distri, unsigned char model_index, std::vector<uint32_t> & bits);

		static unsigned char get_model_index(bool use_sqroot, size_t num_bits);
		static void get_model_type(unsigned char model_index, bool & use_sqroot, size_t & num_bits);

		static void counts_to_distribution(const std::vector<size_t> & counts, std::vector<double> & distri);
		static double count_cost_eval(const std::vector<size_t> & counts, const std::vector<uint32_t> & enc_widths);
		static void model_to_costs(const std::vector<uint32_t> bits, unsigned char model_index, uint32_t target_total, std::vector<double> & costs);
		static void width_correct(std::vector<uint32_t> & widths, uint32_t target_total);

		static double context_selector(const std::map< std::string, std::vector<size_t> > & context_counts, const std::vector<unsigned char> & alphabet, uint32_t target_total, 
			std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & cpref_bits, size_t max_length );
		static void context_subtree_cost(const std::string & cur_context, const std::map< std::string, std::vector<size_t> > & context_counts, const std::vector<unsigned char> & alphabet, size_t max_length, 
			std::map< std::string, std::vector<size_t> > & sum_counts);
		static void model_to_widths(const std::vector<uint32_t> bits, unsigned char model_index, uint32_t target_total, std::vector<uint32_t> & widths);
		static void model_bits_to_full_high_values(const std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & bitmodel, uint32_t target_total, size_t target_pattern_length, 
			const std::vector<unsigned char> & alphabet, std::map< std::string, std::vector<uint32_t> > & hv_results, std::map<std::string, std::vector<double> > & cost_results);

		static void model_widths_to_highs(const std::vector<uint32_t> & widths, uint32_t target_total, std::vector<uint32_t> & high_values);

		static void alignment_counts_add(std::map<std::string, std::vector<size_t> > & countsmap, const std::map<std::string, std::vector<size_t> > & newcounts);
		
		static void bits_write(const std::vector<bool> & bits, std::ofstream & out);
		static void bits_read(std::ifstream & in, std::vector<bool> & bits, size_t num_bits_to_read);

		static void ints_to_bits(const std::vector<uint32_t> & ints, std::vector<bool> & bits, size_t bits_per_int);
		static void bits_to_ints(const std::vector<bool> & bits, std::vector<uint32_t> & ints, size_t num_ints, size_t bits_per_int);

		static void read_paramters(std::ifstream & in, 	std::vector<std::string> & acc_names, 
			std::vector<std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > & sequence_models,
			std::vector<uint32_t> & acc_sequence_all_bases_total,
			std::vector<uint32_t> & acc_sequence_alignment_flag,
			std::vector<std::vector<uint32_t> > & alignments_all_modification_total,
			std::vector<std::vector< std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > > > & alignment_models);

		static void check_read_error(std::ifstream & in);

		static void get_sequence_alphabet(std::vector<unsigned char> & alphabet);
		static void get_alignment_alphabet(std::vector<unsigned char> & alphabet);


		void check_parameters_in_file(const std::string & file) const;
		static void mapmodel_write(std::ofstream & outs, const std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & mapmodel, const std::vector<unsigned char> & alphabet, size_t max_level);
		static void mapmodel_read(std::ifstream & in,  std::map< std::string, std::pair<unsigned char, std::vector<uint32_t> > > & mapmodel, const std::vector<unsigned char> & alphabet, size_t max_level);

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
