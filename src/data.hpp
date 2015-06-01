#ifndef DATA_HPP
#define DATA_HPP

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include<map>
#include <math.h> 
#include <algorithm>

#include "pw_alignment.hpp"
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"


using namespace std;

#include <boost/iostreams/stream.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > btokenizer;
void strsep(string str, const char * sep, vector<string> & parts);

class dnastring {
	public:
	dnastring(string str);
	dnastring();
	dnastring(const dnastring &d);
	~dnastring();

	char at(size_t pos) const;
	size_t length() const;

	static bool found_iupac_ambiguity;
	static char complement(char c);
	static size_t base_to_index(char base);
	static char index_to_base(size_t index);


	private:
	vector<bool> bits;

	static char base_translate_back(bool bit1, bool bit2, bool bit3);

};



class all_data {

	public:
		all_data(string fasta_all_sequences, string maf_all_alignments);
		// no copy constructor, never copy all data
		~all_data();


		const dnastring & getSequence(size_t index) const;
		const pw_alignment & getAlignment(size_t index) const;
	//	const multimap<size_t, size_t> & getAlOnRefMap(size_t seq_idx) const;

		size_t numSequences() const;
		size_t numAlignments() const;
		const size_t numAcc() const;
		bool alignment_fits_ref(const pw_alignment * al) const;
		void print_ref(const pw_alignment * al)const;
		const vector<size_t> & getAcc(size_t acc)const;//return the vector of all sequences for a certain acc.
		size_t accNumber(size_t sequence_id) const;
		const string get_acc(size_t acc)const; //get the accession number and return its name.
		const size_t get_acc_id(string acc)const;
		string get_seq_name(size_t s) const;
		size_t get_seq_size(size_t s) const;
		void set_accession(const string & acc);
		size_t numOfAcc() const;

	
	private:
		// data
		vector<dnastring> sequences;
		vector<string> sequence_names;
		vector<pw_alignment> alignments;
		// fast access indices
		map< string, vector< size_t> > acc_sequences; // acc name -> sequences of that acc
		vector<size_t> sequence_to_accession; // sequence id -> accession id
		map<string, size_t> accession_name; // accession name -> accession id

		map< string, size_t> longname2seqidx; // long sequence name ("Accession:sequence name") -> sequence index

		

		void insert_sequence(const string & acc, const string & seq_name, const string & dna);
		static void name_split(const string & longname, string & acc, string & name);


};


class compute_cc;

class overlap{
public:
	overlap(const all_data&);
	overlap(const overlap & o);
	~overlap();
	// insert p into overlap data structure, return pointer to inserted alignment
	pw_alignment * insert_without_partial_overlap(const pw_alignment & p);
	void insert_without_partial_overlap_p(pw_alignment * p);
	// This function removes an alignment with adress identity to remove from overlap, then deletes remove
	void remove_alignment(const pw_alignment  *remove);
	void remove_alignment_nodelete(const pw_alignment * remove);

	void test_all() const;
	void test_all_part()const;// Change it later to check all those pieces with gap in one sample. Put them all in a set and check if the only missed parts of coverage are those parts.
	void test_overlap()const;
	void print_all_alignment() const;
	const pw_alignment * get_al_at_left_end(size_t ref1, size_t ref2, size_t left1, size_t left2) const;
	multimap<size_t, pw_alignment*>& get_als_on_reference(size_t sequence) ;
	const multimap<size_t, pw_alignment*>& get_als_on_reference_const(size_t sequence) const ;
	void test_multimaps()  ;
	bool checkAlignments(pw_alignment* const p)const;

	const set<pw_alignment*, compare_pw_alignment> & get_all() const;

	void test_partial_overlap() const;
	static void test_partial_overlap_set(set< const pw_alignment *, compare_pw_alignment> & als);
	static void test_partial_overlap_vec(vector< const pw_alignment *> & als);
	static bool check_po(size_t l1, size_t r1, size_t l2, size_t r2);


	size_t size() const;
private:
	const all_data & data;
	set<pw_alignment*, compare_pw_alignment> alignments;
	vector< multimap< size_t, pw_alignment *> > als_on_reference; // sequence index -> pos on that sequence -> alignment reference


	
};



class splitpoints {
	public:
	splitpoints(const pw_alignment & , const overlap &, const all_data &);
	~splitpoints();
	void find_initial_split_points(size_t sequence, size_t left, size_t right);
	void find_initial_split_points_nonrecursive(size_t sequence, size_t left, size_t right);
	void recursive_splits();//initial split points + all induced split points everywhere
	void nonrecursive_splits(); // initial split points only 
	void insert_split_point(size_t sequence, size_t position);//recursively find the other split points
	void insert_split_point_nonrecursive(size_t sequence, size_t position);//insert without recursion
	void split_all(set<const pw_alignment*, compare_pw_alignment> & remove_alignments, vector<pw_alignment> & insert_alignments);
	void splits(const pw_alignment * p,  vector<pw_alignment> & insert_alignments);
	bool onlyGapSample(const pw_alignment* p);
	 vector<pw_alignment>  get_insert () const;

	private:
	const overlap & overl;
	const pw_alignment & newal;
	const all_data & data;
	vector<set<size_t> > split_points;
	vector<pw_alignment> insert_alignments;	
/*
	initial points in sets (method with and without sets)
	compute remove and insert alignments
	cost change level 1, abort if no gain
	if there is gain recurse


*/
	
};


class model{
public:
	model(all_data &);
	~model();
	void acc_base_frequency();
	void alignment_modification();
	void cost_function(pw_alignment& p) const;
	void cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2) const ;
	void gain_function(const pw_alignment& p, double & g1, double & g2) const;

	void train();

private:
	all_data & data;
//	vector<vector<size_t> > transform;
	vector<vector<double> > cost_on_acc;
	vector<vector<vector<vector<double> > > >modification;
};
#define Alignment_level 1
#define Sequence_level 2
#define NUM_DELETE 5
#define NUM_KEEP 10

class abstract_context_functor {
	public:
	abstract_context_functor();
	virtual void see_context(size_t acc1, size_t acc2,const pw_alignment& p, size_t pos, string context, char last_char,ofstream &);
	virtual void see_entire_context(size_t acc1, size_t acc2, string entireContext);
	

};


class counting_functor : public abstract_context_functor {
	public:
	counting_functor(all_data &);
	virtual void see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, string context, char last_char,ofstream&);
	const map<string, vector<double> > & get_context(size_t acc1, size_t acc2)const;
	void total_context();
	double get_total(size_t acc1, size_t acc2, string context)const;
	void create_context(size_t acc1, size_t acc2, string context);
	private:
	all_data & data;
	vector<vector<map<string, vector<double> > > >successive_modification;
	vector<vector<map <string, double > > > total;	




};
class mc_model;

class cost_functor : public abstract_context_functor {
	public:
	cost_functor(all_data &, const vector<vector<map<string, vector<double> > > >&);
	virtual void see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, string context, char last_char,ofstream& );
	double get_modify(const pw_alignment & p, size_t acc1, size_t acc2)const;
	private:
	vector<vector< map<string, vector<double> > > >  modification; 
	all_data & data;
	double modify1;
	double modify2;	

};


class adding_functor : public abstract_context_functor {
	public:

};
class encoding_functor : public abstract_context_functor {
	public:
	encoding_functor(all_data& , mc_model*, wrapper &, dlib::entropy_encoder_kernel_1 &);
	virtual void see_context(size_t acc1, size_t acc2,const pw_alignment& p, size_t pos, string context, char last_char,ofstream&);	
	virtual void see_entire_context(size_t acc1, size_t acc2, string entireContext);
	const map<string, vector<double> > & get_alignment_context()const;
//	vector<string> & get_alignment_context(pw_alignment& p)const;
	private:
	all_data & data;
	mc_model * model;	
	wrapper& wrappers;
	dlib::entropy_encoder_kernel_1 & enc;
	string alignment_pattern;//shayadam behtar bashe ye vector of string tarif konam
	map<string, vector<double> > alignment_context;

};

class clustering_functor : public abstract_context_functor{
	public:
	virtual void see_context(size_t acc1, size_t acc2, const pw_alignment& p, size_t pos, string context, char last_char, ofstream &);//computing_modification_oneToTwo is used to fill in the map of modification between center and its associated member.
	//fek konam hamoon encoding_functor ok bashe, lazem nist ino benvisim

	private:
	map<string, vector<double> >modification;


};
class decoding_functor : public abstract_context_functor {


};

class mc_model{
	public:
		mc_model(all_data&);
		~mc_model();	
		void markov_chain();
		void markov_chain_alignment(ofstream&);
		const vector<vector<map <string, vector<double> > > >& get_mod_cost()const;
		void cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2,ofstream&)const ;
		void gain_function(const pw_alignment& p, double & g1, double & g2,ofstream&)const ;
		void train(ofstream &);
		char modification_character(int modify_base, int num_delete, int insert_base, int num_keep)const;
		void modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep)const;
		void computing_modification_oneToTwo(const pw_alignment & p, abstract_context_functor & functor,ofstream&)const;
		void computing_modification_twoToOne(const pw_alignment & p, abstract_context_functor & functor, ofstream&)const;
		void cost_function( pw_alignment& p, ofstream&) const;
		string print_modification_character(char enc)const;
		const map<string, vector<double> > & getPattern(size_t acc)const;
		const vector<double> & get_create_cost(size_t acc) const;
		const vector<map<string, vector<double> > > & model_parameters()const;
		void write_parameters(ofstream &);
		void write_alignments_pattern(ofstream&);
		vector<unsigned int> get_high_at_position(size_t seq_index, size_t position) const;
		vector<unsigned int> get_center_high_at_position(size_t cent_ref, size_t cent_left, size_t position)const;
		vector<unsigned int> get_reverse_center_high_at_position(size_t cent_ref, size_t cent_left, size_t position)const;
		const map<string, vector<unsigned int> > & get_high(size_t acc)const;
		void make_all_the_patterns();
		void make_all_alignments_patterns();
		void set_patterns(ifstream&);
		void set_alignment_pattern(ifstream&);
		string get_context(size_t position, size_t seq_id)const;
		vector<size_t> get_powerOfTwo()const;
		string get_firstPattern()const;
		string get_firstAlignmentPattern() const;
		const map<string, vector<double> > & get_alignment_context(size_t al_id, size_t seq_id, encoding_functor & functor)const;
		const map<string, vector<unsigned int> >& get_highValue(size_t acc1, size_t acc2)const;
		void computing_modification_in_cluster(string center, string member)const;
		const map<string, vector<double> > & get_cluster_member_context(pw_alignment & al, size_t center_id, encoding_functor & functor)const;
	//	void print_modification(char enc)const;
		size_t modification_length(char mod)const;
		void get_encoded_member(pw_alignment & al, size_t center_ref, size_t center_left, encoding_functor & functor,ofstream&)const;
	private:
	all_data & data;
	vector<map<string, vector<double> > >sequence_successive_bases;
	vector<vector<double> > create_cost;
	vector<size_t> powersOfTwo;
	vector<vector<map<string, vector<double> > > >mod_cost; // alignment modificaton information cost
	vector<map<string, vector<unsigned int> > > high;//sequences patterns
	map<string,vector<double> > all_the_patterns; // TODO vector<double> part is wrong (independent of accession)

	set<string> all_alignment_patterns; // all possible alignment patterns
	vector<vector<map<string , vector<unsigned int> > > >highValue;//alignments patterns
		

};

class encoding{
	public:
	encoding(all_data&);
	~encoding();
	void calculate_bounds();
	private:
	all_data & data;	




};
	
	









#endif
