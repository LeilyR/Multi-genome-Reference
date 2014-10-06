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
		size_t numAcc() const;
		bool alignment_fits_ref(const pw_alignment * al) const;
		void print_ref(const pw_alignment * al)const;
		vector<size_t> getAcc(string accName)const;//return the vector of all sequences for a certain acc.
		size_t accNumber(size_t sequence_id);

	
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
	void split_partial_overlap(const pw_alignment * new_alignment, set<const pw_alignment*, compare_pw_alignment> & remove_alignments, vector<pw_alignment> & insert_alignments, size_t level) const;
	void insert_without_partial_overlap(const pw_alignment & p);
	void remove_alignment(const pw_alignment  *remove);

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
	bool check_po(size_t l1, size_t r1, size_t l2, size_t r2) const;


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
	void find_initial_split_point();//initial split points
	void insert_split_point(size_t sequence, size_t position);//recursively find the other split points
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

#define NUM_DELETE 5
#define NUM_KEEP 10

class abstract_context_functor {
	public:
	abstract_context_functor();
	virtual void see_context(const pw_alignment * p, size_t pos, string context, char last_char);
	

};


class counting_functor : public abstract_context_functor {
	public:
	counting_functor(all_data &);
	virtual void see_context(const pw_alignment * p, size_t pos, string context, char last_char);
	const map<string, vector<double> > & get_context(size_t acc1, size_t acc2)const;
	void total_context();
	double get_total(size_t acc1, size_t acc2, string context)const;
	private:
	all_data & data;
	vector<vector<map<string, vector<double> > > >successive_modification;
	vector<vector<map <string, double > > > total;			




};

class adding_functor : public abstract_context_functor {
	public:

};
class encoding_functor : public abstract_context_functor {


};
class decoding_functor : public abstract_context_functor {


};

class mc_model{
	public:
		mc_model(all_data&);
		~mc_model();	
		void markov_chain();
		void markov_chain_alignment();
		void cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2)const ;
		void gain_function(const pw_alignment& p, double & g1, double & g2)const ;
		void train();
		char modification_character(int modify_base, int num_delete, int insert_base, int num_keep)const;
		void modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep);
		void computing_modification_oneToTwo(const pw_alignment & p, abstract_context_functor & functor)const;
		void computing_modification_twoToOne(const pw_alignment & p, abstract_context_functor & functor)const;
		void cost_function( pw_alignment& p) const;
		string print_modification_character(char enc);
	private:
	all_data & data;
	vector<map<string, vector<double> > >sequence_successive_bases;
	vector<vector<double> > create_cost;
	vector<size_t> powersOfTwo;
	vector<vector<map<string, vector<double> > > >mod_cost;

};
	
	









#endif
