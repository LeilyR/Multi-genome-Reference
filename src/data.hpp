#ifndef DATA_HPP
#define DATA_HPP

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>


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
		const multimap<size_t, size_t> & getAlOnRefMap(size_t seq_idx) const;

	private:
		// data
		vector<dnastring> sequences;
		vector<string> sequence_names;
		vector<pw_alignment> alignments;
		// fast access indices
		map< string, vector< size_t> > acc_sequences; // acc name -> sequences of that acc
		map< string, size_t> longname2seqidx; // long sequence name ("Accession:sequence name") -> sequence index
		vector< multimap< size_t, size_t> > als_on_reference; // sequence index -> pos on that sequence -> alignment index


		void insert_sequence(const string & acc, const string & seq_name, const string & dna);
		static void name_split(const string & longname, string & acc, string & name);


};

















#endif
