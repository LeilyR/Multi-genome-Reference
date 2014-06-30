#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "pw_alignment.hpp"
#include "data.hpp"

using namespace std;

#define TEST 0


#if TEST

int main(int argc, char * argv[]) {
	cout << " hello " << endl;

	pw_alignment p(string("ATT----TTCTT"), string("AGTGATAT----"), 12, 15, 23, 26,1,1);
	pw_alignment s(string("ATT----TTCTT"), string("ACTGATG---AC"),13, 18, 24, 29,2,1);

	for(size_t i=0; i<p.alignment_length(); ++i) {
		char s1;
		char s2;
		p.alignment_col(i, s1, s2);
		cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << endl;
	  
	}

	pw_alignment p1;
	pw_alignment p2;
	p.split(false,4, p1, p2);
		for(size_t i=0; i<p2.alignment_length(); ++i) {
		char s1;
		char s2;
		p2.alignment_col(i, s1, s2);
		cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << endl;
	 
	}
	overlap o;
	o.add_alignment(s);
	typedef multimap <size_t, pw_alignment> mmap;
	mmap alignment_map;
	alignment_map.insert(make_pair(12,p));
	alignment_map.insert(make_pair(23,p));
	alignment_map  alignments_on_reference1;
	mmap::const_iterator alignit = alignments_on_reference1.begin();
	for(; alignit != alignments_on_reference1.end(); ++alignit) {
		char s1;
		char s2;
		for(size_t i=0;i<p.alignment_length();++i){
		s.alignment_col(i, s1, s2);
		cout << "pos " << i<< " s1 " << s1 << " s2 " << s2 << endl;
}

}
	return 0;

}



#endif
void usage() {
	cerr << "Graph program" << endl;
	cerr << "Usage:" << endl;
}

int do_fasta_prepare(int argc, char * argv[]) {

	if(argc < 4) {
		usage();
		cerr << "Program fasta_prepare" << endl; 
		cerr << "Combine fasta files from different species/accessions into a common file" << endl;
		cerr << "Parameters: " << endl;
		cerr << "output file -- fasta file" << endl;
		cerr << "<list of fasta files> -- input file, one file per species/accession, no colon in file name" << endl;
		return 1;
	}


	std::ofstream outf(argv[2]);
	if(outf) {
		vector<ifstream*> inputs;
		vector<string> names;
		for(int i=3; i<argc; ++i) {
			string istr(argv[i]);
			vector<string> slashparts;
			strsep(istr, "/", slashparts);
			string afterslash = slashparts.at(slashparts.size()-1);
			vector<string> periodparts;
			strsep(afterslash, ".", periodparts);
			string name = periodparts.at(0);
			if(periodparts.size()>1) {
				stringstream sstr;
				for(size_t j=0; j<periodparts.size()-1; j++) {
					if(j!=0) sstr << ".";
					sstr << periodparts.at(j);
				}
				name = sstr.str();
			}
			for(size_t j=0; j<name.length(); ++j) {
				if(':' == name.at(j)) {
					cerr << "Error: " << argv[i] << " file name contains a colon" << endl;
					exit(1);
				}
			}
			ifstream * iin = new ifstream(argv[i]);
			if(!(*iin)) {
				cerr << "Error: cannot read: " << argv[i] << endl;
				exit(1);
			}
			inputs.push_back(iin);
			names.push_back(name);
		}

		for(size_t i=0; i<inputs.size(); ++i) {


			string str;
			while(getline((*(inputs.at(i))), str)) {
				if(str.length()>0) { // ignore empty lines
				if(str.at(0)=='>') {
					string nname = str.substr(1);
					stringstream nhead;
					nhead << ">" << names.at(i) << ":" << nname;
					outf << nhead.str() << endl;
				} else {
					outf << str << endl;
				}
				}
			}
			inputs.at(i)->close();
			delete inputs.at(i);
		}
		outf.close();
	} else {
		cerr << "Error: cannot write: " << argv[2] << endl;
		exit(1);
	}

	outf.close();

	return 0;
}
 


int do_model(int argc, char * argv[]) {

	if(argc < 4) {
		usage();
		cerr << "Program: model" << endl;
		cerr << "Parameters:" << endl;
		cerr << "* fasta file from fasta_prepare" << endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << endl;
	}

	string fastafile(argv[2]);
	string maffile(argv[3]);

	all_data data(fastafile, maffile);



	return 0;
}


#if !TEST
	
int main(int argc, char * argv[]) {

	if(argc < 2) {
		usage();
		exit(1);
	}
	string program(argv[1]);

	if(0==program.compare("fasta_prepare")) {
		return do_fasta_prepare(argc, argv);
	} else if(0==program.compare("model")) {
		return do_model(argc, argv);
	}

	else {
		usage();
	}


	return 1;
}



#endif























