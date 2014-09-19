#include "data.hpp"


void strsep(string str, const char * sep, vector<string> & parts) {
        boost::char_separator<char> tsep(sep);
        btokenizer tok(str, tsep);
        for(btokenizer::iterator tit=tok.begin(); tit!=tok.end(); ++tit) {
                parts.push_back(*tit);
        }
        
}

bool dnastring::found_iupac_ambiguity = false;

char dnastring::complement(char c) {
	switch(c) {
		case 'a':
		case 'A':
		return 'T';
		break;
		case 't':
		case 'T':
		return 'A';
		break;
		case 'c':
		case 'C':
		return 'G';
		break;
		case 'g':
		case 'G':
		return 'C';
		break;
	
			case 'R':
			case 'r':
			case 'Y':
			case 'y':
			case 'M':
			case 'm':
			case 'K':
			case 'k':
			case 'W':
			case 'w':
			case 'S':
			case 's':
			case 'B':
			case 'b':
			case 'D':
			case 'd':
			case 'H':
			case 'h':
			case 'V':
			case 'v':
			case 'N':	
			case 'n':
		return 'N';
		break;
		case '.':
		case '-':
		return '-';
		break;

	}
	return 'X';
}


dnastring::dnastring(string str):bits(str.length()*3) {
	for(size_t i=0; i<str.length(); ++i) {
		char c = str.at(i);
		bool bit1 = false;
		bool bit2 = false;
		bool bit3 = false;
		switch(c) {
			case 'A':
			case 'a':
				bit1 = false;
				bit2 = false;
				bit3 = false;

			break;
			case 'C':
			case 'c':
				bit1 = true;
				bit2 = false;
				bit3 = false;

			break;
			case 'T':
			case 't':
				bit1 = false;
				bit2 = true;
				bit3 = false;

			break;
			case 'G':
			case 'g':
				bit1 = true;
				bit2 = true;
				bit3 = false;

			break;
		
			case 'R':
			case 'r':
			case 'Y':
			case 'y':
			case 'M':
			case 'm':
			case 'K':
			case 'k':
			case 'W':
			case 'w':
			case 'S':
			case 's':
			case 'B':
			case 'b':
			case 'D':
			case 'd':
			case 'H':
			case 'h':
			case 'V':
			case 'v':
				found_iupac_ambiguity = true;
			case 'N':	
			case 'n':
				bit1 = false;
				bit2 = false;
				bit3 = true;
			break;

			default:
				cerr << "Error: Illegal character in DNA sequence: " << c << endl;
				exit(1);
			break;
		}

		bits.at(3*i) = bit1;
		bits.at(3*i+1) = bit2;
		bits.at(3*i+2) = bit3;
	}


}

dnastring::dnastring() {}

dnastring::dnastring(const dnastring & d): bits(d.bits) {
}

dnastring::~dnastring() {}

char dnastring::at(size_t pos) const {
	bool bit1 = bits.at(3*pos);
	bool bit2 = bits.at(3*pos+1);
	bool bit3 = bits.at(3*pos+2);
	return base_translate_back(bit1, bit2, bit3);
}

size_t dnastring::length() const {
	return bits.size() / 3;
}

char dnastring::base_translate_back(bool bit1, bool bit2, bool bit3) {

	if(bit1) {
		if(bit2) {
			if(bit3) {
				cerr << "Error: wrong bit vector 111" << endl;
				exit(1);
				return 'X';
			} else {
				return 'G';
			}
		
		} else {
			if(bit3) {
				cerr << "Error wrong bit vector 101" << endl;
				exit(1);
				return 'X';
			} else {
				return 'C';
			}

		}	
	
	} else {
		if(bit2) {
			if(bit3) {
				cerr << "Error: wrong bit vector 011" << endl;
				exit(1);
				return 'X';
			} else {
				return 'T';
			}

		} else {
			if(bit3) {
				return 'N';
			} else {
				return 'A';
			}

		}
	}
}

size_t dnastring::base_to_index(char base) {
	switch(base){
		case 'A':	
		return 0;
		break;
		case 'T':
		return 1;
		break;
		case 'C':
		return 2;
		break;
		case 'G':
		return 3;
		break;
		case 'N':		
		return 4;
		break;
		case '-':
		return 5;
		break;
		}
	return 100;
}
char dnastring::index_to_base(size_t index) {
	switch(index){
		case 0:	
		return 'A';
		break;
		case 1:
		return 'T';
		break;
		case 2:
		return 'C';
		break;
		case 3:
		return 'G';
		break;
		case 4:		
		return 'N';
		break;
		case 5:
		return '-';
		break;
		}
	return 'X';
}


all_data::all_data(string fasta_all_sequences, string maf_all_alignments) {
	ifstream fastain(fasta_all_sequences.c_str());
	if(fastain) {
		string str;
		stringstream curseq;
		string curname("");
		string curacc("");
		while(getline(fastain, str)) {
			if(str.at(0)=='>') {
				if(0!=curname.compare("")) { // store previous sequence
					insert_sequence(curacc, curname, curseq.str());
					curname = "";
					curacc = "";
					curseq.str("");
				}

				// read next header
				string fhead = str.substr(1);
				name_split(fhead, curacc, curname);
			} else {
				curseq << str;
			}
		}
		// store last sequence
		insert_sequence(curacc, curname, curseq.str());
		//cout << "R " << curseq.str().substr(62750, 15) << endl;
		curname = "";
		curacc = "";
		curseq.str("");
		fastain.close();
	} else {
		cerr << "Error: cannot read: " << fasta_all_sequences << endl;
		exit(1);
	}

	cout << "loaded " << sequences.size() << " sequences" << endl;

	vector< multimap< size_t, size_t> > als_on_reference; // sequence index -> begin pos on that sequence -> alignment index
	// a new multimap for each ref sequence
	als_on_reference.resize(sequences.size());
	for(size_t i=0; i<sequences.size(); ++i) {
		als_on_reference.at(i) = multimap<size_t, size_t>();
	}



	ifstream mafin(maf_all_alignments.c_str());
	size_t skip_self = 0;
	size_t skip_alt = 0;
	if(mafin) {
		string str;
		while(getline(mafin, str)) {
			if(str.at(0)!='#') { // skip headers
				if(str.at(0)=='a') { // next alignment
					string aline1;
					string aline2;
					getline(mafin, aline1);
					getline(mafin, aline2);
					getline(mafin, str); // empty line at end of alignment
					if(aline1.at(0)!='s' || aline2.at(0)!='s') {
						cerr << "Error: exspected maf sequence lines, found: " << endl << aline1 << endl<< aline2 << endl;
						exit(1);
					}

					vector<string> parts;
					strsep(aline1, " ", parts);
					if(parts.size()!=7) {
						cerr << "Error: exspected 7 fields in sequence line: " << aline1 << endl;
						exit(1);
					}
					string acc1;
					string name1;
					name_split(parts.at(1), acc1, name1);
					size_t start1 = atoi(parts.at(2).c_str());
					size_t size1 = atoi(parts.at(3).c_str());
					char strand1 = parts.at(4).at(0);
					size_t seqlength1 = atoi(parts.at(5).c_str());
					string al1 = parts.at(6);
					map<string, size_t>::iterator findseq1 = longname2seqidx.find(parts.at(1));
					if(findseq1==longname2seqidx.end()) {
						cerr << "Error: unknown sequence: " << parts.at(1) << endl;
						exit(1);
					}
					size_t idx1 = findseq1->second;
					// maf uses coordinates in recom sequence for - strand alignments, transform to end < start coordinates
					size_t incl_end1 = start1 + size1 - 1;
					if(strand1!='+') {
						size_t tmp = start1;
						start1 = seqlength1 - start1 -1;
						incl_end1 = seqlength1 - tmp - size1;
						for(size_t j=0; j<al1.length(); ++j) {
							al1.at(j) = dnastring::complement(al1.at(j));
						}
					}


					parts.clear();
					strsep(aline2, " ", parts);
					string acc2;
					string name2;
					name_split(parts.at(1), acc2, name2);
					size_t start2 = atoi(parts.at(2).c_str());
					size_t size2 = atoi(parts.at(3).c_str());
					size_t seqlength2 = atoi(parts.at(5).c_str());
					char strand2 = parts.at(4).at(0);
					string al2 = parts.at(6);
					map<string, size_t>::iterator findseq2 = longname2seqidx.find(parts.at(1));
					if(findseq2==longname2seqidx.end()) {
						cerr << "Error: unknown sequence: " << parts.at(1) << endl;
						exit(1);
					}
					size_t idx2 = findseq2->second;
					size_t incl_end2 = start2 + size2 - 1;
					if(strand2!='+') {
						size_t tmp = start2;
						start2 = seqlength2 - start2 -1;
						incl_end2 = seqlength2 - tmp - size2;
						for(size_t j=0; j<al1.length(); ++j) {
							al2.at(j) = dnastring::complement(al2.at(j));
						}
					}
					
					// both al parts not identical	
					if(idx1 != idx2 || !(start1==start2 && incl_end1 == incl_end2    )) {



						bool skip = false;

						// check if we already have an alignment with same coordinates
						// because we are looking for identical alignment it suffices to go over only one multimap
						pair<multimap<size_t, size_t>::iterator, multimap<size_t, size_t>::iterator> eqr =
						als_on_reference.at(idx1).equal_range(start1);
						for(multimap<size_t, size_t>::iterator it = eqr.first; it!=eqr.second; ++it) {
							pw_alignment & al = alignments.at(it->second);
							size_t same_ref = 2;
							if(al.getreference(0)==idx1) {
								same_ref = 0;
							} else if(al.getreference(1)==idx1) {
								same_ref = 1;
							}
							assert(same_ref < 2);
							size_t other_ref = 1 - same_ref;
							if(al.getreference(other_ref) == idx2) { // al and new alignment on same sequences
								if(start1 == al.getbegin(same_ref) && start2 == al.getbegin(other_ref) && 
										incl_end1 == al.getend(same_ref) && incl_end2 == al.getend(other_ref)) {
									skip =true;
									break;
								}
							}
						 }
					
						if(skip) {
							skip_alt++;
						} else {



							pw_alignment al(al1, al2, start1, start2, incl_end1, incl_end2, idx1, idx2);
							size_t alidx = alignments.size();
							alignments.push_back(al);

							als_on_reference.at(idx1).insert(make_pair(start1, alidx));
							als_on_reference.at(idx2).insert(make_pair(start2, alidx));


						}

					} else {
						skip_self++;
						// cerr << "Warning: Skip self alignment in seq " << idx1 << " from " << start1 << " to " << incl_end1 << endl;
					
					}

					
				
				}
			
			
			}	
		}

		mafin.close();

		if(dnastring::found_iupac_ambiguity) {
			cerr << "Warning: IUPAC DNA ambiguity characters were replaced by N" << endl;
		}

		cout << "Loaded: " << sequences.size() << " sequences and " << alignments.size() << " pairwise alignments " << endl;
		cout << skip_self << " self alignments were skipped" << endl;
		cout << skip_alt << " alternative alignments of identical regions were skipped" << endl;
	
	} else {
		cerr << "Error: cannot read: " << maf_all_alignments << endl;
		exit(1);
	}


}

all_data::~all_data() {

}


	const dnastring & all_data::getSequence(size_t index) const {
		return sequences.at(index);
	}


	const pw_alignment & all_data::getAlignment(size_t index) const {
		return alignments.at(index);
	}
	vector<size_t> all_data::getAcc(string accName)const{
		return acc_sequences.at(accName);
	}
	size_t all_data::accNumber(size_t sequence_id){
		return sequence_to_accession.at(sequence_id);
	}


//const multimap<size_t, size_t> & all_data::getAlOnRefMap(size_t seq_idx) const {
//	return als_on_reference.at(seq_idx);
//}

	size_t all_data::numSequences() const {
		return sequences.size();
	}

	size_t all_data::numAlignments() const {
 		return alignments.size();
	}
	size_t all_data::numAcc() const {
		return  acc_sequences.size();
	}


	void all_data::insert_sequence(const string & acc, const string & seq_name, const string & dna) {
	size_t seq_id = sequences.size();
	sequences.push_back(dnastring(dna));
	sequence_names.push_back(seq_name);
	map<string, vector<size_t> >::iterator findacc = acc_sequences.find(acc);
	size_t acc_id;
	if(findacc == acc_sequences.end()) {
		acc_sequences.insert(make_pair(acc, vector<size_t>(0)));
		findacc = acc_sequences.find(acc);
		acc_id = accession_name.size();	
		accession_name.insert(make_pair(acc, acc_id));	
	} else {
	map<string, size_t>::iterator find_acc_id = accession_name.find(acc);
	acc_id = find_acc_id -> second;
	}
	findacc->second.push_back(seq_id);
	sequence_to_accession.push_back(acc_id);
	stringstream longname;
	longname << acc << ":" << seq_name;
	longname2seqidx.insert(make_pair(longname.str(), seq_id));
	}


/**
	Transform 
	"Col0:scaffold_0 comment"
	to 
	"Col0", "scaffold_0"
	
**/		
	void all_data::name_split(const string & longname, string & acc, string & name) {
	vector<string> wparts; // remove all after first space
	strsep(longname, " ", wparts);
	string wname = wparts.at(0);
	vector<string> fparts;
	strsep(wname, ":", fparts);
	if(fparts.size()< 2) {
		cerr << "Error: sequence " << longname << " does not contain an Accession name separated by a colon" << endl;
		exit(1);
	}
	acc = fparts.at(0);
	stringstream namewriter;
	for(size_t i=1; i<fparts.size();++i) {
		if(i!=1) namewriter << ":";
		namewriter << fparts.at(i);
	}
	name = namewriter.str();


}


	bool all_data::alignment_fits_ref(const pw_alignment * al) const {
	const dnastring & ref1 = getSequence(al->getreference1());
	const dnastring & ref2 = getSequence(al->getreference2());
	bool gapOnSample1 = true;
	bool gapOnSample2 = true;
	bool al1fw = true;
	bool al2fw = true;
	size_t al1at = al->getbegin1();
	size_t al2at = al->getbegin2();
	if(al->getbegin1() > al->getend1()) {
		al1fw = false;
	}

	if(al->getbegin2() > al->getend2()) {
		al2fw = false;
	}
//al->print();
//	cout << " directions " << al1fw << al2fw << endl;

	for(size_t col=0; col < al->alignment_length(); ++col) {
		char al1char='X';
		char al2char='X';
		al->alignment_col(col, al1char, al2char);	
		
	//	cout << "c " << col << " " << al1char << " " << al2char << endl;	
		if(al1char!='-') { 
			char rchar = ref1.at(al1at);
		//	cout << " al 1 at " << al1at;
			if(al1fw) {
				al1at++;
			} else {
				al1at--;
			}
		//	cout << " is " << rchar << endl;
			if(rchar != al1char) {
			//	al->print();
				cerr << "Warning: alignment test failed. Alignment has " << al1char << " at " << col <<" where reference sequence has " << rchar << endl;
			//	print_ref(al);
				return false;
			}
		}
		if(al2char!='-') {
			char rchar = ref2.at(al2at);
		//	cout << " al 2 at " << al2at ;
			if(al2fw) {
				al2at++;
			} else {
				al2at--;
			}
		//	cout << " is " << rchar << endl;
			if(rchar != al2char) {
			//	al->print();
				cerr << " at " << al2at <<endl;
				cerr << "Warning: alignment test failed. Alignment has " << al2char << " where reference sequence has " << rchar << endl;
				al->print();
			//	print_ref(al);
				return false;
			}
		}
	}
	for(size_t i = 0;i<al->alignment_length();i++){
		char al1char = 'X';
		char al2char = 'X';
		al->alignment_col(i, al1char, al2char);	
		if(al1char !='-' ) gapOnSample1 = false;
		if(al2char !='-' ) gapOnSample2 = false;		
	}
		if(gapOnSample1 || gapOnSample2 ) { }
		else{

			if(al1fw) {
				size_t al1e = al1at-1;
				if(al1e!=al->getend1()) {
					cerr << "Warning: alignment end wrong. Should be " << al1e << " but is " << al->getend1() << endl;
					//print_ref(al);
					return false;
				} 
			
			} else {
				size_t al1b = al1at+1;
				if(al1b!=al->getend1()) {
					print_ref(al);

					cerr << "Warning: rev alignment end wrong. Should be " << al->getend1() << " but is " << al1b << endl;

					return false;
				}
			}
		
			if(al2fw) {
			size_t al2e = al2at-1;
			if(al2e!=al->getend2()) {
				print_ref(al);
				cerr << "Warning: alignment end wrong. Should be " << al->getend2()<< " but is " << al2e << endl;
				return false;
			}
		 	} else {
				size_t al2b = al2at+1;
				if(al2b!=al->getend2()) {
			//	print_ref(al);
				cerr << "Warning: rev alignment end wrong. Should be " << al2b  << " but is " << al->getend2() << endl;
				print_ref(al);
			
				return false;
				}
			}
		}
	
	return true;
		
	
}


	void all_data::print_ref(const pw_alignment * al)const{
		const dnastring & ref1 = getSequence(al->getreference1());
		const dnastring & ref2 = getSequence(al->getreference2());
		size_t al1at = al->getbegin1();
		size_t al2at = al->getbegin2();

		if(al->getbegin1() < al->getend1())
			cout << " direction1: forward "<< endl;
		else 	cout << " direction1: reverse "<< endl;
		if(al->getbegin2() < al->getend2())
			cout << " direction2: forward "<< endl;
		else 	cout << " direction2: reverse "<< endl;
		for(size_t col=0; col < al->alignment_length(); ++col) {
			char al1char='X';
			char al2char='X';
			al->alignment_col(col, al1char, al2char);	
		
			cout << "c " << col << " " << al1char << " " << al2char << endl;
				if(al1char!='-') { 
				char rchar = ref1.at(al1at);
				cout << " al 1 at " << al1at;
					if(al->getbegin1() < al->getend1()) {
						al1at++;
					} else {
						al1at--;
					}
				cout << " is " << rchar << endl;
				}
				if(al2char!='-') {
				char rchar = ref2.at(al2at);
				cout << " al 2 at " << al2at ;
					if(al->getbegin2() < al->getend2()) {
						al2at++;
					} else {
						al2at--;
					}
				cout << " is " << rchar << endl;
				}

		}	
		al -> print();
}


	overlap::overlap(const all_data & d): data(d), als_on_reference(d.numSequences()){

}

	overlap::overlap(const overlap & o): data(o.data), als_on_reference(o.data.numSequences()){}


	
overlap::~overlap(){}

void overlap::remove_alignment(const pw_alignment * remove){
	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqrb1 =
	als_on_reference.at(remove->getreference1()).equal_range(remove->getbegin1());
	for(multimap<size_t, pw_alignment*>::iterator it = eqrb1.first; it!=eqrb1.second; ++it) {
		if ( it->second == remove ) {
			als_on_reference.at(remove->getreference1()).erase(it);
			break;
		}		
	}
	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqre1 =
	als_on_reference.at(remove->getreference1()).equal_range(remove->getend1());
	for(multimap<size_t, pw_alignment*>::iterator it = eqre1.first; it!=eqre1.second; ++it) {
		if ( (it->second) == remove ) {
			als_on_reference.at(remove->getreference1()).erase(it);
			break;
		}		
	}

	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqrb2 =
	als_on_reference.at(remove->getreference2()).equal_range(remove->getbegin2());
	for(multimap<size_t, pw_alignment*>::iterator it = eqrb2.first; it!=eqrb2.second; ++it) {
		if ( (it->second) == remove ) {
			als_on_reference.at(remove->getreference2()).erase(it);
			break;
		}		
	}

	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqre2 =
	als_on_reference.at(remove->getreference2()).equal_range(remove->getend2());
	for(multimap<size_t, pw_alignment*>::iterator it = eqre2.first; it!=eqre2.second; ++it) {
		if ( (it->second) == remove ) {
			als_on_reference.at(remove->getreference2()).erase(it);
			break;
		}		
	}

	pw_alignment * rem = const_cast<pw_alignment *>(remove);
	set<pw_alignment*, compare_pw_alignment>::iterator findr = alignments.find(rem);
//	remove->print();
	assert(findr!=alignments.end());
	alignments.erase(findr);
	//delete remove;
}

	void overlap::test_all() const {

	for(set<pw_alignment*, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		pw_alignment * al = *it;
		bool alok = data.alignment_fits_ref(al);
		if(!alok) {
			cout << " test all failed " << endl;
			al->print();
			exit(1);
		}
		
	
	}


}




	void overlap::print_all_alignment() const {
		for(set<pw_alignment*, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it ) {
			const pw_alignment * al = * it;
			//cout << " x " << al << endl;
			al->print();
		}

}

	void overlap::test_multimaps() {

		for(set<pw_alignment*, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it ) {
			const pw_alignment * remove = * it;
			pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqrb1 =
			als_on_reference.at(remove->getreference1()).equal_range(remove->getbegin1());
			size_t found =0;
			for(multimap<size_t, pw_alignment*>::iterator it = eqrb1.first; it!=eqrb1.second; ++it) {
			if ( it->second == remove ) {
				found++;
				break;
			}		
			}
			pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqre1 =
			als_on_reference.at(remove->getreference1()).equal_range(remove->getend1());
			for(multimap<size_t, pw_alignment*>::iterator it = eqre1.first; it!=eqre1.second; ++it) {
				if ( (it->second) == remove ) {
				found++;
				break;
			}		
			}

			pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqrb2 =
			als_on_reference.at(remove->getreference2()).equal_range(remove->getbegin2());
			for(multimap<size_t, pw_alignment*>::iterator it = eqrb2.first; it!=eqrb2.second; ++it) {
			if ( (it->second) == remove ) {
				found++;
				break;
			}		
			}

			pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqre2 =
			als_on_reference.at(remove->getreference2()).equal_range(remove->getend2());
			for(multimap<size_t, pw_alignment*>::iterator it = eqre2.first; it!=eqre2.second; ++it) {
			if ( (it->second) == remove ) {
				found++;
				break;
			}		
			}
			assert(found==4);

		}
		for(size_t s=0; s<data.numSequences(); s++) {
			for(multimap< size_t, pw_alignment *>::const_iterator it = als_on_reference.at(s).begin(); it!=als_on_reference.at(s).end(); ++it) {
				cout << it->second << endl;
				pw_alignment * p = const_cast<pw_alignment *>(it->second);
				//cout << p << endl;
				set<pw_alignment*, compare_pw_alignment>::iterator findp = alignments.find(p);
				if(findp==alignments.end()){
					cout<<"wrong one: "<<p<<endl;
					p->print();
				}
				assert(findp!=alignments.end());
			}

		}
		
	}



	
	void overlap::insert_without_partial_overlap(const pw_alignment & p){
	pw_alignment * np = new pw_alignment(p);
	//cout << " insert " << np << endl;
	assert(alignments.find(np) == alignments.end());
	alignments.insert(np);

 	
	multimap<size_t , pw_alignment *> & alignment_on_reference1 = als_on_reference.at(np->getreference1());

	
//	if(np->getreference2()!=np->getreference1()) {
		multimap<size_t , pw_alignment *> & alignment_on_reference2 = als_on_reference.at(np->getreference2());
		pair<size_t,  pw_alignment *> begin2(np->getbegin2(),np);
		alignment_on_reference2.insert(begin2);
		pair<size_t, pw_alignment *> end2(np->getend2(),np);
		alignment_on_reference2.insert(end2);
//	}

	pair<size_t,  pw_alignment *> begin1(np->getbegin1(),np);
		alignment_on_reference1.insert(begin1);
	pair<size_t,  pw_alignment *> end1(np->getend1(),np);
		alignment_on_reference1.insert(end1);
}


	const pw_alignment * overlap::get_al_at_left_end(size_t ref1, size_t ref2, size_t left1, size_t left2) const {
		const multimap<size_t, pw_alignment *> & r1map = als_on_reference.at(ref1);
		pair<multimap<size_t, pw_alignment *>::const_iterator,multimap<size_t, pw_alignment *>::const_iterator> eqr = r1map.equal_range(left1);
		for( multimap<size_t, pw_alignment *>::const_iterator it = eqr.first; it!= eqr.second; ++it) {
			const pw_alignment * cal = it->second;
			if(cal->getreference1()==ref1) {
				size_t cal_left1 = cal->getbegin1();
				if(cal->getend1() < cal_left1) cal_left1 = cal->getend1();
				if(cal_left1!=left1) continue;
				if(ref2!=cal->getreference2()) continue;
				size_t cal_left2 = cal->getbegin2();
				if(cal->getend2() < cal_left2) cal_left2 = cal->getend2();
				if(cal_left2!=left2) continue;


				return cal;
			} 
			else { // ref1 is reference2 of cal
				size_t cal_left1 = cal->getbegin1();
				if(cal->getend1() < cal_left1) cal_left1 = cal->getend1();
				if(cal_left1!=left2) continue;
				if(ref2!=cal->getreference1()) continue;	
				size_t cal_left2 = cal->getbegin2();
				if(cal->getend2() < cal_left2) cal_left2 = cal->getend2();
				if(cal_left2!=left1) continue;
				return cal;
			}
		}
		return NULL;
	}

	multimap<size_t, pw_alignment*> &  overlap::get_als_on_reference(size_t sequence)  {
		return als_on_reference.at(sequence);
	}
	const multimap<size_t, pw_alignment*> &  overlap::get_als_on_reference_const(size_t sequence) const {
		return als_on_reference.at(sequence);
	}

	void overlap::test_all_part()const{
	
	/*set < const pw_alignment*, compare_pw_alignment> alignment_parts;
		
		for(size_t i =0 ; i< data.numAlignments(); ++i){
		const pw_alignment & al = data.getAlignment(i);
		size_t al_right1 = al.getend1();
		size_t al_right2 = al.getend2();
		size_t al_left1 = al.getbegin1();
		size_t al_left2 = al.getbegin2();

		if(al_right1 < al_left1){al_right1 = al.getbegin1();
					 al_left1 = al.getend1();
					}
		if(al_right2 < al_left2){al_right2 = al.getbegin2();
					 al_left2= al.getend2();
					}
	//	size_t cal_right1 = 0;
		size_t cal_left1 = al_left1;
		size_t cal_left2 = al_left2;
		do {
			const pw_alignment * cal = get_al_at_left_end(al.getreference1(),al.getreference2(),cal_left1, cal_left2);
			if(cal!=NULL) {
			size_t cal_right1 = cal -> getend1();
			size_t cal_right2 = cal -> getend2();
			cal_left1 = cal -> getbegin1();
			cal_left2 = cal -> getbegin2();
				if(cal_right1 < cal_left1){cal_right1 = cal -> getbegin1();
					                   cal_left1 = cal -> getend1();
				 		}
				if(cal_right2 < cal_left2){cal_right2 = cal -> getbegin2();
			 			 	   cal_left2= cal -> getend2();
						}
			
			cal_left1 = cal_right1 +1;
			cal_left2 = cal_right2 +1;
			alignment_parts.insert(cal);

			if(cal_right1 > al_right1 || cal_right2 > al_right2) {
				cout << "Small alignment part is to long" << endl;
				exit(1);
			}
			}
			 else { 
			cout<< "There is no alignment!"<<endl;
			exit(1);
			}


		} while(cal_left1 <= al_right1);
		
	}
		if(alignment_parts.size()!=alignments.size()) {
		cout<< " Small parts don't cover the original one!"<<endl;
		exit(1);
		}*/
	}
	void overlap::test_overlap()const{
		for(set<pw_alignment*, compare_pw_alignment>::iterator it1 = alignments.begin(); it1!=alignments.end(); ++it1){	
			pw_alignment * al1 = *it1;
			size_t l1ref1 = al1->getbegin1();
			size_t r1ref1 = al1->getend1();
			size_t l1ref2 = al1->getbegin2();
			size_t r1ref2 = al1->getend2();
			if(l1ref1>r1ref1){
			l1ref1 = al1->getend1();
			r1ref1 = al1->getbegin1();		
			}
			if(l1ref2>r1ref2){
			l1ref2 = al1->getend2();
			r1ref2 = al1->getbegin2();		
			}
			for(set<pw_alignment*, compare_pw_alignment>::iterator it2 = alignments.begin(); it2!=alignments.end(); ++it2){
				pw_alignment * al2 = *it2;
				size_t l2ref1 = al2->getbegin1();
				size_t r2ref1 = al2->getend1();
				size_t l2ref2 = al2->getbegin2();
				size_t r2ref2 = al2->getend2();
				if(l2ref1>r2ref1){
				l2ref1 = al2->getend1();
				r2ref1 = al2->getbegin1();		
				}
				if(l2ref2>r2ref2){
				l2ref2 = al2->getend2();
				r2ref2 = al2->getbegin2();		
				}
			if ((l1ref1 < l2ref1 && l2ref1 < r1ref1) || (l1ref1 < l2ref1 && r2ref1 < r2ref1)){
				cout<< "There are some overlap!"<<endl;
				exit(1);
			}
			else if((l1ref2 < l2ref2 && l2ref2 < r1ref2) || (l1ref2 < l2ref2 && r2ref2 < r2ref2)){
				cout<< "There are some overlap!"<<endl;
				exit(1);
			}
			else{
				cout<< "There is no overlap" << endl;
				exit(0);
			}

		}
	}	
}
	bool overlap::checkAlignments(pw_alignment* const p)const{
		set<pw_alignment*,compare_pw_alignment>::iterator it=alignments.find(p);
		if(it !=alignments.end()) return true;	
		else return false;	
	}
	



size_t overlap::size() const {

//	for(size_t i=0; i<als_on_reference.size(); i++) {
//		cout << " ref " << i << " al start/end points: " << als_on_reference.at(i).size() << endl;
//	}
	return alignments.size();
}
		
	splitpoints::splitpoints(const pw_alignment & p, const overlap & o, const all_data & d):overl(o), newal(p),data(d), split_points(d.numSequences()) {

		}

	splitpoints::~splitpoints(){}

	void splitpoints::find_initial_split_points(size_t sequence, size_t left, size_t right) {
		const multimap<size_t , pw_alignment *> & alignments_on_reference = overl.get_als_on_reference_const(sequence);
		for( multimap<size_t, pw_alignment *>::const_iterator it=alignments_on_reference.lower_bound(left);it!=alignments_on_reference.end(); ++it){
			const pw_alignment * al = it->second;
			bool canbreak = true;
			
			if(al->getreference1() == sequence) {
				size_t alleft;
				size_t alright;
				al->get_lr1(alleft, alright);
				if(alleft > right) canbreak = false;
				if(alleft < left && left <= alright) {
	//				cout<<"here1"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
					insert_split_point(sequence, left);
				}
				if(alleft <= right && right < alright) {
					cout<<"here2"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
					insert_split_point(sequence, right+1);

				}
				if(left < alleft && alleft < right) {
	//				cout<<"here3"<<endl;
	//				cout<<"al: "<<endl;
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						cout <<col <<"\t"<< c1<<"\t"<<c2<<endl;
					}*/

				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point(sequence, alleft);

				}
				if(left < alright && alright < right) {
					cout<<"here4"<<endl;
				//	cout<<"al: "<<endl;
					al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point(sequence, alright+1);
					
				}

			}
		
			if(al->getreference2() == sequence) {
				size_t alleft;
				size_t alright;
				al->get_lr2(alleft, alright);
				if(alleft > right) canbreak = false;
				if(alleft < left && left <= alright) {
				//	cout<<"here5"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();

					insert_split_point(sequence, left);
				
				}
				if(alleft <= right && right < alright) {
					cout<<"here6"<<endl;
				//	al->print();
					insert_split_point(sequence, right+1);
				}
				if(left < alleft && alleft < right) {
					cout<<"here7"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point(sequence, alleft);
				}
				if(left < alright && alright < right) {
					cout<<"here8"<<endl;
				//	cout<<"al: "<<endl;
					al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point(sequence, alright+1);
				}

			
			}

			if(canbreak) break;
		}

	}

	void splitpoints::find_initial_split_point(){




		size_t left1;
		size_t right1;
		size_t left2 ;
		size_t right2;
		newal.get_lr1(left1, right1);
		newal.get_lr2(left2, right2);
		find_initial_split_points(newal.getreference1(), left1, right1);
		find_initial_split_points(newal.getreference2(), left2, right2);
	
		if(newal.getreference1()==newal.getreference2()) {
			if(right1>=left2 && left1 < left2){
				insert_split_point(newal.getreference2(),left2);
				insert_split_point(newal.getreference2(),right1+1);
			}
			if(right2>=left1 && left2 < left1){
				insert_split_point(newal.getreference2(),left1);
				insert_split_point(newal.getreference2(),right2+1);
			}
		}
	
}

	void splitpoints::insert_split_point(size_t sequence, size_t position) {
		pair<set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
		if(insertes.second) {
		//	cout << " new split point " << sequence << " at " << position << endl;
		
			size_t left1;
			size_t right1;
			newal.get_lr1(left1, right1);
			size_t left2;
			size_t right2;
			newal.get_lr2(left2, right2);

			if(newal.getreference1()==sequence) {
				if(left1<position && right1>=position) {
					pw_alignment fp;
					pw_alignment sp;
					size_t spleft2;
					size_t spright2;
					newal.split(true, position, fp, sp);
					sp.get_lr2(spleft2, spright2);
					insert_split_point(newal.getreference2(), spleft2);
				}
			}
			if(newal.getreference2()==sequence) {
				if(left2<position && right2>=position) {
				//	cout<<"HERE!!"<<endl;
					pw_alignment fp;
					pw_alignment sp;
					size_t spleft1;
					size_t spright1;
					newal.split(false, position, fp, sp);
					sp.get_lr1(spleft1, spright1);
					insert_split_point(newal.getreference1(), spleft1);
				}
			}


			const multimap<size_t, pw_alignment *> & alonref = overl.get_als_on_reference_const(sequence);
			for(multimap<size_t, pw_alignment *>::const_iterator it = alonref.lower_bound(position); it!=alonref.end(); ++it) {
				const pw_alignment * al = it-> second;
				if(al->getreference1()!=al->getreference2()) {
					size_t alleft;
					size_t alright;
					al->get_lr_on_reference(sequence, alleft, alright);
					if(alright < position) {
						break;
					}
					if(alleft == position) {
						break;
					} 
					if(alright+1== position) {
						break;
					}
					if(alleft>position){
						break;
					}
			//		al->print(); 
				//	cout<<"position"<<position<<endl;
				//	cout<<"alleft"<<alleft<<endl;
				//	cout<<"alright"<<alright<<endl;
					assert(alleft < position && position <= alright);
					pw_alignment fp;
					pw_alignment sp;
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						cout <<col <<"\t"<< c1<<"\t"<<c2<<endl;
					}*/
				//	cout<<"al length: "<< al->alignment_length() <<endl;
				//	if(al->alignment_length()>1){
						al->split_on_reference(sequence, position, fp, sp);
						if(sp.getreference1()==sequence) {
							insert_split_point(sp.getreference2(), sp.getbegin2());
						} else {
						//	cout<<"Heya!!!"<<endl;
							insert_split_point(sp.getreference1(), sp.getbegin1());
						}
					//}
				} else {
					
					size_t alleft1;
					size_t alright1;
					al->get_lr1(alleft1, alright1);
					size_t alleft2;
					size_t alright2;
					al->get_lr2(alleft2, alright2);
					if(position > alleft1 && position <= alright1) {
						pw_alignment fp;
						pw_alignment sp;
						al->split(true, position, fp, sp);
						size_t spleft2;
						size_t spright2;
						newal.get_lr1(spleft2, spright2);
						insert_split_point(sp.getreference2(), spleft2);
					}
											
					if(position > alleft2 && position <= alright2) {
						pw_alignment fp;
						pw_alignment sp;
					/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						cout <<col <<"\t"<< c1<<"\t"<<c2<<endl;
					}*/
					//	cout<<"al length: "<< al->alignment_length() <<endl;
					//	if(!onlyGapSample(al)){
					//	cout<<"al: "<<endl;
					//	cout<<"........"<<endl;
					//	al->print();
					//	cout<<"position: "<< position<<endl;
							al->split(false, position, fp, sp);
							data.alignment_fits_ref(&fp);
							data.alignment_fits_ref(&sp);
							size_t spleft1;
							size_t spright1;
							newal.get_lr1(spleft1, spright1);
							insert_split_point(sp.getreference1(), spleft1);
					//	}
					}
				}
			}
			
		}
		

	}


	void splitpoints::split_all(set<const pw_alignment*,compare_pw_alignment> & remove_alignments, vector<pw_alignment> & insert_alignments ){
		pw_alignment p1;
		pw_alignment p2;
		for(size_t i = 0; i <data.numSequences();i++){
			const multimap<size_t, pw_alignment *> & als_on_ref = overl.get_als_on_reference_const(i);
			for(set<size_t>::iterator split = split_points.at(i).begin(); split!= split_points.at(i).end(); ++split){
				for(multimap<size_t, pw_alignment*>::const_iterator it =als_on_ref.lower_bound(*split); it!=als_on_ref.end(); ++it){
					pw_alignment * p = it-> second;
					size_t pleft1;
					size_t pright1;
					size_t pleft2;
					size_t pright2;
					p->get_lr1(pleft1,pright1);
					p->get_lr1(pleft2,pright2);
					if( pleft1<*split && pright1>*split){
						remove_alignments.insert(p);
					//	cout<<"remove alignment1: "<<endl;
					//	p->print();
					}
					if( pleft2<*split && pright2>*split){
					//	pw_alignment *np = new pw_alignment(*p);
						remove_alignments.insert(p);
					//	cout<<"remove alignment2: "<<endl;
				//		p->print();
					}
				}		
			}
		}

	//	cout << " in split all "<< remove_alignments.size() << " remove_alignments " << endl;
		vector<pw_alignment>  split_pieces;
		splits(&newal, split_pieces);	
		for(set<const pw_alignment*,compare_pw_alignment>::iterator removed = remove_alignments.begin(); removed != remove_alignments.end(); ++removed){
			splits(*removed, split_pieces);			
		}
	//	cout<<"split_pieces: "<<endl;
	//	for(size_t i = 0 ; i < split_pieces.size();i++){split_pieces.at(i).print();}	
//		cout<<"size of split pieces"<<split_pieces.size()<<endl;
		set<pw_alignment*,compare_pw_alignment> inserted_pieces;
		for(size_t i = 0; i<split_pieces.size();i++) {
				if(!data.alignment_fits_ref(&split_pieces.at(i))) {
					exit(1);
				}
			set<pw_alignment*,compare_pw_alignment>::const_iterator it = inserted_pieces.find(&(split_pieces.at(i)));
			if (overl.checkAlignments(&split_pieces.at(i))){}
			else if ( it != inserted_pieces.end()){}
			else {
				if(!onlyGapSample(&split_pieces.at(i))){	
					insert_alignments.push_back(split_pieces.at(i));
					inserted_pieces.insert(&split_pieces.at(i));
				}
			}
		}
			
			
	}
	void splitpoints::splits(const pw_alignment * p,  vector<pw_alignment> & insert_alignments){
		pw_alignment p1;
		pw_alignment p2;
		size_t left1;
		size_t right1;	
		size_t left2;
		size_t right2;
		p->get_lr1(left1,right1);
		p->get_lr1(left2,right2);
		for(set<size_t>::iterator splitp = split_points.at(p->getreference1()).upper_bound(left1); splitp!= split_points.at(p->getreference1()).end(); splitp++){
			if(right1>=*splitp){
				p->split(true,*splitp,p1,p2);
				p = &p2;
				if(!onlyGapSample(&p1)){	
					insert_alignments.push_back(p1);
				}
			}
			else break;
		}
		if(!onlyGapSample(p)){	
			insert_alignments.push_back(*p);		
		}
	}
	bool splitpoints::onlyGapSample(const pw_alignment* p){
		bool gapOnSample1 = true;
		bool gapOnSample2 = true;
		for(size_t i = 0; i< p->alignment_length();i++){
			char p1char = 'X';
			char p2char = 'X';
			p->alignment_col(i, p1char, p2char);	
			if(p1char !='-' ) gapOnSample1 = false;
			if(p2char !='-' ) gapOnSample2 = false;					
		}
		return gapOnSample1 || gapOnSample2;
		}
	vector<pw_alignment>  splitpoints::get_insert()const{
		return insert_alignments;
	}
		
	model::model(all_data & d): data(d),cost_on_acc (5, vector<double>(data.numAcc(),1)),modification(data.numAcc(), vector<vector<vector<double> > >(data.numAcc(),vector<vector<double> > (6, vector<double>(6,1) ))) {}

	model::~model(){}

	void model::acc_base_frequency(){
//	vector<vector<double> > cost_on_acc (5, vector<double>(data.numAcc(),0));
	vector<size_t> total_base_number_on_acc (data.numAcc(),0);
	for(size_t k = 0; k < data.numSequences(); k++){
		vector<size_t> number(5,0);
		total_base_number_on_acc.at(data.accNumber(k)) += data.getSequence(k).length();
		for(size_t i = 0 ; i< data.getSequence(k).length(); i++ ){
			size_t base = dnastring::base_to_index(data.getSequence(k).at(i));
			number.at(base) ++ ;
				}
		for(size_t j = 0; j< 5; j++){	
			cost_on_acc.at(j).at(data.accNumber(k)) += number.at(j);
				
					}
			}
			for(size_t m=0; m<data.numAcc();m++){
				double numbases = 0;
				for(size_t t=0; t<5;t++){
					numbases+= cost_on_acc.at(t).at(m);
				}
		for(size_t t=0; t<5;t++){
				cost_on_acc.at(t).at(m)= cost_on_acc.at(t).at(m)/numbases;
			//	cout<< "cost of base "<< t <<" on acc "<< m << " is "<< cost_on_acc.at(t).at(m)<<endl;
				}
			}
		}

		void model::alignment_modification(){
	//	vector<vector<vector<vector<double> > > >modification(data.numAcc(), vector<vector<vector<double> > >(data.numAcc(),vector<vector<double> > (6, vector<double>(6,0) )));
		// count:
		for(size_t i = 0 ; i< data.numAlignments() ; i++){
			const pw_alignment & p = data.getAlignment(i);
			size_t acc1 = data.accNumber(p.getreference1());	
			size_t acc2 = data.accNumber(p.getreference2());
			for(size_t j = 0 ; j< p.getsample1().size()/3 ; j++){
				char s1ch;
				char s2ch;
				p.alignment_col(j, s1ch, s2ch);
				size_t s1 = dnastring::base_to_index(s1ch);
				size_t s2 = dnastring::base_to_index(s2ch);
				modification.at(acc1).at(acc2).at(s1).at(s2)+=1.0;
				modification.at(acc2).at(acc1).at(s2).at(s1)+=1.0;
			}
		}

		// sum up and divide to get frequencies
		for(size_t acc1=0; acc1 < data.numAcc(); ++acc1) {
			for(size_t acc2=0; acc2 < data.numAcc(); ++acc2) {
				for(size_t j=0; j<6; ++j) {
					double sum = 0;
					for(size_t k=0; k<6; ++k) {

						cout << "k " << k << " num " << modification.at(acc1).at(acc2).at(j).at(k) << endl;
						sum+=modification.at(acc1).at(acc2).at(j).at(k);
					}	
					for(size_t k=0; k<6; ++k) {
						modification.at(acc1).at(acc2).at(j).at(k)/=sum;
						cout << " simple model cost: acc " << acc1 << " to " << acc2 << " modify " << j << " to " << k << " costs " << -log2(modification.at(acc1).at(acc2).at(j).at(k)) << " bits " << endl; 
					
					}
				}

			}
		}
	}

	void model::cost_function( pw_alignment& p) const {
//		cout << " cf2 " << endl;
		vector<double> cost_on_sample(2);
		vector<double> modify_cost(2);
		double c1;
		double c2;
		double m1;
		double m2;
		cost_function(p, c1, c2, m1, m2);
		cost_on_sample.at(0) = c1;
		cost_on_sample.at(1) = c2;
		modify_cost.at(0) = m1;
		modify_cost.at(1) = m2;
		p.set_cost(cost_on_sample, modify_cost);
	}

	void model::cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2) const {
		cout << " cf " << endl;
		p.print(); // TODO 

		vector<double> cost_on_sample(2,0);
		vector<double> modify_cost(2,0);
		vector<vector<size_t> >al_base_number(2,vector<size_t>(6,0));
		vector<vector<double> >sequence_cost(2,vector<double>(5,0));
		vector<vector<double> >create_cost(2,vector<double>(5,0));		
	//	vector<vector<vector<double> > > modification_cost(2,vector<vector<double> >(6,vector<double>(6,0)));
		vector<vector<double> > modification_number(6, vector<double>(6,0)); // in (i,j): store modifications from 1 to 2
	//	vector<double> cost_on_sample (2,0);
		for(size_t i = 0 ; i < 5; i++){
		 sequence_cost.at(0).at(i)= -log2(cost_on_acc.at(i).at(data.accNumber(p.getreference1())));
		 sequence_cost.at(1).at(i)= -log2(cost_on_acc.at(i).at(data.accNumber(p.getreference2())));
		}
//		cout << " alsample size " << p.getsample1().size() << endl;
		for(size_t i = 0 ; i< p.getsample1().size()/3; i++ ){
			char s1ch;
			char s2ch;
			p.alignment_col(i, s1ch, s2ch);
			size_t s1 = dnastring::base_to_index(s1ch);
			size_t s2 = dnastring::base_to_index(s2ch);
			al_base_number.at(0).at(s1) ++ ;
			al_base_number.at(1).at(s2) ++ ;

			modification_number.at(s1).at(s2)+=1.0;
		}
		for(size_t j=0; j<5; j++){	
			create_cost.at(0).at(j)= al_base_number.at(0).at(j)*sequence_cost.at(0).at(j);
			cost_on_sample.at(0)=cost_on_sample.at(0)+ create_cost.at(0).at(j);
//			cout<<"creat cost of "<< dnastring::index_to_base(j)<<"  on reference1 is "<<  (al_base_number.at(0).at(j))*(sequence_cost.at(0).at(j))<<endl;
			create_cost.at(1).at(j)= al_base_number.at(1).at(j)*sequence_cost.at(1).at(j);
			cost_on_sample.at(1)=cost_on_sample.at(1)+ create_cost.at(1).at(j);			
//			cout<<"creat cost of "<< dnastring::index_to_base(j)<<"  on reference2 is "<<  al_base_number.at(1).at(j)*sequence_cost.at(1).at(j)<<endl;
		}
//			cout<<"creat cost of the alignment on sample 1: "<< cost_on_sample.at(0);
//			cout<<"creat cost of the alignment on sample 2: "<< cost_on_sample.at(1);	
		
		for(size_t j=0; j<6; j++){	
			for (size_t k= 0; k<6; k++){
			//	modification_cost.at(0).at(j).at(k) = -log2(modification.at(data.accNumber(p.getreference1())).at(data.accNumber(p.getreference2())).at(j).at(k));
				modify_cost.at(0) +=  -log2(modification.at(data.accNumber(p.getreference1())).at(data.accNumber(p.getreference2())).at(j).at(k)) * modification_number.at(j).at(k);
			//	modification_cost.at(1).at(j).at(k) = -log2(modification.at(data.accNumber(p.getreference2())).at(data.accNumber(p.getreference1())).at(j).at(k));
				modify_cost.at(1) +=  -log2(modification.at(data.accNumber(p.getreference2())).at(data.accNumber(p.getreference1())).at(j).at(k)) * modification_number.at(k).at(j);
			}
		}	
	
		c1 = cost_on_sample.at(0);
		c2 = cost_on_sample.at(1);
		m1 = modify_cost.at(0);
		m2 = modify_cost.at(1);
//		cout << " create 2 " << c2 << " m1 " << m1 << endl; 

	}
	void model::gain_function(const pw_alignment& p, double & g1, double & g2) const {

		double c1;
		double c2;
		double m1;
		double m2;
		cost_function(p, c1, c2, m1, m2);
		/*
				computing the gain:
				without using the alignment: we pay c1 + c2
				using the alignment: we pay c1 + m1 (or c2 + m2)
				the gain of using the alignment is: c2 - m1 

		*/
		cout << "in gain: create 2 " << c2 << " m1 " << m1 << endl; 

		g1 = c2 - m1;
		g2 = c1 - m2;

	}

	void model::train() {
		acc_base_frequency();
		alignment_modification();
	}

	mc_model::mc_model(all_data & d):data(d), sequence_successive_bases(data.numAcc()), create_cost(data.numAcc(),vector<double>(5,1)),mod_cost(data.numAcc(),vector<map<string, vector<double> > >(data.numAcc())){
		size_t numberOfPowers = NUM_DELETE;
		if(NUM_KEEP>numberOfPowers) {
			numberOfPowers = NUM_KEEP;
		}
		powersOfTwo = vector<size_t>(numberOfPowers, 1);
		for(size_t i=1; i< numberOfPowers; ++i) {
			powersOfTwo.at(i) = powersOfTwo.at(i-1)*2;
		
		}


	}
	mc_model::~mc_model(){}
	void mc_model::markov_chain(){
		size_t level = 2;
	//	vector<vector<double> >create_cost(data.numAcc(),vector<double>(5,1));		
		for(size_t k = 0; k < data.numSequences(); k++){
			size_t acc = data.accNumber(k);
			for(size_t i = 0 ; i< data.getSequence(k).length(); i++ ){
				char schr = data.getSequence(k).at(i);
				size_t s = dnastring::base_to_index(schr);
				stringstream context;		
				for(size_t j = level; j>0; j--){
					char rchr;
					signed long temp = i - level;
					if(temp >= 0){
						rchr = data.getSequence(k).at(i-j);
					}else{
						rchr = 'A';

					}
					context << rchr;
				}
				string seq;
				context>>seq;
				map <string, vector<double> >::iterator it= sequence_successive_bases.at(acc).find(seq);
				if(it==sequence_successive_bases.at(acc).end()) {
					sequence_successive_bases.at(acc).insert(make_pair(seq, vector<double>(5,1)));
					it= sequence_successive_bases.at(acc).find(seq);
				}
				it->second.at(s)++;
				
			}
					}
		for(size_t i =0 ; i< data.numAcc(); i++){
			int total;
			for(map <string, vector<double> >::iterator it= sequence_successive_bases.at(i).begin();it!=sequence_successive_bases.at(i).end();it++){
				string seq = it->first;
				vector<double> & base = sequence_successive_bases.at(i).at(seq);
				for(size_t j = 0; j<5;j++){
					total += base.at(j);
			//	cout<<"base: "<<base.at(0)<<endl;
				}
		//	cout<<"totalnumber of bases for each seq:  "<< total<<endl;
				for(size_t k=0; k<5;k++){
				base.at(k) = -log2(base.at(k)/total);
				create_cost.at(i).at(k) += base.at(k);
				}
			//	cout<<"cr_cost: "<<base.at(0) <<endl;		
			}
		}
			


	}
	
	void mc_model::markov_chain_alignment(){
		counting_functor functor(data);
		for(size_t k = 0; k < data.numAlignments(); k++){
			const pw_alignment & p = data.getAlignment(k);
			computing_modification_oneToTwo(p,functor);
			computing_modification_twoToOne(p,functor);
			}
		functor.total_context();
		for(size_t i = 0; i< data.numAcc();i++){
			for(size_t j = 0; j < data.numAcc();j++){
				for(map <string, vector<double> >::const_iterator it= functor.get_context(i,j).begin();it!=functor.get_context(i,j).end();it++){
					string seq1 = it->first;
					const vector<double> & base = functor.get_context(i,j).at(seq1);
				/*	cout<<"base is: "<<endl;
					for(size_t a= 0; a< base.size();a++){
						cout<< "base at "<<  print_modification_character(a)<<" is "<<base.at(a)<<endl;
					}*/
				/*	cout<< "context is: "<<endl;
					for(size_t m = 0 ; m < seq1.size() ; m++){
						cout<< int(seq1.at(m))<<endl;
					}*/
				//	cout<<"the total number of happening the above context between "<<i<<" and "<<j<<" is "<< functor.get_total(i,j,seq1) <<endl;
					for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++){
				//		cout<<"The number of happening "<< print_modification_character(k) << " between acc " <<i<< " and acc " << j << " after a certain context is "<< base.at(k)<<endl;
						cout<<"In MC model, the cost value of modifying a pattern on acc " << i << " to " <<print_modification_character(k) << " on acc " << j << " is " << -log2(base.at(k)/functor.get_total(i,j,seq1)) << " bits " << endl; 
						map <string, vector<double> >::iterator it1= mod_cost.at(i).at(j).find(seq1);
						if(it1==mod_cost.at(i).at(j).end()) {
							mod_cost.at(i).at(j).insert(make_pair(seq1, vector<double>((NUM_DELETE+NUM_KEEP+10),1)));
							it1= mod_cost.at(i).at(j).find(seq1);
						}
						it1->second.at(k)=-log2(base.at(k)/functor.get_total(i,j,seq1));
					}
				}		
			}			
		}
	}

	void mc_model::cost_function( pw_alignment& p) const {
		vector<double> cost_on_sample(2);
		vector<double> modify_cost(2);
		double c1;
		double c2;
		double m1;
		double m2;
		cost_function(p, c1, c2, m1, m2);
		cost_on_sample.at(0) = c1;
		cost_on_sample.at(1) = c2;
		modify_cost.at(0) = m1;
		modify_cost.at(1) = m2;
		p.set_cost(cost_on_sample, modify_cost);
	}

	void mc_model::cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2)const {
		counting_functor  f(data);
		size_t level = 2;
		computing_modification_oneToTwo(p,f);
		computing_modification_twoToOne(p,f);
		f.total_context();
		vector<double> cost_on_sample(2,0);
		vector<double> modify_cost(2,0);
		size_t acc1 = data.accNumber(p.getreference1());
		size_t acc2 = data.accNumber(p.getreference2());
		vector <map < string, vector<double> > >alignment_successive_bases(2);
		char s1chr;
		char s2chr;
		for(size_t i = level; i<p.alignment_length(); i++){
			p.alignment_col(i, s1chr, s2chr);				
			size_t s1 = dnastring::base_to_index(s1chr);
			size_t s2 = dnastring::base_to_index(s2chr);
			stringstream context1;	
			stringstream context2;	
				for(size_t j = level; j>0; j--){
					char r1chr;
					char r2chr;
					p.alignment_col(i-j, r1chr, r2chr);				
					context1 << r1chr;
					context2 << r2chr;
			}
			string seq1;
			string seq2;
			context1>>seq1;
			context2>>seq2;
			map <string, vector<double> >::iterator it= alignment_successive_bases.at(0).find(seq1);
				if(it==alignment_successive_bases.at(0).end()) {
					alignment_successive_bases.at(0).insert(make_pair(seq1, vector<double>(6,1)));
					it= alignment_successive_bases.at(0).find(seq1);
				}
				it->second.at(s1)++;
			map <string, vector<double> >::iterator it1= alignment_successive_bases.at(1).find(seq2);
				if(it1==alignment_successive_bases.at(1).end()) {
					alignment_successive_bases.at(1).insert(make_pair(seq2, vector<double>(6,1)));
					it1= alignment_successive_bases.at(1).find(seq2);
				}
				it1->second.at(s2)++;			
		}	
		for(map <string, vector<double> >::iterator it= alignment_successive_bases.at(0).begin();it!=alignment_successive_bases.at(0).end();it++){
			string seq = it->first;
			int total1 = 1;
			vector<double> & base = alignment_successive_bases.at(0).at(seq);
			for(size_t j = 0; j<6;j++){
				total1 += base.at(j);
			}
			for(size_t m = 0; m<6;m++){
				cost_on_sample.at(0) += -log2(base.at(m)/total1);
			}

		}
		for(map <string, vector<double> >::iterator it= alignment_successive_bases.at(1).begin();it!=alignment_successive_bases.at(1).end();it++){
			string seq = it->first;
			int total2 = 1;
			vector<double> & base = alignment_successive_bases.at(1).at(seq);
			for(size_t j = 0; j<6;j++){
				total2 += base.at(j);
			}
			for(size_t m = 0; m<6;m++){
				cost_on_sample.at(1) += -log2(base.at(m)/total2);
			}
		}
		for(map <string, vector<double> >::const_iterator it= f.get_context(acc1,acc2).begin();it!=f.get_context(acc1,acc2).end();it++){
					string seq1 = it->first;
					const vector<double> & base = f.get_context(acc1,acc2).at(seq1);
					for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++){
						modify_cost.at(0) +=-log2(base.at(k)/f.get_total(acc1,acc2,seq1));
					}
				}

		for(map <string, vector<double> >::const_iterator it= f.get_context(acc2,acc1).begin();it!=f.get_context(acc2,acc1).end();it++){
					string seq1 = it->first;
					const vector<double> & base = f.get_context(acc2,acc1).at(seq1);
					for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++){
						modify_cost.at(1) +=-log2(base.at(k)/f.get_total(acc2,acc1,seq1));
					}
				}

		c1 = cost_on_sample.at(0);
		c2 = cost_on_sample.at(1);
		m1 = modify_cost.at(0);
		m2 = modify_cost.at(1);
		
	}
	void mc_model::gain_function(const pw_alignment& p, double & g1, double & g2)const {
		double c1;
		double c2;
		double m1;
		double m2;
		cost_function(p, c1, c2, m1,m2);

		cout << " gain function c2 " << c2 << " m1 " << m1 << " gain " << g1 << endl; 
		g1 = c2 - m1;
		g2 = c1 - m2;

	}
	
	void mc_model::train(){
		markov_chain();
		markov_chain_alignment();
		
	}


	char mc_model::modification_character(int modify_base, int num_delete, int insert_base, int num_keep)const {
		//return the enc
		if(num_delete != -1){
		//	cout<< "there is a delete of length "<< num_delete <<endl; 
			return NUM_DELETE+num_delete;
		}
		if(num_keep != -1){
		//	cout<< "there is a keep of length" << num_keep << endl;
			return NUM_KEEP+num_keep;
		}
		if(modify_base != -1) {
		//	cout<< "there is a modification at " << dnastring::index_to_base(modify_base) << endl;
			return modify_base;
		}
		if(insert_base != -1){
		//	cout<< " there is a insertion at " << dnastring::index_to_base(insert_base) <<endl;
			return insert_base + NUM_KEEP + NUM_DELETE + 5;
		}
		assert (false);
		return -1;
	}
		string mc_model::print_modification_character(char enc){
			int modify_base = -1;
			int num_delete =-1;
			int insert_base = -1;
			int num_keep = -1;
			modification(enc, modify_base, num_delete, insert_base, num_keep);

			stringstream s;
			if(num_delete != -1){
				s<< " a delete of length " <<  num_delete; 
			}
			if(num_keep != -1){
				s<< " a keep of length" << num_keep;
			}
			if(modify_base != -1) {
				s<< " a modification at " << dnastring::index_to_base(modify_base);
			}
			if(insert_base != -1){
				s<< " a insertion at " << dnastring::index_to_base(insert_base);
			}
			return s.str();
		}

	void mc_model::modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep) {
		modify_base = -1;
		num_delete =-1;
		insert_base = -1;
		num_keep = -1;

		if(enc < 5) {
			modify_base = enc;
			return;
		}

		if(enc < 5 + NUM_DELETE) {
			num_delete = powersOfTwo.at(enc - 5);
			return;
		}
		if (enc < 5 + NUM_DELETE + NUM_KEEP){
			num_keep = powersOfTwo.at(enc-NUM_DELETE-5);
			return;
		}
		if(enc< 5+ 5 + NUM_KEEP + NUM_DELETE){
			insert_base = enc- 5 - 5 - NUM_KEEP - NUM_DELETE;
		}
	}
	void mc_model::computing_modification_oneToTwo(const pw_alignment & p, abstract_context_functor & functor)const{
		size_t level =2;
		string seq = "" ;
		for(size_t j = 0; j < level; j++){
			seq+=modification_character(1,-1,-1,-1);
		}
		for (size_t i = 0; i< p.alignment_length(); i++){
			size_t n = 0;
			int modify_base =-1;
			int num_delete=-1;
			int insert_base=-1;
			int num_keep=-1;
			string seq1(" ",level+1);
			char seq2;
			for(size_t w = level; w>0 ;w--){	
				seq1.at(level-w)=seq.at(seq.size()-w);
			}
			char s1chr;
			char s2chr;
			size_t s1;
			size_t s2;
			p.alignment_col(i, s1chr, s2chr);				
			s1 = dnastring::base_to_index(s1chr);
			s2 = dnastring::base_to_index(s2chr);
		//	cout<<"context at position zero is: "<< int(seq1.at(0))<<endl;
			seq1.at(level)=s1;
			if(s1 == s2){
				size_t klength = 0;
				size_t l1 = 0;													
				for(size_t j = i; j<p.alignment_length(); j++){
					char q1chr;
					char q2chr;
					p.alignment_col(j, q1chr, q2chr);				
					size_t q1 = dnastring::base_to_index(q1chr);
					size_t q2 = dnastring::base_to_index(q2chr);
					if(q1 == q2){
						klength +=1;
					}else {	
						break;
					}
				}
			/*	if(i<51){
					cout<<"keep length at "<< i << " is  "<< klength<<endl;
				}*/
				n = klength-1;
				for (size_t m = 0; m< powersOfTwo.size();m++){
					if((klength & powersOfTwo.at(m)) != 0){
						num_keep=m;
						seq += modification_character(modify_base,num_delete,insert_base,num_keep);
						l1++;
					}
				}
				seq2 = seq.at(seq.size()-l1);
			//	cout<<"last pattern: "<< int(seq2)<<endl;
				functor. see_context(&p,i,seq1,seq2);
			}else{
				if((s1!=5) & (s2!=5)){
					modify_base = s1;
				/*	if(i<51){
						cout<<"modification at "<< i << " is  "<<s1 <<endl;
					}*/
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(&p,i,seq1,seq2);
				}
				if(s1 == 5){
					insert_base = s2;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);						
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(&p,i,seq1,seq2);						
				}
				if(s2 == 5){
					size_t dlength = 0;
					size_t l1 = 0;												
					for(size_t j = i; j < p.alignment_length(); j++){
						char q1chr;
						char q2chr;
						p.alignment_col(j, q1chr, q2chr);				
						size_t q1 = dnastring::base_to_index(q1chr);
						size_t q2 = dnastring::base_to_index(q2chr);
						if(q2 == 5 ){
							dlength +=1;
						}else {
							break;
						}
					}
					n = dlength-1;
					for (size_t m = 0; m< powersOfTwo.size();m++){
						if((dlength & powersOfTwo.at(m)) != 0){
							num_delete = m;
							seq += modification_character(modify_base,num_delete,insert_base,num_keep);		
							l1++;							
						}
					}
					seq2 = seq.at(seq.size()-l1);
					functor. see_context(&p,i,seq1,seq2);	
				}

			}
			i=i+n;
			//	uint32_t initial = 1;
			//	for (uint32_t m = 0; m <32 ; m++)
			//		uint32_t power_of_two = initial << m; 
			//		if((klength & power_of_two) !=0)
		}
	/*	cout<<"The alignment is: "<<endl;
		p.print();
		cout<<"encoded sequence from one to two is: "<<endl;
		for(size_t m = 0; m < seq.size(); m ++){
			cout<< int(seq.at(m))<<endl;
			if (m > 50) break;
		}*/
	}

	
	void mc_model::computing_modification_twoToOne(const pw_alignment & p, abstract_context_functor & functor)const{
		size_t level =2;
		string seq = "" ;
		for(size_t j = 0; j < level; j++){
			seq+=modification_character(1,-1,-1,-1);
		}
		for (size_t i = 0; i< p.alignment_length(); i++){
			int modify_base =-1;
			int num_delete=-1;
			int insert_base=-1;
			int num_keep=-1;
			size_t n = 0;			
			string seq1(" ",level+1);
			char seq2;
			for(size_t w = level; w>0 ;w--){
				seq1.at(level-w)=seq.at(seq.size()-w);
			}
			char s1chr;
			char s2chr;
			size_t s1;
			size_t s2;
			p.alignment_col(i, s1chr, s2chr);				
			s1 = dnastring::base_to_index(s1chr);
			s2 = dnastring::base_to_index(s2chr);
			seq1.at(level)=s1;
			if(s1 == s2){
				size_t klength = 0;
				size_t l1 = 0;													
				for(size_t j = i; j<p.alignment_length(); j++){
					char q1chr;
					char q2chr;
					p.alignment_col(j, q1chr, q2chr);				
					size_t q1 = dnastring::base_to_index(q1chr);
					size_t q2 = dnastring::base_to_index(q2chr);
					if(q1 == q2){
						klength +=1;
					}else {	
						break;
					}
				}
				n = klength-1;
			/*	if(i<51){
					cout<<"keep length at "<< i << " is  "<< klength<<endl;
				}*/
				for (size_t m = 0; m< powersOfTwo.size();m++){
					if((klength & powersOfTwo.at(m)) != 0){
						num_keep=m;
						seq += modification_character(modify_base,num_delete,insert_base,num_keep);
						l1++;
					}
				}
				seq2 = seq.at(seq.size()-l1);
				functor. see_context(&p,i,seq1,seq2);				
			}else{
				if((s1!=5) & (s2!=5)){
					modify_base = s2;
				/*	if(i<51){
						cout<<"modification at "<< i << " is  "<<s2 <<endl;
					}*/
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(&p,i,seq1,seq2);
				}
				if(s2 == 5){
					insert_base = s1;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);						
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(&p,i,seq1,seq2);						
				}
				if(s1 == 5){
					size_t dlength = 0;
					size_t l1 = 0;												
					for(size_t j = i; j < p.alignment_length(); j++){
						char q1chr;
						char q2chr;
						p.alignment_col(j, q1chr, q2chr);				
						size_t q1 = dnastring::base_to_index(q1chr);
						size_t q2 = dnastring::base_to_index(q2chr);
						if(q1 == 5 ){
							dlength +=1;
						}else {
							break;
						}
					}
					n = dlength-1;
					for (size_t m = 0; m< powersOfTwo.size();m++){
						if((dlength & powersOfTwo.at(m)) != 0){
							num_delete = m;
							seq += modification_character(modify_base,num_delete,insert_base,num_keep);		
							l1++;							
						}
					}
					seq2 = seq.at(seq.size()-l1);
					functor. see_context(&p,i,seq1,seq2);	
				}
			}
			i=i+n;
		}		
	/*	cout<<"The alignment is: "<<endl;
		p.print();
		cout<<"encoded sequence from two to one is: "<<endl;
		for(size_t m = 0; m < seq.size(); m ++){
			cout<< int(seq.at(m))<<endl;
			if (m > 50) break;
		}*/

	}
	abstract_context_functor::abstract_context_functor(){
	
	}
	void abstract_context_functor::see_context(const pw_alignment * p, size_t pos, string context, char last_char){
		
	}
	counting_functor::counting_functor(all_data & d):data(d), successive_modification(data.numAcc(),vector<map<string, vector<double> > >(data.numAcc())),total(data.numAcc(),vector<map<string, double > >(data.numAcc())) {}

	void counting_functor::see_context(const pw_alignment * p, size_t pos, string context, char last_char){
		size_t acc1 = data.accNumber(p->getreference1());
		size_t acc2 = data.accNumber(p->getreference2());
		map <string, vector<double> >::iterator it1= successive_modification.at(acc1).at(acc2).find(context);
			if(it1==successive_modification.at(acc1).at(acc2).end()) {
				successive_modification.at(acc1).at(acc2).insert(make_pair(context, vector<double>((NUM_DELETE+NUM_KEEP+10),1)));
				it1= successive_modification.at(acc1).at(acc2).find(context);
			}
			it1->second.at(last_char)++;
	}
	void counting_functor::total_context(){
		for(size_t i = 0; i < data.numAcc(); i++){
			for(size_t j = 0; j<data.numAcc(); j++){
				for(map<string, vector<double> >::iterator it = successive_modification.at(i).at(j).begin(); it!= successive_modification.at(i).at(j).end();it++){
					string context = it->first;
					map<string, double >::iterator it1=total.at(i).at(j).find(context);
					if(it1 == total.at(i).at(j).end()){
						total.at(i).at(j).insert(make_pair(context,0));
						it1=total.at(i).at(j).find(context);
					}
					for(size_t k = 0; k < NUM_DELETE+NUM_KEEP+10; k++){
					it1->second += it->second.at(k);	
					}
				}
			}
		}


	}

	double counting_functor::get_total(size_t acc1, size_t acc2, string context)const{
		return total.at(acc1).at(acc2).at(context);
	}

			
	
	
	
	const map<string, vector<double> > & counting_functor::get_context(size_t acc1, size_t acc2)const{
		return successive_modification.at(acc1).at(acc2);
	}

