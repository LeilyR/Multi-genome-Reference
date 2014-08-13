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
		case 'A':
		return 'T';
		break;
		case 'T':
		return 'A';
		break;
		case 'C':
		return 'G';
		break;
		case 'G':
		return 'C';
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
						start1 = incl_end1;
						incl_end1 = tmp;
					}


					parts.clear();
					strsep(aline2, " ", parts);
					string acc2;
					string name2;
					name_split(parts.at(1), acc2, name2);
					size_t start2 = atoi(parts.at(2).c_str());
					size_t size2 = atoi(parts.at(3).c_str());
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
						start2 = incl_end2;
						incl_end2 = tmp;
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

//const multimap<size_t, size_t> & all_data::getAlOnRefMap(size_t seq_idx) const {
//	return als_on_reference.at(seq_idx);
//}

size_t all_data::numSequences() const {
	return sequences.size();
}

size_t all_data::numAlignments() const {
	return alignments.size();
}



void all_data::insert_sequence(const string & acc, const string & seq_name, const string & dna) {
	size_t seq_id = sequences.size();
	sequences.push_back(dnastring(dna));
	sequence_names.push_back(seq_name);
	map<string, vector<size_t> >::iterator findacc = acc_sequences.find(acc);
	if(findacc == acc_sequences.end()) {
		acc_sequences.insert(make_pair(acc, vector<size_t>(0)));
		findacc = acc_sequences.find(acc);
	}
	findacc->second.push_back(seq_id);
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

	bool al1fw = true;
	bool al2fw = true;
	size_t al1at = al->getbegin1();
	size_t al2at = al->getbegin1();
	if(al->getbegin1() > al->getend1()) {
		al1fw = false;
		al1at = al->getend1();
	}

	if(al->getbegin2() > al->getend2()) {
		al2fw = false;
		al2at = al->getend1();
	}



	for(size_t col=0; col < al->alignment_length(); ++col) {
		char al1char='X';
		char al2char='X';
		al->alignment_col(col, al1char, al2char);	
		
		
		if(al1char!='-') { 
			char rchar = ref1.at(al1at);
			if(al1fw) {
				al1at++;
			} else {
				al1at--;
				rchar = dnastring::complement(rchar);
			}
			if(rchar != al1char) {
				cerr << "Warning: alignment test failed. Alignment has " << al1char << " where reference sequence has " << rchar << endl;
				return false;
			}
		}
		if(al2char!='-') {
			char rchar = ref1.at(al2at);
			if(al2fw) {
				al2at++;
			} else {
				al2at--;
				rchar = dnastring::complement(rchar);
			}
			if(rchar != al2char) {
				cerr << "Warning: alignment test failed. Alignment has " << al2char << " where reference sequence has " << rchar << endl;
				return false;
			}
		}
	}

	if(al1fw) {
		size_t al1e = al1at-1;
		if(al1e!=al->getend1()) {
			cerr << "Warning: alignment end wrong. Should be " << al1e << " but is " << al->getend1() << endl;
			return false;
		} 
	
	} else {
		size_t al1b = al1at+1;
		if(al1b!=al->getbegin1()) {
			cerr << "Warning: rev alignment begin wrong. Should be " << al1b << " but is " << al->getbegin1() << endl;
			return false;
		}
	}
	if(al2fw) {
		size_t al2e = al2at-1;
		if(al2e!=al->getend2()) {
			cerr << "Warning: alignment end wrong. Should be " << al2e << " but is " << al->getend2() << endl;
			return false;
		} 
	
	} else {
		size_t al2b = al2at+1;
		if(al2b!=al->getbegin2()) {
			cerr << "Warning: rev alignment begin wrong. Should be " << al2b << " but is " << al->getbegin2() << endl;
			return false;
		}
	}

	return true;
}



overlap::overlap(all_data & d): data(d), als_on_reference(d.numSequences()){

}

void overlap::split_partial_overlap(const pw_alignment * new_alignment, set<pw_alignment*, compare_pw_alignment> & remove_alignments, vector<pw_alignment> & insert_alignments) const {

	pw_alignment p1;
	pw_alignment p2;
	pw_alignment p3;
	pw_alignment p4;

	const multimap<size_t , pw_alignment *> & alignments_on_reference1 = als_on_reference.at(new_alignment->getreference1());

	const multimap<size_t , pw_alignment *> & alignments_on_reference2 = als_on_reference.at(new_alignment->getreference2());



	size_t leftinsert1 = new_alignment->getbegin1();
	size_t rightinsert1 = new_alignment->getend1();
	if(new_alignment->getend1() < leftinsert1) {
		leftinsert1 = new_alignment->getend1();
		rightinsert1 = new_alignment->getbegin1();
	}
	multimap<size_t, pw_alignment*>::const_iterator allalignit1 = alignments_on_reference1.lower_bound(leftinsert1);


	if(allalignit1 == alignments_on_reference1.end()) allalignit1 = alignments_on_reference1.begin();

	for(; allalignit1 != alignments_on_reference1.end(); ++allalignit1) {
		const pw_alignment * al1 = allalignit1->second;
		size_t alleft1;
		size_t alright1;
		alleft1 = al1->getbegin1();
		alright1 = al1->getend1();
		if(alleft1 > al1->getend1()) {
		alleft1 = al1->getend1();
		alright1 = al1->getbegin1();
		}
			
		if( allalignit1->second->getreference1() == new_alignment->getreference1()){

		if(leftinsert1 < alright1 && rightinsert1 > alleft1) {
			if (leftinsert1 > alleft1 && rightinsert1 > alright1){
				al1->split (true,leftinsert1,p1,p2);
				new_alignment->split(true,alright1,p3,p4);

				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p2,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}
			if (leftinsert1 < alleft1 && rightinsert1 < alright1){
				new_alignment->split (true,alleft1,p1,p2);
				al1->split(true,rightinsert1,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p2,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);

			}

			if (leftinsert1 > alleft1 && rightinsert1 < alright1){
				al1->split (true,leftinsert1,p1,p2);
				p2.split(true,rightinsert1,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);

			}
			if (leftinsert1 < alleft1 && rightinsert1 > alright1){
				new_alignment->split (true,alleft1,p1,p2);
				p2.split(true,alright1,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);

			}
		
		}

		else break;	
		}
		else {

		if(leftinsert1 < alright1 && rightinsert1 > alleft1) {
			if (leftinsert1 > alleft1 && rightinsert1 > alright1){
				al1->split (true,leftinsert1,p1,p2);
				new_alignment -> split(false,alright1,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p2,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert1 < alleft1 && rightinsert1 < alright1){
				new_alignment->split (true,alleft1,p1,p2);
				al1 -> split(false,rightinsert1,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p2,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert1 > alleft1 && rightinsert1 < alright1){
				al1->split (true,leftinsert1,p1,p2);
				p2.split(false,rightinsert1,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert1 < alleft1 && rightinsert1 > alright1){
				new_alignment->split (true,alleft1,p1,p2);
				p2.split(false,alright1,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

	
		}



		
	}
}


	
	size_t leftinsert2 = new_alignment->getbegin2();
	size_t rightinsert2 = new_alignment->getend2();
	if(new_alignment->getend2() < leftinsert2) {
		leftinsert2 = new_alignment->getend2();
		rightinsert2 = new_alignment->getbegin2();
	}
	multimap<size_t, pw_alignment*>::const_iterator allalignit2 = alignments_on_reference2.lower_bound(leftinsert2);


	if(allalignit2 == alignments_on_reference2.end()) allalignit2 = alignments_on_reference2.begin();

	for(; allalignit2 != alignments_on_reference2.end(); ++allalignit2) {
		const pw_alignment * al2 = allalignit2->second;
		size_t alleft2;
		size_t alright2;
		alleft2 = al2->getbegin2();
		alright2 = al2->getend2();
		if(alleft2 > al2->getend2()) {
		alleft2 = al2->getend2();
		alright2 = al2->getbegin2();
		}
			
		if( allalignit2->second->getreference2() == new_alignment->getreference2()){

		if(leftinsert2 < alright2 && rightinsert2 > alleft2) {
			if (leftinsert2 > alleft2 && rightinsert2 > alright2){
			al2->split (true,leftinsert2,p1,p2);
			new_alignment-> split(true,alright2,p3,p4);
			split_partial_overlap(&p1,remove_alignments,insert_alignments);
			split_partial_overlap(&p2,remove_alignments,insert_alignments);
			split_partial_overlap(&p3,remove_alignments,insert_alignments);
			split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert2 < alleft2 && rightinsert2 < alright2){
			new_alignment->split (true,alleft2,p1,p2);
			al2 -> split(true,rightinsert2,p3,p4);
			split_partial_overlap(&p1,remove_alignments,insert_alignments);
			split_partial_overlap(&p2,remove_alignments,insert_alignments);
			split_partial_overlap(&p3,remove_alignments,insert_alignments);
			split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert2 > alleft2 && rightinsert2 < alright2){
			al2->split (true,leftinsert2,p1,p2);
			p2.split(true,rightinsert2,p3,p4);
			split_partial_overlap(&p1,remove_alignments,insert_alignments);
			split_partial_overlap(&p3,remove_alignments,insert_alignments);
			split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert2 < alleft2 && rightinsert2 > alright2){
			new_alignment->split (true,alleft2,p1,p2);
			p2.split(true,alright2,p3,p4);
			split_partial_overlap(&p1,remove_alignments,insert_alignments);
			split_partial_overlap(&p3,remove_alignments,insert_alignments);
			split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

	
		}

		else break;	
	}
		else {

		if(leftinsert2 < alright2 && rightinsert2 > alleft2) {
			if (leftinsert2 > alleft2 && rightinsert2 > alright2){
				al2->split (true,leftinsert2,p1,p2);
				new_alignment -> split(false,alright2,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p2,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert2 < alleft2 && rightinsert2 < alright2){
				new_alignment-> split (true,alleft2,p1,p2);
				al2 -> split(false,rightinsert2,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p2,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert2 > alleft2 && rightinsert2 < alright2){
				al2->split (true,leftinsert2,p1,p2);
				p2.split(false,rightinsert2,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

			if (leftinsert2 < alleft2 && rightinsert2 > alright2){
				new_alignment->split (true,alleft2,p1,p2);
				p2.split(false,alright2,p3,p4);
				split_partial_overlap(&p1,remove_alignments,insert_alignments);
				split_partial_overlap(&p2,remove_alignments,insert_alignments);
				split_partial_overlap(&p3,remove_alignments,insert_alignments);
				split_partial_overlap(&p4,remove_alignments,insert_alignments);
			}

	
		}



		
	}
  }

	
}

	
overlap::~overlap(){}

void overlap::remove_alignment(pw_alignment *remove){
	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqrb1 =
	als_on_reference.at(remove->getreference1()).equal_range(remove->getbegin1());
	for(multimap<size_t, pw_alignment*>::iterator it = eqrb1.first; it!=eqrb1.second; ++it) {
		if ( &(it->second) == &remove ) {
			als_on_reference.at(remove->getreference1()).erase(it);
			break;
		}		
	}
	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqre1 =
	als_on_reference.at(remove->getreference1()).equal_range(remove->getend1());
	for(multimap<size_t, pw_alignment*>::iterator it = eqre1.first; it!=eqre1.second; ++it) {
		if ( &(it->second) == &remove ) {
			als_on_reference.at(remove->getreference1()).erase(it);
			break;
		}		
	}

	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqrb2 =
	als_on_reference.at(remove->getreference2()).equal_range(remove->getbegin2());
	for(multimap<size_t, pw_alignment*>::iterator it = eqrb2.first; it!=eqrb2.second; ++it) {
		if ( &(it->second) == &remove ) {
			als_on_reference.at(remove->getreference1()).erase(it);
			break;
		}		
	}

	pair< multimap<size_t, pw_alignment*>::iterator, multimap<size_t, pw_alignment*>::iterator > eqre2 =
	als_on_reference.at(remove->getreference2()).equal_range(remove->getend2());
	for(multimap<size_t, pw_alignment*>::iterator it = eqre2.first; it!=eqre2.second; ++it) {
		if ( &(it->second) == &remove ) {
			als_on_reference.at(remove->getreference2()).erase(it);
			break;
		}		
	}
	alignments.erase(remove);
}

	void overlap::test_all() const {

	for(set<pw_alignment*, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		pw_alignment * al = *it;
		bool alok = data.alignment_fits_ref(al);
		if(!alok) {
			exit(1);
		}
		
	
	}


}









	
	void overlap::insert_without_partial_overlap(const pw_alignment & p){
	pw_alignment * np = new pw_alignment(p);
	alignments.insert(np);

 	
	multimap<size_t , pw_alignment *> & alignment_on_reference1 = als_on_reference.at(np->getreference1());
	multimap<size_t , pw_alignment *> & alignment_on_reference2 = als_on_reference.at(np->getreference2());

	pair<size_t,  pw_alignment *> begin1(np->getbegin1(),np);
		alignment_on_reference1.insert(begin1);
	pair<size_t,  pw_alignment *> begin2(np->getbegin2(),np);
		alignment_on_reference2.insert(begin2);
	pair<size_t,  pw_alignment *> end1(np->getend1(),np);
		alignment_on_reference1.insert(end1);
	pair<size_t, pw_alignment *> end2(np->getend2(),np);
		alignment_on_reference2.insert(end2);
}
