#include "data.hpp"

// Modification of the file
// Reading of Sam files

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


all_data::all_data(string fasta_all_sequences, string samFile) {
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

	// READ SAM FILE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	bool verbose = true;
	  if(verbose) cout << "readAlignment sam file"<<endl;

	  SamFile samIn;
	  samIn.OpenForRead(samFile.c_str());

	  // Read the sam header.
	  SamFileHeader samHeader;
	  samIn.ReadHeader(samHeader);

	  SamRecord samRecord;
	  Cigar* tmpCigar;
	  string alRefSeq, alReadSeq;
	  int skip_alt = 0;
	  int skip_self = 0;
	  //GenomeSequence myRefSeq("/ebio/abt6_projects7/small_projects/mdubarry/Documents/SampleProgram/bin/output/tmpOut.fasta");
	  GenomeSequence myRefSeq(fasta_all_sequences.c_str());
	  while(samIn.ReadRecord(samHeader, samRecord))
	  {
		// Achtung : sam file considere the file reference as an uniq sequence even if there are differentes sequences in reference, so the position are concatenate.
		//e.g. if I take the first base of the second sequence from my reference.fasta, it will not return 1 but lengthOfTheSequence1 + 1 !

		// For each Record do :
		tmpCigar = samRecord.getCigarInfo(); //Pointer to Cigar object
		const char* currentSeq = samRecord.getSequence();
		int myStartOfReadOnRefIndex = myRefSeq.getGenomePosition(samRecord.getReferenceName(),samRecord.get1BasedPosition()); // Select the start on the good sequence in the fasta File
		if(myStartOfReadOnRefIndex == -1){
			std::cerr << "Error with the reference "<< samRecord.getReferenceName() << std::endl;
			exit(1);
		}
		if (verbose) cout<<myStartOfReadOnRefIndex<<" getReferenceName "<<samRecord.getReferenceName()<<" get1BasedPosition "<< samRecord.get1BasedPosition() <<" myRefSeq.sequenceLength " <<myRefSeq.sequenceLength()
		  << " getAlignmentLength "<< samRecord.getAlignmentLength() << " getReadLength " << samRecord.getReadLength()<< " getNumBeginClips "<< tmpCigar->getNumBeginClips()<< endl;
		if(verbose) std::cout <<"flag "<<samRecord.getFlag() <<" reverse ? " << SamFlag::isReverse(samRecord.getFlag()) << std::endl;

		// Loop through the read and determine the difference with the reference
		int prevRef = -1;
		for(int index = 0; index < samRecord.getReadLength(); index++)
		{
			int refOffset = tmpCigar->getRefOffset(index);
			if(refOffset == Cigar::INDEX_NA)
			{
				// No reference base, meaning it is not in the reference, so add missing
				alRefSeq += '-';
				alReadSeq += currentSeq[index];
			}
			else
			{
				// While the reference offset is not 1 more than the previous, it means
				// there is a deletion/N, so add missing to the read, and add the reference bases.
				while(refOffset != ++prevRef)
				{
					alReadSeq += '-';
					alRefSeq += myRefSeq[myStartOfReadOnRefIndex + prevRef];

				}
				//now we are at a spot in both the reference and the read, so add it.
				alReadSeq += currentSeq[index];
				alRefSeq += myRefSeq[myStartOfReadOnRefIndex + refOffset];
			}
		}
		if(verbose)  std::cout << "\n" << alRefSeq <<"\n" << alReadSeq  <<std::endl;
		int myEndReadOnRefIndex =  myRefSeq.getGenomePosition(samRecord.getReferenceName(),samRecord.get1BasedAlignmentEnd()) + 1 ;
		if(verbose) cout << "reference start "<< myStartOfReadOnRefIndex << " end " << myEndReadOnRefIndex <<endl;

		// Reference
		string acc1;
		string name1;
		name_split(samRecord.getReferenceName(), acc1, name1);
		size_t start1 = samRecord.get1BasedPosition() - 1 ;
		size_t incl_end1 = start1 + tmpCigar->getExpectedQueryBaseCount() - 1 ; // change to sequence partial length
		map<string, size_t>::iterator findseq1 = longname2seqidx.find(samRecord.getReferenceName());
		if(findseq1==longname2seqidx.end()) {
			cerr << "Error: unknown sequence: " << samRecord.getReferenceName() << endl;
			exit(1);
		}
		size_t idx1 = findseq1->second;

		// Read
		string acc2;
		string name2;
		name_split(samRecord.getReadName(), acc2, name2);
		size_t start2;
		if(tmpCigar->getNumBeginClips()==0)
			start2 = 0;
		else
			start2 = tmpCigar->getNumBeginClips()+1; // getNumbeginClips give the number of clip and the read start at the next one ( so +1)

		size_t incl_end2 = start2 + samRecord.getReadLength() -1;//TODO -1 ?

		// In Sam file : Reverse and Forward ??
		string tmpString = samRecord.getReadName();
		map<string, size_t>::iterator findseq2 = longname2seqidx.find((tmpString));
		if(findseq2==longname2seqidx.end()) {
			cerr << "Error: unknown sequence: " << samRecord.getReadName() << endl;
			exit(1);
		}
		size_t idx2 = findseq2->second;


    
	/*ifstream mafin(maf_all_alignments.c_str());

	
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
					if(0!=str.compare("")) {
						cerr << "Error: exspected empty line after alignment. Seen: " << str << endl;
						exit(1);
					}
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
					}


					parts.clear();
					strsep(aline2, " ", parts);
					if(parts.size()!=7) {
						cerr << "Error: exspected 7 fields in sequence line: " << aline2 << endl;
						exit(1);
					}

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
					}
	*/
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
			if(al.getreference(0)==idx1)
				same_ref = 0;
			else if(al.getreference(1)==idx1) 
				same_ref = 1;
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
	
		if(skip)
			skip_alt++;
		else {
			pw_alignment al(alRefSeq, alReadSeq, start1, start2, incl_end1, incl_end2, idx1, idx2);
			size_t alidx = alignments.size();
			alignments.push_back(al);
			als_on_reference.at(idx1).insert(make_pair(start1, alidx));
			als_on_reference.at(idx2).insert(make_pair(start2, alidx));
		}

	}
	else 
		skip_self++;
		// cerr << "Warning: Skip self alignment in seq " << idx1 << " from " << start1 << " to " << incl_end1 << endl;
// Clean sequences before reloop
    alReadSeq = "";
    alRefSeq = "";
}

	if(dnastring::found_iupac_ambiguity) {
		cerr << "Warning: IUPAC DNA ambiguity characters were replaced by N" << endl;
	}

	cout << "Loaded: " << sequences.size() << " sequences and " << alignments.size() << " pairwise alignments " << endl;
	cout << skip_self << " self alignments were skipped" << endl;
	cout << skip_alt << " alternative alignments of identical regions were skipped" << endl;
}

all_data::~all_data() {

}


	const dnastring & all_data::getSequence(size_t index) const {
		return sequences.at(index);
	}


	const pw_alignment & all_data::getAlignment(size_t index) const {
		return alignments.at(index);
	}
	const vector<pw_alignment>& all_data::getAlignments()const{
		return alignments;
	}
	const vector<size_t> & all_data::getAcc(size_t acc)const{
		string accName;
		for(map<string,size_t>::const_iterator it = accession_name.begin();it != accession_name.end(); it++){
			if(it->second == acc){
				accName = it->first;
			}else continue;
		}
		return acc_sequences.at(accName);		
	}
	size_t all_data::accNumber(size_t sequence_id) const {
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
	const size_t all_data::numAcc() const {
		return  acc_sequences.size();
	}
	size_t all_data::numOfAcc()const{
		return accession_name.size();
	}
	void all_data::set_accession(const string & acc){
		accession_name.insert(make_pair(acc, accession_name.size()));
	//	cout<<"accession name size: " << accession_name.size() <<endl;
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
				rchar =  dnastring::complement(rchar);

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
				rchar =  dnastring::complement(rchar);

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


const string all_data::get_acc(size_t acc)const{	
	string Accession;
	for(map<string, size_t>::const_iterator it = accession_name.begin() ; it != accession_name.end();it++){
		if(it->second == acc){
			Accession = it -> first;
		}else continue;
	}
	return Accession;
}

const size_t all_data::get_acc_id(string acc)const{
	map<string, size_t>::const_iterator it = accession_name.find(acc);
	return it->second;
}


string all_data::get_seq_name(size_t s) const {
	return sequence_names.at(s);
}

size_t all_data::get_seq_size(size_t s) const {
	return sequences.at(s).length();
}
/*
overlap::overlap(const all_data & d): data(d), als_on_reference(d.numSequences()){

	}

overlap::overlap(const overlap & o): data(o.data), als_on_reference(o.data.numSequences()){}


	
	overlap::~overlap(){
		for(set<pw_alignment*, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
			delete *it;
		}
	}

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
		delete remove;
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
		cout << endl;
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



	
pw_alignment * overlap::insert_without_partial_overlap(const pw_alignment & p){
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

	return np;
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
	
	set < const pw_alignment*, compare_pw_alignment> alignment_parts;
		
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
		}
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
		
const set<pw_alignment*, compare_pw_alignment> & overlap::get_all() const {
	return alignments;
}


// true if true partial overlap
bool overlap::check_po(size_t l1, size_t r1, size_t l2, size_t r2) {
	if(r2 >= l1 && l2 <= r1) {
		if(l1!=l2 || r1!=r2) {
			return true;	
		}
	}
	return false;
}

void overlap::test_partial_overlap_set(set< const pw_alignment *, compare_pw_alignment> & als) {
	
	vector<const pw_alignment *> all;
	for(set<const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
		const pw_alignment * al = *it;
		all.push_back(al);
	}
	overlap::test_partial_overlap_vec(all);
}


void overlap::test_partial_overlap_vec(vector< const pw_alignment *> & all) {
	for(size_t i=0; i<all.size(); ++i) {
		for(size_t j=i+1; j<all.size(); ++j) {
			const pw_alignment * a = all.at(i);
			const pw_alignment * b = all.at(j);
			size_t al1, ar1, al2, ar2;
			size_t bl1, br1, bl2, br2;
			a->get_lr1(al1, ar1);
			a->get_lr2(al2, ar2);
			b->get_lr1(bl1, br1);
			b->get_lr2(bl2, br2);
			size_t aref1, aref2, bref1, bref2;
			aref1 = a->getreference1();
			aref2 = a->getreference2();
			bref1 = b->getreference1();
			bref2 = b->getreference2();
			if(aref1 == bref1) {
				if(check_po(al1, ar1, bl1, br1)) {
					cout << "partial overlap error 1: " << endl;
					a->print();
					b->print();
					exit(1);
				}
			}
			
			if(aref1 == bref2) {
				if(check_po(al1, ar1, bl2, br2)) {
					cout << "partial overlap error 2: " << endl;
					a->print();
					b->print();
					exit(1);
				}
			}

			if(aref2 == bref1) {
				if(check_po(al2, ar2, bl1, br1)) {
					cout << "partial overlap error 3: " << endl;
					a->print();
					b->print();
					exit(1);
				}
			}

			if(aref2 == bref2) {
				if(check_po(al2, ar2, bl2, br2)) {
					cout << "partial overlap error 4: " << endl;
					a->print();
					b->print();
					exit(1);
				}
			}

		
		}
	}




}
void overlap::test_partial_overlap() const {

	vector<const pw_alignment *> all;
	for(set<pw_alignment *, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment * al = *it;
		all.push_back(al);
	}

	overlap::test_partial_overlap_vec(all);

}

#define SPLITPRINT 1

splitpoints::splitpoints(const pw_alignment & p, const overlap & o, const all_data & d):overl(o), newal(p),data(d), split_points(d.numSequences()) {

}

splitpoints::~splitpoints() {}

void splitpoints::find_initial_split_points_nonrecursive(size_t sequence, size_t left, size_t right) {
		const multimap<size_t , pw_alignment *> & alignments_on_reference = overl.get_als_on_reference_const(sequence);

#if SPLITPRINT		
		cout << " seach for initial split points on " << sequence << " from " << left << endl;
		size_t count = 0;
#endif

		for( multimap<size_t, pw_alignment *>::const_iterator it=alignments_on_reference.lower_bound(left);it!=alignments_on_reference.end(); ++it){
			const pw_alignment * al = it->second;
			
#if SPLITPRINT
			cout << count++ <<" See " << endl;
			al->print();
			cout << endl;
#endif

			// loop break condition. This special treatment is needed to avoid to-early break in case of both parts of al being on the same reference
			size_t al_leftmost_leftbound = (size_t)-1;

			if(al->getreference1() == sequence) {
				size_t alleft;
				size_t alright;
				al->get_lr1(alleft, alright);
				if(alleft < al_leftmost_leftbound) {
					al_leftmost_leftbound = alleft;
				}



				if(alleft < left && left <= alright) {
	//				cout<<"here1"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
					insert_split_point_nonrecursive(sequence, left);
				}
				if(alleft <= right && right < alright) {
	//				cout<<"here2"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
					insert_split_point_nonrecursive(sequence, right+1);

				}
				if(left < alleft && alleft < right) {
	//				cout<<"here3"<<endl;
	//				cout<<"al: "<<endl;
	*/
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						cout <<col <<"\t"<< c1<<"\t"<<c2<<endl;
					}*/

				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
	/*
					insert_split_point_nonrecursive(sequence, alleft);

				}
				if(left < alright && alright < right) {
		//			cout<<"here4"<<endl;
				//	cout<<"al: "<<endl;
		//			al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point_nonrecursive(sequence, alright+1);
					
				}

			}
		
			if(al->getreference2() == sequence) {
				size_t alleft;
				size_t alright;
				al->get_lr2(alleft, alright);
				if(alleft < al_leftmost_leftbound) {
					al_leftmost_leftbound = alleft;
				}



				if(alleft < left && left <= alright) {
				//	cout<<"here5"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();

					insert_split_point_nonrecursive(sequence, left);
				
				}
				if(alleft <= right && right < alright) {
		//			cout<<"here6"<<endl;
				//	al->print();
					insert_split_point_nonrecursive(sequence, right+1);
				}
				if(left < alleft && alleft < right) {
		//			cout<<"here7"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point_nonrecursive(sequence, alleft);
				}
				if(left < alright && alright < right) {
		//			cout<<"here8"<<endl;
				//	cout<<"al: "<<endl;
		//			al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point_nonrecursive(sequence, alright+1);
				}

			
			}




			if(al_leftmost_leftbound > right)  {
#if SPLITPRINT
				cout << "Break after See " << count << " al rightmost leftbound is " << al_leftmost_leftbound<<  endl;	
#endif
				break;
			}
		}


}

void splitpoints::find_initial_split_points(size_t sequence, size_t left, size_t right) {
		const multimap<size_t , pw_alignment *> & alignments_on_reference = overl.get_als_on_reference_const(sequence);

#if SPLITPRINT		
		cout << " seach for initial split points on " << sequence << " from " << left << endl;
		size_t count = 0;
#endif

		for( multimap<size_t, pw_alignment *>::const_iterator it=alignments_on_reference.lower_bound(left);it!=alignments_on_reference.end(); ++it){
			const pw_alignment * al = it->second;
			
#if SPLITPRINT
			cout << count++ <<" See " << endl;
			al->print();
			cout << endl;
#endif

			// loop break condition. This special treatment is needed to avoid to-early break in case of both parts of al being on the same reference
			size_t al_leftmost_leftbound = (size_t)-1;

			if(al->getreference1() == sequence) {
				size_t alleft;
				size_t alright;
				al->get_lr1(alleft, alright);
				if(alleft < al_leftmost_leftbound) {
					al_leftmost_leftbound = alleft;
				}



				if(alleft < left && left <= alright) {
	//				cout<<"here1"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
					insert_split_point(sequence, left);
				}
				if(alleft <= right && right < alright) {
	//				cout<<"here2"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
					insert_split_point(sequence, right+1);

				}
				if(left < alleft && alleft < right) {
	//				cout<<"here3"<<endl;
	//				cout<<"al: "<<endl;
	*/

				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						cout <<col <<"\t"<< c1<<"\t"<<c2<<endl;
					}*/

				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
	/*
					insert_split_point(sequence, alleft);

				}
				if(left < alright && alright < right) {
		//			cout<<"here4"<<endl;
				//	cout<<"al: "<<endl;
		//			al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point(sequence, alright+1);
					
				}

			}
		
			if(al->getreference2() == sequence) {
				size_t alleft;
				size_t alright;
				al->get_lr2(alleft, alright);
				if(alleft < al_leftmost_leftbound) {
					al_leftmost_leftbound = alleft;
				}



				if(alleft < left && left <= alright) {
				//	cout<<"here5"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();

					insert_split_point(sequence, left);
				
				}
				if(alleft <= right && right < alright) {
		//			cout<<"here6"<<endl;
				//	al->print();
					insert_split_point(sequence, right+1);
				}
				if(left < alleft && alleft < right) {
		//			cout<<"here7"<<endl;
				//	cout<<"al: "<<endl;
				//	al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point(sequence, alleft);
				}
				if(left < alright && alright < right) {
		//			cout<<"here8"<<endl;
				//	cout<<"al: "<<endl;
		//			al->print();
				//	cout<<"newal: "<<endl;
				//	newal.print();
					insert_split_point(sequence, alright+1);
				}

			
			}




			if(al_leftmost_leftbound > right)  {
#if SPLITPRINT
				cout << "Break after See " << count << " al rightmost leftbound is " << al_leftmost_leftbound<<  endl;	
#endif
				break;
			}
		}

	}

	void splitpoints::recursive_splits(){




		size_t left1;
		size_t right1;
		size_t left2 ;
		size_t right2;
		newal.get_lr1(left1, right1);
		newal.get_lr2(left2, right2);

#if SPLITPRINT
		cout << " Check ref 1 overlaps " << endl;
#endif

		find_initial_split_points(newal.getreference1(), left1, right1);

#if SPLITPRINT
		cout << " Check ref 2 overlaps " << endl;
#endif
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

		insert_split_point(newal.getreference1(), left1);
		insert_split_point(newal.getreference1(), right1+1);
		insert_split_point(newal.getreference2(), left2);
		insert_split_point(newal.getreference2(), right2+1);
	
}

void splitpoints::nonrecursive_splits(){


		size_t left1;
		size_t right1;
		size_t left2 ;
		size_t right2;
		newal.get_lr1(left1, right1);
		newal.get_lr2(left2, right2);

#if SPLITPRINT
		cout << " Check ref 1 overlaps " << endl;
#endif

		find_initial_split_points_nonrecursive(newal.getreference1(), left1, right1);

#if SPLITPRINT
		cout << " Check ref 2 overlaps " << endl;
#endif
		find_initial_split_points_nonrecursive(newal.getreference2(), left2, right2);
	
		if(newal.getreference1()==newal.getreference2()) {
			if(right1>=left2 && left1 < left2){
				insert_split_point_nonrecursive(newal.getreference2(),left2);
				insert_split_point_nonrecursive(newal.getreference2(),right1+1);
			}
			if(right2>=left1 && left2 < left1){
				insert_split_point_nonrecursive(newal.getreference2(),left1);
				insert_split_point_nonrecursive(newal.getreference2(),right2+1);
			}
		}

		insert_split_point_nonrecursive(newal.getreference1(), left1);
		insert_split_point_nonrecursive(newal.getreference1(), right1+1);
		insert_split_point_nonrecursive(newal.getreference2(), left2);
		insert_split_point_nonrecursive(newal.getreference2(), right2+1);
	


}


void splitpoints::insert_split_point_nonrecursive(size_t sequence, size_t position) {
	pair<set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
	if(insertes.second) {
#if SPLITPRINT
		cout << " new split point " << sequence << " at " << position << endl;
#endif	
	}
}




void splitpoints::insert_split_point(size_t sequence, size_t position) {
	pair<set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
	if(insertes.second) {
#if SPLITPRINT
		cout << " new split point " << sequence << " at " << position << endl;
#endif	
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
				newal.split(true, position, fp, sp);

				pw_alignment fpe;
				pw_alignment spe;


				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(&fp);
				bool sgaps = onlyGapSample(&sp);
					
				// find split part alignment ends on reference after removing end gaps
				size_t fpeleft2= (size_t) -1;
				size_t fperight2 = (size_t) -1;
				size_t speleft2 = (size_t) -1;
				size_t speright2 = (size_t) -1;
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
					fpe.get_lr2(fpeleft2, fperight2);
#if SPLITPRINT
					cout << "newal fpe " << endl;
					fpe.print();
					cout  << endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
					spe.get_lr2(speleft2, speright2);
#if SPLITPRINT
					cout << "newal spe " << endl;
					spe.print();
					cout << endl;
#endif
				}


				if(newal.getbegin2() < newal.getend2()) {
			//			cout << " try ins " << newal.getreference2() << " : " << spleft2 << endl;
					if(!sgaps) {
						insert_split_point(newal.getreference2(), speleft2);
					}
					// additional split point if gaps at split point
					if(!fgaps) {
						insert_split_point(newal.getreference2(), fperight2+1);
					}
				} else {
			//			cout << " try ins " << newal.getreference2() << " : " << fpleft2 << endl;
					if(!fgaps) {
						insert_split_point(newal.getreference2(), fpeleft2);
					}
					if(!sgaps) {
						insert_split_point(newal.getreference2(), speright2+1);
					}

				}
			}
		}
		if(newal.getreference2()==sequence) {
			if(left2<position && right2>=position) {
			//	cout<<"HERE!!"<<endl;
				pw_alignment fp;
				pw_alignment sp;
				newal.split(false, position, fp, sp);
				pw_alignment fpe;
				pw_alignment spe;
				
				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(&fp);
				bool sgaps = onlyGapSample(&sp);
					
				// find split part alignment ends on reference after removing end gaps
				size_t fpeleft1 = (size_t) -1;
				size_t fperight1 = (size_t) -1;
				size_t speleft1 = (size_t) -1;
				size_t speright1 = (size_t) -1;
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
					fpe.get_lr1(fpeleft1, fperight1);
#if SPLITPRINT
					cout << "newal fpe " << endl;
					fpe.print();
					cout  << endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
					spe.get_lr1(speleft1, speright1);
#if SPLITPRINT
					cout << "newal spe " << endl;
					spe.print();
					cout << endl;
#endif
				}

				
				if(newal.getbegin1() < newal.getend1()) {	
		//			cout << " try ins " << newal.getreference1() << " : " << spleft1 << endl;
					if(!sgaps) {
						insert_split_point(newal.getreference1(), speleft1);
					}
					if(!fgaps) {
						insert_split_point(newal.getreference1(), fperight1+1);
					}
				} else {
			//			cout << " try ins " << newal.getreference1() << " : " << fpleft1 << endl;
					if(!fgaps) {
						insert_split_point(newal.getreference1(), fpeleft1);
					}
					if(!sgaps) {
						insert_split_point(newal.getreference1(), speright1+1);
					}
				}
			}
		}





		const multimap<size_t, pw_alignment *> & alonref = overl.get_als_on_reference_const(sequence);
		for(multimap<size_t, pw_alignment *>::const_iterator it = alonref.lower_bound(position); it!=alonref.end(); ++it) {
			const pw_alignment * al = it-> second;
				
				
	//			cout << " in ins " << sequence << " : " << position << " see: " << endl;
	//			al->print();
	//			cout << endl;


			if(al->getreference1()!=al->getreference2()) {
				size_t alleft;
				size_t alright;
				al->get_lr_on_reference(sequence, alleft, alright);
				if(alright < position) {
			//			cout << " break " << endl;
					break;
				}
				if(alleft == position) {
			//			cout << " break " << endl;
					break;
				} 
				if(alright+1== position) {
			//			cout << " break " << endl;
					break;
				}
				if(alleft>position){
			//			cout << " break " << endl;
					break;
				}
			//		al->print(); 
				//	cout<<"position"<<position<<endl;
				//	cout<<"alleft"<<alleft<<endl;
				//	cout<<"alright"<<alright<<endl;
				assert(alleft < position && position <= alright);
				pw_alignment fp;
				pw_alignment sp;
				pw_alignment fpe;
				pw_alignment spe;
	*/
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						cout <<col <<"\t"<< c1<<"\t"<<c2<<endl;
					}*/
				//	cout<<"al length: "<< al->alignment_length() <<endl;
				//	if(al->alignment_length()>1){
/*
				al->split_on_reference(sequence, position, fp, sp);
				


				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(&fp);
				bool sgaps = onlyGapSample(&sp);
					
				// find split part alignment ends on reference after removing end gaps
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
#if SPLITPRINT
					cout << "al fpe " << endl;
					fpe.print();
					cout  << endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
#if SPLITPRINT
					cout << "al spe " << endl;
					spe.print();
					cout << endl;
#endif
				}




					
				if(sp.getreference1()==sequence) {


					
						if(al->getbegin2() < al->getend2()) {
							if(!sgaps) {
								insert_split_point(sp.getreference2(), spe.getbegin2());
							} 
							if(!sgaps) {
								insert_split_point(sp.getreference2(), fpe.getend2()+1);
							}


						} else {
							if(!fgaps) {
								insert_split_point(sp.getreference2(), fpe.getend2());
							}
							if(!sgaps) {
								insert_split_point(sp.getreference2(), spe.getbegin2()+1);
							}
						}
					} else {




						//	cout<<"Heya!!!"<<endl;
						if(al->getbegin2() < al->getend2()) {
							if(!sgaps) {
								insert_split_point(sp.getreference1(), spe.getbegin1());
							}
							if(!fgaps) {
								insert_split_point(sp.getreference1(), fpe.getend1()+1);
							}
						} else {
							if(!fgaps) {
								insert_split_point(sp.getreference1(), fpe.getend1());
							}
							if(!sgaps) {
								insert_split_point(sp.getreference1(), spe.getbegin1()+1);
							}
						}
					}
				} else {
					
					size_t alleft1;
					size_t alright1;
					al->get_lr1(alleft1, alright1);
					size_t alleft2;
					size_t alright2;
					al->get_lr2(alleft2, alright2);
					if(position > alleft1 && position <= alright1) {
		//				cout << "inone" << endl;
						pw_alignment fp;
						pw_alignment sp;
						al->split(true, position, fp, sp);
					pw_alignment fpe;
					pw_alignment spe;
					// do splitted parts have only gaps one of their reference sequences
					bool fgaps = onlyGapSample(&fp);
					bool sgaps = onlyGapSample(&sp);
					
					// find split part alignment ends on reference after removing end gaps
					size_t fpeleft2 = (size_t) -1;
					size_t fperight2 = (size_t) -1;
					size_t speleft2 = (size_t) -1;
					size_t speright2 = (size_t) -1;
					if(!fgaps) {
						fp.remove_end_gaps(fpe);
						fpe.get_lr2(fpeleft2, fperight2);
#if SPLITPRINT
						cout << "al fpe " << endl;
						fpe.print();
						cout  << endl;
#endif
					}
					if(!sgaps) {
						sp.remove_end_gaps(spe);
						spe.get_lr2(speleft2, speright2);
#if SPLITPRINT
						cout << "al spe " << endl;
						spe.print();
						cout << endl;
#endif
					}


						
							
						if(al->getbegin2() < al->getend2()) {
			//				cout << " spleft2 " << spleft2 << endl;
							if(!sgaps) {
								insert_split_point(sp.getreference2(), speleft2);
							}
							if(!fgaps) {
								insert_split_point(sp.getreference2(), fperight2+1);
							}
						} else {
			//				cout << " fpgetend2 " << fp.getend2() << endl;
							if(!fgaps) {
								insert_split_point(fp.getreference2(), fpe.getend2());
							}
							if(!sgaps) {
								insert_split_point(fp.getreference2(), spe.getbegin2()+1);
							}

						}
					}
											
					if(position > alleft2 && position <= alright2) {
			//			cout << "intwo" << endl;
						pw_alignment fp;
						pw_alignment sp;
						al->split(false, position, fp, sp);
					pw_alignment fpe;
					pw_alignment spe;
					
			// do splitted parts have only gaps one of their reference sequences
					bool fgaps = onlyGapSample(&fp);
					bool sgaps = onlyGapSample(&sp);
					
					// find split part alignment ends on reference after removing end gaps
					size_t fpeleft1 = (size_t) -1;
					size_t fperight1 = (size_t) -1;
					size_t speleft1 = (size_t) -1;
					size_t speright1 = (size_t) -1;
					if(!fgaps) {
						fp.remove_end_gaps(fpe);
						fpe.get_lr1(fpeleft1, fperight1);
#if SPLITPRINT
						cout << "al fpe " << endl;
						fpe.print();
						cout  << endl;
#endif
					}
					if(!sgaps) {
						sp.remove_end_gaps(spe);
						spe.get_lr1(speleft1, speright1);
#if SPLITPRINT
						cout << "al spe " << endl;
						spe.print();
						cout << endl;
#endif
					}




//							data.alignment_fits_ref(&fp);
//							data.alignment_fits_ref(&sp);
						
						if(al->getbegin1() < al->getend1()) {
							if(!sgaps) {
								insert_split_point(sp.getreference1(), speleft1);
							}
							if(!fgaps){
								insert_split_point(sp.getreference1(), fperight1+1);
							}
						} else {
							if(!fgaps) {
								insert_split_point(sp.getreference1(), fpe.getend1());
							}
							if(!sgaps) {
								insert_split_point(sp.getreference1(), spe.getbegin1()+1);
							}
						}
					//	}
					}
				}
			}
			
		}
		

	}


	void splitpoints::split_all(set<const pw_alignment*,compare_pw_alignment> & remove_alignments, vector<pw_alignment> & insert_alignments ){

#if SPLITPRINT		
		cout << " All split points " << endl;
		for(size_t i=0; i<data.numSequences(); ++i) {
			cout << "ref " << i << " ";
			for(set<size_t>::iterator it = split_points.at(i).begin(); it!=split_points.at(i).end(); ++it) {
				cout << " " << *it;
			}
			cout << endl;
		}
#endif

		pw_alignment p1;
		pw_alignment p2;
		for(size_t i = 0; i <data.numSequences();i++){

#if SPLITPRINT			
			cout << "REFERENCE " << i << endl;
#endif
			const multimap<size_t, pw_alignment *> & als_on_ref = overl.get_als_on_reference_const(i);
			for(set<size_t>::iterator split = split_points.at(i).begin(); split!= split_points.at(i).end(); ++split){
#if SPLITPRINT
				cout << "SPLIT " << *split << endl;
#endif
				for(multimap<size_t, pw_alignment*>::const_iterator it =als_on_ref.lower_bound(*split); it!=als_on_ref.end(); ++it){
					pw_alignment * p = it-> second;


//					p->print();
//					cout << endl;

					if(p->getreference1() == i) {
						size_t pleft1;
						size_t pright1;
						p->get_lr1(pleft1,pright1);
						if( pleft1<*split && pright1>=*split){
							remove_alignments.insert(p);
#if SPLITPRINT
							cout<<"remove alignment1: "<<endl;
							p->print();
#endif

						}
					}
					if(p->getreference2()==i) {
						size_t pleft2;
						size_t pright2;
						p->get_lr2(pleft2,pright2);
						if( pleft2<*split && pright2>=*split){
							//	pw_alignment *np = new pw_alignment(*p);
							remove_alignments.insert(p);
#if SPLITPRINT
							cout<<"remove alignment2: "<<endl;
									p->print();
#endif
						}
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
#ifndef NDEBUG
				if(!data.alignment_fits_ref(&split_pieces.at(i))) {
				//	cout<<"fails here!"<<endl;
					exit(1);
				}
#endif
			if(!onlyGapSample(&split_pieces.at(i))){	
				pw_alignment noendgaps;
				split_pieces.at(i).remove_end_gaps(noendgaps);
#if SPLITPRINT
				cout << "ENDGAPS" << endl;
				split_pieces.at(i).print();
				cout << "TO " << endl;
				noendgaps.print();
				cout << endl;
#endif
				set<pw_alignment*,compare_pw_alignment>::const_iterator it = inserted_pieces.find(&noendgaps);
				if (overl.checkAlignments(&noendgaps)) {}
				else if ( it != inserted_pieces.end()){}
				else {
					insert_alignments.push_back(noendgaps);
				
					pw_alignment * ipp = new pw_alignment(noendgaps); //&(insert_alignments.at(insert_alignments.size()-1));

					inserted_pieces.insert(ipp);
				}
			}
		}
		
		for(set<pw_alignment*,compare_pw_alignment>::iterator it = inserted_pieces.begin(); it!=inserted_pieces.end(); ++it) {
			delete *it;
		}
			
	}
	void splitpoints::splits(const pw_alignment * p,  vector<pw_alignment> & insert_alignments){
#if SPLITPRINT	
		cout << "SPL" << endl;
		p->print();
		cout << endl;

		cout << " split on " << p->getreference1() << endl;
#endif

		pw_alignment p1;
		pw_alignment p2;
		size_t left1;
		size_t right1;	
//		size_t left2;
//		size_t right2;
		p->get_lr1(left1,right1);
//		p->get_lr2(left2,right2);
		for(set<size_t>::iterator splitp = split_points.at(p->getreference1()).upper_bound(left1); splitp!= split_points.at(p->getreference1()).end(); splitp++){
			if(right1>=*splitp){
#if SPLITPRINT
				cout << "sp " << *splitp << endl;
#endif
				p->split(true,*splitp,p1,p2);
				if(p->getbegin1() < p->getend1()) {
					p = &p2;

#if SPLITPRINT
					p1.print();
#endif
					if(!onlyGapSample(&p1) && p1.alignment_length() >1 &&p1.getbegin1()!=p1.getend1() && p1.getbegin2()!=p1.getend2()) {
						insert_alignments.push_back(p1);
					}
				} else if(p->getbegin1() > p->getend1()) {
					p = &p1;

#if SPLITPRINT
					p2.print();
#endif

					if(!onlyGapSample(&p2) && p2.alignment_length()>1 && p2.getbegin1()!=p2.getend1() && p2.getbegin2()!=p2.getend2() ){	
						insert_alignments.push_back(p2);
					}
					
				}
			}
			else break;
		}

#if SPLITPRINT
		cout << " last part" << endl;
		p->print();
		cout << endl;
#endif

		if(!onlyGapSample(p)&& p->alignment_length()>1 && p->getbegin1()!=p->getend1() && p->getbegin2()!=p->getend2() ){	
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
			if(! gapOnSample1 && !gapOnSample2) {
				return false;
			}		
		}
		return gapOnSample1 || gapOnSample2;
	}


	vector<pw_alignment>  splitpoints::get_insert()const{
		return insert_alignments;
	}
		*/
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

					//	cout << "k " << k << " num " << modification.at(acc1).at(acc2).at(j).at(k) << endl;
						sum+=modification.at(acc1).at(acc2).at(j).at(k);
					}	
					for(size_t k=0; k<6; ++k) {
						modification.at(acc1).at(acc2).at(j).at(k)/=sum;
					//	cout << " simple model cost: acc " << acc1 << " to " << acc2 << " modify " << j << " to " << k << " costs " << -log2(modification.at(acc1).at(acc2).at(j).at(k)) << " bits " << endl; 
					
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
		cout << "gain_function "<<endl;
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
		//cout << "in gain: create 2 " << c2 << " m1 " << m1 << endl; 

		g1 = c2 - m1;
		g2 = c1 - m2;

	}

	void model::train() {
		acc_base_frequency();
		alignment_modification();
	}

mc_model::mc_model(all_data & d):data(d), sequence_successive_bases(d.numAcc()), create_cost(d.numAcc(),vector<double>(5,1)),mod_cost(d.numAcc(),vector<map<string, vector<double> > >(d.numAcc())),high(d.numAcc()),highValue(d.numAcc(),vector<map<string, vector<unsigned int> > >(d.numAcc())){
	cout << "mc_model "<<endl;
	size_t numberOfPowers = 32;
//	size_t numberOfPowers = NUM_DELETE;
	if(NUM_KEEP>numberOfPowers) {
		numberOfPowers = NUM_KEEP;
	}
	powersOfTwo = vector<size_t>(numberOfPowers, 1);
	for(size_t i=1; i< numberOfPowers; ++i) {
		powersOfTwo.at(i) = powersOfTwo.at(i-1)*2;
		
	}
}


mc_model::~mc_model(){}

/*
   train sequence model

*/
void mc_model::markov_chain(){
		for(size_t k = 0; k < data.numSequences(); k++){
			size_t acc = data.accNumber(k);
			for(size_t i = 0 ; i< data.getSequence(k).length(); i++ ){
				char schr = data.getSequence(k).at(i);
				size_t s = dnastring::base_to_index(schr);
				stringstream context;		
				for(size_t j = Sequence_level; j>0; j--){
					char rchr;
					signed long temp = i - j;
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
			for(map <string, vector<double> >::iterator it= sequence_successive_bases.at(i).begin();it!=sequence_successive_bases.at(i).end();it++){
				string seq = it->first;
				int total = 0;
				vector<double> & base = sequence_successive_bases.at(i).at(seq);
				for(size_t j = 0; j<5;j++){
					total += base.at(j);
			//	cout<<"base: "<<base.at(j)<<endl;
				}
		//	cout<<"totalnumber of bases for each seq:  "<< total<<endl;
				for(size_t k=0; k<5;k++){
					base.at(k) = -log2(base.at(k)/total);
					create_cost.at(i).at(k) += base.at(k);
				//	cout<<"cr_cost: "<<base.at(k) <<endl;		
				}
			}
		}
			


	}



/*
	train alignment/modification model
 TODO does this write the alignment training paramters to outs?
*/
	
void mc_model::markov_chain_alignment(ofstream& outs){
	// zero content counting functor
	counting_functor functor(data);
	// make all possible patterns in this class (all_alignment_patterns)
	make_all_alignments_patterns();

	// set all counts to 1 for all context/accession pairs
	for(size_t i = 0; i< data.numAcc();i++){
		for(size_t j = 0; j < data.numAcc();j++){
			for(set<string>::iterator it= all_alignment_patterns.begin(); it != all_alignment_patterns.end() ; it++){
				string pattern = *it;	
				functor.create_context(i, j, pattern);			
			}
		}
	}


	// count all edit operations contained in the alignments (per context)
	for(size_t k = 0; k < data.numAlignments(); k++){
		const pw_alignment & p = data.getAlignment(k);
		computing_modification_oneToTwo(p,functor,outs);
		computing_modification_twoToOne(p,functor,outs);
	}
	functor.total_context();
	for(size_t i = 0; i< data.numAcc();i++){
		for(size_t j = 0; j < data.numAcc();j++){



			for(map <string, vector<double> >::const_iterator it= functor.get_context(i,j).begin();it!=functor.get_context(i,j).end();it++){
				string seq1 = it->first;
				const vector<double> & base = functor.get_context(i,j).at(seq1);
			//	cout<<"base is: "<<endl;
			//	for(size_t a= 0; a< base.size();a++){
				//		cout<< "base at "<< a<< " which is  " << print_modification_character(a)<<" is "<<base.at(a)<<endl;
			//	}
			//	cout<< "context is: "<<endl;
				//	for(size_t m = 0 ; m < seq1.size() ; m++){
				//		cout<< int(seq1.at(m))<<endl;
				//	}
				//	cout<<"the total number of happening the above context between "<<i<<" and "<<j<<" is "<< functor.get_total(i,j,seq1) <<endl;
				for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++) {
				//		cout<<"The number of happening "<< print_modification_character(k) << " between acc " <<i<< " and acc " << j << " after a certain context is "<< base.at(k)<<endl;
				//		cout<<"In MC model, the cost value of" <<print_modification_character(k) <<  " after above pattern between acc " << i << " and acc " << j << " is " << -log2(base.at(k)/functor.get_total(i,j,seq1)) << " bits " << endl; 
					map <string, vector<double> >::iterator it1= mod_cost.at(i).at(j).find(seq1);
					if(it1==mod_cost.at(i).at(j).end()) {
						mod_cost.at(i).at(j).insert(make_pair(seq1, vector<double>((NUM_DELETE+NUM_KEEP+10),1)));
						it1= mod_cost.at(i).at(j).find(seq1);
					}
					it1->second.at(k)=-log2(base.at(k)/functor.get_total(i,j,seq1));
				}
			}
			
			// Compute low and high values TODO 			
			for(set<string>::iterator it= all_alignment_patterns.begin(); it != all_alignment_patterns.end() ; it++){	
				vector<double> num(NUM_DELETE+NUM_KEEP+10,0);
				vector<unsigned int> low(NUM_DELETE+NUM_KEEP+10,0);
				vector<unsigned int> high_value(NUM_DELETE+NUM_KEEP+10,0);
				unsigned int l = 0;
			//	unsigned int total = 0;
				size_t bit = 12; // number of bits to use for encoding event width
				string current_pattern	= *it;
				// it1: high values for current pattern/accession pair
                               	highValue.at(i).at(j).insert(make_pair(current_pattern,vector<unsigned int>(NUM_DELETE+NUM_KEEP+10,0)));
				map<string, vector<unsigned int> >::iterator it1=highValue.at(i).at(j).find(current_pattern);
				assert(it1 != highValue.at(i).at(j).end());
				// it3: get counts for current pattern/accession pair
				map <string, vector<double> >::const_iterator it3= functor.get_context(i,j).find(current_pattern);
				double total =  functor.get_total(i,j,current_pattern);

				for (size_t f=0; f < NUM_DELETE+NUM_KEEP+10;f++){
					low.at(f) = l;
					num.at(f)=it3->second.at(f);
					size_t rescaledNum = (num.at(f)/total)*(powersOfTwo.at(bit) - NUM_DELETE - NUM_KEEP - 11) + 1;
					assert(rescaledNum >= 1);
					assert(rescaledNum < powersOfTwo.at(bit));
				//	cout << "rescled num: "<< rescaledNum << "num: " << num.at(j) << endl;
					high_value.at(f) = l + rescaledNum;
					l = high_value.at(f);
				}
				/*
				for(size_t f = 0; f < NUM_DELETE+NUM_KEEP+10 ; f++){
					if(high_value.at(f)==low.at(f)){
						for(size_t m = 0; m < NUM_DELETE+NUM_KEEP+10-1 ; m++){
							double t = functor.get_total(i,j,current_pattern);
							int rescaledNum = (num.at(f)/t)*(powersOfTwo.at(bit)-5);
							high_value.at(m)= low.at(m)+rescaledNum+1;
							low.at(m+1)=  high_value.at(m);
						}
						high_value.at(NUM_DELETE+NUM_KEEP+9) = low.at(NUM_DELETE+NUM_KEEP+9)+(num.at(f)/t)*(powersOfTwo.at(bit)-5)+1;
						break;
					}
				}
				*/

				// store high values
				for(size_t f = 0; f < NUM_DELETE+NUM_KEEP+10; f++){
					it1->second.at(f)=high_value.at(f);
				}
				/*	for(size_t k =0; k < NUM_DELETE+NUM_KEEP+10; k++){
						cout<<"high value: "<< it1 ->second.at(k)<<endl;
					}*/
			}

		}			
	}
}


const map<string, vector<unsigned int> > & mc_model::get_highValue(size_t acc1, size_t acc2)const{
	return highValue.at(acc1).at(acc2);
}



	void mc_model::cost_function( pw_alignment& p,ofstream & outs) const {
		vector<double> cost_on_sample(2);
		vector<double> modify_cost(2);
		double c1;
		double c2;
		double m1;
		double m2;
		cost_function(p, c1, c2, m1, m2,outs);
		cost_on_sample.at(0) = c1;
		cost_on_sample.at(1) = c2;
		modify_cost.at(0) = m1;
		modify_cost.at(1) = m2;

		p.set_cost(cost_on_sample, modify_cost);
	}
	const vector<vector< map<string, vector<double> > > > &mc_model::get_mod_cost()const{
		return mod_cost;

	}

	void mc_model::cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2,ofstream & outs)const {
	//	p.print();
	//	cout<<"data address in cost function: "<< &data <<endl;
	//	data.numAcc();
		cost_functor f(data,mod_cost);
	//	p.print();
	//	size_t length = p.alignment_length();
	//	cout<<"length: "<< length<<endl;
		computing_modification_oneToTwo(p,f,outs);
		computing_modification_twoToOne(p,f,outs);
		vector<double> cost_on_sample(2,0);
		vector<double> modify_cost(2,0);
		size_t acc1 = data.accNumber(p.getreference1());
		size_t acc2 = data.accNumber(p.getreference2());
		char s1chr;
		char s2chr;
		size_t left1;
		size_t right1;	
		size_t left2;
		size_t right2;
		p.get_lr1(left1,right1);
		p.get_lr2(left2,right2);
		cout<<"left: "<<left1<<"right: "<<right1<<endl;
		for(size_t i = left1; i< right1; i++){
			s1chr = data.getSequence(p.getreference1()).at(i);
			cout<<"check point!"<<endl;
			size_t s1 = dnastring::base_to_index(s1chr);
			stringstream context1;
			for (size_t j = Sequence_level; j>0; j--){
				if(i<j){
					char r1chr = 'A';
					context1 << r1chr;					
				}else{
					char r1chr = data.getSequence(p.getreference1()).at(i-j);
					context1 << r1chr;
					cout<<"rchr: "<<r1chr<< i <<endl;
				}
			}
			string seq1;
			context1>>seq1;
			cout<<"seq1: " << seq1<<endl;
			map <string, vector<double> >::const_iterator it= sequence_successive_bases.at(acc1).find(seq1);
			assert(it != sequence_successive_bases.at(acc1).end());
			cout<<"sequence successive at "<< s1 << " is "<<it->second.at(s1)<<endl;
			cost_on_sample.at(0) += it->second.at(s1);
		}	
	//	cout<<"cost: "<<cost_on_sample.at(0)<<endl;	
		for(size_t i = left2; i<right2; i++){
			s2chr = data.getSequence(p.getreference2()).at(i);
			size_t s2 = dnastring::base_to_index(s2chr);
			stringstream context2;
			for(size_t j = Sequence_level; j>0; j--){
				if(i<j){
						char r2chr = 'A';
						context2 << r2chr;	
				}else{
						char r2chr = data.getSequence(p.getreference2()).at(i-j);
						context2 << r2chr;
				}
			}
			string seq2;
			context2>>seq2;
			map <string, vector<double> >::const_iterator it1= sequence_successive_bases.at(acc2).find(seq2);
			assert(it1 != sequence_successive_bases.at(acc2).end());
		//	cout<<"sequence successive at "<< s2 << " is "<<it1->second.at(s2)<<endl;
			cost_on_sample.at(1) += it1->second.at(s2);
		}	
	//	cout<<"cost: "<<cost_on_sample.at(1)<<endl;										
	/*	for(map <string, vector<double> >::const_iterator it= f.get_context(acc1,acc2).begin();it!=f.get_context(acc1,acc2).end();it++){
			string seq1 = it->first;
			const vector<double> & base = f.get_context(acc1,acc2).at(seq1);
			map <string, vector<double> >::const_iterator it1= mod_cost.at(acc1).at(acc2).find(seq1);
			// cout << " pattern " << seq1 << endl;
//			cout << " acc " << acc1 << " " << acc2 << endl;
//			for(size_t y=0; y<seq1.length(); ++y) {
//				size_t c = seq1.at(y);
//				cout << " p " << y << " = " << c << endl;
//				cout << print_modification_character(c);
//				cout << endl;
//			}
			assert(it1!=mod_cost.at(acc1).at(acc2).end());
			for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++){
				modify_cost.at(0) +=(base.at(k)-1)*(it1->second.at(k));
			}
			cout<< " context is "<<endl;
			for(size_t i = 0 ; i < seq1.size(); i++){
				cout<< int(seq1.at(i))<<endl;
			}
			cout<<"base is: "<<endl;
			for(size_t a= 0; a< base.size();a++){
				cout<< "base at "<< a << " which is " << print_modification_character(a)<<" is "<<base.at(a)<<endl;
				cout<< "modification cost of it is "<< -log2(base.at(a)/f.get_total(acc1,acc2,seq1))<<endl;
			}	

		}
		for(map <string, vector<double> >::const_iterator it= f.get_context(acc2,acc1).begin();it!=f.get_context(acc2,acc1).end();it++){
			string seq1 = it->first;
			cout<< " context2 is "<<endl;
			for(size_t i = 0 ; i < seq1.size(); i++){
				cout<< int(seq1.at(i))<<endl;
			}
			const vector<double> & base = f.get_context(acc2,acc1).at(seq1);
			map <string, vector<double> >::const_iterator it1= mod_cost.at(acc2).at(acc1).find(seq1);
			assert(it1!=mod_cost.at(acc2).at(acc1).end());
			for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++){
				modify_cost.at(1) +=(base.at(k)-1)*(it1->second.at(k));
			}
		//	cout<<"Modification cost on the second ref: " << modify_cost.at(1) << endl;			
		}*/

		c1 = cost_on_sample.at(0);
		c2 = cost_on_sample.at(1);
//		m1 = modify_cost.at(0);
//		m2 = modify_cost.at(1);
		m1 = f.get_modify(p,acc1,acc2);
		m2 = f.get_modify(p,acc2,acc1);
	//	cout<< "length: " << length<<endl;
	//	cout<< "c1: " << c1 << " c2: "<< c2 << " m1: "<< m1<< " m2: "<< m2 <<endl;
	}
	void mc_model::gain_function(const pw_alignment& p, double & g1, double & g2,ofstream & outs)const {
		double c1;
		double c2;
		double m1;
		double m2;
		cost_function(p, c1, c2, m1,m2,outs);

	
		g1 = c2 - m1;
		g2 = c1 - m2;

		//cout << " gain function c2 " << c2 << " m1 " << m1 << " gain1 " << g1 << endl; 
		//cout << " gain function c1 " << c1 << " m2 " << m2 << " gain2 " << g2 << endl; 


	}
	void mc_model::write_parameters(ofstream & outs){
		make_all_the_patterns();
	//	ofstream outs("encode",std::ofstream::binary);
	//	if(outs.is_open()){
			for(size_t i = 0 ; i < data.numAcc(); i++){
				outs << data.get_acc(i);
				cout<<"acc: "<<data.get_acc(i)<<endl;
				outs<< (char)0;
				for(map<string, vector<double> >::iterator it= all_the_patterns.begin(); it != all_the_patterns.end() ; it++){
					string pattern = it ->first;
					map<string,vector<double> >::iterator it1=sequence_successive_bases.at(i).find(pattern);
					if(it1 != sequence_successive_bases.at(i).end()){
						for(size_t n = 0; n <5; n ++){
							it->second.at(n) = it1 ->second.at(n);
						}
					}else{
						for(size_t n =0; n < 5 ; n++)
							it -> second.at(n) = -log2(0.2);
					}
				}
				map <string, vector<unsigned int> >lower_bound;
				for(map<string, vector<double> >::iterator it= all_the_patterns.begin(); it != all_the_patterns.end() ; it++){	
					vector<double> num(5,0);
					vector<bool> bit_to_byte(0);
					vector<unsigned int> low(5,0);
					vector<unsigned int> high_value(5,0);
					unsigned int l = 0;
			//		unsigned int total = 0;
					size_t bit = 12;
					string current_pattern	= it ->first;
                                 	high.at(i).insert(make_pair(current_pattern,vector<unsigned int>(5,0)));
					map<string, vector<unsigned int> >::iterator it1=high.at(i).find(current_pattern);
					assert(it1 != high.at(i).end());
					for (size_t j=0; j < 5;j++){
						low.at(j) = l;
					//	cout<< " low: " << low.at(j)<<endl;
						num.at(j)=it->second.at(j);
						//cout<< "num at " << j << " is " << exp((-num.at(j))*log(2))<<endl;
						int power_of_two = exp((-num.at(j))*log(2))*powersOfTwo.at(bit);
						high_value.at(j) = l + power_of_two;
						l = high_value.at(j);
					//	cout<<"high: "<< high_value.at(j)<<endl;
					}
					for(size_t j = 0; j < 5 ; j++){
						if(high_value.at(j)==low.at(j)){
						//	cout << "low = high" <<endl;
							for(size_t m = 0; m < 4 ; m++){
								int power_of_two = exp((-num.at(m))*log(2))*(powersOfTwo.at(bit)-5);
								high_value.at(m)= low.at(m)+power_of_two+1;
								low.at(m+1)=  high_value.at(m);
							}
							high_value.at(4) = low.at(4)+exp((-num.at(4))*log(2))*(powersOfTwo.at(bit)-5)+1;
							break;
						}
					}
					for(size_t j = 0; j < 5; j++){
						it1->second.at(j)=high_value.at(j);
						int h = high_value.at(j);
					//	cout<< "low1: "<<low.at(j)<<endl;
					/*	if(it->first == "AC"){
							cout<< "high: "<<h<<endl;
						}*/
						for(size_t m = 0; m < bit; m++){
							bit_to_byte.push_back(h%2);
				//			cout<< "h%2: "<< h%2;
							h = h/2;
						}
				//		cout<< " "<<endl;
						assert(high_value.at(j)!=low.at(j));
					}
				//	total = it1->second.at(4);
				//	cout<< "total in stream: "<< it1->second.at(4)<<endl;
				//	cout<<bit_to_byte.size()<<endl;					
					for(size_t n =0; n < bit_to_byte.size()-8; n++){
						unsigned char a = 0;
						for(size_t m = n; m <n+8; m++){
							a+= powersOfTwo.at(m-n)* bit_to_byte.at(m);
						}
						n= n+7;
						outs<< a;
				//		cout<< "eight bits of high: "<<int(a); 
					}
				//	cout << " " << endl;
				}
			}
			outs<<(char)8;
		//}
	//	outs.close();
	}
	void mc_model::write_alignments_pattern(ofstream & outs){
	//	ofstream outs("encode",std::ofstream::binary);
	//	ofstream outs("encode",std::ofstream::binary|std::ofstream::app);
		size_t bit = 12;
	//	if(outs.is_open()){//Per accession!
			for(size_t i = 0 ; i < data.numAcc(); i++){
				for(size_t j =0; j < data.numAcc(); j++){
					outs << data.get_acc(i);
					outs<< (char) 0;
					outs<< data.get_acc(j);
					outs<< (char) 0;
			//		cout<<"acc1 : " << data.get_acc(i) << " acc2: " << data.get_acc(j) << endl;
					for(map<string, vector<unsigned int> >::iterator it= highValue.at(i).at(j).begin(); it != highValue.at(i).at(j).end(); it++){	
						vector<bool> bit_to_byte(0);
						for(size_t j = 0; j < NUM_DELETE+NUM_KEEP+10; j++){
							int h =it->second.at(j);
							for(size_t m = 0; m < bit; m++){
								bit_to_byte.push_back(h%2);
								h = h/2;
							}
						}
					//	size_t counter =0;
					//	cout<< "bit to byte: "<< bit_to_byte.size()<<endl;
						for(size_t n =0; n < bit_to_byte.size()-8; n++){
							unsigned char a = 0;
							for(size_t m = n; m <n+8; m++){
								a+= powersOfTwo.at(m-n)* bit_to_byte.at(m);
							}
							n= n+7;
					//		counter = counter+1;
							outs<< a;
						}	
					//	cout<< "counter: "<< counter << endl;
					}
				//	outs << (char) 0;
				}				
			}
			outs<<(char) 8;
	//	}
	//	outs.close();
	}
/*
	Modification instructions:
	5 - single base modification (incl N)
	NUM_DELETE - delete 2^nd bases
	NUM_KEEP - keep 2^nd bases
	5 - insert a base
	

*/   
void mc_model::make_all_alignments_patterns(){
	string context; 
	set<string>pattern;
	// Alignment_level is makov chain level for alignments
	for(size_t i = 0 ; i < Alignment_level ; i++){
		context += (char)0;
	}
	pattern.insert(context);
	// this will create about (ND+NK+10)^Alignment_length patterns:
	for(size_t i =0; i < Alignment_level; i++) {
		set<string> intermediate_pattern;
		for(size_t j = 0; j <NUM_DELETE+NUM_KEEP+10; j++){
			// For each current pattern: modify position i to character j
			for(set<string>::iterator it = pattern.begin(); it!= pattern.end();it++){
				string seq = *it;
				seq.at(i)=j;	
				
				// throw away some patterns to enforce decreasing order of binary encoding for num keep/delete
				size_t keepthispattern  = 1;
				if(i>0) {
					int mod, del, ins, keep;
					int mod_prev, del_prev, ins_prev, keep_prev;
					modification(seq.at(i), mod, del, ins, keep);
					modification(seq.at(i-1), mod_prev, del_prev, ins_prev, keep_prev);
					if(del > -1 && del_prev > -1) {
						if(del <= del_prev) { // enforce increasing order binary encoding
							keepthispattern = 0;
						}
					
					} else if (keep >-1 && keep_prev > -1) {
						if(keep <= keep_prev) {
							keepthispattern = 0;
						}
					}
				}
				if(keepthispattern)
					intermediate_pattern.insert(seq);
			//	cout<< "pattern: " << seq << endl;
			}
		}
		pattern.clear();
		for(set<string>::iterator it1 = intermediate_pattern.begin(); it1 != intermediate_pattern.end();++it1){
			string seq1 = *it1;
			pattern.insert(seq1);
		}
	} // for Alignment_level
	set<string> intermediate1_pattern;
	for(set<string>::iterator it = pattern.begin(); it != pattern.end();++it){
		for(size_t i = 0 ; i < 6 ; i++){
			string seq = *it;
			char c = i;
			string seq1 = seq + c;
		//	cout << "seq1: " << seq1 <<endl;
			intermediate1_pattern.insert(seq1);
		}
	}
	pattern.clear();
	for(set<string>::iterator it1 = intermediate1_pattern.begin(); it1 != intermediate1_pattern.end();++it1){
			string seq1 = *it1;
			pattern.insert(seq1);
	}
	for(set<string>::iterator it = pattern.begin();it !=pattern.end();it++){
		string seq = *it;
		all_alignment_patterns.insert(seq);
	}
	/*
		size_t number =0;
		for(set<string>::iterator it = pattern.begin(); it != pattern.end();++it){
			string seq = *it;
			number ++ ;
			string str = print_modification_character(seq.at(0));
//			string str1 = print_modification_character(seq.at(1));
			cout<< "" << str << "" <<  "" << int(seq.at(1)) << endl;			
		}
		cout<< "number: "<< number<<endl; 
 	 */
}

// TODO do we want to make markov chain levels dependent on input sequence length?
	void mc_model::train(ofstream & outs){
		make_all_the_patterns();
		markov_chain();
		markov_chain_alignment(outs);
		write_parameters(outs);
		write_alignments_pattern(outs);
	}
	
	void mc_model::make_all_the_patterns(){
		string seq = "";
		set<string> pattern;
		for(size_t i = 0 ; i < Sequence_level ; i++){
			seq += dnastring::index_to_base(0);
		}
		pattern.insert(seq);
		for(size_t i = 0; i < Sequence_level;i++){	
			set<string> intermediate_pattern;
			for(size_t j = 0; j <5; j++){
				for(set<string>::iterator pat = pattern.begin();pat != pattern.end();++pat){
					string seq1 = *pat;
					seq1.at(Sequence_level-1-i)=dnastring::index_to_base(j);	
					intermediate_pattern.insert(seq1);
				}
			}
			for(set<string>::iterator pat1 = intermediate_pattern.begin(); pat1 != intermediate_pattern.end();++pat1){
				string seq2 = *pat1;
				pattern.insert(seq2);
			}
		}	
		for(set<string>::iterator it = pattern.begin();it !=pattern.end();it++){
			string seq3 = *it;
			all_the_patterns.insert(make_pair(seq3,vector<double>(5,0)));
		}	

	}

	void mc_model::set_patterns(ifstream& in){
		make_all_the_patterns();
		size_t bit = 12;
//		ifstream in("encode", std::ifstream::binary);
		char c;
		char h;
		c= in.get();	
		while(c != 8){
		//	cout << " here2"<<endl;
			size_t accession = 0;
			string acc;
			stringstream s;
			while(c != 0){
		//	cout<< "c: " << int(c)<<endl;
		//	cout << " here3"<<endl;
				s << c;
				c = in.get();
			}
			s >> acc;
			data.set_accession(acc);//since in decoding we have no access to our fasta file we need to set accession names in data class
			accession = data.get_acc_id(acc);
		//	cout<<"acc name: "<< acc <<endl;
			for(map<string,vector<double> >::const_iterator it= all_the_patterns.begin(); it!= all_the_patterns.end();it++){
				string pattern = it ->first;
				high.at(accession).insert(make_pair(pattern, vector<unsigned int>(5,0)));
			}
			for(map<string,vector<unsigned int> >::iterator it= high.at(accession).begin(); it!= high.at(accession).end();it++){
				vector<bool> binary_high_value(0);
				size_t bound = (5*bit)/8;
			//	cout<< "bound: " << bound << endl;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					size_t H = size_t(h); // should be a better solution!
			//		cout<<"h: "<< size_t(h)<<endl;
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(H%2);
				//		cout<< " "<< H%2 ;	
						H = H/2;
					}
				//		cout << " " << endl;
				}
				size_t counter = 0;
			//	cout<< "binary size= "<< binary_high_value.size()<<endl;
			//	cout<< "size- 8 "<<  binary_high_value.size() -( binary_high_value.size()%bit)<<endl;
				for(size_t i = 0; i < binary_high_value.size()-bit;i++){
					unsigned int high_value = 0;					
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*powersOfTwo.at(j-i);
					}
					i=i+bit-1;
					it -> second.at(counter)=high_value;
				//	if(it->first == "AC"){
				//		cout<< "high value of AC in set pattern: ";
				//		cout<< high_value<<endl;
				//	}
			//		cout<<"high value in model class: "<< high_value << " at " << counter << " i " << i <<endl;
			//		cout<< " "<<endl;
					counter = counter + 1;
				}
			//	cout<< "counter: "<< counter << endl;
			/*	unsigned int high_value = 0;
				for(size_t j =  binary_high_value.size() -( binary_high_value.size()%bit) ; j <  binary_high_value.size(); j++){
						high_value += binary_high_value.at(j)*powersOfTwo.at(j-binary_high_value.size()+( binary_high_value.size()%bound));					
				}*/
			//	it -> second.at(4)=powersOfTwo.at(bit); 
			//	cout << "high at 4 "<< it->second.at(4)<<endl;
			/*	for(size_t j =0; j < 5 ; j++){
					cout<< "high from stream: " << it->second.at(j)<<endl;
				}*/
			}
			c= in.get();	
		}
	//	cout<<"last c: "<< int(c)<<endl;
	}
	void mc_model::set_alignment_pattern(ifstream & in){
		make_all_alignments_patterns();
		size_t bit = 12;
		char c;
		char h;
		c= in.get();	
		while(c != 8){
			size_t accession1 = 0;
			string acc1;
			stringstream s1;
			while(c != 0){
				s1 << c;
				c = in.get();
			}
			s1 >> acc1;
			data.set_accession(acc1);//since in decoding we have no access to our fasta file we need to set accession names in data class
			accession1 = data.get_acc_id(acc1);
			size_t accession2 = 0;
			string acc2;
			stringstream s2;
			c= in.get();
			while(c != 0){
				s2 << c;
				c = in.get();
			}
			s2 >> acc2;
			data.set_accession(acc2);
			accession2 = data.get_acc_id(acc2);	
			for(set<string>::const_iterator it= all_alignment_patterns.begin(); it!= all_alignment_patterns.end();it++){
				string pattern = *it;
			//	cout<<"acc1: "<<accession1 << " acc2 " << accession2<<endl;
				highValue.at(accession1).at(accession2).insert(make_pair(pattern, vector<unsigned int>(NUM_KEEP+NUM_DELETE+10,0)));
			}
			for(map<string,vector<unsigned int> >::iterator it= highValue.at(accession1).at(accession2).begin(); it!= highValue.at(accession1).at(accession2).end();it++){
				vector<bool> binary_high_value(0);
				size_t bound = ((NUM_KEEP+NUM_DELETE+10)*bit)/8;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					size_t H = size_t(h);
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(H%2);
						H = H/2;	
					}
				}
				size_t counter =0;
				//cout<< "bit to byte 1: "<< binary_high_value.size() << endl;
				for(size_t i = 0; i < (binary_high_value.size())-bit;i++){
					unsigned int high_value = 0;					
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*powersOfTwo.at(j-i);
					}
					i=i+bit-1;
					it -> second.at(counter)= high_value;
					counter = counter +1;
				}
		//		cout<< "counter: "<< counter<<endl;
				it -> second.at(NUM_KEEP+NUM_DELETE+9)=powersOfTwo.at(bit); 
		//		cout<< "size of high value: "<< it ->second.size() << endl;
			}
			c= in.get();	
		}


	}
	vector<unsigned int> mc_model::get_high_at_position(size_t seq_index, size_t position)const{
		const dnastring & sequence = data.getSequence(seq_index);
	//	cout << " sequence " << seq_index << " length " << sequence.length() << endl;
		size_t accession = data.accNumber(seq_index);
		stringstream context;
		for(size_t j = Sequence_level; j>0; j--){
			if(position < j){
				char chr = 'A';
				context<<chr;
			}else{
				char chr = sequence.at(position-j);
				context<<chr;
			}
		}
		string current_pattern;
		context >> current_pattern;
	/*	cout<<"current pattern in seq  "<< seq_index << " of accession " << accession <<endl;
		for(size_t j=0; j< current_pattern.length(); j++){
			cout<<current_pattern.at(j)<<endl;
		}*/
		map<string, vector<unsigned int> >::const_iterator it=high.at(accession).find(current_pattern);
		assert(it!=high.at(accession).end());
	//	cout << " sequence " << seq_index << " length " << sequence.length() << " " << current_pattern<<  endl;
	//	cout<<"accession: "<< accession << endl;
	//	for(size_t k =0; k< 5; k++){
	//		cout<< "high at " << k << " is "<< it->second.at(k)<<endl;
	//	}
	/*	if(current_pattern == "AC"){
			cout<< "high at AC: ";
			for(size_t j = 0; j < 5 ; j++){
				cout<< it->second.at(j) <<endl;
			}
		}*/
		return it->second;
	}
	vector<unsigned int> mc_model::get_center_high_at_position(size_t cent_ref, size_t cent_left, size_t position)const{
		const dnastring & sequence = data.getSequence(cent_ref);
	//	cout << " sequence " << seq_index << " length " << sequence.length() << endl;
		size_t accession = data.accNumber(cent_ref);
		stringstream context;
		for(size_t j = Sequence_level; j>0; j--){
			if(position < cent_left+j){
				char chr = 'A';
				context<<chr;
			}else{
				char chr = sequence.at(position-j);
				context<<chr;
			}
		}
		string current_pattern;
		context >> current_pattern;
	//	cout<<"current pattern: "<< current_pattern<<endl;					
	/*	cout<<"current pattern in seq  "<< seq_index << " of accession " << accession <<endl;
		for(size_t j=0; j< current_pattern.length(); j++){
			cout<<current_pattern.at(j)<<endl;
		}*/
		map<string, vector<unsigned int> >::const_iterator it=high.at(accession).find(current_pattern);
		assert(it!=high.at(accession).end());
	//	cout << " sequence " << seq_index << " length " << sequence.length() << " " << current_pattern<<  endl;
	//	cout<<"accession: "<< accession << endl;
	//	for(size_t k =0; k< 5; k++){
	//		cout<< "high at " << k << " is "<< it->second.at(k)<<endl;
	//	}
		return it->second;
	}
	vector<unsigned int> mc_model::get_reverse_center_high_at_position(size_t cent_ref, size_t cent_right, size_t position)const{
		const dnastring & sequence = data.getSequence(cent_ref);
	//	cout << " sequence " << seq_index << " length " << sequence.length() << endl;
		size_t accession = data.accNumber(cent_ref);
		stringstream context;
		for(size_t j = Sequence_level; j>0; j--){
			if(position > cent_right-j){
				char chr = 'A';
				context<<chr;
			}else{
				char chr= dnastring::complement(sequence.at(position+j));				
			//	char chr = sequence.at(position+j);
				context<<chr;
			}
		}
		string current_pattern;
		context >> current_pattern;
	//	cout<<"current pattern: "<< current_pattern<<endl;					
	/*	cout<<"current pattern in seq  "<< seq_index << " of accession " << accession <<endl;
		for(size_t j=0; j< current_pattern.length(); j++){
			cout<<current_pattern.at(j)<<endl;
		}*/
		map<string, vector<unsigned int> >::const_iterator it=high.at(accession).find(current_pattern);
		assert(it!=high.at(accession).end());
	//	cout << " sequence " << seq_index << " length " << sequence.length() << " " << current_pattern<<  endl;
	//	cout<<"accession: "<< accession << endl;
	//	for(size_t k =0; k< 5; k++){
	//		cout<< "high at " << k << " is "<< it->second.at(k)<<endl;
	//	}
		return it->second;
	}

	vector<size_t> mc_model::get_powerOfTwo()const{
		return powersOfTwo;
	}		
	const map<string, vector<unsigned int> >&  mc_model::get_high(size_t acc)const{
		return high.at(acc);
	}
	string mc_model::get_context(size_t position, size_t seq_id)const{
		const dnastring & sequence = data.getSequence(seq_id);
		stringstream context;
		for(size_t j = Sequence_level; j > 0; j--){
			if(position<j){
				char chr = 'A';
				context << chr;
			}else{
				char chr = sequence.at(position-j);
				context << chr;
			}
		}
		string current_pattern;
		context >> current_pattern;	
		return current_pattern;
	}
	string mc_model::get_firstPattern()const{
		string first_pattern;
		for(size_t i = 0; i< Sequence_level; i++){
			first_pattern += "A";
		}
		return first_pattern;
	}
	string mc_model::get_firstAlignmentPattern()const{
		string first_pattern;
		size_t firstPatterns = Alignment_level;
		for(size_t j = 0; j < Alignment_level; j++){
			firstPatterns --;
			first_pattern+=modification_character(-1,-1,-1,firstPatterns);
		}
		return first_pattern;
	}	
	char mc_model::modification_character(int modify_base, int num_delete, int insert_base, int num_keep)const {
		//return the enc
		if(num_delete != -1){
		//	cout<< "there is a delete of length "<< num_delete <<endl; 
		//	return NUM_DELETE+num_delete; TODO this was really wrong
			return 5 + num_delete;
		}
		if(num_keep != -1){
		//	cout<< "there is a keep of length" << num_keep << endl;
		//	return NUM_KEEP+num_keep;
			return 5 + NUM_DELETE + num_keep;
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
	string mc_model::print_modification_character(char enc)const{
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
			s<< " a modification to " << dnastring::index_to_base(modify_base);
		}
		if(insert_base != -1){
			s<< " a insertion at " << dnastring::index_to_base(insert_base);
		}
		return s.str();
	}

// TODO what does this do?
// It is wrong because modification returns 2^num_delete
/*
	void mc_model::print_modification(char enc)const{
		int modify_base = -1;
		int num_delete =-1;
		int insert_base = -1;
		int num_keep = -1;
		modification(enc, modify_base, num_delete, insert_base, num_keep);
		size_t modification_type = -1;
		if(num_delete != -1){
			 modification_type = 5+num_delete;
			 cout << " delete " << num_delete << endl;
		}
		if(num_keep != -1){
			modification_type = 5+NUM_DELETE+num_keep;
			 cout << " keep " << num_keep << endl;
		}
		if(modify_base != -1) {
			modification_type = modify_base; 
			cout << " replace to " << modify_base << endl;
		}
		if(insert_base != -1){
			cout << " insert: " << insert_base << endl;
			modification_type = insert_base + NUM_KEEP + NUM_DELETE + 5;

		}
		assert(enc == modification_type);
	//	return modification_type;

	}
#
*/
	size_t mc_model::modification_length(char mod)const{
		int modify_base = -1;
		int num_delete =-1;
		int insert_base = -1;
		int num_keep = -1;
		size_t length = 0;
		modification(mod, modify_base, num_delete, insert_base, num_keep);
		if(modify_base != -1){
			length = 1;
		}
		if(insert_base != -1){
			length = 1;
		}
		if(num_delete != -1){
			length = num_delete;
		}
		if(num_keep != -1){
			length = num_keep;
		}
		return length;
	}
	void mc_model::modification(char enc, int & modify_base, int & num_delete, int & insert_base, int & num_keep)const {
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
			insert_base = enc- 5 - NUM_KEEP - NUM_DELETE;
			return;
		}
		assert(0);

	}

// TODO 
// this function should get a better name. It applies the functor on each position in the alignment
	void mc_model::computing_modification_oneToTwo(const pw_alignment & p, abstract_context_functor & functor,ofstream & outs)const{
		// TODO new entropy encoder makes no sense here
		string seq = "";
	//	cout<<"data ad in computing mod: "<< & data << endl;
	//	p.print();
		size_t acc1 = data.accNumber(p.getreference1());
		size_t acc2 = data.accNumber(p.getreference2());
		size_t first_patterns = Alignment_level;
		size_t power = powersOfTwo.at(NUM_KEEP-1);
	//	cout<< "power: "<<powersOfTwo.at(0)<<endl;	
		for(size_t j = 0; j < Alignment_level; j++){// two keeps of length 2^1 and 2^0 has been created at the begining of each alignment.
			first_patterns--;
			seq+=modification_character(-1,-1,-1,first_patterns);
		}
	//	cout<<p.alignment_length()<<endl;
		for (size_t i = 0; i< p.alignment_length(); i++){
			size_t n = 0;
			int modify_base =-1;
			int num_delete=-1;
			int insert_base=-1;
			int num_keep=-1;
			string seq1(" ",Alignment_level+1);
			char seq2;
			for(size_t w = Alignment_level; w>0 ;w--){
				seq1.at(Alignment_level-w)=seq.at(seq.size()-w);
			}
			char s1chr;
			char s2chr;
			size_t s1;
			size_t s2;
			p.alignment_col(i, s1chr, s2chr);				
			s1 = dnastring::base_to_index(s1chr);
			s2 = dnastring::base_to_index(s2chr);
			seq1.at(Alignment_level)=s1;
			if(s1 == s2){
				size_t klength = 0;
				for(size_t j = i; j<p.alignment_length(); j++){
					char q1chr;
					char q2chr;
					p.alignment_col(j, q1chr, q2chr);				
					size_t q1 = dnastring::base_to_index(q1chr);
					size_t q2 = dnastring::base_to_index(q2chr);
					if(q1 == q2){
						klength +=1;
					}else{	
						break;
					}
				}
			/*	if(i<51){
					cout<<"keep length at "<< i << " is  "<< klength<<endl;
				}*/
			//	cout<<"keep length: "<<klength<<endl;
			//	cout<<"NUM KEEP: "<< NUM_KEEP-1<<endl;
				if(klength > power){
					num_keep = NUM_KEEP-1;
					n=powersOfTwo.at(num_keep)-1;
				//	cout<<"Long keep"<<endl;
					seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
				}else{
					for (size_t m = (NUM_KEEP-1); m > 0; m--){
				//	for (size_t m = powersOfTwo.size()-1; m >= 0; m--)
						if((klength & powersOfTwo.at(m)) != 0){
							num_keep=m;
						//	cout<<"m: "<<m <<endl;
							n= powersOfTwo.at(num_keep)-1;
							seq += modification_character(modify_base,num_delete,insert_base,num_keep);
							break;
						}
					}
				}
				seq2 = seq.at(seq.size()-1);
		//		cout<< "recorded keep length: " << n+1 << endl;
		//		cout<< "size of seq: "<< seq.size()<<endl;
		//		cout << "seq 2 " << int(seq2) << endl;
		//		cout<< "recorded keep length: " << n+1 << endl;
				functor. see_context(acc1,acc2,p,i,seq1,seq2, outs);
		//		cout<<"n: "<< n << endl;
			}else{
				if((s1!=5) & (s2!=5)){
					modify_base = s2;
				/*	if(i<51){
						cout<<"modification at "<< i << " is  "<<s1 <<endl;
					}*/
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc1,acc2,p,i,seq1,seq2,outs);
				//	cout<< "seq1" << seq1 <<endl;
				}
				if(s1 == 5){
					insert_base = s2;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);						
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc1,acc2,p,i,seq1,seq2,outs);
			//		cout<< "seq1" << seq1 <<endl;						
				}
				if(s2 == 5){
					size_t dlength = 0;
					for(size_t j = i; j < p.alignment_length(); j++){
						char q1chr;
						char q2chr;
						p.alignment_col(j, q1chr, q2chr);				
					//	size_t q1 = dnastring::base_to_index(q1chr);
						size_t q2 = dnastring::base_to_index(q2chr);
						if(q2 == 5 ){
							dlength +=1;
						}else {
							break;
						}
					}
					if(dlength > powersOfTwo.at(NUM_DELETE-1)){
//					if(dlength > powersOfTwo.at(powersOfTwo.size()-1))
//						num_delete= powersOfTwo.size()-1;
						num_delete = NUM_DELETE-1;
						n= powersOfTwo.at(num_delete)-1;
						seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
					}else{
						for(size_t m = (NUM_DELETE-1); m > 0; m--){
//						for (size_t m = powersOfTwo.size()-1;m>=0;m--)
							if((dlength & powersOfTwo.at(m)) != 0){
								num_delete = m;
								n= powersOfTwo.at(num_delete)-1;
								seq += modification_character(modify_base,num_delete,insert_base,num_keep);		
								break;
							}
						}
					}
					seq2 = seq.at(seq.size()-1);
					functor. see_context(acc1,acc2,p,i,seq1,seq2,outs);
				}

			}
			/*cout<< " context1 is "<<endl;
			for(size_t h = 0 ; h < seq1.size(); h++){
				cout<< int(seq1.at(h))<<endl;
			}*/
			//cout<<"last char: "<<int(seq2)<<endl;
			//cout<<"i+n: "<< i+n<<endl;
//			cout<< " context1 is "<<endl;
//			for(size_t h = 0 ; h < seq1.size(); h++){
//				cout<< int(seq1.at(h))<<endl;
//			}
	//		cout<<"last char in mod function: "<<int(seq2)<<endl;
//			cout<<"i+n: "<< i+n<<endl;
			i=i+n;
			//cout<<"n"<< n << endl;
			//	uint32_t initial = 1;
			//	for (uint32_t m = 0; m <32 ; m++)
			//		uint32_t power_of_two = initial << m; 
			//		if((klength & power_of_two) !=0)
		}

	//	cout<<"The alignment is: "<<endl;
	//	p.print();
	//	cout<<"encoded sequence from one to two is: "<<endl;
	//	for(size_t m = 0; m < seq.size(); m ++){
	//		cout<< int(seq.at(m))<<endl;
	//	}
	//	cout<< "one to two was done! " << endl;
		functor.see_entire_context(acc1,acc2,seq);
	}

	
	void mc_model::computing_modification_twoToOne(const pw_alignment & p, abstract_context_functor & functor,ofstream & outs)const{
		string seq = "";
		size_t acc1 = data.accNumber(p.getreference1());
		size_t acc2 = data.accNumber(p.getreference2());
		size_t first_patterns = Alignment_level;
		for(size_t j = 0; j < Alignment_level; j++){
			first_patterns --;
			seq+=modification_character(-1,-1,-1,first_patterns);
		}
		for (size_t i = 0; i< p.alignment_length(); i++){
			int modify_base =-1;
			int num_delete=-1;
			int insert_base=-1;
			int num_keep=-1;
			size_t n = 0;			
			string seq1(" ",Alignment_level+1);
			char seq2;
			for(size_t w = Alignment_level; w>0 ;w--){
				seq1.at(Alignment_level-w)=seq.at(seq.size()-w);
			//	cout<< "seq at size - w: " << seq.at(seq.size()-w)<<endl;
			}
			char s1chr;
			char s2chr;
			size_t s1;
			size_t s2;
			p.alignment_col(i, s1chr, s2chr);				
			s1 = dnastring::base_to_index(s1chr);
			s2 = dnastring::base_to_index(s2chr);
			seq1.at(Alignment_level)=s2;
			if(s1 == s2){
				size_t klength = 0;
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
				if(klength > powersOfTwo.at(NUM_KEEP-1)){
//				if(klength > powersOfTwo.at(powersOfTwo.size()-1))
//					num_keep= powersOfTwo.size()-1;
					num_keep = NUM_KEEP-1;
					n=powersOfTwo.at(num_keep)-1;
					seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
				}else{
					for(size_t m = NUM_KEEP-1; m > 0; m--){
//					for (size_t m =powersOfTwo.size()-1; m>= 0; m--)
						if((klength & powersOfTwo.at(m)) != 0){
							num_keep=m;
							n=powersOfTwo.at(num_keep)-1;						
							seq += modification_character(modify_base,num_delete,insert_base,num_keep);
							break;
						}
					}
				}
				seq2 = seq.at(seq.size()-1);
				functor. see_context(acc2,acc1,p,i,seq1,seq2,outs);
			}else{
				if((s1!=5) & (s2!=5)){
					modify_base = s1;
				/*	if(i<51){
						cout<<"modification at "<< i << " is  "<<s2 <<endl;
					}*/
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc2,acc1,p,i,seq1,seq2,outs);
				}
				if(s2 == 5){
					insert_base = s1;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);					
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc2,acc1,p,i,seq1,seq2,outs);						
				}
				if(s1 == 5){
					size_t dlength = 0;
					for(size_t j = i; j < p.alignment_length(); j++){
						char q1chr;
						char q2chr;
						p.alignment_col(j, q1chr, q2chr);				
						size_t q1 = dnastring::base_to_index(q1chr);
					//	size_t q2 = dnastring::base_to_index(q2chr);
						if(q1 == 5 ){
							dlength +=1;
						}else {
							break;
						}
					}
					if(dlength > powersOfTwo.at(NUM_DELETE-1)){
//					if(dlength > powersOfTwo.at(powersOfTwo.size()-1))
//						num_delete= powersOfTwo.size()-1;
						num_delete = NUM_DELETE-1;
						n = powersOfTwo.at(num_delete)-1;
						seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
					}else{
						for(size_t m = NUM_DELETE -1 ; m > 0; m--){
//						for (size_t m = powersOfTwo.size()-1; m>0; m--)
							if((dlength & powersOfTwo.at(m)) != 0){
								num_delete = m;
								n = powersOfTwo.at(num_delete)-1;
								seq += modification_character(modify_base,num_delete,insert_base,num_keep);			
							}
						}
					}
					seq2 = seq.at(seq.size()-1);
					functor. see_context(acc2,acc1,p,i,seq1,seq2,outs);
				}
			}
	/*		cout<< " context2 is "<<endl;
			for(size_t h = 0 ; h < seq1.size(); h++){
				cout<< int(seq1.at(h))<<endl;
			}*/
		//	cout<<"last char in mod function: "<<int(seq2)<<endl;
			i=i+n;
		//	cout<< "i in modification : " << i << endl;
		}		
	//	cout<<"encoded sequence from two to one is: "<<endl;
	/*	for(size_t m = 0; m < seq.size(); m ++){
			cout<< int(seq.at(m))<<endl;
		}*/
		functor.see_entire_context(acc2,acc1,seq);
	}
	const map<string, vector<double> > & mc_model::getPattern(size_t acc)const{
		return sequence_successive_bases.at(acc);

	}
	
	const vector<double> & mc_model::get_create_cost(size_t acc)const{
		return create_cost.at(acc);		

	}
	
	const vector<map<string, vector<double> > > & mc_model::model_parameters()const{
		return sequence_successive_bases;
//			cout<< " context2 is "<<endl;
//			for(size_t h = 0 ; h < seq1.size(); h++){
//				cout<< int(seq1.at(h))<<endl;
//		}
//		cout<<"last char: "<<int(seq2)<<endl;	
//		cout<<"encoded sequence from two to one is: "<<endl;
//		for(size_t m = 0; m < seq.size(); m ++){
//			cout<< int(seq.at(m))<<endl;
//		}

	}
//	const map<string, vector<double> > & mc_model::get_alignment_context(size_t al_id, size_t seq_id, encoding_functor & functor)const{
	/*	const pw_alignment & al = data.getAlignment(al_id);
		size_t acc1 = data.accNumber(al.getreference1());
	//	size_t acc2 = data.accNumber(al.getreference2());
		size_t accession = data.accNumber(seq_id);
		if(accession == acc1){
			computing_modification_oneToTwo(al, functor);			
			return functor.get_alignment_context();
		}else{
			computing_modification_twoToOne(al, functor);
			return functor.get_alignment_context();
		}*/
//		const map<string, vector<double> > & res = *((const map<string, vector<double> >*) NULL);
//		return res;
//	}
//	const map<string, vector<double> > &mc_model::get_cluster_member_context(pw_alignment & al, size_t center_id, encoding_functor & functor)const{
	/*	size_t acc1 = data.accNumber(al.getreference1());
	//	size_t acc2 = data.accNumber(al.getreference2());
		size_t accession = data.accNumber(center_id);
		if(accession == acc1){
			computing_modification_oneToTwo(al, functor);	
			return functor.get_alignment_context();	
		}else{
			computing_modification_twoToOne(al, functor);
			return functor.get_alignment_context();
		}*/
//		const map<string, vector<double> > & res = *((const map<string, vector<double> >*) NULL);
//		return res;
//	}
/*	void mc_model ::get_encoded_member(pw_alignment & al, size_t center_id, encoding_functor & functor,ofstream& outs)const{
		size_t acc1 = data.accNumber(al.getreference1());
		size_t accession = data.accNumber(center_id);
		if(accession == acc1){
			cout<< "center is on acc1"<<endl;
			computing_modification_oneToTwo(al, functor,outs);	
		}else{
			computing_modification_twoToOne(al, functor,outs);
		}
	}
*/
	void mc_model::computing_modification_in_cluster(string center, string member)const{
	// we should read their bases one by one find the modification pattern between them, information be gained from clustering class
	// in main function call it over global results, but then it is not part of training! :|
	// different sequences in a cluster may have different length, but we always need to check the modification from the center. Thus just go through the center size.
		//alignment between center and member	
		for(size_t i = 0; i < center.size(); i ++){
			size_t c = dnastring::base_to_index(center.at(i));
			size_t m = dnastring::base_to_index(member.at(i));
			if(c==m){
				size_t klength = 0;
				for(size_t j = i; j<center.size(); j++){
					size_t c1 = dnastring::base_to_index(center.at(j));
					size_t m1 = dnastring::base_to_index(member.at(j));
					if(c1 == m1){
						klength +=1;
					}else {	
						break;
					}
				}


			}else {


			}
		}

	}
	abstract_context_functor::abstract_context_functor(){
	
	}
	void abstract_context_functor::see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, string context, char last_char, ofstream& outs){
		
	}
	void abstract_context_functor::see_entire_context(size_t acc1,size_t acc2, string entireContext){

	}
	counting_functor::counting_functor(all_data & d):data(d), successive_modification(d.numAcc(),vector<map<string, vector<double> > >(d.numAcc())),total(d.numAcc(),vector<map<string, double > >(d.numAcc())) {}
// TODO why do we use double for counting?
	void counting_functor::see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, string context, char last_char, ofstream & outs){
	//	cout<< "accession 1: " << acc1 << " accession 2: " << acc2 << " size: " << pos << " last char: " << dnastring::base_to_index(last_char) << " " << int(last_char)<<endl;
	//	cout<< "context is: "<< endl;
	/*	for(size_t i = 0 ; i < context.size(); i++){
			cout<< int(context.at(i))<<endl;
		}*/
		map <string, vector<double> >::iterator it1= successive_modification.at(acc1).at(acc2).find(context);
		if(it1==successive_modification.at(acc1).at(acc2).end()) {
			successive_modification.at(acc1).at(acc2).insert(make_pair(context, vector<double>((NUM_DELETE+NUM_KEEP+10),1)));
			it1= successive_modification.at(acc1).at(acc2).find(context);
			assert(it1 != successive_modification.at(acc1).at(acc2).end());
		}
			it1->second.at(last_char)++;
		//	cout<< "context is: "<< context.length() << endl;
			//cout<< context <<endl;			
		//	cout<<"number of happening "<<int(last_char)<< " after above context is "<< it1->second.at(last_char)<<endl;
	}

/*
	compute total (sum over all counts)

*/
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
				// we assume to find everything at least twice in the training to avoid events with zero information cost
				if(it1->second < 2) {
					it1->second = 2;
				}
			}
		}
	}
}

double counting_functor::get_total(size_t acc1, size_t acc2, string context)const{
	map<string, double>::const_iterator it = total.at(acc1).at(acc2).find(context);
//	cout << " a " << acc1 << " " << acc2 << " : " << total.at(acc1).at(acc2).size() << endl;
	assert(it!=total.at(acc1).at(acc2).end());
	return it->second;
}

			
/*
	initialize context counting with 1
*/   	
	void counting_functor::create_context(size_t acc1, size_t acc2, string context) {
		map<string, vector<double> >::iterator it = successive_modification.at(acc1).at(acc2).find(context);
		if(it==successive_modification.at(acc1).at(acc2).end()) {
			successive_modification.at(acc1).at(acc2).insert(make_pair(context, vector<double>(NUM_DELETE+NUM_KEEP+10, 1)));
		}
	}
	
	const  map<string, vector<double> > & counting_functor::get_context(size_t acc1, size_t acc2)const{
		return successive_modification.at(acc1).at(acc2);
	}

	cost_functor::cost_functor(all_data & d, const vector<vector<map<string, vector<double> > > > & mod_cost):data(d){
		modify1 = 0;
		modify2 = 0;		
		modification = mod_cost;
	}
	void cost_functor::see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, string context, char last_char, ofstream & outs){
		size_t ref1 = p.getreference1();
		size_t ref2 = p.getreference2();
		size_t accession1 = data.accNumber(ref1);
		size_t accession2 = data.accNumber(ref2);
		if(acc1 == accession1){//if acc1 is the first accession
			map<string, vector<double> >::const_iterator it = modification.at(acc1).at(acc2).find(context);
		//	cout<< "modification cost at "<< pos << " is "<< it->second.at(last_char)<<endl;
			modify1 +=it->second.at(last_char);
		}
		if(acc1 == accession2){
			map<string, vector<double> >::const_iterator it = modification.at(acc1).at(acc2).find(context);
		//	cout<< "modification cost at "<< pos << " is "<< it->second.at(last_char)<<endl;
			modify2 +=it->second.at(last_char);
		}
	}
	double cost_functor::get_modify(const pw_alignment & p,size_t acc1, size_t acc2)const{
		double modify;
		size_t ref1 = p.getreference1();
		size_t ref2 = p.getreference2();
		size_t accession1 = data.accNumber(ref1);
		size_t accession2 = data.accNumber(ref2);
		if(acc1 == accession1){
			modify = modify1; 
		}
		if(acc1 == accession2){
			modify = modify2;
		}
		return modify;
	}

	/*
	encoding_functor::encoding_functor(all_data & d, mc_model * m, wrapper & wrap):data(d),model(m),wrappers(wrap){
	}
	
	void encoding_functor::see_context(size_t acc1, size_t acc2,const pw_alignment & p, size_t pos, string context, char last_char, ofstream& outs){//last_char is infact a pattern!
		size_t bit = 13;
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		enc->set_stream(outs);
		unsigned int total = model->get_powerOfTwo().at(bit)+20;
		map<string, vector<unsigned int> >::const_iterator it1 = model->get_highValue(acc1,acc2).find(context);// if modification is from acc2 to acc1 the order is already exchanged. So this is true
		vector<unsigned int> low(NUM_DELETE+NUM_KEEP+10,0);
		vector<unsigned int> high(NUM_DELETE+NUM_KEEP+10,0);
		for(size_t m = 0; m < it1->second.size(); m++){
			if(m ==0){
				low.at(m) = 0;
			}else{
			 	low.at(m) = high.at(m-1);
			}
			high.at(m) = it1->second.at(m);
			if(m == NUM_DELETE+NUM_KEEP+10-1){
				high.at(m)=  model->get_powerOfTwo().at(bit);
			}
		}
		cout<< "context: ";
		for(size_t i =0; i < context.size(); i++){
			cout<< int(context.at(i));
		}
		cout<< " " <<endl;
		cout << " ended char in al: "<< int(last_char)<<endl;
		enc->encode(low.at(last_char),high.at(last_char),total);
		wrappers.encode(low.at(last_char),high.at(last_char),total);
		delete enc;
	}



	void encoding_functor::see_entire_context(size_t acc1, size_t acc2, string entireContext){
		alignment_pattern = entireContext;
	}
	const map<string, vector<double> > & encoding_functor::get_alignment_context()const{
		for(map<string, vector<double> >::const_iterator it = alignment_context.begin(); it !=alignment_context.end(); it++){
*/
		/*	for(size_t n =0; n < NUM_DELETE+NUM_KEEP+10; n++){
				cout << "counts: "<< it->second.at(n)<<endl;
			}*/
/*
		}
		return alignment_context;
	}
*/
/*
	vector<string>  &  encoding_functor::get_alignment_context(pw_alignment & p)const{
		vector<string>longName;
*/
	/*	pw_alignment * copy_p = new pw_alignment(p);
		map<pw_alignment*, vector<string> >::const_iterator it = pattern.find(copy_p);
		assert(it != pattern.end());
		delete copy_p;
		return it->second;*
		for(map<string , vector<double> >::const_iterator it = alignment_context.begin(); it!=alignment_context.end(); it++){
			string context = it->first;
			for(size_t i = 0; i < it->second.size();i++){
				stringstream longname;
				if(it->second.at(i)!=0){
					longname << context << ":" << it->second.at(i);
					longName.push_back(longname.str());
				}
				
			}
		}
		cout<<"new longname: "<<endl;
		for(size_t j =0; j < longName.size(); j++){
			cout<< longName.at(j);
		}
		cout<< " " <<endl;
		return longName;
	}

	*/

