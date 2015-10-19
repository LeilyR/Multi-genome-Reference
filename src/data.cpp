#include "data.hpp"


void strsep(std::string str, const char * sep, std::vector<std::string> & parts) {
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


dnastring::dnastring(std::string str):bits(str.length()*3) {
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
				std::cerr << "Error: Illegal character in DNA sequence: " << c << std::endl;
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

std::string dnastring::str() const {
	std::string res("");
	for(size_t i=0; i<length(); ++i) {
		res.append(1, at(i)); // TODO faster?
	}
	return res;
}


char dnastring::base_translate_back(bool bit1, bool bit2, bool bit3) {

	if(bit1) {
		if(bit2) {
			if(bit3) {
				std::cerr << "Error: wrong bit std::vector 111" << std::endl;
				exit(1);
				return 'X';
			} else {
				return 'G';
			}
		
		} else {
			if(bit3) {
				std::cerr << "Error wrong bit std::vector 101" << std::endl;
				exit(1);
				return 'X';
			} else {
				return 'C';
			}

		}	
	
	} else {
		if(bit2) {
			if(bit3) {
				std::cerr << "Error: wrong bit std::vector 011" << std::endl;
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
all_data::all_data() {}

void all_data::read_fasta_sam(std::string fasta_all_sequences, std::string sam_all_alignments) {
	std::cout << "read sam File "<< std::endl;
	std::ifstream fastain(fasta_all_sequences.c_str());
	if(fastain) {
		std::string str;
		std::stringstream curseq;
		std::string curname("");
		std::string curacc("noAcc"); //TODO put default acc !
		while(getline(fastain, str)) {
			if(str.at(0)=='>') {
				if(0!=curname.compare("")) { // store previous sequence
					insert_sequence(curacc, curname, curseq.str());
					curname = "";
					curacc = "noAcc";
					curseq.str("");
				}

				// read next header
				//std::string fhead = str.substr(1);
				curname = str.substr(1);
			//	name_split(fhead, curacc, curname);
			} else {
				curseq << str;
			}
		}
		// store last sequence
		insert_sequence(curacc, curname, curseq.str());
		//cout << "R " << curseq.str().substr(62750, 15) << endl;
		curname = "";
		curacc = "noAcc";
		curseq.str("");

		fastain.close();
	} else {
		std::cerr << "Error: cannot read: " << fasta_all_sequences << std::endl;
		exit(1);
	}

	std::vector< std::multimap< size_t, size_t> > als_on_reference; // sequence index -> begin pos on that sequence -> alignment index
	// a new multimap for each ref sequence
	als_on_reference.resize(sequences.size());
	for(size_t i=0; i<sequences.size(); ++i) {
		als_on_reference.at(i) = std::multimap<size_t, size_t>();
	}

	// READ SAM FILE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	/* Achtung
	 * Library libStatGen creates an output -bs.umfa : if already exist can have some conflict, need to remove it before rerunning this code
	 */
	bool verbose = false;
	  std::cout << "readAlignment sam file"<<std::endl;

	  SamFile samIn;
	  samIn.OpenForRead(sam_all_alignments.c_str());

	  // Read the sam header.
	  SamFileHeader samHeader;
	  samIn.ReadHeader(samHeader);

	  SamRecord samRecord;
	  Cigar* tmpCigar;
	  std::string alRefSeq, alReadSeq;
	  int skip_alt = 0;
	  int skip_self = 0;
	  GenomeSequence myRefSeq(fasta_all_sequences.c_str());
	  myRefSeq.setDebugFlag(1);
	  while(samIn.ReadRecord(samHeader, samRecord))
	  {
		// Achtung : sam file considere the file reference as an uniq sequence even if there are differentes sequences in reference, so the position are concatenate.
		//e.g. if I take the first base of the second sequence from my reference.fasta, it will not return 1 but lengthOfTheSequence1 + 1 !

		// For each Record do :
		tmpCigar = samRecord.getCigarInfo(); //Pointer to Cigar object
		const char* currentSeq = samRecord.getSequence();

		int myStartOfReadOnRefIndex = myRefSeq.getGenomePosition(samRecord.getReferenceName(),samRecord.get1BasedPosition()); // Select the start on the good sequence in the fasta File
		if(myStartOfReadOnRefIndex == -1){
			std::cerr << "Error with the reference " <<samRecord.get1BasedPosition()<<" " <<samRecord.getReferenceName() << std::endl;
			exit(1);
		}

		if (verbose) std::cout<<myStartOfReadOnRefIndex<<" getReferenceName "<<samRecord.getReferenceName()<<" get1BasedPosition "<< samRecord.get1BasedPosition() <<" myRefSeq.sequenceLength " <<myRefSeq.sequenceLength()
		  << " getAlignmentLength "<< samRecord.getAlignmentLength() << " getReadLength " << samRecord.getReadLength()<< " getNumBeginClips "<< tmpCigar->getNumBeginClips()<< std::endl;
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
		//int myEndReadOnRefIndex =  myRefSeq.getGenomePosition(samRecord.getReferenceName(),samRecord.get1BasedAlignmentEnd()) + 1 ;
		//if(verbose) cout << "reference start "<< myStartOfReadOnRefIndex << " end " << myEndReadOnRefIndex <<endl;

		// Reference
		std::string acc1("noAcc");
		std::string name1;
		//name_split(samRecord.getReferenceName(), acc1, name1);
		size_t start1 = samRecord.get1BasedPosition() - 1 ;
		size_t incl_end1 = start1 + tmpCigar->getExpectedReferenceBaseCount() - 1 ;
		std::string tmpString = acc1 + ":" + samRecord.getReferenceName();
		//std::cout << "tmpString " << tmpString << " long " << longname2seqidx.begin()->first << std::endl;
		std::map<std::string, size_t>::iterator findseq1 = longname2seqidx.find(tmpString);
		if(findseq1==longname2seqidx.end()) {
			std::cerr << "Error: unknown sequence in sam file : " << samRecord.getReferenceName() << std::endl;
			exit(1);
		}
		size_t idx1 = findseq1->second;

		// Read
		std::string acc2("noAcc");
		std::string name2;
		//name_split(samRecord.getReadName(), acc2, name2);
		size_t start2;
		if(tmpCigar->getNumBeginClips()==0)
			start2 = 0;
		else
			start2 = tmpCigar->getNumBeginClips()+1; // getNumbeginClips give the number of clip and the read start at the next one ( so +1)

		size_t incl_end2 = start2 + samRecord.getReadLength() -1;//TODO -1 ?

		// In Sam file : Reverse and Forward ??
		std::string tmpString2 = acc2 + ":" + samRecord.getReadName();

		std::map<std::string, size_t>::iterator findseq2 = longname2seqidx.find((tmpString2));
		if(findseq2==longname2seqidx.end()) {
			std::cerr << "Error: unknown sequence in sam file : " << samRecord.getReadName() << std::endl;
			exit(1);
		}
		size_t idx2 = findseq2->second;


    	// both al parts not identical	
	if(idx1 != idx2 || !(start1==start2 && incl_end1 == incl_end2    )) {
		bool skip = false;
		// check if we already have an alignment with same coordinates
		// because we are looking for identical alignment it suffices to go over only one multimap
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> eqr =
		als_on_reference.at(idx1).equal_range(start1);
		for(std::multimap<size_t, size_t>::iterator it = eqr.first; it!=eqr.second; ++it) {
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
			als_on_reference.at(idx1).insert(std::make_pair(start1, alidx));
			als_on_reference.at(idx2).insert(std::make_pair(start2, alidx));
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
		std::cerr << "Warning: IUPAC DNA ambiguity characters were replaced by N" << std::endl;
	}

	std::cout << "Loaded: " << sequences.size() << " sequences and " << alignments.size() << " pairwise alignments " << std::endl;
	std::cout << skip_self << " self alignments were skipped" << std::endl;
	std::cout << skip_alt << " alternative alignments of identical regions were skipped" << std::endl;




}


void all_data::read_fasta_maf(std::string fasta_all_sequences, std::string maf_all_alignments) {
	std::ifstream fastain(fasta_all_sequences.c_str());
	if(fastain) {
		std::string str;
		std::stringstream curseq;
		std::string curname("");
		std::string curacc("");
		while(getline(fastain, str)) {
			if(str.at(0)=='>') {
				if(0!=curname.compare("")) { // store previous sequence
					insert_sequence(curacc, curname, curseq.str());
					curname = "";
					curacc = "";
					curseq.str("");
				}

				// read next header
				std::string fhead = str.substr(1);
				name_split(fhead, curacc, curname);
			} else {
				curseq << str;
			}
		}
		// store last sequence
		insert_sequence(curacc, curname, curseq.str());
		//std::cout << "R " << curseq.str().substr(62750, 15) << std::endl;
		curname = "";
		curacc = "";
		curseq.str("");
		fastain.close();
	} else {
		std::cerr << "Error: cannot read: " << fasta_all_sequences << std::endl;
		exit(1);
	}

	std::cout << "loaded " << sequences.size() << " sequences" << std::endl;

	std::vector< std::multimap< size_t, size_t> > als_on_reference; // sequence index -> begin pos on that sequence -> alignment index
	// a new multimap for each ref sequence
	als_on_reference.resize(sequences.size());
	for(size_t i=0; i<sequences.size(); ++i) {
		als_on_reference.at(i) = std::multimap<size_t, size_t>();
	}



	std::ifstream mafin(maf_all_alignments.c_str());
	size_t skip_self = 0;
	size_t skip_alt = 0;
	if(mafin) {
		std::string str;
		while(getline(mafin, str)) {
			if(str.at(0)!='#') { // skip headers
				if(str.at(0)=='a') { // next alignment
					std::string aline1;
					std::string aline2;
					getline(mafin, aline1);
					getline(mafin, aline2);
					getline(mafin, str); // empty line at end of alignment
					if(0!=str.compare("")) {
						std::cerr << "Error: exspected empty line after alignment. Seen: " << str << std::endl;
						exit(1);
					}
					if(aline1.at(0)!='s' || aline2.at(0)!='s') {
						std::cerr << "Error: exspected maf sequence lines, found: " << std::endl << aline1 << std::endl<< aline2 << std::endl;
						exit(1);
					}

					std::vector<std::string> parts;
					strsep(aline1, " ", parts);
					if(parts.size()!=7) {
						std::cerr << "Error: exspected 7 fields in sequence line: " << aline1 << std::endl;
						exit(1);
					}
					std::string acc1;
					std::string name1;
					name_split(parts.at(1), acc1, name1);
					size_t start1 = atoi(parts.at(2).c_str());
					size_t size1 = atoi(parts.at(3).c_str());
					char strand1 = parts.at(4).at(0);
					size_t seqlength1 = atoi(parts.at(5).c_str());
					std::string al1 = parts.at(6);
					std::map<std::string, size_t>::iterator findseq1 = longname2seqidx.find(parts.at(1));
					if(findseq1==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence: " << parts.at(1) << std::endl;
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
						std::cerr << "Error: exspected 7 fields in sequence line: " << aline2 << std::endl;
						exit(1);
					}

					std::string acc2;
					std::string name2;
					name_split(parts.at(1), acc2, name2);
					size_t start2 = atoi(parts.at(2).c_str());
					size_t size2 = atoi(parts.at(3).c_str());
					size_t seqlength2 = atoi(parts.at(5).c_str());
					char strand2 = parts.at(4).at(0);
					std::string al2 = parts.at(6);
					std::map<std::string, size_t>::iterator findseq2 = longname2seqidx.find(parts.at(1));
					if(findseq2==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence: " << parts.at(1) << std::endl;
						exit(1);
					}
					size_t idx2 = findseq2->second;
					size_t incl_end2 = start2 + size2 - 1;
					if(strand2!='+') {
						size_t tmp = start2;
						start2 = seqlength2 - start2 -1;
						incl_end2 = seqlength2 - tmp - size2;
					}
					
					// both al parts not identical	
					if(idx1 != idx2 || !(start1==start2 && incl_end1 == incl_end2    )) {



						bool skip = false;

						// check if we already have an alignment with same coordinates
						// because we are looking for identical alignment it suffices to go over only one multimap
						std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> eqr =
						als_on_reference.at(idx1).equal_range(start1);
						for(std::multimap<size_t, size_t>::iterator it = eqr.first; it!=eqr.second; ++it) {
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

							als_on_reference.at(idx1).insert(std::make_pair(start1, alidx));
							als_on_reference.at(idx2).insert(std::make_pair(start2, alidx));


						}

					} else {
						skip_self++;
						// std::cerr << "Warning: Skip self alignment in seq " << idx1 << " from " << start1 << " to " << incl_end1 << std::endl;
					
					}

				} else {
					std::cerr << "Error: next alignment should start with 'a': " << str << std::endl;
					exit(1);
				}
			
			
			}	
		}

		mafin.close();

		if(dnastring::found_iupac_ambiguity) {
			std::cerr << "Warning: IUPAC DNA ambiguity characters were replaced by N" << std::endl;
		}

		std::cout << "Loaded: " << sequences.size() << " sequences and " << alignments.size() << " pairwise alignments " << std::endl;
		std::cout << skip_self << " self alignments were skipped" << std::endl;
		std::cout << skip_alt << " alternative alignments of identical regions were skipped" << std::endl;
	
	} else {
		std::cerr << "Error: cannot read: " << maf_all_alignments << std::endl;
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
	const std::vector<pw_alignment>& all_data::getAlignments()const{
		return alignments;
	}
	void all_data::add_pw_alignment(const pw_alignment& p){
		alignments.push_back(p);
	}
	const std::vector<size_t> & all_data::getAcc(size_t acc)const{
		std::string accName;
		for(std::map<std::string,size_t>::const_iterator it = accession_name.begin();it != accession_name.end(); it++){
			if(it->second == acc){
				accName = it->first;
			}else continue;
		}
		return acc_sequences.at(accName);		
	}
	size_t all_data::accNumber(size_t sequence_id) const {
		return sequence_to_accession.at(sequence_id);
	}


//const std::multimap<size_t, size_t> & all_data::getAlOnRefMap(size_t seq_idx) const {
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
	const std::map< std::string, size_t>& all_data::getLongname2seqidx()const{
		return longname2seqidx;
	}

	int all_data::findIdAlignment(std::string name,int start, int end){ // declare in graph or data ?
		for(std::vector<pw_alignment>::iterator it = alignments.begin(); it != alignments.end(); ++it){
			size_t idReference = it->getreference1();
			std::string nameReference = get_seq_name(idReference);
			//std::cout << name << " " << start << " " << end << " data " << nameReference << " " << it->getbegin2() << " " <<  it->getend2() << std::endl;
			if(start == it->getbegin2() && end == it->getend2() && name == nameReference ){
				return std::distance(alignments.begin(), it);
			}

		}
		std::cerr << "Alignment not found in findIdAlignment " << std::endl;
		return 0;
	}

	void all_data::set_accession(const std::string & acc){
		accession_name.insert(std::make_pair(acc, accession_name.size()));
		std::cout<< "acc "<< acc <<" accession name size: " << accession_name.size() <<std::endl;
	}
	void all_data::insert_sequence(const std::string & acc, const std::string & seq_name, const std::string & dna) {
	size_t seq_id = sequences.size();
	sequences.push_back(dnastring(dna));
	sequence_names.push_back(seq_name);
	std::map<std::string, std::vector<size_t> >::iterator findacc = acc_sequences.find(acc);
	size_t acc_id;
	if(findacc == acc_sequences.end()) {
		acc_sequences.insert(std::make_pair(acc, std::vector<size_t>(0)));
		findacc = acc_sequences.find(acc);
		acc_id = accession_name.size();	
		accession_name.insert(std::make_pair(acc, acc_id));	
	} else {
	std::map<std::string, size_t>::iterator find_acc_id = accession_name.find(acc);
	acc_id = find_acc_id -> second;
	}
	findacc->second.push_back(seq_id);
	sequence_to_accession.push_back(acc_id);
	std::stringstream longname;
	longname << acc << ":" << seq_name;
	longname2seqidx.insert(std::make_pair(longname.str(), seq_id));
	}


/**
	Transform 
	"Col0:scaffold_0 comment"
	to 
	"Col0", "scaffold_0"
	
**/		
	void all_data::name_split(const std::string & longname, std::string & acc, std::string & name) {
	std::vector<std::string> wparts; // remove all after first space
	strsep(longname, " ", wparts);
	std::string wname = wparts.at(0);
	std::vector<std::string> fparts;
	strsep(wname, ":", fparts);
	if(fparts.size()< 2) {
		std::cerr << "Error: sequence " << longname << " does not contain an Accession name separated by a colon" << std::endl;
		exit(1);
	}
	acc = fparts.at(0);
	std::stringstream namewriter;
	for(size_t i=1; i<fparts.size();++i) {
		if(i!=1) namewriter << ":";
		namewriter << fparts.at(i);
	}
	name = namewriter.str();


}


bool all_data::alignment_fits_ref(const pw_alignment & alr) const {
	const pw_alignment * al = &alr;
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
//	std::cout << " directions " << al1fw << al2fw << std::endl;

	for(size_t col=0; col < al->alignment_length(); ++col) {
		char al1char='X';
		char al2char='X';
		al->alignment_col(col, al1char, al2char);	
		
	//	std::cout << "c " << col << " " << al1char << " " << al2char << std::endl;	
		if(al1char!='-') { 
			char rchar = ref1.at(al1at);
		//	std::cout << " al 1 at " << al1at;
			if(al1fw) {
				al1at++;
			} else {
				al1at--;
				rchar =  dnastring::complement(rchar);

			}
		//	std::cout << " is " << rchar << std::endl;
			if(rchar != al1char) {
			//	al->print();
				std::cerr << "Warning: alignment test failed. Alignment has " << al1char << " at " << col <<" where reference sequence has " << rchar << std::endl;
			//	print_ref(al);
				return false;
			}
		}
		if(al2char!='-') {
			char rchar = ref2.at(al2at);
		//	std::cout << " al 2 at " << al2at ;
			if(al2fw) {
				al2at++;
			} else {
				al2at--;
				rchar =  dnastring::complement(rchar);

			}
		//	std::cout << " is " << rchar << std::endl;
			if(rchar != al2char) {
			//	al->print();
				std::cerr << " at " << al2at <<std::endl;
				std::cerr << "Warning: alignment test failed. Alignment has " << al2char << " where reference sequence has " << rchar << std::endl;
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
					std::cerr << "Warning: alignment end1 wrong. Should be " << al1e << " but is " << al->getend1() << std::endl;
					//print_ref(al);
					return false;
				} 
			
			} else {
				size_t al1b = al1at+1;
				if(al1b!=al->getend1()) {
					print_ref(al);

					std::cerr << "Warning: rev alignment end wrong. Should be " << al->getend1() << " but is " << al1b << std::endl;

					return false;
				}
			}
		
			if(al2fw) {
			size_t al2e = al2at-1;
			if(al2e!=al->getend2()) {
				print_ref(al);
				std::cerr << "Warning: alignment end2 wrong. Should be " << al->getend2()<< " but is " << al2e << std::endl;
				return false;
			}
		 	} else {
				size_t al2b = al2at+1;
				if(al2b!=al->getend2()) {
			//	print_ref(al);
				std::cerr << "Warning: rev alignment end wrong. Should be " << al2b  << " but is " << al->getend2() << std::endl;
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
			std::cout << " direction1: forward "<< std::endl;
		else 	std::cout << " direction1: reverse "<< std::endl;
		if(al->getbegin2() < al->getend2())
			std::cout << " direction2: forward "<< std::endl;
		else 	std::cout << " direction2: reverse "<< std::endl;
		for(size_t col=0; col < al->alignment_length(); ++col) {
			char al1char='X';
			char al2char='X';
			al->alignment_col(col, al1char, al2char);	
		
			std::cout << "c " << col << " " << al1char << " " << al2char << std::endl;
				if(al1char!='-') { 
				char rchar = ref1.at(al1at);
				std::cout << " al 1 at " << al1at;
					if(al->getbegin1() < al->getend1()) {
						al1at++;
					} else {
						al1at--;
					}
				std::cout << " is " << rchar << std::endl;
				}
				if(al2char!='-') {
				char rchar = ref2.at(al2at);
				std::cout << " al 2 at " << al2at ;
					if(al->getbegin2() < al->getend2()) {
						al2at++;
					} else {
						al2at--;
					}
				std::cout << " is " << rchar << std::endl;
				}

		}	
		al -> print();
	}


const std::string all_data::get_acc(size_t acc)const{	
	std::string Accession;
	for(std::map<std::string, size_t>::const_iterator it = accession_name.begin() ; it != accession_name.end();it++){
		if(it->second == acc){
			Accession = it -> first;
		}else continue;
	}
	return Accession;
}

const size_t all_data::get_acc_id(std::string acc)const{
	std::map<std::string, size_t>::const_iterator it = accession_name.find(acc);
	return it->second;
}


std::string all_data::get_seq_name(size_t s) const {
	return sequence_names.at(s);
}

size_t all_data::get_seq_size(size_t s) const {
	return sequences.at(s).length();
}
void all_data::compare_seq_with_decoding(std::ifstream & in){
	char c;
	for(size_t i = 0 ; i < numSequences();i++){
		unsigned int n;
		in >> n;
		assert( n == i);
		dnastring sequence = sequences.at(i);
		std::cout<< "length: "<<sequence.length()<<std::endl;
		for(size_t j =0 ; j < sequence.length(); j++){
			c = in.get();
			std::cout<< "seq "<< i << " at " << j << " is " << sequence.at(j) << " and c is " << c <<std::endl;
			assert(sequence.at(j)== c);
		}
	}
}

overlap::overlap(const all_data & d): data(d), als_on_reference(d.numSequences()){

	}

overlap::overlap(const overlap & o): data(o.data), als_on_reference(o.data.numSequences()){}


	
overlap::~overlap(){
}

void overlap::remove_alignment(const pw_alignment & removeR){
	const pw_alignment * remove = &removeR;

	std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqrb1 =
	als_on_reference.at(remove->getreference1()).equal_range(remove->getbegin1());
	for(std::multimap<size_t, pw_alignment>::iterator it = eqrb1.first; it!=eqrb1.second; ++it) {
		if ( removeR.equals(it->second) ) {
			als_on_reference.at(remove->getreference1()).erase(it);
			break;
		}		
	}
	std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqre1 =
	als_on_reference.at(remove->getreference1()).equal_range(remove->getend1());
	for(std::multimap<size_t, pw_alignment>::iterator it = eqre1.first; it!=eqre1.second; ++it) {
		if ( removeR.equals(it->second) ) {
			als_on_reference.at(remove->getreference1()).erase(it);
			break;
		}		
	}
	std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqrb2 =
	als_on_reference.at(remove->getreference2()).equal_range(remove->getbegin2());
	for(std::multimap<size_t, pw_alignment>::iterator it = eqrb2.first; it!=eqrb2.second; ++it) {
		if ( removeR.equals(it->second) ) {
			als_on_reference.at(remove->getreference2()).erase(it);
			break;
		}		
	}
	std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqre2 =
	als_on_reference.at(remove->getreference2()).equal_range(remove->getend2());
	for(std::multimap<size_t, pw_alignment>::iterator it = eqre2.first; it!=eqre2.second; ++it) {
		if ( removeR.equals(it->second) ) {
			als_on_reference.at(remove->getreference2()).erase(it);
			break;
		}		
	}
	std::set<pw_alignment, compare_pw_alignment>::iterator findr = alignments.find(removeR);
	//	remove->print();
//	assert(findr!=alignments.end()); // TODO check again
	alignments.erase(removeR);

}


void overlap::test_all() const {

	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment & al = *it;
		bool alok = data.alignment_fits_ref(al);
		if(!alok) {
			std::cout << " test all failed " << std::endl;
			al.print();
			exit(1);
		}
		
	
	}
}




void overlap::print_all_alignment() const {
	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it ) {
		const pw_alignment * al = &(*it);
		//std::cout << " x " << al << std::endl;
		al->print();
		std::cout << std::endl;
	}

}

void overlap::test_multimaps() {

	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it ) {
		const pw_alignment * remove = &(*it);
		std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqrb1 =
		als_on_reference.at(remove->getreference1()).equal_range(remove->getbegin1());
		size_t found =0;
		for(std::multimap<size_t, pw_alignment>::iterator it = eqrb1.first; it!=eqrb1.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}
		std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqre1 =
		als_on_reference.at(remove->getreference1()).equal_range(remove->getend1());
		for(std::multimap<size_t, pw_alignment>::iterator it = eqre1.first; it!=eqre1.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}

		std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqrb2 =
		als_on_reference.at(remove->getreference2()).equal_range(remove->getbegin2());
		for(std::multimap<size_t, pw_alignment>::iterator it = eqrb2.first; it!=eqrb2.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}

		std::pair< std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > eqre2 =
		als_on_reference.at(remove->getreference2()).equal_range(remove->getend2());
		for(std::multimap<size_t, pw_alignment>::iterator it = eqre2.first; it!=eqre2.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}
		assert(found==4);

	}
	for(size_t s=0; s<data.numSequences(); s++) {
		for(std::multimap< size_t, pw_alignment>::const_iterator it = als_on_reference.at(s).begin(); it!=als_on_reference.at(s).end(); ++it) {
			const pw_alignment & p = it->second;
	
			//std::cout << p << std::endl;
			std::set<pw_alignment, compare_pw_alignment>::iterator findp = alignments.find(p);
			if(findp==alignments.end()){
				std::cout<<"wrong one: "<<std::endl;
				p.print();
			}
			assert(findp!=alignments.end());
		}

	}
		
}



void overlap::insert_without_partial_overlap(const pw_alignment & p){
	p.print();
	std::cout<<" "<<std::endl;
	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = alignments.begin(); it != alignments.end(); it++){
		const pw_alignment & al = *it;
		al.print();
	}
	std::cout << " insert " << std::endl;
//	assert(alignments.find(p) == alignments.end());
	std::pair<std::set<pw_alignment, compare_pw_alignment>::iterator, bool > npp = alignments.insert(p);
	std::set<pw_alignment, compare_pw_alignment>::iterator npi = npp.first;
	const pw_alignment & np = *(npi);
 	
	std::multimap<size_t , pw_alignment> & alignment_on_reference1 = als_on_reference.at(np.getreference1());
	std::multimap<size_t , pw_alignment> & alignment_on_reference2 = als_on_reference.at(np.getreference2());

	std::pair<size_t, pw_alignment> begin2(np.getbegin2(),np);
	alignment_on_reference2.insert(begin2);
	std::pair<size_t, pw_alignment> end2(np.getend2(),np);
	alignment_on_reference2.insert(end2);

	std::pair<size_t, pw_alignment > begin1(np.getbegin1(),np);
	alignment_on_reference1.insert(begin1);
	std::pair<size_t, pw_alignment > end1(np.getend1(),np);
	alignment_on_reference1.insert(end1);

}

const pw_alignment * overlap::get_al_at_left_end(size_t ref1, size_t ref2, size_t left1, size_t left2) const {
		const std::multimap<size_t, pw_alignment> & r1map = als_on_reference.at(ref1);
		std::pair<std::multimap<size_t, pw_alignment >::const_iterator,std::multimap<size_t, pw_alignment>::const_iterator> eqr = r1map.equal_range(left1);
		for( std::multimap<size_t, pw_alignment>::const_iterator it = eqr.first; it!= eqr.second; ++it) {
			const pw_alignment & calr = it->second;
			const pw_alignment * cal = & calr;
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

std::multimap<size_t, pw_alignment> &  overlap::get_als_on_reference(size_t sequence)  {
	return als_on_reference.at(sequence);
}
const std::multimap<size_t, pw_alignment> &  overlap::get_als_on_reference_const(size_t sequence) const {
	return als_on_reference.at(sequence);
}



// TODO why not used
void overlap::test_all_part()const{
	
	/*std::set < const pw_alignment*, compare_pw_alignment> alignment_parts;
		
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
				std::cout << "Small alignment part is to long" << std::endl;
				exit(1);
			}
			}
			 else { 
			std::cout<< "There is no alignment!"<<std::endl;
			exit(1);
			}


		} while(cal_left1 <= al_right1);
		
	}
		if(alignment_parts.size()!=alignments.size()) {
		std::cout<< " Small parts don't cover the original one!"<<std::endl;
		exit(1);
		}*/
}

void overlap::test_overlap()const{
	for(std::set<pw_alignment, compare_pw_alignment>::iterator it1 = alignments.begin(); it1!=alignments.end(); ++it1){	
		const pw_alignment & alr = *it1;
		const pw_alignment * al1 = &alr;
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
		for(std::set<pw_alignment, compare_pw_alignment>::iterator it2 = alignments.begin(); it2!=alignments.end(); ++it2){
			const pw_alignment & alr2 = *it2;
			const pw_alignment * al2 = &alr2;
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
				std::cout<< "There are some overlap!"<<std::endl;
				exit(1);
			}
			else if((l1ref2 < l2ref2 && l2ref2 < r1ref2) || (l1ref2 < l2ref2 && r2ref2 < r2ref2)){
				std::cout<< "There are some overlap!"<<std::endl;
				exit(1);
			}
			else{
				std::cout<< "There is no overlap" << std::endl;
				exit(0);
			}

		}
	}	
}
	
bool overlap::checkAlignments(const pw_alignment & p)const{
	std::set<pw_alignment, compare_pw_alignment>::iterator it=alignments.find(p);
	if(it !=alignments.end()) return true;	
	else return false;	
}
	



size_t overlap::size() const {

//	for(size_t i=0; i<als_on_reference.size(); i++) {
//		std::cout << " ref " << i << " al start/end points: " << als_on_reference.at(i).size() << std::endl;
//	}
	return alignments.size();
}
		
const std::set<pw_alignment, compare_pw_alignment> & overlap::get_all() const {
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

void overlap::test_partial_overlap_set(std::set< const pw_alignment *, compare_pw_alignment> & als) {
	
	std::vector<const pw_alignment *> all;
	for(std::set<const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
		const pw_alignment * al = *it;
		all.push_back(al);
	}
	overlap::test_partial_overlap_vec(all);
}

// TODO ???
void overlap::test_partial_overlap_vec(std::vector< const pw_alignment *> & all) {
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
					std::cout << "partial overlap error 1: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}
			
			if(aref1 == bref2) {
				if(check_po(al1, ar1, bl2, br2)) {
					std::cout << "partial overlap error 2: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}

			if(aref2 == bref1) {
				if(check_po(al2, ar2, bl1, br1)) {
					std::cout << "partial overlap error 3: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}

			if(aref2 == bref2) {
				if(check_po(al2, ar2, bl2, br2)) {
					std::cout << "partial overlap error 4: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}

		
		}
	}




}
void overlap::test_partial_overlap() const {

	std::vector<const pw_alignment *> all;
	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment * al = &(*it);
		all.push_back(al);
	}

	overlap::test_partial_overlap_vec(all);

}

#define SPLITPRINT 0

splitpoints::splitpoints(const pw_alignment & p, const overlap & o, const all_data & d):overl(o), newal(p),data(d), split_points(d.numSequences()) {

#if SPLITPRINT
	std::cout << "RUN SPLITPOINTS ON " << std::endl;
	o.print_all_alignment();

	std::cout << " split with " << std::endl;
	p.print();
	std::cout << " end " << std::endl;

#endif

}

splitpoints::~splitpoints() {}

void splitpoints::find_initial_split_points_nonrecursive(size_t sequence, size_t left, size_t right) {
		const std::multimap<size_t , pw_alignment > & alignments_on_reference = overl.get_als_on_reference_const(sequence);

#if SPLITPRINT		
		std::cout << " seach for initial split points on " << sequence << " from " << left << std::endl;
		size_t count = 0;
#endif

		for( std::multimap<size_t, pw_alignment>::const_iterator it=alignments_on_reference.lower_bound(left);it!=alignments_on_reference.end(); ++it){
			const pw_alignment & alr = it->second;
			const pw_alignment * al = &alr;
			
#if SPLITPRINT
			std::cout << count++ <<" See " << std::endl;
			alr.print();
			std::cout << std::endl;
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
	//				std::cout<<"here1"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
					insert_split_point_nonrecursive(sequence, left);
				}
				if(alleft <= right && right < alright) {
	//				std::cout<<"here2"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
					insert_split_point_nonrecursive(sequence, right+1);

				}
				if(left < alleft && alleft < right) {
	//				std::cout<<"here3"<<std::endl;
	//				std::cout<<"al: "<<std::endl;
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
					}*/

				//	al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();
					insert_split_point_nonrecursive(sequence, alleft);

				}
				if(left < alright && alright < right) {
		//			std::cout<<"here4"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
		//			al->print();
				//	std::cout<<"newal: "<<std::endl;
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
				//	std::cout<<"here5"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();

					insert_split_point_nonrecursive(sequence, left);
				
				}
				if(alleft <= right && right < alright) {
		//			std::cout<<"here6"<<std::endl;
				//	al->print();
					insert_split_point_nonrecursive(sequence, right+1);
				}
				if(left < alleft && alleft < right) {
		//			std::cout<<"here7"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();
					insert_split_point_nonrecursive(sequence, alleft);
				}
				if(left < alright && alright < right) {
		//			std::cout<<"here8"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
		//			al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();
					insert_split_point_nonrecursive(sequence, alright+1);
				}

			
			}




			if(al_leftmost_leftbound > right)  {
#if SPLITPRINT
				std::cout << "Break after See " << count << " al rightmost leftbound is " << al_leftmost_leftbound<<  std::endl;	
#endif
				break;
			}
		}


}

void splitpoints::find_initial_split_points(size_t sequence, size_t left, size_t right) {
		const std::multimap<size_t , pw_alignment > & alignments_on_reference = overl.get_als_on_reference_const(sequence);

#if SPLITPRINT		
		std::cout << " seach for initial split points on " << sequence << " from " << left << std::endl;
		size_t count = 0;
#endif

		for( std::multimap<size_t, pw_alignment >::const_iterator it=alignments_on_reference.lower_bound(left);it!=alignments_on_reference.end(); ++it){
			const pw_alignment & alr = it->second;
			const pw_alignment * al = &alr;
			
#if SPLITPRINT
			std::cout << count++ <<" See " << std::endl;
			al->print();
			std::cout << std::endl;
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
	//				std::cout<<"here1"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
					insert_split_point(sequence, left);
				}
				if(alleft <= right && right < alright) {
	//				std::cout<<"here2"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
					insert_split_point(sequence, right+1);

				}
				if(left < alleft && alleft < right) {
	//				std::cout<<"here3"<<std::endl;
	//				std::cout<<"al: "<<std::endl;
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
					}*/

				//	al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();
					insert_split_point(sequence, alleft);

				}
				if(left < alright && alright < right) {
		//			std::cout<<"here4"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
		//			al->print();
				//	std::cout<<"newal: "<<std::endl;
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
				//	std::cout<<"here5"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();

					insert_split_point(sequence, left);
				
				}
				if(alleft <= right && right < alright) {
		//			std::cout<<"here6"<<std::endl;
				//	al->print();
					insert_split_point(sequence, right+1);
				}
				if(left < alleft && alleft < right) {
		//			std::cout<<"here7"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
				//	al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();
					insert_split_point(sequence, alleft);
				}
				if(left < alright && alright < right) {
		//			std::cout<<"here8"<<std::endl;
				//	std::cout<<"al: "<<std::endl;
		//			al->print();
				//	std::cout<<"newal: "<<std::endl;
				//	newal.print();
					insert_split_point(sequence, alright+1);
				}

			
			}




			if(al_leftmost_leftbound > right)  {
#if SPLITPRINT
				std::cout << "Break after See " << count << " al rightmost leftbound is " << al_leftmost_leftbound<<  std::endl;	
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
		std::cout << " Check ref 1 overlaps " << std::endl;
#endif

		find_initial_split_points(newal.getreference1(), left1, right1);

#if SPLITPRINT
		std::cout << " Check ref 2 overlaps " << std::endl;
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
		std::cout << " Check ref 1 overlaps " << std::endl;
#endif

		find_initial_split_points_nonrecursive(newal.getreference1(), left1, right1);

#if SPLITPRINT
		std::cout << " Check ref 2 overlaps " << std::endl;
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
	std::pair<std::set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
	if(insertes.second) {
#if SPLITPRINT
		std::cout << " new split point " << sequence << " at " << position << std::endl;
#endif	
	}
}




void splitpoints::insert_split_point(size_t sequence, size_t position) {
	std::pair<std::set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
	if(insertes.second) {
#if SPLITPRINT
		std::cout << " new split point " << sequence << " at " << position << std::endl;
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
				bool fgaps = onlyGapSample(fp);
				bool sgaps = onlyGapSample(sp);
					
				// find split part alignment ends on reference after removing end gaps
				size_t fpeleft2= (size_t) -1;
				size_t fperight2 = (size_t) -1;
				size_t speleft2 = (size_t) -1;
				size_t speright2 = (size_t) -1;
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
					fpe.get_lr2(fpeleft2, fperight2);
#if SPLITPRINT
					std::cout << "newal fpe " << std::endl;
					fpe.print();
					std::cout  << std::endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
					spe.get_lr2(speleft2, speright2);
#if SPLITPRINT
					std::cout << "newal spe " << std::endl;
					spe.print();
					std::cout << std::endl;
#endif
				}


				if(newal.getbegin2() < newal.getend2()) {
			//			std::cout << " try ins " << newal.getreference2() << " : " << spleft2 << std::endl;
					if(!sgaps) {
						insert_split_point(newal.getreference2(), speleft2);
					}
					// additional split point if gaps at split point
					if(!fgaps) {
						insert_split_point(newal.getreference2(), fperight2+1);
					}
				} else {
			//			std::cout << " try ins " << newal.getreference2() << " : " << fpleft2 << std::endl;
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
			//	std::cout<<"HERE!!"<<std::endl;
				pw_alignment fp;
				pw_alignment sp;
				newal.split(false, position, fp, sp);
				pw_alignment fpe;
				pw_alignment spe;
				
				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(fp);
				bool sgaps = onlyGapSample(sp);
					
				// find split part alignment ends on reference after removing end gaps
				size_t fpeleft1 = (size_t) -1;
				size_t fperight1 = (size_t) -1;
				size_t speleft1 = (size_t) -1;
				size_t speright1 = (size_t) -1;
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
					fpe.get_lr1(fpeleft1, fperight1);
#if SPLITPRINT
					std::cout << "newal fpe " << std::endl;
					fpe.print();
					std::cout  << std::endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
					spe.get_lr1(speleft1, speright1);
#if SPLITPRINT
					std::cout << "newal spe " << std::endl;
					spe.print();
					std::cout << std::endl;
#endif
				}

				
				if(newal.getbegin1() < newal.getend1()) {	
		//			std::cout << " try ins " << newal.getreference1() << " : " << spleft1 << std::endl;
					if(!sgaps) {
						insert_split_point(newal.getreference1(), speleft1);
					}
					if(!fgaps) {
						insert_split_point(newal.getreference1(), fperight1+1);
					}
				} else {
			//			std::cout << " try ins " << newal.getreference1() << " : " << fpleft1 << std::endl;
					if(!fgaps) {
						insert_split_point(newal.getreference1(), fpeleft1);
					}
					if(!sgaps) {
						insert_split_point(newal.getreference1(), speright1+1);
					}
				}
			}
		}





		const std::multimap<size_t, pw_alignment > & alonref = overl.get_als_on_reference_const(sequence);
		for(std::multimap<size_t, pw_alignment >::const_iterator it = alonref.lower_bound(position); it!=alonref.end(); ++it) {
			const pw_alignment & alr = it-> second;
			const pw_alignment * al = &alr;
				
				
	//			std::cout << " in ins " << sequence << " : " << position << " see: " << std::endl;
	//			al->print();
	//			std::cout << std::endl;


			if(al->getreference1()!=al->getreference2()) {
				size_t alleft;
				size_t alright;
				al->get_lr_on_reference(sequence, alleft, alright);
				if(alright < position) {
			//			std::cout << " break " << std::endl;
					break;
				}
				if(alleft == position) {
			//			std::cout << " break " << std::endl;
					break;
				} 
				if(alright+1== position) {
			//			std::cout << " break " << std::endl;
					break;
				}
				if(alleft>position){
			//			std::cout << " break " << std::endl;
					break;
				}
			//		al->print(); 
				//	std::cout<<"position"<<position<<std::endl;
				//	std::cout<<"alleft"<<alleft<<std::endl;
				//	std::cout<<"alright"<<alright<<std::endl;
				assert(alleft < position && position <= alright);
				pw_alignment fp;
				pw_alignment sp;
				pw_alignment fpe;
				pw_alignment spe;
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
					}*/
				//	std::cout<<"al length: "<< al->alignment_length() <<std::endl;
				//	if(al->alignment_length()>1){

				al->split_on_reference(sequence, position, fp, sp);
				


				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(fp);
				bool sgaps = onlyGapSample(sp);
					
				// find split part alignment ends on reference after removing end gaps
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
#if SPLITPRINT
					std::cout << "al fpe " << std::endl;
					fpe.print();
					std::cout  << std::endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
#if SPLITPRINT
					std::cout << "al spe " << std::endl;
					spe.print();
					std::cout << std::endl;
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




						//	std::cout<<"Heya!!!"<<std::endl;
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
		//				std::cout << "inone" << std::endl;
						pw_alignment fp;
						pw_alignment sp;
						al->split(true, position, fp, sp);
					pw_alignment fpe;
					pw_alignment spe;
					// do splitted parts have only gaps one of their reference sequences
					bool fgaps = onlyGapSample(fp);
					bool sgaps = onlyGapSample(sp);
					
					// find split part alignment ends on reference after removing end gaps
					size_t fpeleft2 = (size_t) -1;
					size_t fperight2 = (size_t) -1;
					size_t speleft2 = (size_t) -1;
					size_t speright2 = (size_t) -1;
					if(!fgaps) {
						fp.remove_end_gaps(fpe);
						fpe.get_lr2(fpeleft2, fperight2);
#if SPLITPRINT
						std::cout << "al fpe " << std::endl;
						fpe.print();
						std::cout  << std::endl;
#endif
					}
					if(!sgaps) {
						sp.remove_end_gaps(spe);
						spe.get_lr2(speleft2, speright2);
#if SPLITPRINT
						std::cout << "al spe " << std::endl;
						spe.print();
						std::cout << std::endl;
#endif
					}


						
							
						if(al->getbegin2() < al->getend2()) {
			//				std::cout << " spleft2 " << spleft2 << std::endl;
							if(!sgaps) {
								insert_split_point(sp.getreference2(), speleft2);
							}
							if(!fgaps) {
								insert_split_point(sp.getreference2(), fperight2+1);
							}
						} else {
			//				std::cout << " fpgetend2 " << fp.getend2() << std::endl;
							if(!fgaps) {
								insert_split_point(fp.getreference2(), fpe.getend2());
							}
							if(!sgaps) {
								insert_split_point(fp.getreference2(), spe.getbegin2()+1);
							}

						}
					}
											
					if(position > alleft2 && position <= alright2) {
			//			std::cout << "intwo" << std::endl;
						pw_alignment fp;
						pw_alignment sp;
						al->split(false, position, fp, sp);
					pw_alignment fpe;
					pw_alignment spe;
					
			// do splitted parts have only gaps one of their reference sequences
					bool fgaps = onlyGapSample(fp);
					bool sgaps = onlyGapSample(sp);
					
					// find split part alignment ends on reference after removing end gaps
					size_t fpeleft1 = (size_t) -1;
					size_t fperight1 = (size_t) -1;
					size_t speleft1 = (size_t) -1;
					size_t speright1 = (size_t) -1;
					if(!fgaps) {
						fp.remove_end_gaps(fpe);
						fpe.get_lr1(fpeleft1, fperight1);
#if SPLITPRINT
						std::cout << "al fpe " << std::endl;
						fpe.print();
						std::cout  << std::endl;
#endif
					}
					if(!sgaps) {
						sp.remove_end_gaps(spe);
						spe.get_lr1(speleft1, speright1);
#if SPLITPRINT
						std::cout << "al spe " << std::endl;
						spe.print();
						std::cout << std::endl;
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


	void splitpoints::split_all(std::set<pw_alignment, compare_pw_alignment> & remove_alignments, std::vector<pw_alignment> & insert_alignments ){

#if SPLITPRINT		
		std::cout << " All split points " << std::endl;
		for(size_t i=0; i<data.numSequences(); ++i) {
			std::cout << "ref " << i << " ";
			for(std::set<size_t>::iterator it = split_points.at(i).begin(); it!=split_points.at(i).end(); ++it) {
				std::cout << " " << *it;
			}
			std::cout << std::endl;
		}
#endif

		pw_alignment p1;
		pw_alignment p2;
		for(size_t i = 0; i <data.numSequences();i++){

#if SPLITPRINT			
			std::cout << "REFERENCE " << i << std::endl;
#endif
			const std::multimap<size_t, pw_alignment> & als_on_ref = overl.get_als_on_reference_const(i);
			for(std::set<size_t>::iterator split = split_points.at(i).begin(); split!= split_points.at(i).end(); ++split){
#if SPLITPRINT
				std::cout << "SPLIT " << *split << std::endl;
#endif
				std::cout<< "split all!" <<std::endl;
				for(std::multimap<size_t, pw_alignment>::const_iterator it =als_on_ref.lower_bound(*split); it!=als_on_ref.end(); ++it){
					const pw_alignment & pr = it-> second;
					const pw_alignment * p = & pr;
					p->print();
					std::cout << std::endl;

					if(p->getreference1() == i) {
						size_t pleft1;
						size_t pright1;
						p->get_lr1(pleft1,pright1);
						if( pleft1<*split && pright1>=*split){
							remove_alignments.insert(pr);
#if SPLITPRINT
							std::cout<<"remove alignment1: "<<std::endl;
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
							remove_alignments.insert(pr);
#if SPLITPRINT
							std::cout<<"remove alignment2: "<<std::endl;
									p->print();
#endif
						}
					}
				}		
			}
		}

	//	std::cout << " in split all "<< remove_alignments.size() << " remove_alignments " << std::endl;
		std::vector<pw_alignment>  split_pieces;
		splits(newal, split_pieces);
		std::cout << "new al "	<< std::endl;
		newal.print();
		for(std::set<pw_alignment,compare_pw_alignment>::iterator removed = remove_alignments.begin(); removed != remove_alignments.end(); ++removed){
			splits(*removed, split_pieces);			
		}
	//	std::cout<<"split_pieces: "<<std::endl;
	//	for(size_t i = 0 ; i < split_pieces.size();i++){split_pieces.at(i).print();}	
		std::cout<<"size of split pieces"<<split_pieces.size()<<std::endl;
		std::set<pw_alignment,compare_pw_alignment> inserted_pieces;
		for(size_t i = 0; i<split_pieces.size();i++) {
#ifndef NDEBUG
				if(!data.alignment_fits_ref(split_pieces.at(i))) {
				//	std::cout<<"fails here!"<<std::endl;
					exit(1);
				}
#endif
			if(!onlyGapSample(split_pieces.at(i))){	
				pw_alignment noendgaps;
				split_pieces.at(i).remove_end_gaps(noendgaps);
#if SPLITPRINT
				std::cout << "ENDGAPS" << std::endl;
				split_pieces.at(i).print();
				std::cout << "TO " << std::endl;
				noendgaps.print();
				std::cout << std::endl;
#endif
				std::set<pw_alignment,compare_pw_alignment>::const_iterator it = inserted_pieces.find(noendgaps);
				if (overl.checkAlignments(noendgaps)) {}
				else if ( it != inserted_pieces.end()){}
				else {
					insert_alignments.push_back(noendgaps);
					inserted_pieces.insert(noendgaps);
				}
			}
		}
		if(insert_alignments.size()==0 && remove_alignments.size()==0){
			insert_alignments.push_back(newal);
			std::cout << "new_al is pushed back" <<std::endl;
		}
		
					
	}
	void splitpoints::splits(const pw_alignment & pr,  std::vector<pw_alignment> & insert_alignments){
		pw_alignment p = pr;
#if SPLITPRINT	
		std::cout << "SPL" << std::endl;
		p.print();
		std::cout << std::endl;

		std::cout << " split on " << p.getreference1() << std::endl;
#endif

		pw_alignment p1;
		pw_alignment p2;
		size_t left1;
		size_t right1;	
//		size_t left2;
//		size_t right2;
		p.get_lr1(left1,right1);
//		p->get_lr2(left2,right2);
		for(std::set<size_t>::iterator splitp = split_points.at(p.getreference1()).upper_bound(left1); splitp!= split_points.at(p.getreference1()).end(); splitp++){
			if(right1>=*splitp){
#if SPLITPRINT
				std::cout << "sp " << *splitp << std::endl;
#endif
				p.split(true,*splitp,p1,p2);
				if(p.getbegin1() < p.getend1()) {
					p = p2;

#if SPLITPRINT
					p1.print();
#endif
					if(!onlyGapSample(p1) && p1.alignment_length() >1 &&p1.getbegin1()!=p1.getend1() && p1.getbegin2()!=p1.getend2()) {
						insert_alignments.push_back(p1);
					}
				} else if(p.getbegin1() > p.getend1()) {
					p = p1;

#if SPLITPRINT
					p2.print();
#endif

					if(!onlyGapSample(p2) && p2.alignment_length()>1 && p2.getbegin1()!=p2.getend1() && p2.getbegin2()!=p2.getend2() ){	
						insert_alignments.push_back(p2);
					}
					
				}
			}
			else break;
		}

#if SPLITPRINT
		std::cout << " last part" << std::endl;
		p.print();
		std::cout << std::endl;
#endif

		if(!onlyGapSample(p)&& p.alignment_length()>1 && p.getbegin1()!=p.getend1() && p.getbegin2()!=p.getend2() ){	
			insert_alignments.push_back(p);		
		}
	}

	bool splitpoints::onlyGapSample(const pw_alignment & p) const {
		bool gapOnSample1 = true;
		bool gapOnSample2 = true;
		for(size_t i = 0; i< p.alignment_length();i++){
			char p1char = 'X';
			char p2char = 'X';
			p.alignment_col(i, p1char, p2char);	
			if(p1char !='-' ) gapOnSample1 = false;
			if(p2char !='-' ) gapOnSample2 = false;			
			if(! gapOnSample1 && !gapOnSample2) {
				return false;
			}		
		}
		return gapOnSample1 || gapOnSample2;
	}

/*
	std::vector<pw_alignment>  splitpoints::get_insert()const{
		return insert_alignments;
	}
*/




	model::model(all_data & d): data(d),cost_on_acc (5, std::vector<double>(data.numAcc(),1)),modification(data.numAcc(), vector<vector<vector<double> > >(data.numAcc(),vector<vector<double> > (6, vector<double>(6,1) ))) {}

	model::~model(){}

	void model::acc_base_frequency(){
	//	std::vector<vector<double> > cost_on_acc (5, vector<double>(data.numAcc(),0));
		std::vector<size_t> total_base_number_on_acc (data.numAcc(),0);
		for(size_t k = 0; k < data.numSequences(); k++){
			std::vector<size_t> number(5,0);
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
			//	std::cout<< "cost of base "<< t <<" on acc "<< m << " is "<< cost_on_acc.at(t).at(m)<<std::endl;
			}
		}
	}
	void model::alignment_modification(){
	//	std::vector<vector<vector<vector<double> > > >modification(data.numAcc(), vector<vector<vector<double> > >(data.numAcc(),vector<vector<double> > (6, vector<double>(6,0) )));
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

					//	std::cout << "k " << k << " num " << modification.at(acc1).at(acc2).at(j).at(k) << std::endl;
						sum+=modification.at(acc1).at(acc2).at(j).at(k);
					}	
					for(size_t k=0; k<6; ++k) {
						modification.at(acc1).at(acc2).at(j).at(k)/=sum;
					//	std::cout << " simple model cost: acc " << acc1 << " to " << acc2 << " modify " << j << " to " << k << " costs " << -log2(modification.at(acc1).at(acc2).at(j).at(k)) << " bits " << std::endl; 
					
					}
				}

			}
		}
	}

	void model::cost_function( pw_alignment& p) const {
//		std::cout << " cf2 " << std::endl;
		std::vector<double> cost_on_sample(2);
		std::vector<double> modify_cost(2);
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

		std::vector<double> cost_on_sample(2,0);
		std::vector<double> modify_cost(2,0);
		std::vector<vector<size_t> >al_base_number(2,vector<size_t>(6,0));
		std::vector<vector<double> >sequence_cost(2,vector<double>(5,0));
		std::vector<vector<double> >create_cost(2,vector<double>(5,0));		
	//	std::vector<vector<vector<double> > > modification_cost(2,vector<vector<double> >(6,vector<double>(6,0)));
		std::vector<vector<double> > modification_number(6, vector<double>(6,0)); // in (i,j): store modifications from 1 to 2
	//	std::vector<double> cost_on_sample (2,0);
		for(size_t i = 0 ; i < 5; i++){
		 sequence_cost.at(0).at(i)= -log2(cost_on_acc.at(i).at(data.accNumber(p.getreference1())));
		 sequence_cost.at(1).at(i)= -log2(cost_on_acc.at(i).at(data.accNumber(p.getreference2())));
		}
//		std::cout << " alsample size " << p.getsample1().size() << std::endl;
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
//			std::cout<<"creat cost of "<< dnastring::index_to_base(j)<<"  on reference1 is "<<  (al_base_number.at(0).at(j))*(sequence_cost.at(0).at(j))<<std::endl;
			create_cost.at(1).at(j)= al_base_number.at(1).at(j)*sequence_cost.at(1).at(j);
			cost_on_sample.at(1)=cost_on_sample.at(1)+ create_cost.at(1).at(j);			
//			std::cout<<"creat cost of "<< dnastring::index_to_base(j)<<"  on reference2 is "<<  al_base_number.at(1).at(j)*sequence_cost.at(1).at(j)<<std::endl;
		}
		//	std::cout<<"creat cost of the alignment on sample 1: "<< cost_on_sample.at(0);
		//	std::cout<<"creat cost of the alignment on sample 2: "<< cost_on_sample.at(1);	
		
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
		std::cout << " create 2 " << c2 << " m1 " << m1 << std::endl; 
		std::cout << " create 1 " << c1 << " m2 " << m2 << std::endl; 


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
				the gain of using the alignment is: c2 - m1 (c1+c2-(c1+m1))

		*/
//		std::cout << "in gain: create 2 " << c2 << " m1 " << m1 << std::endl; 

		g1 = c2 - m1;
		g2 = c1 - m2;

	}
	void model::total_information(size_t & information){
		information = 0;
		size_t number_of_bases =0;
		for(size_t i =0; i < data.numAcc(); i++){
			for(size_t k = 0; k < data.getAcc(i).size() ; k++){
				const dnastring & sequence = data.getSequence(data.getAcc(i).at(k));
			/*	if(data.getAcc(i).at(k)==21){
					for(size_t j =212100; j < 212175; j++){
						size_t base = dnastring::base_to_index(sequence.at(j));
						std::cout << j << " " << base << std::endl;
					}
				}*/
				for(size_t j =0; j < sequence.length(); j++){
					number_of_bases ++;
					size_t base = dnastring::base_to_index(sequence.at(j));
					information += -log2(cost_on_acc.at(base).at(i));
				}
			}

		}
	//	std::cout<< "number of bases: "<< number_of_bases << std::endl;

	}
	void model::train() {
		acc_base_frequency();
		alignment_modification();
	}

mc_model::mc_model(all_data & d):data(d), sequence_successive_bases(d.numAcc()), create_cost(d.numAcc(),std::vector<double>(5,1)),mod_cost(d.numAcc(),vector<std::map<std::string, vector<double> > >(d.numAcc())),high(d.numOfAcc()),highValue(d.numAcc(),std::vector<std::map<std::string, std::vector<unsigned int> > >(d.numAcc())){
	size_t numberOfPowers = 32;
//	size_t numberOfPowers = NUM_DELETE;
	if(NUM_KEEP>numberOfPowers) {
		numberOfPowers = NUM_KEEP;
	}
	powersOfTwo = std::vector<size_t>(numberOfPowers, 1);
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
				std::stringstream context;		
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
				std::string seq;
				context>>seq;
				std::map <std::string, std::vector<double> >::iterator it= sequence_successive_bases.at(acc).find(seq);
				if(it==sequence_successive_bases.at(acc).end()) {
					sequence_successive_bases.at(acc).insert(std::make_pair(seq, std::vector<double>(5,1)));
					it= sequence_successive_bases.at(acc).find(seq);
				}
				it->second.at(s)++;	
			}
		}
		for(size_t i =0 ; i< data.numAcc(); i++){
			for(std::map <std::string, std::vector<double> >::iterator it= sequence_successive_bases.at(i).begin();it!=sequence_successive_bases.at(i).end();it++){
				std::string seq = it->first;
				int total = 0;
				std::vector<double> & base = sequence_successive_bases.at(i).at(seq);
				for(size_t j = 0; j<5;j++){
					total += base.at(j);
			//	std::cout<<"base: "<<base.at(j)<<std::endl;
				}
		//	std::cout<<"totalnumber of bases for each seq:  "<< total<<std::endl;
				for(size_t k=0; k<5;k++){
					base.at(k) = -log2(base.at(k)/total);
					create_cost.at(i).at(k) += base.at(k);
				//	std::cout<<"cr_cost: "<<base.at(k) <<std::endl;		
				}
			}
		}
			


	}

/*
	train alignment/modification model
*/
void mc_model::markov_chain_alignment(){
	
//void mc_model::markov_chain_alignment(std::ofstream& outs){
	// zero content counting functor
	counting_functor functor(data);
	// make all possible patterns in this class (all_alignment_patterns)
	make_all_alignments_patterns();

	// Set all counts to 1 for all context/accession pairs
	for(size_t i = 0; i< data.numAcc();i++){
		for(size_t j = 0; j < data.numAcc();j++){
			for(std::set<std::string>::iterator it= all_alignment_patterns.begin(); it != all_alignment_patterns.end() ; it++){
				std::string pattern = *it;	
				functor.create_context(i, j, pattern);			
			}
		}
	}


	// count all edit operations contained in the alignments (per context)
	for(size_t k = 0; k < data.numAlignments(); k++){
		const pw_alignment & p = data.getAlignment(k);
		computing_modification_oneToTwo(p,functor);
		computing_modification_twoToOne(p,functor);
	}
	functor.total_context();
	for(size_t i = 0; i< data.numAcc();i++){
		for(size_t j = 0; j < data.numAcc();j++){



			for(std::map <std::string, std::vector<size_t> >::const_iterator it= functor.get_context(i,j).begin();it!=functor.get_context(i,j).end();it++){
				std::string seq1 = it->first;
				const std::vector<size_t> & base = functor.get_context(i,j).at(seq1);
			//	std::cout<<"base is: "<<std::endl;
			//	for(size_t a= 0; a< base.size();a++){
				//		std::cout<< "base at "<< a<< " which is  " << print_modification_character(a)<<" is "<<base.at(a)<<std::endl;
			//	}
			//	std::cout<< "context is: "<<std::endl;
				//	for(size_t m = 0 ; m < seq1.size() ; m++){
				//		std::cout<< int(seq1.at(m))<<std::endl;
				//	}
				//	std::cout<<"the total number of happening the above context between "<<i<<" and "<<j<<" is "<< functor.get_total(i,j,seq1) <<std::endl;
				for(size_t k = 0; k< (NUM_DELETE+NUM_KEEP+10);k++) {
				//		std::cout<<"The number of happening "<< print_modification_character(k) << " between acc " <<i<< " and acc " << j << " after a certain context is "<< base.at(k)<<std::endl;
				//		std::cout<<"In MC model, the cost value of" <<print_modification_character(k) <<  " after above pattern between acc " << i << " and acc " << j << " is " << -log2(base.at(k)/functor.get_total(i,j,seq1)) << " bits " << std::endl; 
					std::map <std::string, std::vector<double> >::iterator it1= mod_cost.at(i).at(j).find(seq1);
					if(it1==mod_cost.at(i).at(j).end()) {
						mod_cost.at(i).at(j).insert(std::make_pair(seq1, std::vector<double>((NUM_DELETE+NUM_KEEP+10),1)));
						it1= mod_cost.at(i).at(j).find(seq1);
					}
					it1->second.at(k)=-log2(base.at(k)/functor.get_total(i,j,seq1));
				}
			}
			
			// Compute low and high values TODO 			
			for(std::set<std::string>::iterator it= all_alignment_patterns.begin(); it != all_alignment_patterns.end() ; it++){	
				std::vector<double> num(NUM_DELETE+NUM_KEEP+10,0);
				std::vector<unsigned int> low(NUM_DELETE+NUM_KEEP+10,0);
				std::vector<unsigned int> high_value(NUM_DELETE+NUM_KEEP+10,0);
				unsigned int l = 0;
			//	unsigned int total = 0;
				size_t bit = 12; // number of bits to use for encoding event width
				std::string current_pattern = *it;
				// it1: high values for current pattern/accession pair
                               	highValue.at(i).at(j).insert(std::make_pair(current_pattern,std::vector<unsigned int>(NUM_DELETE+NUM_KEEP+10,0)));
				std::map<std::string, std::vector<unsigned int> >::iterator it1=highValue.at(i).at(j).find(current_pattern);
				assert(it1 != highValue.at(i).at(j).end());
				// it3: get counts for current pattern/accession pair
				std::map <std::string, std::vector<size_t> >::const_iterator it3= functor.get_context(i,j).find(current_pattern);
				double total =  functor.get_total(i,j,current_pattern);

				for (size_t f=0; f < NUM_DELETE+NUM_KEEP+10;f++){
					low.at(f) = l;
					num.at(f)=it3->second.at(f);
					size_t rescaledNum = (num.at(f)/total)*(powersOfTwo.at(bit) - NUM_DELETE - NUM_KEEP - 11) + 1;
					assert(rescaledNum >= 1);
					assert(rescaledNum < powersOfTwo.at(bit));
				//	std::cout << "rescled num: "<< rescaledNum << "num: " << num.at(j) << std::endl;
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
						std::cout<<"high value: "<< it1 ->second.at(k)<<std::endl;
					}*/
			}

		}			
	}
}


const std::map<std::string, std::vector<unsigned int> > & mc_model::get_highValue(size_t acc1, size_t acc2)const{
	return highValue.at(acc1).at(acc2);
}



void mc_model::cost_function( pw_alignment& p) const {
	std::vector<double> cost_on_sample(2);
	std::vector<double> modify_cost(2);
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

const std::vector<vector< std::map<std::string, vector<double> > > > &mc_model::get_mod_cost()const{
	return mod_cost;
}


void mc_model::cost_function(const pw_alignment& p, double & c1, double & c2, double & m1, double & m2)const {
	if(p.is_cost_cached()) {
//		std::cout << " from cache " << &p << std::endl;
		c1 = p.get_create1();
		c2 = p.get_create2();
		m1 = p.get_modify1();
		m2 = p.get_modify2();
//			std::cout << " c2 " << c2 << " m1 " << m1 << std::endl; 
//		std::cout << " c1 " << c1 << " m2 " << m2 << std::endl; 

		return; 
	}



//	p.print();
//	std::cout<<"data address in cost function: "<< &data <<std::endl;
//	data.numAcc();
	cost_functor f(data,mod_cost);
//	p.print();
//	size_t length = p.alignment_length();
//	std::cout<<"length: "<< length<<std::endl;
	computing_modification_oneToTwo(p,f);
	computing_modification_twoToOne(p,f);
	std::vector<double> cost_on_sample(2,0);
	std::vector<double> modify_cost(2,0);
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
//	std::cout<<"left: "<<left1<<"right: "<<right1<<std::endl;
	for(size_t i = left1; i< right1; i++){
		s1chr = data.getSequence(p.getreference1()).at(i);
		size_t s1 = dnastring::base_to_index(s1chr);
		std::stringstream context1;
		for (size_t j = Sequence_level; j>0; j--){
			
			if(i<j){
				char r1chr = 'A';
				context1 << r1chr;					
			}else{
				char r1chr = data.getSequence(p.getreference1()).at(i-j);
				context1 << r1chr;
			//	std::cout<<"rchr: "<<r1chr<<std::endl;	
			}
		}
		std::string seq1;
		context1>>seq1;
	//	std::cout<<"seq1: " << seq1<<std::endl;
		std::map <std::string, std::vector<double> >::const_iterator it= sequence_successive_bases.at(acc1).find(seq1);
		assert(it != sequence_successive_bases.at(acc1).end());
	//	std::cout<<"sequence successive at "<< s1 << " is "<<it->second.at(s1)<<std::endl;
		cost_on_sample.at(0) += it->second.at(s1);
	}	
//	std::cout<<"cost: "<<cost_on_sample.at(0)<<std::endl;	
	for(size_t i = left2; i<right2; i++){
		s2chr = data.getSequence(p.getreference2()).at(i);
		size_t s2 = dnastring::base_to_index(s2chr);
		std::stringstream context2;
		for(size_t j = Sequence_level; j>0; j--){
			if(i<j){
				char r2chr = 'A';
				context2 << r2chr;	
			}else{
				char r2chr = data.getSequence(p.getreference2()).at(i-j);
				context2 << r2chr;
			}
		}
		std::string seq2;
		context2>>seq2;
		std::map <std::string, std::vector<double> >::const_iterator it1= sequence_successive_bases.at(acc2).find(seq2);
		assert(it1 != sequence_successive_bases.at(acc2).end());
	//	std::cout<<"sequence successive at "<< s2 << " is "<<it1->second.at(s2)<<std::endl;
		cost_on_sample.at(1) += it1->second.at(s2);
	}



	c1 = cost_on_sample.at(0);
	c2 = cost_on_sample.at(1);
//		m1 = modify_cost.at(0);
//		m2 = modify_cost.at(1);
	m1 = f.get_modify(p,acc1,acc2);
	m2 = f.get_modify(p,acc2,acc1);
	modify_cost.at(0) = m1;
	modify_cost.at(1) = m2;
	p.set_cost(cost_on_sample, modify_cost);
//		std::cout << " c2 " << c2 << " m1 " << m1 << std::endl; 
//		std::cout << " c1 " << c1 << " m2 " << m2 << std::endl; 
	assert(p.is_cost_cached());
	assert(m1==p.get_modify1());
	assert(m2==p.get_modify2());
	assert(c1==p.get_create1());
	assert(c2==p.get_create2());

	//	std::cout<< "length: " << length<<std::endl;
	//	std::cout<< "c1: " << c1 << " c2: "<< c2 << " m1: "<< m1<< " m2: "<< m2 <<std::endl;
}


void mc_model::gain_function(const pw_alignment& p, double & g1, double & g2)const {
	double c1;
	double c2;
	double m1;
	double m2;
	cost_function(p, c1, c2, m1,m2);

	
	g1 = c2 - m1;
	g2 = c1 - m2;

//		std::cout << " gain function c2 " << c2 << " m1 " << m1 << " gain1 " << g1 << std::endl; 
//		std::cout << " gain function c1 " << c1 << " m2 " << m2 << " gain2 " << g2 << std::endl; 


}

void mc_model::write_parameters(std::ofstream & outs){
	for(size_t i = 0 ; i < data.numAcc(); i++){
		outs << data.get_acc(i);
		std::cout<<"acc: "<<data.get_acc(i)<<std::endl;
		outs<< (char)0;
	}
	outs<<(char)8;
	make_all_the_patterns();
	for(size_t i = 0 ; i < data.numAcc(); i++){
		outs<< (char)0;
		outs << i;
		outs << (char)0;
		std::cout<<"acc: "<<i<<std::endl;
		for(std::map<std::string, std::vector<double> >::iterator it= all_the_patterns.begin(); it != all_the_patterns.end() ; it++){
			std::string pattern = it ->first;
			std::cout << "pattern: " << pattern << std::endl;
			std::map<std::string,std::vector<double> >::iterator it1=sequence_successive_bases.at(i).find(pattern);
			if(it1 != sequence_successive_bases.at(i).end()){
				for(size_t n = 0; n <5; n ++){
					it->second.at(n) = it1 ->second.at(n);
				}
			}else{
				for(size_t n =0; n < 5 ; n++){
					it -> second.at(n) = -log2(0.2);
				}
			}
		}
			std::map <std::string, std::vector<unsigned int> >lower_bound;
			for(std::map<std::string, std::vector<double> >::iterator it= all_the_patterns.begin(); it != all_the_patterns.end() ; it++){	
					std::vector<double> num(5,0);
					std::vector<bool> bit_to_byte(0);
					std::vector<unsigned int> low(5,0);
					std::vector<unsigned int> high_value(5,0);
					unsigned int l = 0;
			//		unsigned int total = 0;
					size_t bit = 12;
					std::string current_pattern= it ->first;
                                 	high.at(i).insert(std::make_pair(current_pattern,std::vector<unsigned int>(5,0)));
					std::map<std::string, std::vector<unsigned int> >::iterator it1=high.at(i).find(current_pattern);
					assert(it1 != high.at(i).end());
					for (size_t j=0; j < 5;j++){
						low.at(j) = l;
					//	std::cout<< " low: " << low.at(j)<<std::endl;
						num.at(j)=it->second.at(j);
						//std::cout<< "num at " << j << " is " << exp((-num.at(j))*log(2))<<std::endl;
						int power_of_two = exp((-num.at(j))*log(2))*powersOfTwo.at(bit);
						high_value.at(j) = l + power_of_two;
						l = high_value.at(j);
					//	std::cout<<"high: "<< high_value.at(j)<<std::endl;
					}
					for(size_t j = 0; j < 5 ; j++){
						if(high_value.at(j)==low.at(j)){
						//	std::cout << "low = high" <<std::endl;
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
					//	std::cout<< "low1: "<<low.at(j)<<std::endl;
					/*	if(it->first == "AC"){
							std::cout<< "high: "<<h<<std::endl;
						}*/
						for(size_t m = 0; m < bit; m++){
							bit_to_byte.push_back(h%2);
				//			std::cout<< "h%2: "<< h%2;
							h = h/2;
						}
				//		std::cout<< " "<<std::endl;
						assert(high_value.at(j)!=low.at(j));
					}
				//	total = it1->second.at(4);
				//	std::cout<< "total in stream: "<< it1->second.at(4)<<std::endl;
				//	std::cout<<bit_to_byte.size()<<std::endl;					
					for(size_t n =0; n < bit_to_byte.size()-8; n++){
						unsigned char a = 0;
						for(size_t m = n; m <n+8; m++){
							a+= powersOfTwo.at(m-n)* bit_to_byte.at(m);
						}
						n= n+7;
						outs<< a;
				//		std::cout<< "eight bits of high: "<<int(a); 
					}
				//	std::cout << " " << std::endl;
				}
			}
			outs<<(char)8;
		//}
	//	outs.close();
	}
	void mc_model::write_alignments_pattern(std::ofstream & outs){
	//	std::ofstream outs("encode",std::ofstream::binary);
	//	std::ofstream outs("encode",std::ofstream::binary|std::ofstream::app);
		size_t bit = 12;
	//	if(outs.is_open()){//Per accession!
			for(size_t i = 0 ; i < data.numAcc(); i++){
				for(size_t j =0; j < data.numAcc(); j++){
					outs<< (char) 0;
					outs<< i;
					outs<< (char) 0;
					outs<< j;
				//	outs << data.get_acc(i);
					outs<< (char) 0;
				//	outs<< data.get_acc(j);
				//	outs<< (char) 0;
			//		std::cout<<"acc1 : " << data.get_acc(i) << " acc2: " << data.get_acc(j) << std::endl;
					for(std::map<std::string, std::vector<unsigned int> >::iterator it= highValue.at(i).at(j).begin(); it != highValue.at(i).at(j).end(); it++){	
						std::vector<bool> bit_to_byte(0);
						for(size_t j = 0; j < NUM_DELETE+NUM_KEEP+10; j++){
							int h =it->second.at(j);
							for(size_t m = 0; m < bit; m++){
								bit_to_byte.push_back(h%2);
								h = h/2;
							}
						}
					//	size_t counter =0;
					//	std::cout<< "bit to byte: "<< bit_to_byte.size()<<std::endl;
						for(size_t n =0; n < bit_to_byte.size()-8; n++){
							unsigned char a = 0;
							for(size_t m = n; m <n+8; m++){
								a+= powersOfTwo.at(m-n)* bit_to_byte.at(m);
							}
							n= n+7;
					//		counter = counter+1;
							outs<< a;
						}	
					//	std::cout<< "counter: "<< counter << std::endl;
					}
				//	outs << (char) 0;
				}				
			}
			outs<<(char) 8;
	//	}
	//	outs.close();
		std::cout << "all the alignment patterns are written! " <<std::endl;
	}
/*
	Modification instructions:
	5 - single base modification (incl N)
	NUM_DELETE - delete 2^nd bases
	NUM_KEEP - keep 2^nd bases
	5 - insert a base
	


*/   
void mc_model::make_all_alignments_patterns(){
	std::string context; 
	std::set<std::string>pattern;
	// Alignment_level is makov chain level for alignments
	for(size_t i = 0 ; i < Alignment_level ; i++){
		context += (char)0;
	}
	pattern.insert(context);
	// this will create about (ND+NK+10)^Alignment_length patterns:
	for(size_t i =0; i < Alignment_level; i++) {
		std::set<std::string> intermediate_pattern;
		for(size_t j = 0; j <NUM_DELETE+NUM_KEEP+10; j++){
			// For each current pattern: modify position i to character j
			for(std::set<std::string>::iterator it = pattern.begin(); it!= pattern.end();it++){
				std::string seq = *it;
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
			//	std::cout<< "pattern: " << seq << std::endl;
			}
		}
		pattern.clear();
		for(std::set<std::string>::iterator it1 = intermediate_pattern.begin(); it1 != intermediate_pattern.end();++it1){
			std::string seq1 = *it1;
			pattern.insert(seq1);
		}
	} // for Alignment_level
	std::set<std::string> intermediate1_pattern;
	for(std::set<std::string>::iterator it = pattern.begin(); it != pattern.end();++it){
		for(size_t i = 0 ; i < 6 ; i++){
			std::string seq = *it;
			char c = i;
			std::string seq1 = seq + c;
		//	std::cout << "seq1: " << seq1 <<std::endl;
			intermediate1_pattern.insert(seq1);
		}
	}
	pattern.clear();
	for(std::set<std::string>::iterator it1 = intermediate1_pattern.begin(); it1 != intermediate1_pattern.end();++it1){
			std::string seq1 = *it1;
			pattern.insert(seq1);
	}
	for(std::set<std::string>::iterator it = pattern.begin();it !=pattern.end();it++){
		std::string seq = *it;
		all_alignment_patterns.insert(seq);
	}
	/*
		size_t number =0;
		for(std::set<std::string>::iterator it = pattern.begin(); it != pattern.end();++it){
			std::string seq = *it;
			number ++ ;
			std::string str = print_modification_character(seq.at(0));
//			std::string str1 = print_modification_character(seq.at(1));
			std::cout<< "" << str << "" <<  "" << int(seq.at(1)) << std::endl;			
		}
		std::cout<< "number: "<< number<<std::endl; 
 	 */
}

// TODO do we want to make markov chain levels dependent on input sequence length?
	void mc_model::train(std::ofstream & outs){
		make_all_the_patterns();
		markov_chain();
		markov_chain_alignment();
		write_parameters(outs);
		write_alignments_pattern(outs);

	}
	
	void mc_model::make_all_the_patterns(){
		std::string seq = "";
		std::set<std::string> pattern;
		for(size_t i = 0 ; i < Sequence_level ; i++){
			seq += dnastring::index_to_base(0);
		}
		pattern.insert(seq);
		for(size_t i = 0; i < Sequence_level;i++){	
			std::set<std::string> intermediate_pattern;
			for(size_t j = 0; j <5; j++){
				for(std::set<std::string>::iterator pat = pattern.begin();pat != pattern.end();++pat){
					std::string seq1 = *pat;
					seq1.at(Sequence_level-1-i)=dnastring::index_to_base(j);	
					intermediate_pattern.insert(seq1);
				}
			}
			for(std::set<std::string>::iterator pat1 = intermediate_pattern.begin(); pat1 != intermediate_pattern.end();++pat1){
				std::string seq2 = *pat1;
				pattern.insert(seq2);
			}
		}	
		for(std::set<std::string>::iterator it = pattern.begin();it !=pattern.end();it++){
			std::string seq3 = *it;
			all_the_patterns.insert(std::make_pair(seq3,std::vector<double>(5,0)));
		//	std::cout << " sequence pattern length " <<seq3.length() << std::endl;
		}
		std::cout << "All the sequence patterns are made! "<<std::endl;
	}
	void mc_model::set_patterns(std::ifstream& in){
		make_all_the_patterns();
		size_t bit = 12;
//		std::ifstream in("encode", std::ifstream::binary);
		std::map<std::string,std::vector<unsigned int> > values;
		char c;
		c = in.get();
		while(c != 8){
			std::string acc;
			std::stringstream s;
			while(c != 0){
				s << c;
				c = in.get();
			}
			s >> acc;
			std::cout << " acc " << acc <<std::endl;
			data.set_accession(acc);//since in decoding we have no access to our fasta file we need to set accession names in data class
			high.push_back(values);//initializing the high vector
			c = in.get();
		}
		std::vector<std::map<std::string , std::vector<unsigned int> > > al_values (data.numOfAcc());
		for(size_t i =0; i < data.numOfAcc(); i++){
			highValue.push_back(al_values);
		}
		char h;
		c = in.get();
		while(c != 8){
			unsigned int accession = 0;
			in >> accession;
			c = in.get();
			for(std::map<std::string,std::vector<double> >::const_iterator it= all_the_patterns.begin(); it!= all_the_patterns.end();it++){
				std::string pattern = it ->first;
			//	std::cout << "pattern: "<<pattern <<std::endl;
			//	std::cout<< "accession " << accession << "high size: " << high.size() << " numOfAcc "<< data.numOfAcc() << std::endl;
				high.at(accession).insert(std::make_pair(pattern, std::vector<unsigned int>(5,0)));
			}
			for(std::map<std::string,std::vector<unsigned int> >::iterator it= high.at(accession).begin(); it!= high.at(accession).end();it++){
				std::vector<bool> binary_high_value(0);
				size_t bound = (5*bit)/8;
			//	std::cout<< "bound: " << bound << std::endl;
				for(size_t j = 0 ; j < bound ; j++){ 
					h=in.get();
					size_t H = size_t(h); // should be a better solution!
			//		std::cout<<"h: "<< size_t(h)<<std::endl;
					for(size_t k = 0; k < 8 ; k++){
						binary_high_value.push_back(H%2);
				//		std::cout<< " "<< H%2 ;	
						H = H/2;
					}
				//		std::cout << " " << std::endl;
				}
				size_t counter = 0;
			//	std::cout<< "binary size= "<< binary_high_value.size()<<std::endl;
			//	std::cout<< "size- 8 "<<  binary_high_value.size() -( binary_high_value.size()%bit)<<std::endl;
				for(size_t i = 0; i < binary_high_value.size()-bit;i++){
					unsigned int high_value = 0;					
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*powersOfTwo.at(j-i);
					}
					i=i+bit-1;
					it -> second.at(counter)=high_value;
				//	if(it->first == "AC"){
				//		std::cout<< "high value of AC in std::set pattern: ";
				//		std::cout<< high_value<<std::endl;
				//	}
			//		std::cout<<"high value in model class: "<< high_value << " at " << counter << " i " << i <<std::endl;
			//		std::cout<< " "<<std::endl;
					counter = counter + 1;
				}
			//	std::cout<< "counter: "<< counter << std::endl;
			/*	unsigned int high_value = 0;
				for(size_t j =  binary_high_value.size() -( binary_high_value.size()%bit) ; j <  binary_high_value.size(); j++){
						high_value += binary_high_value.at(j)*powersOfTwo.at(j-binary_high_value.size()+( binary_high_value.size()%bound));					
				}*/
			//	it -> second.at(4)=powersOfTwo.at(bit); 
			//	std::cout << "high at 4 "<< it->second.at(4)<<std::endl;
			/*	for(size_t j =0; j < 5 ; j++){
					std::cout<< "high from stream: " << it->second.at(j)<<std::endl;
				}*/
			}
			c= in.get();	
		}
	//	std::cout<<"last c: "<< int(c)<<std::endl;
	}
	void mc_model::set_alignment_pattern(std::ifstream & in){
		make_all_alignments_patterns();
		std::cout << "all alignment patterns have been made! "<<std::endl;
		size_t bit = 12;
		char c;
		char h;
		c= in.get();	
		while(c != 8){
			unsigned int accession1 = 0;
			in >> accession1;
		//	std::string acc1;
		//	std::stringstream s1;
		//	while(c != 0){
		//		s1 << c;
				c = in.get();
		//	}
		//	s1 >> acc1;
		//	data.set_accession(acc1);//since in decoding we have no access to our fasta file we need to set accession names in data class
		//	accession1 = data.get_acc_id(acc1);
			unsigned int accession2 = 0;
			in >> accession2;
		//	std::string acc2;
		//	std::stringstream s2;
			c= in.get();
		//	while(c != 0){
		//		s2 << c;
		//		c = in.get();
		//	}
		//	s2 >> acc2;
		//	data.set_accession(acc2);
		//	accession2 = data.get_acc_id(acc2);	
			std::cout << "accession 1 "<< accession1 << " accession 2 "<<accession2 << std::endl;
			for(std::set<std::string>::const_iterator it= all_alignment_patterns.begin(); it!= all_alignment_patterns.end();it++){
				std::string pattern = *it;
		//		std::cout<<"acc1: "<<accession1 << " acc2 " << accession2<<std::endl;
				highValue.at(accession1).at(accession2).insert(std::make_pair(pattern, std::vector<unsigned int>(NUM_KEEP+NUM_DELETE+10,0)));
			}
			for(std::map<std::string,std::vector<unsigned int> >::iterator it= highValue.at(accession1).at(accession2).begin(); it!= highValue.at(accession1).at(accession2).end();it++){
				std::vector<bool> binary_high_value(0);
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
			//	std::cout<< "bit to byte 1: "<< binary_high_value.size() << std::endl;
				for(size_t i = 0; i < (binary_high_value.size())-bit;i++){
					unsigned int high_value = 0;					
					for(size_t j =i; j < i+bit; j++){
						high_value += binary_high_value.at(j)*powersOfTwo.at(j-i);
					}
					i=i+bit-1;
					it -> second.at(counter)= high_value;
					counter = counter +1;
				}
		//		std::cout<< "counter: "<< counter<<std::endl;
				it -> second.at(NUM_KEEP+NUM_DELETE+9)=powersOfTwo.at(bit); 
		//		std::cout<< "size of high value: "<< it ->second.size() << std::endl;
			}
			c= in.get();	
		}


	}
	std::vector<unsigned int> mc_model::get_high_at_position(size_t seq_index, size_t position)const{
		const dnastring & sequence = data.getSequence(seq_index);
	//	std::cout << " sequence " << seq_index << " length " << sequence.length() << std::endl;
		size_t accession = data.accNumber(seq_index);
		std::stringstream context;
		for(size_t j = Sequence_level; j>0; j--){
			if(position < j){
				char chr = 'A';
				context<<chr;
			}else{
				char chr = sequence.at(position-j);
				context<<chr;
			}
		}
		std::string current_pattern;
		context >> current_pattern;
/*		if(position ==112060){
			std::cout<<"current pattern in seq  "<< seq_index << " of accession " << accession <<std::endl;
			for(size_t j=0; j< current_pattern.length(); j++){
				std::cout<<current_pattern.at(j)<<std::endl;
			}
		}*/
		std::map<std::string, std::vector<unsigned int> >::const_iterator it=high.at(accession).find(current_pattern);
		assert(it!=high.at(accession).end());
	//	std::cout << " sequence " << seq_index << " length " << sequence.length() << " " << current_pattern<<  std::endl;
	//	std::cout<<"accession: "<< accession << std::endl;
	//	for(size_t k =0; k< 5; k++){
	//		std::cout<< "high at " << k << " is "<< it->second.at(k)<<std::endl;
	//	}
	//	if(current_pattern == "GA"){
	//		std::cout<< "high at GA: ";
	//		for(size_t j = 0; j < 5 ; j++){
	//			std::cout<< it->second.at(j) <<std::endl;
	//		}
	//	}
		return it->second;
	}
	std::vector<unsigned int> mc_model::get_center_high_at_position(size_t cent_ref, size_t cent_left, size_t position)const{
		const dnastring & sequence = data.getSequence(cent_ref);
	//	std::cout << " sequence " << seq_index << " length " << sequence.length() << std::endl;
		size_t accession = data.accNumber(cent_ref);
		std::stringstream context;
		for(size_t j = Sequence_level; j>0; j--){
			if(position < cent_left+j){
				char chr = 'A';
				context<<chr;
			}else{
				char chr = sequence.at(position-j);
				context<<chr;
			}
		}
		std::string current_pattern;
		context >> current_pattern;
	//	std::cout<<"current pattern: "<< current_pattern<<std::endl;					
	/*	std::cout<<"current pattern in seq  "<< seq_index << " of accession " << accession <<std::endl;
		for(size_t j=0; j< current_pattern.length(); j++){
			std::cout<<current_pattern.at(j)<<std::endl;
		}*/
		std::map<std::string, std::vector<unsigned int> >::const_iterator it=high.at(accession).find(current_pattern);
		assert(it!=high.at(accession).end());
	//	std::cout << " sequence " << seq_index << " length " << sequence.length() << " " << current_pattern<<  std::endl;
	//	std::cout<<"accession: "<< accession << std::endl;
	//	for(size_t k =0; k< 5; k++){
	//		std::cout<< "high at " << k << " is "<< it->second.at(k)<<std::endl;
	//	}
		return it->second;
	}
	std::vector<unsigned int> mc_model::get_reverse_center_high_at_position(size_t cent_ref, size_t cent_left, size_t position)const{
		const dnastring & sequence = data.getSequence(cent_ref);
	//	std::cout << " sequence " << seq_index << " length " << sequence.length() << std::endl;
		size_t accession = data.accNumber(cent_ref);
		std::stringstream context;
		for(size_t j = Sequence_level; j>0; j--){
		//	if(position > cent_right-j){
			if(position < cent_left+j){
				char chr = 'A';
				context<<chr;
			}else{
			//	char chr= dnastring::complement(sequence.at(position+j));				
				char chr = sequence.at(position-j);
				context<<chr;
			}
		}
		std::string current_pattern;
		context >> current_pattern;
	//	std::cout<<"current pattern: "<< current_pattern<<std::endl;					
	/*	std::cout<<"current pattern in seq  "<< seq_index << " of accession " << accession <<std::endl;
		for(size_t j=0; j< current_pattern.length(); j++){
			std::cout<<current_pattern.at(j)<<std::endl;
		}*/
		std::map<std::string, std::vector<unsigned int> >::const_iterator it=high.at(accession).find(current_pattern);
		assert(it!=high.at(accession).end());
	//	std::cout << " sequence " << seq_index << " length " << sequence.length() << " " << current_pattern<<  std::endl;
	//	std::cout<<"accession: "<< accession << std::endl;
	//	for(size_t k =0; k< 5; k++){
	//		std::cout<< "high at " << k << " is "<< it->second.at(k)<<std::endl;
	//	}
		return it->second;
	}

	std::vector<size_t> mc_model::get_powerOfTwo()const{
		return powersOfTwo;
	}		
	const std::map<std::string, std::vector<unsigned int> >&  mc_model::get_high(size_t acc)const{
		return high.at(acc);
	}
	std::string mc_model::get_context(size_t position, size_t seq_id)const{
		const dnastring & sequence = data.getSequence(seq_id);
		std::stringstream context;
		for(size_t j = Sequence_level; j > 0; j--){
			if(position<j){
				char chr = 'A';
				context << chr;
			}else{
				char chr = sequence.at(position-j);
				context << chr;
			}
		}
		std::string current_pattern;
		context >> current_pattern;	
		return current_pattern;
	}
	std::string mc_model::get_firstPattern()const{
		std::string first_pattern;
		for(size_t i = 0; i< Sequence_level; i++){
			first_pattern += "A";
		}
		return first_pattern;
	}
	std::string mc_model::get_firstAlignmentPattern()const{
		std::string first_pattern;
		size_t firstPatterns = Alignment_level;
		for(size_t j = 0; j < Alignment_level; j++){
			firstPatterns --;
			first_pattern += modification_character(-1,-1,-1,firstPatterns);
		}
		return first_pattern;
	}	
	char mc_model::modification_character(int modify_base, int num_delete, int insert_base, int num_keep)const {
		//return the enc
		if(num_delete != -1){
		//	std::cout<< "there is a delete of length "<< num_delete <<std::endl; 
		//	return NUM_DELETE+num_delete;
			return 5 + num_delete;
		}
		if(num_keep != -1){
		//	std::cout<< "there is a keep of length" << num_keep << std::endl;
		//	return NUM_KEEP+num_keep;
			return 5 + NUM_DELETE + num_keep;
		}
		if(modify_base != -1) {
		//	std::cout<< "there is a modification at " << dnastring::index_to_base(modify_base) << std::endl;
			return modify_base;
		}
		if(insert_base != -1){
		//	std::cout<< " there is a insertion at " << dnastring::index_to_base(insert_base) <<std::endl;
			return insert_base + NUM_KEEP + NUM_DELETE + 5;
		}
		assert (false);
		return -1;
	}
	std::string mc_model::print_modification_character(char enc)const{
		int modify_base = -1;
		int num_delete =-1;
		int insert_base = -1;
		int num_keep = -1;
		modification(enc, modify_base, num_delete, insert_base, num_keep);
		std::stringstream s;
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
			length = 0;
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
	void mc_model::computing_modification_oneToTwo(const pw_alignment & p, abstract_context_functor & functor)const{
		std::string seq = "";
	//	std::cout<<"data ad in computing mod: "<< & data << std::endl;
	//	p.print();
		size_t left_1; 
	//	size_t left_2;
		size_t right_1;
	//	size_t right_2;
		p.get_lr1(left_1,right_1);
		size_t acc1 = data.accNumber(p.getreference1());
		size_t acc2 = data.accNumber(p.getreference2());
		size_t first_patterns = Alignment_level;
		size_t power = powersOfTwo.at(NUM_KEEP-1);
	//	std::cout<< "power: "<<powersOfTwo.at(0)<<std::endl;	
		for(size_t j = 0; j < Alignment_level; j++){// two keeps of length 2^1 and 2^0 has been created at the begining of each alignment (if level is 2).
			first_patterns--;
			seq+=modification_character(-1,-1,-1,first_patterns);
		}
	//	std::cout<<p.alignment_length()<<std::endl;
		for (size_t i = 0; i< p.alignment_length(); i++){
			size_t n = 0; // number of alignment columns we look at
			int modify_base =-1;
			int num_delete=-1;
			int insert_base=-1;
			int num_keep=-1;
			std::string seq1(" ",Alignment_level+1);
			char seq2;
			for(size_t w = Alignment_level; w>0 ;w--){
				seq1.at(Alignment_level-w)=seq.at(seq.size()-w);
			}
			char s1chr;
			char s2chr;
			size_t s1;
			size_t s2;
			char s1nchr;
			char s2nchr;
			size_t s1n;
			size_t s2n;
			p.alignment_col(i, s1chr, s2chr);				
			s1 = dnastring::base_to_index(s1chr);
			s2 = dnastring::base_to_index(s2chr);
			for(size_t j =i; j < p.alignment_length(); j++){
				p.alignment_col(j, s1nchr, s2nchr);				
				s1n = dnastring::base_to_index(s1nchr);
				s2n = dnastring::base_to_index(s2nchr);
				if(s1n == 5){
				
				}else break;	
			}
		//	if(p.getreference1() == 28 &&left_1 == 176557){
		//		std::cout<<"using s1n! " << s1n <<std::endl;
		//	}
			seq1.at(Alignment_level)=s1n;
			if(s1 == s2){
			//	std::cout<< "keep at "<< i <<std::endl;
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
					std::cout<<"keep length at "<< i << " is  "<< klength<<std::endl;
				}*/
			//	std::cout<<"keep length: "<<klength<<std::endl;
			//	std::cout<<"NUM KEEP: "<< NUM_KEEP-1<<std::endl;
				if(klength > power){
					num_keep = NUM_KEEP-1;
					n=powersOfTwo.at(num_keep)-1;
				//	std::cout<<"Long keep"<<std::endl;
					seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
				}else{
					for (size_t m = NUM_KEEP; m > 0; m--){
				//	for (size_t m = powersOfTwo.size()-1; m >= 0; m--)
						if((klength & powersOfTwo.at(m-1)) != 0){
							num_keep=m-1;
						//	std::cout<<"m: "<<m <<std::endl;
							n= powersOfTwo.at(num_keep)-1;
							seq += modification_character(modify_base,num_delete,insert_base,num_keep);
							break;
						}
					}
				}
				seq2 = seq.at(seq.size()-1);
		//		std::cout<< "recorded keep length: " << n+1 << std::endl;
		//		std::cout<< "size of seq: "<< seq.size()<<std::endl;
		//		std::cout << "seq 2 " << int(seq2) << std::endl;
		//		std::cout<< "recorded keep length: " << n+1 << std::endl;
		//		std::cout<< "a keep is observed!" <<std::endl;
				functor. see_context(acc1,acc2,p,i,seq1,seq2);
		//		std::cout<<"n: "<< n << std::endl;
			}else{
				if((s1!=5) & (s2!=5)){
					modify_base = s2;
				/*	if(i<51){
						std::cout<<"modification at "<< i << " is  "<<s1 <<std::endl;
					}*/
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc1,acc2,p,i,seq1,seq2);
				//	std::cout<< "seq1" << seq1 <<std::endl;
				}
				if(s1 == 5){
					insert_base = s2;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);						
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc1,acc2,p,i,seq1,seq2);
			//		std::cout<< "seq1" << seq1 <<std::endl;						
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
						for(size_t m = NUM_DELETE; m > 0; m--){
					//	std::cout<< "m in delete1: "<< m << std::endl;
//						for (size_t m = powersOfTwo.size()-1;m>=0;m--)
							if((dlength & powersOfTwo.at(m-1)) != 0){
								num_delete = m-1;
								n= powersOfTwo.at(num_delete)-1;
								seq += modification_character(modify_base,num_delete,insert_base,num_keep);		
								break;
							}
						}
					}
					seq2 = seq.at(seq.size()-1);
					functor. see_context(acc1,acc2,p,i,seq1,seq2);
				}

			}
			/*std::cout<< " context1 is "<<std::endl;
			for(size_t h = 0 ; h < seq1.size(); h++){
				std::cout<< int(seq1.at(h))<<std::endl;
			}*/
			//std::cout<<"last char: "<<int(seq2)<<std::endl;
			//std::cout<<"i+n: "<< i+n<<std::endl;
//			std::cout<< " context1 is "<<std::endl;
//			for(size_t h = 0 ; h < seq1.size(); h++){
//				std::cout<< int(seq1.at(h))<<std::endl;
//			}
	//		std::cout<<"last char in mod function: "<<int(seq2)<<std::endl;
//			std::cout<<"i+n: "<< i+n<<std::endl;
			i=i+n;
			//std::cout<<"n"<< n << std::endl;
			//	uint32_t initial = 1;
			//	for (uint32_t m = 0; m <32 ; m++)
			//		uint32_t power_of_two = initial << m; 
			//		if((klength & power_of_two) !=0)
		}

	//	std::cout<<"The alignment is: "<<std::endl;
	//	p.print();
	//	std::cout<<"encoded sequence from one to two is: "<<std::endl;
	//	for(size_t m = 0; m < seq.size(); m ++){
	//		std::cout<< int(seq.at(m))<<std::endl;
	//	}
	//	std::cout<< "one to two was done! " << std::endl;
		functor.see_entire_context(acc1,acc2,seq);
	}

	
	void mc_model::computing_modification_twoToOne(const pw_alignment & p, abstract_context_functor & functor)const{
		std::string seq = "";
		size_t acc1 = data.accNumber(p.getreference1());
		size_t acc2 = data.accNumber(p.getreference2());
		size_t first_patterns = Alignment_level;
	/*	if(p.getreference2()==28&&p.getreference1()==0 && p.getbegin2()==176557){
			std::cout << "forward"<<std::endl;
			p.print();
		}
		if(p.getreference2()==28&&p.getreference1()==0 && p.getend2()==176557){
			std::cout << "reverse"<<std::endl;
			p.print();
		}*/
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
			std::string seq1(" ",Alignment_level+1);
			char seq2;
			for(size_t w = Alignment_level; w>0 ;w--){
				seq1.at(Alignment_level-w)=seq.at(seq.size()-w);
			//	std::cout<< "seq at size - w: " << seq.at(seq.size()-w)<<std::endl;
			}
			char s1chr;
			char s2chr;
			size_t s1;
			size_t s2;
			p.alignment_col(i, s1chr, s2chr);				
			s1 = dnastring::base_to_index(s1chr);
			s2 = dnastring::base_to_index(s2chr);
			char s1nchr;
			char s2nchr;
			size_t s2n;
			if(s2 == 5){
			//	if(p.getreference2()==28&&p.getreference1()==0){std::cout<< "i " << i << " s2 " << s2 << std::endl;}
				for(size_t j =i; j < p.alignment_length(); j++){
					p.alignment_col(j, s1nchr, s2nchr);				
					s2n = dnastring::base_to_index(s2nchr);
					if(s2n != 5){
						seq1.at(Alignment_level)=s2n;
						break;
					}else continue;
				}
			}else seq1.at(Alignment_level)=s2;
	//		std::cout<<"using s2n! " << s2n <<std::endl;
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
					std::cout<<"keep length at "<< i << " is  "<< klength<<std::endl;
				}*/
				if(klength > powersOfTwo.at(NUM_KEEP-1)){
//				if(klength > powersOfTwo.at(powersOfTwo.size()-1))
//					num_keep= powersOfTwo.size()-1;
					num_keep = NUM_KEEP-1;
					n=powersOfTwo.at(num_keep)-1;
					seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
				}else{
					for(size_t m = NUM_KEEP; m > 0; m--){
				//	std::cout<<"m: "<<m <<std::endl;
//					for (size_t m =powersOfTwo.size()-1; m>= 0; m--)
						if((klength & powersOfTwo.at(m-1)) != 0){
							num_keep=m-1;
							n=powersOfTwo.at(num_keep)-1;						
							seq += modification_character(modify_base,num_delete,insert_base,num_keep);
							break;
						}
					}
				}
				seq2 = seq.at(seq.size()-1);
				functor. see_context(acc2,acc1,p,i,seq1,seq2);
			}else{
				if((s1!=5) & (s2!=5)){
					modify_base = s1;
				/*	if(i<51){
						std::cout<<"modification at "<< i << " is  "<<s2 <<std::endl;
					}*/
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc2,acc1,p,i,seq1,seq2);
				}
				if(s2 == 5){
					insert_base = s1;
					seq += modification_character(modify_base,num_delete,insert_base,num_keep);					
					seq2 = modification_character(modify_base,num_delete,insert_base,num_keep);
					functor. see_context(acc2,acc1,p,i,seq1,seq2);						
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
				//	std::cout<< "dlength: " << dlength << std::endl;
					if(dlength > powersOfTwo.at(NUM_DELETE-1)){
//					if(dlength > powersOfTwo.at(powersOfTwo.size()-1))
//						num_delete= powersOfTwo.size()-1;
						num_delete = NUM_DELETE-1;
						n = powersOfTwo.at(num_delete)-1;
						seq+=modification_character(modify_base,num_delete,insert_base,num_keep);
					}else{
						for(size_t m = NUM_DELETE ; m > 0; m--){
//						for (size_t m = powersOfTwo.size()-1; m>0; m--)
						//	std::cout<< "m in delete2: "<< m << std::endl;
							if((dlength & powersOfTwo.at(m-1)) != 0){
								num_delete = m-1;
								n = powersOfTwo.at(num_delete)-1;
								seq += modification_character(modify_base,num_delete,insert_base,num_keep);
								break;			
							}
						}
					}
					seq2 = seq.at(seq.size()-1);
					functor. see_context(acc2,acc1,p,i,seq1,seq2);
				}
			}
	/*		std::cout<< " context2 is "<<std::endl;
			for(size_t h = 0 ; h < seq1.size(); h++){
				std::cout<< int(seq1.at(h))<<std::endl;
			}*/
		//	std::cout<<"last char in mod function: "<<int(seq2)<<std::endl;
			i=i+n;
		//	std::cout<< "i in modification : " << i << std::endl;
		}		
	//	std::cout<<"encoded sequence from two to one is: "<<std::endl;
	/*	for(size_t m = 0; m < seq.size(); m ++){
			std::cout<< int(seq.at(m))<<std::endl;
		}*/
		functor.see_entire_context(acc2,acc1,seq);
	}
	const std::map<std::string, std::vector<double> > & mc_model::getPattern(size_t acc)const{
		return sequence_successive_bases.at(acc);

	}
	
	const std::vector<double> & mc_model::get_create_cost(size_t acc)const{
		return create_cost.at(acc);		

	}
	
	const std::vector<std::map<std::string, vector<double> > > & mc_model::model_parameters()const{
		return sequence_successive_bases;
//			std::cout<< " context2 is "<<std::endl;
//			for(size_t h = 0 ; h < seq1.size(); h++){
//				std::cout<< int(seq1.at(h))<<std::endl;
//		}
//		std::cout<<"last char: "<<int(seq2)<<std::endl;	
//		std::cout<<"encoded sequence from two to one is: "<<std::endl;
//		for(size_t m = 0; m < seq.size(); m ++){
//			std::cout<< int(seq.at(m))<<std::endl;
//		}

	}
	void mc_model ::get_encoded_member(pw_alignment & al,size_t center_ref,size_t center_left, encoding_functor & functor,std::ofstream& outs)const{
		size_t acc1 = data.accNumber(al.getreference1());
		size_t accession = data.accNumber(center_ref);
		size_t left_1; 
		size_t left_2;
		size_t right_1;
		size_t right_2;
		std::cout<< "al length: "<< al.alignment_length()<<std::endl;
		al.get_lr1(left_1,right_1);
		al.get_lr2(left_2,right_2);
		std::cout<< " left1: "<< left_1 << " left2: "<<left_2<< "center left: "<< center_left << std::endl;
		std::cout << "center accession: "<< accession << " accession1: "<< acc1 << std::endl;
		if(al.getreference1()==center_ref && left_1 == center_left){
			std::cout<< "center is on ref1"<<std::endl;
			computing_modification_oneToTwo(al, functor);	
		}
		else{
			std::cout<< "center is on ref2"<<std::endl;
			computing_modification_twoToOne(al, functor);
		}
	}
	void mc_model::get_insertion_high_between_centers(size_t& seq_id ,char & seq_base, char & last_center_base, unsigned int& center_ref,unsigned int & high, unsigned int & low)const{
		size_t acc1 = data.accNumber(center_ref);		
		size_t acc2 = data.accNumber(seq_id);
		std::string pattern;
		size_t bit = 13;
		for(size_t i =0; i < Alignment_level; i++){
			pattern += modification_character(-1,-1,-1,(Alignment_level-i));
		}
		pattern += last_center_base;
		size_t insertion = dnastring::base_to_index(seq_base);
		char insert = modification_character(-1,-1,insertion,-1);
		std::map<std::string , std::vector<unsigned int> >::const_iterator it = highValue.at(acc1).at(acc2).find(pattern);
		if(insert != NUM_KEEP+NUM_DELETE+10-1){
			high = it->second.at(insert);
		}else{
			high = powersOfTwo.at(bit);
		}
		low = it->second.at(insert -1);
	}
	void mc_model::get_a_keep(unsigned int & cent_left, unsigned int& cent_ref, size_t & length, std::vector<unsigned int> & low, std::vector<unsigned int> & high)const{
		size_t acc1 = data.accNumber(cent_ref);	
		size_t power = powersOfTwo.at(NUM_KEEP-1);	
		const dnastring & sequence = data.getSequence(cent_ref);		
		std::string pattern;
		for(size_t i =0; i < Alignment_level; i++){
			pattern += modification_character(-1,-1,-1,(Alignment_level-i));
		}
		size_t pos = cent_left;
		size_t last_base = dnastring::base_to_index(sequence.at(pos));
		std::string current_pattern = pattern;
		current_pattern += last_base;
		std::cout << "len "<<length << "pow "<< power << std::endl;
		while(length > power){
			char keep = modification_character(-1,-1,-1,NUM_KEEP-1);
			length = length - power;
			std::map<std::string , std::vector<unsigned int> >::const_iterator it = highValue.at(acc1).at(acc1).find(current_pattern);
			std::cout << "current pattern "<<int(current_pattern.at(0)) << " " << int(current_pattern.at(1)) << " keep "<< int(keep)<<std::endl;
			high.push_back(it->second.at(keep));
			low.push_back(it->second.at(keep -1));
			pattern.clear();
			for(size_t i = 1; i < current_pattern.size()-1;i++){
				pattern +=current_pattern.at(i);
			}
			pattern += keep;
			pos += power;
			last_base = dnastring::base_to_index(sequence.at(pos));	
			pattern +=last_base;
			current_pattern = pattern;
		}
		for (size_t m = NUM_KEEP; m > 0; m--){
			if((length & powersOfTwo.at(m-1)) != 0){
				std::map<std::string , std::vector<unsigned int> >::const_iterator it = highValue.at(acc1).at(acc1).find(current_pattern);
				char keep = modification_character(-1,-1,-1,m-1);
				high.push_back(it->second.at(keep));
				low.push_back(it->second.at(keep -1));
				break;
			}
		}
	}

	void mc_model::computing_modification_in_cluster(std::string center, std::string member)const{
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
	void abstract_context_functor::see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, std::string context, char last_char){
		
	}
	void abstract_context_functor::see_entire_context(size_t acc1,size_t acc2, std::string entireContext){

	}
	
	counting_functor::counting_functor(all_data & d):data(d), successive_modification(d.numAcc(),std::vector<std::map<std::string, std::vector<size_t> > >(d.numAcc())),total(d.numAcc(),std::vector<std::map<std::string, double > >(d.numAcc())) {}
	void counting_functor::see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, std::string context, char last_char){
	//	std::cout<< "accession 1: " << acc1 << " accession 2: " << acc2 << " size: " << pos << " last char: " << dnastring::base_to_index(last_char) << " " << int(last_char)<<std::endl;
	//	std::cout<< "context is: "<< std::endl;
	/*	for(size_t i = 0 ; i < context.size(); i++){
			std::cout<< int(context.at(i))<<std::endl;
		}*/
		std::map <std::string, std::vector<size_t> >::iterator it1= successive_modification.at(acc1).at(acc2).find(context);
		if(it1==successive_modification.at(acc1).at(acc2).end()) {
			successive_modification.at(acc1).at(acc2).insert(std::make_pair(context, std::vector<size_t>((NUM_DELETE+NUM_KEEP+10),1)));
			it1= successive_modification.at(acc1).at(acc2).find(context);
			assert(it1 != successive_modification.at(acc1).at(acc2).end());
		}
			it1->second.at(last_char)++;
		//	std::cout<< "context is: "<< context.length() << std::endl;
			//std::cout<< context <<std::endl;			
		//	std::cout<<"number of happening "<<int(last_char)<< " after above context is "<< it1->second.at(last_char)<<std::endl;
	}

/*
	compute total (sum over all counts)

*/
void counting_functor::total_context(){
	for(size_t i = 0; i < data.numAcc(); i++){
		for(size_t j = 0; j<data.numAcc(); j++){
			for(std::map<std::string, std::vector<size_t> >::iterator it = successive_modification.at(i).at(j).begin(); it!= successive_modification.at(i).at(j).end();it++){
				std::string context = it->first;
				std::map<std::string, double >::iterator it1=total.at(i).at(j).find(context);
				if(it1 == total.at(i).at(j).end()){
					total.at(i).at(j).insert(std::make_pair(context,0));
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
void counting_functor::total_context(size_t & acc1, size_t & acc2){
	for(std::map<std::string, std::vector<size_t> >::iterator it = successive_modification.at(acc1).at(acc2).begin(); it!= successive_modification.at(acc1).at(acc2).end();it++){
		std::string context = it->first;
		std::map<std::string, double >::iterator it1=total.at(acc1).at(acc2).find(context);
		if(it1 == total.at(acc1).at(acc2).end()){
			total.at(acc1).at(acc2).insert(std::make_pair(context,0));
			it1=total.at(acc1).at(acc2).find(context);
		}
		for(size_t k = 0; k < NUM_DELETE+NUM_KEEP+10; k++){
			it1->second += it->second.at(k);	
		}
		// we assume to find everything at least twice in the training to avoid events with zero modification cost
		if(it1->second < 2) {
			it1->second = 2;
		}
	}
}
double counting_functor::get_total(size_t acc1, size_t acc2, std::string context)const{
	std::map<std::string, double>::const_iterator it = total.at(acc1).at(acc2).find(context);
//	std::cout << " a " << acc1 << " " << acc2 << " : " << total.at(acc1).at(acc2).size() << std::endl;
	assert(it!=total.at(acc1).at(acc2).end());
	return it->second;
}

			
/*
	initialize context counting with 1
*/   	
	void counting_functor::create_context(size_t acc1, size_t acc2, std::string context) {
		std::map<std::string, std::vector<size_t> >::iterator it = successive_modification.at(acc1).at(acc2).find(context);
		if(it==successive_modification.at(acc1).at(acc2).end()) {
			successive_modification.at(acc1).at(acc2).insert(std::make_pair(context, std::vector<size_t>(NUM_DELETE+NUM_KEEP+10, 1)));
		}
	}
	
	const  std::map<std::string, std::vector<size_t> > & counting_functor::get_context(size_t acc1, size_t acc2)const{
		return successive_modification.at(acc1).at(acc2);
	}

	cost_functor::cost_functor(all_data & d, const std::vector<vector<std::map<std::string, vector<double> > > > & mod_cost):data(d){
		modify1 = 0;
		modify2 = 0;		
		modification = mod_cost;
	}
	void cost_functor::see_context(size_t acc1, size_t acc2, const pw_alignment & p, size_t pos, std::string context, char last_char){
		size_t ref1 = p.getreference1();
		size_t ref2 = p.getreference2();
		size_t accession1 = data.accNumber(ref1);
		size_t accession2 = data.accNumber(ref2);
		if(acc1 == accession1){//if acc1 is the first accession
			std::map<std::string, std::vector<double> >::const_iterator it = modification.at(acc1).at(acc2).find(context);
		//	std::cout<< "modification cost at "<< pos << " is "<< it->second.at(last_char)<<std::endl;
			modify1 +=it->second.at(last_char);
		}
		if(acc1 == accession2){
			std::map<std::string, std::vector<double> >::const_iterator it = modification.at(acc1).at(acc2).find(context);
		//	std::cout<< "modification cost at "<< pos << " is "<< it->second.at(last_char)<<std::endl;
			modify2 +=it->second.at(last_char);
		}
	}
	double cost_functor::get_modify(const pw_alignment & p,size_t acc1, size_t acc2)const{
		double modify= 0.0; // Attention : I change the init, before : double modify;
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
	encoding_functor::encoding_functor(all_data & d, mc_model * a_model, wrapper & wrap, dlib::entropy_encoder_kernel_1 & encode ):data(d),model(a_model),wrappers(wrap),enc(encode){
	}
	void encoding_functor::see_context(size_t acc1, size_t acc2,const pw_alignment & p, size_t pos, std::string context, char last_char){//last_char is infact a encoded pattern!
		size_t bit = 13;
	//	dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
	//	enc->std::set_stream(outs);
		unsigned int total = model->get_powerOfTwo().at(bit)+20;
		std::map<std::string, std::vector<unsigned int> >::const_iterator it1 = model->get_highValue(acc1,acc2).find(context);// if modification is from acc2 to acc1 the order is already exchanged. So this is true
		std::vector<unsigned int> low(NUM_DELETE+NUM_KEEP+10,0);
		std::vector<unsigned int> high(NUM_DELETE+NUM_KEEP+10,0);
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
	//	for(size_t h =0; h < NUM_DELETE+NUM_KEEP+10; h++){
		//	std::cout<< "high values: " << high.at(h)<<std::endl;
	//	}

	//	std::cout<< "context at " << pos << " is: ";
		for(size_t i =0; i < context.size(); i++){
	//		std::cout<<int(context.at(i));
			int con = int(context.at(i));
			wrappers.context(con);
		}
	//	std::cout<< " " <<std::endl;
	//	std::cout<< "actual acc1 " << data.accNumber(p.getreference1())  << " actual acc2 " << data.accNumber(p.getreference2()) <<std::endl;
	//	std::cout << "encoding form acc: " << acc1 << " to acc: "<< acc2 << std::endl;
	//	std::cout<< "center acc should be: " << acc1 << std::endl;
	//	std::cout << " ended char in al: "<< int(last_char)<<std::endl;
	//	std::cout<< "encoded low: "<< low.at(last_char)<<" encoded high: "<< high.at(last_char)<<std::endl;
		enc.encode(low.at(last_char),high.at(last_char),total);
		wrappers.encode(low.at(last_char),high.at(last_char),total);
	//	delete enc;
	}

	void encoding_functor::see_entire_context(size_t acc1, size_t acc2, std::string entireContext){
		alignment_pattern = entireContext;
	}
	const std::map<std::string, std::vector<double> > & encoding_functor::get_alignment_context()const{
		for(std::map<std::string, std::vector<double> >::const_iterator it = alignment_context.begin(); it !=alignment_context.end(); it++){
		//	for(size_t n =0; n < NUM_DELETE+NUM_KEEP+10; n++){
		//		std::cout << "counts: "<< it->second.at(n)<<std::endl;
		//	}
		}
		return alignment_context;
	}

