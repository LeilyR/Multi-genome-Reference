#include <iostream>
#include <cstdlib>
#include <ctime>

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <ostream>
#include <vector>
#include "pw_alignment.hpp"
#include "data.hpp"
#include "dynamic_mc.hpp"
#include "overlap.hpp"
#include "model.hpp"
#include "connected_component.hpp"
#include "encoder.hpp"
#include "dynamic_encoder.hpp"
#include "dynamic_decoder.hpp"
#include "test.hpp"
//suffix tree:
#include "suffix_tree.hpp"
//Merging centers:
#include "merging_centers.hpp"
//interval tree:
#include "intervals.hpp"
#include "alignment_index.hpp"
#include <omp.h>
//Filtering the alignments
#include "filter_alignments.hpp"


//Aligning reads:
#include "needleman_wunsch.hpp"
#include "read_daligner.hpp"
#include "graph.hpp" //alignment graph that is used for mapping read against the reference
#include "map_read.hpp"

#define NO_MAKEFILE
#include "dlib/entropy_encoder/entropy_encoder_kernel_1.h"
#include "dlib/entropy_decoder/entropy_decoder_kernel_1.h"
#undef NO_MAKEFILE

#include "overlap.cpp"

#define VERSION "0.0.1" 
#define NUMVERISON 0
#define MIN_READ_VERSION 0

#define TIMING 1


#define TEST 0



/**

	General TODO

Look at diversity within clusters. Can this diversity be used to differentiate between orthologous and paralogous groups,
should they be treated differently in mapping? Some more stringent approach to place repeat instances? Discuss this with Felix again?

Maybe another round of AP-clustering after all reads are aligned? Maybe include not only reads, but also all origingal sequences in that cluster
maybe do multiple alignments?

**/

#if TEST

int main(int argc, char * argv[]) {
	std::cout << " hello " << std::endl;
	Graph g = Graph();


	//  runNeedleman("CCCGGGGGGTGCA","ATAGTTGCA",2);
	//load data
	all_data data = all_data();

	//Init new Graph //TODO in main or parseData ?
	Graph newGraph = Graph();
	g.parseData(data,newGraph);

	//std::cout << newGraph;
	//std::cout <<" New Graph " <<newGraph << std::endl;




/*
	pw_alignment p(std::string("ATT----TTCTT"), string("AGTGATAT----"), 12, 15, 23, 26,1,1);
	pw_alignment s(std::string("ATT----TTCTT"), string("ACTGATG---AC"),13, 18, 24, 29,2,1);

	for(size_t i=0; i<p.alignment_length(); ++i) {
		char s1;
		char s2;
		p.alignment_col(i, s1, s2);
		std::cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << std::endl;
	  
	}

	pw_alignment p1;
	pw_alignment p2;
	p.split(false,4, p1, p2);
		for(size_t i=0; i<p2.alignment_length(); ++i) {
		char s1;
		char s2;
		p2.alignment_col(i, s1, s2);
		std::cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << std::endl;
	 
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
		std::cout << "pos " << i<< " s1 " << s1 << " s2 " << s2 << std::endl;
}

}
*/
	return 0;

}



#endif
void usage() {
	cerr << "Graph program" << std::endl;
	cerr << "Usage:" << std::endl;
}

int do_fasta_prepare(int argc, char * argv[]) {

	if(argc < 4) {
		usage();
		cerr << "Program fasta_prepare" << std::endl; 
		cerr << "Combine fasta files from different species/accessions into a common file" << std::endl;
		cerr << "Parameters: " << std::endl;
		cerr << "output file -- fasta file" << std::endl;
		cerr << "<list of fasta files> -- input file, one file per species/accession, no colon in file name"<< std::endl;
		return 1;
	}


	std::ofstream outf(argv[2]);
	if(outf) {
		std::vector<std::ifstream*> inputs;
		std::vector<std::string> names;
		for(int i=3; i<argc; ++i) {
			std::string istr(argv[i]);
			std::vector<std::string> slashparts;
			strsep(istr, "/", slashparts);
			std::string afterslash = slashparts.at(slashparts.size()-1);
			std::vector<std::string> periodparts;
			strsep(afterslash, ".", periodparts);
			std::string name = periodparts.at(0);
			if(periodparts.size()>1) {
				std::stringstream sstr;
				for(size_t j=0; j<periodparts.size()-1; j++) {
					if(j!=0) sstr << ".";
					sstr << periodparts.at(j);
				}
				name = sstr.str();
			}
			for(size_t j=0; j<name.length(); ++j) {
				if(':' == name.at(j)) {
					cerr << "Error: " << argv[i] << " file name contains a colon" << std::endl;
					exit(1);
				}
			}
			std::ifstream * iin = new ifstream(argv[i]);
			if(!(*iin)) {
				cerr << "Error: cannot read: " << argv[i] << std::endl;
				exit(1);
			}
			inputs.push_back(iin);
			names.push_back(name);
		}

		for(size_t i=0; i<inputs.size(); ++i) {
			std::string str;
			size_t j = 0;
			while(getline((*(inputs.at(i))), str)) {
				if(str.length()>0) { // ignore empty lines
				if(str.at(0)=='>') {
					std::string nname = str.substr(1);
					std::stringstream nhead;
				//	nhead << ">" << i << ":" << j;
					nhead << ">" << names.at(i) << ":" << nname;
					outf << nhead.str() << std::endl;
				//	outf << nhead.str();
					j = j +1;

						
						
/*
	for(size_t col = 0; col < alignment_length(); col++) {
		char c1;
		char c2;
		alignment_col(col, c1, c2);
		std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
	if(col>50) break;	
	}
*/	

				} else {
				//	outf << str;
					outf << str << std::endl;
				}
				}
			}
			inputs.at(i)->close();
			delete inputs.at(i);
		}
		outf.close();
	} else {
		cerr << "Error: cannot write: " << argv[2] << std::endl;
		exit(1);
	}

	outf.close();

	return 0;
}

/*
    src -- The name of one of the source sequences for the alignment. For sequences that are resident in a browser assembly, the form 'database.chromosome' allows automatic creation of links to other assemblies. Non-browser sequences are typically reference by the species name alone.
    start -- The start of the aligning region in the source sequence. This is a zero-based number. If the strand field is '-' then this is the start relative to the reverse-complemented source sequence.
    size -- The size of the aligning region in the source sequence. This number is equal to the number of non-dash characters in the alignment text field below.
    strand -- Either '+' or '-'. If '-', then the alignment is to the reverse-complemented source.
    srcSize -- The size of the entire source sequence, not just the parts involved in the alignment.
    text -- The nucleotides (or amino acids) in the alignment and any insertions (dashes) as well.
*/

void write_maf_record(ostream & out, const std::string & src, size_t start, size_t size, char strand, size_t srcSize, const string & alignment_part) {
	out << "s " << src << " " << start << " " << size << " " << strand << " " << srcSize << " " << alignment_part << std::endl;
}


/*
	write half a pairwise alignment to a maf file 
	(this write the s-line only)

*/
void write_maf_record(ostream & out, const all_data & data, const pw_alignment & al, size_t reference) {
	assert(reference==0 || reference==1);
	// get the requested half alignment
	size_t print_seq;
	size_t print_start;
	size_t print_end;
	std::string print_al;
	// TODO print size not end
	if(reference==0) {
		print_seq = al.getreference1();
		print_start = al.getbegin1();
		print_end = al.getend1();
		print_al = al.get_al_ref1();
	} else {
		print_seq = al.getreference2();
		print_start = al.getbegin2();
		print_end = al.getend2();
		print_al = al.get_al_ref2();
	}
	// translate reference sequence index number to corresponding long name (accession:contig)
	size_t acc = data.accNumber(print_seq);
	std::string accname = data.get_acc(acc);
	std::string seqname = data.get_seq_name(print_seq);
	std::stringstream write_longname;
	write_longname << accname << ':' << seqname;
	char strand = '+';
	size_t size = data.get_seq_size(print_seq);
	size_t al_on_ref_size = print_end - print_start + 1;
	if(print_end < print_start) {
		al_on_ref_size = print_start - print_end + 1;
		strand = '-';
		print_start = size - print_start -1;


		
	} 
	write_maf_record(out, write_longname.str(), print_start, al_on_ref_size, strand, size, print_al);
	

}

/*
	compute star phylogeny multiple sequence alignment
	modifies input alignments by inserting gaps
*/   

void msa_star_alignment(const std::string & center, std::vector<pw_alignment> & alignments) {
	std::vector<std::string> cparts;
	strsep(center, ":", cparts);
//	size_t center_dir = atoi(cparts.at(0).c_str());
	size_t center_ref = atoi(cparts.at(0).c_str());
	size_t center_left = atoi(cparts.at(1).c_str());
	size_t center_length = 0; // length of the cluster center on its reference
	std::vector<size_t> gaps_before;

	for(size_t i=0; i<alignments.size(); ++i) {
		size_t ref1 = alignments.at(i).getreference1();
		size_t ref2 = alignments.at(i).getreference2();
		size_t left1, right1, left2, right2;
		alignments.at(i).get_lr1(left1, right1);
		alignments.at(i).get_lr2(left2, right2);
		size_t al_dir1;
		size_t al_dir2;
	//	if(alignments.at(i).getbegin1()<alignments.at(i).getend1()){
	//		al_dir1 = 0;
	//	}else{
	//		al_dir1 = 1;
	//	}
	//	if(alignments.at(i).getbegin2()<alignments.at(i).getend2()){
	//		al_dir2 = 0;
	//	}else{
	//		al_dir2 = 1;
	//	}
	 //	std::cout << " al " << i << std::endl;
	//	alignments.at(i).print();
	//	std::cout << std::endl;
		// find cluster center and make sure that each alignment goes to cluster center
		// and all cluster centers are identical on the cluster center reference
		// then find the minimal number of gaps before each cluster center reference position
		if(ref1 == center_ref && left1 == center_left) {
			if(center_length==0) {
				center_length = right1 - left1 + 1;
				gaps_before = std::vector<size_t>(center_length + 1, 0);
			} else {
		//		std::cout << " cl " << center_length << " l1 " << left1 << " r1 " << right1 << std::endl;
				assert(center_length == right1 - left1 +1);
			}
			size_t center_ref_pos=0; // relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0;
			// now go over alignment
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c1 == '-') {
					gapcounter++;
				} else {
					if(gaps_before.at(center_ref_pos) < gapcounter) {
						gaps_before.at(center_ref_pos) = gapcounter;
					}	
					gapcounter = 0;
					center_ref_pos ++;
				}		
			}
			// gaps at end of alignment
			if(gaps_before.at(center_ref_pos) < gapcounter) {
				gaps_before.at(center_ref_pos) = gapcounter;
			}
		
		} else {
			assert(ref2 == center_ref && left2 == center_left);
			if(center_length==0) {
				center_length = right2 - left2 + 1;
				gaps_before = std::vector<size_t>(center_length + 1, 0);
			} else {
			//	std::cout << " cl " << center_length << " l2 " << left2 << " r2 " << right2 << std::endl;

				assert(center_length == right2 - left2 +1);
			}
			size_t center_ref_pos=0; // relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0;
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c2 == '-') {
					gapcounter++;
				} else {
					if(gaps_before.at(center_ref_pos) < gapcounter) {
						gaps_before.at(center_ref_pos) = gapcounter;
					}	
					gapcounter = 0;
					center_ref_pos ++;
				}
			}
			// gaps at end of alignment
			if(gaps_before.at(center_ref_pos) < gapcounter) {
				gaps_before.at(center_ref_pos) = gapcounter;
			}
		} // if ref	
	} // for i 


	// now we have inserted a lot of gaps in to gaps_before
	// all that gaps are inserted into all alignments to create a multiple alignment where all sequences have the same length
	for(size_t i=0; i<alignments.size(); ++i){
		size_t ref1 = alignments.at(i).getreference1();
		size_t ref2 = alignments.at(i).getreference2();
		size_t left1, right1, left2, right2;
		alignments.at(i).get_lr1(left1, right1);
		alignments.at(i).get_lr2(left2, right2);
	//	size_t al_dir1;
	//	size_t al_dir2;
	//	if(alignments.at(i).getbegin1()<alignments.at(i).getend1()){
	//		al_dir1 = 0;
	//	}else{
	//		al_dir1 = 1;
	//	}
	//	if(alignments.at(i).getbegin2()<alignments.at(i).getend2()){
	//		al_dir2 = 0;
	//	}else{
	//		al_dir2 = 1;
	//	}

		// for writing new al std::strings
		std::stringstream centerref_al;
		std::stringstream otherref_al;

		if(ref1 == center_ref && left1 == center_left) {
					
			size_t center_ref_pos=0; // we are before that position relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0; // how many gaps did we already see on cluster center ref
			// now go over alignment
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c1 == '-') {
					gapcounter++;
					centerref_al << c1;
					otherref_al << c2;
				} else {
					// add gaps before next center ref symbol
					for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
						centerref_al << '-';
						otherref_al << '-';
					}
					centerref_al << c1;
					otherref_al << c2;
					gapcounter = 0;
					center_ref_pos ++;
				}
			}
			for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
				centerref_al << '-';
				otherref_al << '-';
			}
			pw_alignment newal(centerref_al.str(), otherref_al.str(), 
					alignments.at(i).getbegin1(), alignments.at(i).getbegin2(),
					alignments.at(i).getend1(), alignments.at(i).getend2(),
					alignments.at(i).getreference1(), alignments.at(i).getreference2());
			alignments.at(i) = newal;

		} else {
			assert(ref2 == center_ref && left2 == center_left);
		
			size_t center_ref_pos=0; // relative to cluster center start (even if cluster center is backwards, we count forwards) 
			size_t gapcounter = 0;
			for(size_t j=0; j<alignments.at(i).alignment_length(); ++j) {
				char c1, c2;
				alignments.at(i).alignment_col(j, c1, c2);
				if(c2 == '-') {
					gapcounter++;
					centerref_al << c2;
					otherref_al << c1;
				} else {
					// add gaps before next center ref symbol
					for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
						centerref_al << '-';
						otherref_al << '-';
					}
					centerref_al << c1;
					otherref_al << c2;
					gapcounter = 0;
					center_ref_pos ++;
				}
			
			}
			for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++){
				centerref_al << '-';
				otherref_al << '-';
			}
			// swap alignment, cluster center is now on reference 1
			pw_alignment newal(centerref_al.str(), otherref_al.str(), 
					alignments.at(i).getbegin2(), alignments.at(i).getbegin1(),
					alignments.at(i).getend2(), alignments.at(i).getend1(),
					alignments.at(i).getreference2(), alignments.at(i).getreference1());
			alignments.at(i) = newal;
		} // else: center on ref 2 of al i
	} // for alignments

// test code
	size_t length = 0;
	for(size_t i=0; i<alignments.size(); ++i) {
		if(length == 0) {
			length = alignments.at(i).alignment_length();
		}
		assert(length == alignments.at(i).alignment_length());
	}

}
void msa_star_al_with_long_centers(const std::vector<std::string> & longCenter, std::vector<pw_alignment> & alignment){//center is always on the second reference! 
	size_t center_length = 0; // length of the cluster center on its reference
	std::vector<size_t> gaps_before;
	std::cout<< "al size "<< alignment.size() <<std::endl;
	for(size_t i=0; i<alignment.size(); ++i) {
		size_t ref1 = alignment.at(i).getreference1();
		size_t ref2 = alignment.at(i).getreference2();
		size_t left1, right1, left2, right2;
		alignment.at(i).get_lr1(left1, right1);
		alignment.at(i).get_lr2(left2, right2);
		std::cout<< "ref 2 "<< ref2<<std::endl;
		pw_alignment p1 = alignment.at(i);
	//	p1.print();
		if(center_length==0) {
			center_length = right2 - left2 + 1;
			gaps_before = std::vector<size_t>(center_length + 1, 0);
		} else {
			std::cout << " cl " << center_length << " l2 " << left2 << " r2 " << right2 << std::endl;
			assert(center_length == right2 - left2 +1);
		}
		size_t center_ref_pos=0; // relative to cluster center start (there are both forwards and backwards centers)
		size_t gapcounter = 0;
		// now go over alignment
		std::cout<< alignment.at(i).alignment_length() <<" "<< center_length+1 << std::endl;
		for(size_t j=0; j<p1.alignment_length(); j++) {
			char c1, c2;
			p1.alignment_col(j, c1, c2);
			if(c2 == '-') {//If there is a gap on the center
				gapcounter++;
			//	std::cout << "gap counter "<< gapcounter <<std::endl;


			} 
			else {
			//	std::cout << center_ref_pos <<std::endl;

				if(gaps_before.at(center_ref_pos) < gapcounter) {
					gaps_before.at(center_ref_pos) = gapcounter;
				}	
				gapcounter = 0;
				center_ref_pos ++;
				
			}		
		}
		std::cout << "done!"<<std::endl;
		// gaps at end of alignment
		if(gaps_before.at(center_ref_pos) < gapcounter) {
			gaps_before.at(center_ref_pos) = gapcounter;
		}
	}
	for(size_t i=0; i<alignment.size(); ++i) {
		size_t ref1 = alignment.at(i).getreference1();
		size_t ref2 = alignment.at(i).getreference2();
		size_t left1, right1, left2, right2;
		alignment.at(i).get_lr1(left1, right1);
		alignment.at(i).get_lr2(left2, right2);
		// for writing new al std::strings
		std::stringstream centerref_al;
		std::stringstream otherref_al;
		size_t center_ref_pos=0; // relative to cluster center start (even if cluster center is backwards, we count forwards) 
		size_t gapcounter = 0;
		for(size_t j=0; j<alignment.at(i).alignment_length(); ++j) {
			char c1, c2;
			alignment.at(i).alignment_col(j, c1, c2);
			if(c2 == '-') {
				gapcounter++;
				centerref_al << c2;
				otherref_al << c1;
			} else {
				// add gaps before next center ref symbol
				for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
					centerref_al << '-';
					otherref_al << '-';
				}
				centerref_al << c2;
				otherref_al << c1;
				gapcounter = 0;
				center_ref_pos ++;
			}	
		}
		for(size_t k=gapcounter; k<gaps_before.at(center_ref_pos); k++) {
			centerref_al << '-';
			otherref_al << '-';
		}	
		pw_alignment newal(otherref_al.str(),centerref_al.str(), 
		alignment.at(i).getbegin1(), alignment.at(i).getbegin2(),
		alignment.at(i).getend1(), alignment.at(i).getend2(),
		alignment.at(i).getreference1(), alignment.at(i).getreference2());
		alignment.at(i) = newal;
	}
}

	/*
		write the graph as multiple alignment data structure	
	*/
	std::map<size_t, std::set<size_t> > adjacencies;
	std::map<int, std::set<int> > adjacencies_with_dir;
void write_graph_fasta(const std::string & graphout, const std::map<string, std::vector<pw_alignment> > & cluster_result_al, all_data & data){
	std::ofstream gout(graphout.c_str());
	std::vector<std::string> centers;
	for(std::map<std::string, std::vector<pw_alignment> >::const_iterator it = cluster_result_al.begin(); it!=cluster_result_al.end(); ++it) {
		std::vector<pw_alignment> als = it->second;
		if(als.size()>0) {
			std::string center = it->first;
			std::vector<std::string> split;
			strsep(center, ":" , split);
			size_t cent_ref = atoi(split.at(0).c_str());
			size_t cent_left = atoi(split.at(1).c_str());
			const pw_alignment & p = it->second.at(0);
			size_t left1;
			size_t right1;
			size_t left2;
			size_t right2;
			p.get_lr1(left1,right1);
			p.get_lr2(left2,right2);
			if(cent_ref == p.getreference1()&& cent_left == left1){
				std::string temp;
				for(size_t j=0; j<p.alignment_length(); ++j) {
					char c1, c2;
					p.alignment_col(j, c1, c2);
					if(c1 != '-'){
						temp+=c1;
					}
				}
				centers.push_back(temp);
			}else{
				assert(cent_ref == p.getreference2()&& cent_left == left2);
				std::string temp;
				for(size_t j=0; j<p.alignment_length(); ++j) {
					char c1, c2;
					p.alignment_col(j, c1, c2);
					if(c2 != '-'){
						temp+=c2;
					}
				}
				centers.push_back(temp);

			}	
		}
	}
	
	data.make_fasta(gout,centers);
}
void write_graph_maf(const std::string & graphout, const std::map<string, std::vector<pw_alignment> > & cluster_result_al, const all_data & data) {	
	std::ofstream gout(graphout.c_str());
	gout << "##maf version=1 scoring=probability" << std::endl;
	gout << "# graph version " << VERSION << std::endl;
	// TODO add more metadata
	gout << "# Coordinates are 0-based.  For - strand matches, coordinates" << std::endl;
	gout << "# in the reverse complement of the 2nd sequence are used." << std::endl;
	gout << "#" << std::endl;
	gout << "# name start alnSize strand seqSize alignment" << std::endl;
	gout << "#" << std::endl;
	size_t cluster_number = 1; 
	for(std::map<std::string, std::vector<pw_alignment> >::const_iterator it = cluster_result_al.begin(); it!=cluster_result_al.end(); ++it) {
		std::string center = it->first;
		std::vector<pw_alignment> als = it->second;
	//	std::cout << "Center: "<< center << " als " << als.size() << std::endl;
		msa_star_alignment(center, als);
		if(als.size()>0) {
			gout << "# cluster " << cluster_number << "( ";
			std::map<int, std::set<int> >::const_iterator adj = adjacencies_with_dir.find(cluster_number);//XXX it is not correct anymore. now the non aligned regions were added and the indeceds are different.
			if(adj != adjacencies_with_dir.end()){
				gout << cluster_number << " -> ";
				for(std::set<int>::iterator it1 = adj->second.begin(); it1 != adj->second.end();it1++){
					gout<< *it1 << ", ";
				}
			}
			//check for the reverse!
			int rev_cluster_number = -1*(cluster_number);
			std::map<int, std::set<int> >::const_iterator rev_adj = adjacencies_with_dir.find(rev_cluster_number);
			if(rev_adj != adjacencies_with_dir.end()){
				gout << rev_cluster_number << " -> ";
				for(std::set<int>::iterator it1 = rev_adj->second.begin(); it1 != rev_adj->second.end();it1++){
					gout<< *it1 << ", ";
				}
			}
			gout << ")" << std::endl;
			gout << "a score=1" << std::endl; 

		

			write_maf_record(gout, data, als.at(0), 0);
			for(size_t i=0; i<als.size(); ++i) {
				write_maf_record(gout, data, als.at(i), 1);
			}	
			cluster_number++;
		}
	}
	gout.close();

}
void write_graph_maf_with_long_centers(const std::string & graphout ,const std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers,const all_data & data){
	std::ofstream gout(graphout.c_str());
	gout << "##maf version=1 scoring=probability" << std::endl;
	gout << "# graph version " << VERSION << std::endl;
	gout << "# Coordinates are 0-based.  For - strand matches, coordinates" << std::endl;
	gout << "# in the reverse complement of the 2nd sequence are used." << std::endl;
	gout << "#" << std::endl;
	gout << "# name start alnSize strand seqSize alignment" << std::endl;
	gout << "#" << std::endl;
	size_t cluster_number = 0;
	for(std::map<std::vector<std::string> , std::vector<pw_alignment> >::const_iterator it = new_centers.begin(); it != new_centers.end(); it++){//It includes all the centers
		std::vector<std::string> longCenter = it->first;
		std::vector<pw_alignment> als = it->second;
	//	size_t sequence_number = data.numSequences(); 
		for(size_t i =0; i < longCenter.size();i++){
			std::cout << longCenter.at(i)<<std::endl;
		}
		std::cout << "cent size " << longCenter.size() << " cluster_number " << cluster_number << std::endl;
		msa_star_al_with_long_centers(longCenter,als);
		if(als.size()>0){
			gout << "# cluster " << cluster_number << std::endl;
			gout << "a score=1" << std::endl;
			//TODO Add the center string itself
			for(size_t i=0; i<als.size(); ++i) {
				write_maf_record(gout, data, als.at(i), 0); 
			}	
			cluster_number++;
		}
	}
	std::cout<< "cluster number "<< cluster_number << " "<< new_centers.size()<< std::endl;
}
std::string get_reverse(const std::string & seq){
	stringstream temp;
	std::vector<std::string> center_parts;
	strsep(seq, ":" , center_parts);
	unsigned int dir = atoi(center_parts.at(0).c_str());
	unsigned int ref = atoi(center_parts.at(1).c_str());
	unsigned int left = atoi(center_parts.at(2).c_str());
	if(dir == 0){
		temp<< 1 <<":"<<ref<<":"<<left;
	}else{
		temp << 0 << ":"<<ref<<":"<<left;
	}
	return temp.str();
}
int do_mc_model(int argc, char * argv[]) {
	typedef mc_model use_model;
	typedef overlap_interval_tree overlap_type;
//	typedef dynamic_mc_model use_model;
//	typedef clustering<use_model> use_clustering;
//	typedef overlap overlap_type;

	typedef initial_alignment_set<use_model,overlap_type> use_ias;
	typedef compute_cc<overlap_type> compute_cc_type;

	typedef affpro_clusters<use_model,overlap_type> use_affpro;
	typedef encoder<use_model> use_encode;
	if(argc < 6) {
		usage();
		cerr << "Program: model" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* output maf file for the graph" << std::endl;
		cerr << "* output binary compressed file (use 'noencode' to skip encoding step)" << std::endl;
		cerr << "* number of threads to use (optional, default 10)" << std::endl;
	}

// TODO read maf or sam

	clock_t read_data_time = clock();
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string graphout(argv[4]);
	size_t num_threads = 1;
	if(argc == 7) {
		num_threads = atoi(argv[6]);
	}
	std::string encoding_out(argv[5]);
	std::ofstream outs(encoding_out.c_str(),std::ofstream::binary);



// Read all data
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	overlap_type ol(data);
	wrapper wrap;
	test_encoder test;
	read_data_time = clock() - read_data_time;


//	size_t all_length =0;
//	for(size_t i=0; i<data.numSequences(); i++) {
//		all_length+= data.getSequence(i).length();
//	}
//	std::cout<< "all sequence bases " << all_length << std::endl;

// Find connected components of alignments with some overlap //Moved it after training. this is the old place of making ccs.
/*	clock_t initial_cc_time = clock();
	compute_cc cccs(data);
	std::vector<std::set< pw_alignment , compare_pw_alignment> > ccs;
	cccs.compute(ccs); //fill in ccs, order by size(Notice that they are just connected to each other if they have overlap!)
	initial_cc_time = clock() - initial_cc_time;
//	std::cout << " Found " << ccs.size() << " connected components" << std::endl;
	for(size_t i=0; i<ccs.size()-1; i++) {//A test function that chckes overlap between components. There shouldn't be any overlap between them.
		std::cout << "Connected component "<< i << " contains " << ccs.at(i).size() << " alignments" << std::endl;
		for(std::set< pw_alignment , compare_pw_alignment>::iterator it = ccs.at(i).begin();it != ccs.at(i).end();it++){
			const pw_alignment & p = *it;
		//	p.print();
			for(size_t j = i+1; j < ccs.size();j++){
				std::set< pw_alignment , compare_pw_alignment> ccs_j = ccs.at(j);
				ol.test_no_overlap_between_ccs(p, ccs_j);
			}
		}
	}*/
// Train the model on all data
	clock_t train_models_time = clock();
//	dynamic_mc_model dy_model(data);
//	dy_model.train(outs);
	use_model m(data);
	m.train(outs);
	size_t total = 0;
	train_models_time = clock() - train_models_time;
//	double sum = 0;
//	double c1;
//	double c2;
//	double m1;
//	double m2;
//	size_t len=0;
/*	for(size_t i =0; i <  data.numSequences(); i++){
		const dnastd::string & sequence = data.getSequence(i);
		size_t length = sequence.length();
		len +=length;
	}*/
// std::cout<<"length: " << len << std::endl;
//	for(size_t i = 0 ; i< data.numAlignments() ; i++){
//		const pw_alignment & p = data.getAlignment(i);
	//	size_t ref1 = p.getreference1();
	//	size_t ref2 = p.getreference2();
	//	const dnastd::string & sequence1 = data.getSequence(ref1);
	//	const dnastd::string & sequence2 = data.getSequence(ref2);
	//	size_t length1 = sequence1.length();
	//	size_t length2 = sequence2.length();
//		std::cout<< "als costs: " <<std::endl;
//		p.print();
//		m.cost_function(p, c1,c2,m1,m2);
//		std::cout << "c1 " << c1 << " c2 " << c2 << " m1 " << m1 << " m2 "<< m2 <<std::endl;
	//	sum +=c1+c2;
	//	len +=length1 + length2;
//	}
//	std::cout<< "sum is "<<sum <<std::endl;
//	std::cout<< "len times two: "<< len*2 <<std::endl;
	// base cost to use an alignment (information need for its adress)
	double cluster_base_cost = log2(data.numAlignments());//TODO +200
	std::cout << " base cost " << cluster_base_cost << std::endl;
// Find connected components of alignments with some overlap // Moved it here after training to be able to remove small als first and then connect them to each other.
	std::set<const pw_alignment*, compare_pointer_pw_alignment> al_with_pos_gain; 
	for(size_t i =0; i < data.numAlignments();i++){
		const pw_alignment & al = data.getAlignment(i);
		double g1 ,g2;
		m.gain_function(al,g1,g2);
		double av_gain = g1+g2/2 - cluster_base_cost;
		if(av_gain > 0){
			al_with_pos_gain.insert(&al);
		//	std::cout << "al " << al <<std::endl;
		}
	}
//	std::cout << "al_with_pos_gain size " << al_with_pos_gain.size() << std::endl;
	clock_t initial_cc_time = clock();
	compute_cc_type cccs(al_with_pos_gain, data.numSequences(),num_threads);
	std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> > ccs; 
	cccs.compute(ccs); //fill in ccs, order by size(Notice that they are just connected to each other if they have overlap!)
	initial_cc_time = clock() - initial_cc_time;
//	std::cout << " Found " << ccs.size() << " connected components" << std::endl;
	for(size_t i=0; i<ccs.size()-1; i++) {//A test function that chckes overlap between components. There shouldn't be any overlap between them.
		std::cout << "Connected component "<< i << " contains " << ccs.at(i).size() << " alignments" << std::endl;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.at(i).begin();it != ccs.at(i).end();it++){
			const pw_alignment * p = *it;
		//	p.print();
		//	for(size_t j = i+1; j < ccs.size();j++){
		//		std::set< const pw_alignment* , compare_pointer_pw_alignment> ccs_j = ccs.at(j);
			//	ol.test_no_overlap_between_ccs(*p, ccs_j);//XXX this is just a test function, comment it the first run!
			//	data.checkAlignmentRange(*p);//XXX this is just a test function, comment it after the first run!
		//	}
		}
	}
	std::vector<overlap_type> cc_overlap(ccs.size(), overlap_type(data));// An overlap class object is needed for each connected component, thus we cannot use ol
// Select an initial alignment std::set for each connected component (in parallel)
	
	std::map<std::string, std::vector<string> > global_results;//for each center returns all its cluster members
	std::map<std::string, std::vector<pw_alignment> > alignments_in_a_cluster;//string ---> center of a cluster, vector ---> alignments with that center
	std::map<vector<std::string>, std::vector<pw_alignment> > new_centers;//Equivalent to alignments_in_a_cluster for long centers.
	std::vector<std::map<size_t, std::vector<std::string> > >centersPositionOnASeq(data.numSequences());//all the long centers of each sequence and their position
	std::vector<std::multimap<size_t, std::string> > centerOnSequence(data.numSequences());//all the centers on a seq and their position(It is used when we use long centers)
	vector<vector<std::string> > long_centers;
	size_t num_clusters = 0;
	size_t num_cluster_members_al = 0;
	size_t num_cluster_seq = 0;
	size_t num_cluster_inputs_al = 0;
#if TIMING
	clock_t ias_time = 0;
	clock_t test_function_time = 0;
	clock_t second_cc_time = 0;
	clock_t ap_time = 0;

#endif

	double all_gain = 0;
	double all_gain_in_ias = 0;
	size_t all_used_input_alignments = 0;
	size_t all_npo_alignments = 0;
	std::cout << "ccs size: " << ccs.size() << std::endl;
	#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for(size_t i=0; i<ccs.size(); ++i) {
#pragma omp critical(print)
{
		std::cout << " on initial CC " << i << " size " << ccs.at(i).size() << std::endl;
}
		clock_t ias_time_local = clock();
		std::set<const pw_alignment*, compare_pointer_pw_alignment> & cc = ccs.at(i);//Before cutting partial overlaps
		use_ias ias(data,cc, m, cluster_base_cost,num_threads);//Bottle neck!!!!
		ias.compute(cc_overlap.at(i));//cc_overlap.at(i) includes all the non_overlapped alignments that retains after cutting
		ias_time_local = clock() - ias_time_local;

		clock_t test_time_local = clock();
//		cc_overlap.at(i).test_partial_overlap(); // TODO remove slow test function
		test_time_local = clock() - test_time_local;
#pragma omp critical(gain_statistics)
{
		all_gain+=ias.get_max_gain();
		all_gain_in_ias+=ias.get_result_gain();
		all_used_input_alignments+=ias.get_used_alignments();
		all_npo_alignments += cc_overlap.at(i).size();
}
		//it is the number of als after cutting partial overlaps:
	//	std::cout << " number of non overlapped alignments " << cc_overlap.at(i).size() << std::endl;
	//	std::set<pw_alignment, compare_pw_alignment> al_in_overlap = cc_overlap.at(i).get_all();
	//	for(std::set<pw_alignment,compare_pw_alignment>::iterator it = al_in_overlap.begin(); it != al_in_overlap.end();it++){
	//		const pw_alignment & p = *it;
		//	data.checkAlignmentRange(p);//TODO this is just a test function, comment it for the next runs!
	
	//	}
		
		// TODO this can be done a lot faster because there is no partial overlap here
		clock_t second_cc_time_local = clock();
		std::vector< std::set<const pw_alignment* , compare_pointer_pw_alignment> > cc_cluster_in; //has shorter length than original sequence pieces
		compute_cc_type occ(cc_overlap.at(i), data.numSequences(), 1);
		occ.compute(cc_cluster_in);//makes different components with related pieces of the alignments with no partial overlap, they are sroted by the size of components

		second_cc_time_local = clock() - second_cc_time_local;
#if TIMING
#pragma omp critical(time)
{
		ias_time += ias_time_local;
		test_function_time += test_time_local;
		second_cc_time += second_cc_time_local;
}
#endif
		std::cout<< "cc_cluster_in size: "<< cc_cluster_in.size()<<std::endl;
	//	std::cout << "cc " << i << ": from " << cc.size() << " original als with total gain " << ias.get_max_gain() << " we made " << cc_overlap.at(i).size() << " pieces with total gain " << ias.get_result_gain() <<  std::endl; 
	//	std::cout << " components for clustering: " << std::endl;
		map<string, vector<pw_alignment> >al_of_a_ccs;//It is filled with the original centers and is used for creating long centers & their als.(string ---> center, vector<pw_al> ---> alignments of a cluster)
		vector<vector<std::string> > local_long_centers;
		for(size_t j=0; j<cc_cluster_in.size(); ++j) {

			std::cout << " run affpro at " << j << " on " << cc_cluster_in.at(j).size()<<std::endl;
		//	data.numAcc();

		
		//	std::cout << " in al std::set size " << cc_cluster_in.at(j).size() << std::endl;
			size_t nums=0;
			for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator ait=cc_cluster_in.at(j).begin(); ait!=cc_cluster_in.at(j).end(); ++ait) {
			//	std::cout << " nums " << nums << std::endl;
				const pw_alignment * cc_al = *ait;
			//	data.checkAlignmentRange(*cc_al);//TODO this is just a test function, comment it after the first run!
				nums++;
		//		cc_al.print();
				std::cout << std::endl;
				double c1;
				double c2;
				double m1;
				double m2;
			//	m.cost_function(*cc_al, c1, c2, m1,m2);
			//	std::cout << "c1 " << c1 << " c2 " << c2 << " m1 " << m1 << " m2 "<< m2 <<std::endl;
			//	if(cc_al.getreference1()== 0 && cc_al.getreference2() == 1){
			//		cc_cluster_in.at(j).erase(ait);
			//	}
			}
		

			clock_t ap_time_local = clock();
			use_affpro uaf(cc_cluster_in.at(j), m, cluster_base_cost);
			std::map<std::string, std::vector<string> >cluster_result;
			std::map<std::string, std::vector<pw_alignment> > local_al_in_a_cluster;

			uaf.run(cluster_result,local_al_in_a_cluster);
			ap_time_local = clock() - ap_time_local;
#if TIMING
#pragma omp critical(aptime) 
{
			ap_time += ap_time_local;
}

#endif



			size_t counter =0;
			for(std::map<std::string,std::vector<string> >::iterator it=cluster_result.begin();it !=cluster_result.end();it++){
				counter++;
				std::string center = it->first;
			//	std::cout << "center is "<<center << std::endl;
				std::map<std::string, std::vector<pw_alignment> >::iterator it1 = local_al_in_a_cluster.find(center);
				if(it1 == local_al_in_a_cluster.end()){
					local_al_in_a_cluster.insert(make_pair(center, std::vector<pw_alignment>()));
					it1 = local_al_in_a_cluster.find(center);
				}
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int dir = atoi(center_parts.at(0).c_str());
				unsigned int ref = atoi(center_parts.at(1).c_str());
				unsigned int left = atoi(center_parts.at(2).c_str());
				for(size_t i =0; i < it->second.size();i++){
				//	std::cout<< "member is "<<it->second.at(i)<<std::endl;
					std::vector<std::string> member_parts;
					strsep(it->second.at(i), ":" , member_parts);
					unsigned int mem_dir = atoi(member_parts.at(0).c_str());
					unsigned int mem_ref = atoi(member_parts.at(1).c_str());
					unsigned int mem_left = atoi(member_parts.at(2).c_str());
					// look at all input alignments and determine which cluster it belongs to
					for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = cc_cluster_in.at(j).begin(); it2!=cc_cluster_in.at(j).end(); ++it2){
						const pw_alignment * al = *it2;
						unsigned int al_dir1;
						if(al->getbegin1() < al->getend1()){
							al_dir1 =0;
						}else{
							al_dir1 = 1;
						}
						unsigned int al_dir2;
						if(al->getbegin2() < al->getend2()){
							al_dir2 =0;
						}else{
							al_dir2 = 1;
						}
						size_t ref1 = al->getreference1();
						size_t ref2 = al->getreference2();
						size_t left1;
						size_t left2;
						size_t right1;
						size_t right2;
						al->get_lr1(left1,right1);
						al->get_lr2(left2,right2);
						if(al_dir1 == dir && ref1 == ref && left1 == left && al_dir2 == mem_dir && ref2 == mem_ref && left2 == mem_left){
							it1->second.push_back(*al);
						}
						if(al_dir2 == dir && ref2 == ref && left2 == left && al_dir1==mem_dir&& ref1 == mem_ref && left1 == mem_left ){
							it1->second.push_back(*al);
						}
					} // for cluster_in set
				}
			}
#pragma omp critical 
{
			num_clusters += cluster_result.size();
			num_cluster_inputs_al += cc_cluster_in.at(j).size();
			std::cout << "cluster result "<< cluster_result.size() << " local al size "<< local_al_in_a_cluster.size() << std::endl;
			//here we check for the reverse center and if that is the case we just keep the one with the higher gain //TODO can it go to the pw_alignment?
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it = local_al_in_a_cluster.begin(); it!=local_al_in_a_cluster.end(); ++it){
				if(it->second.size() != 0){
					al_of_a_ccs.insert(*it);
				}
				std::string reverse_center = get_reverse(it->first);
				std::map<std::string,std::vector<string> >::iterator result_it = cluster_result.find(it->first);
				assert(result_it != cluster_result.end());
				std::map<std::string,std::vector<string> >::iterator rev_result_it = global_results.find(reverse_center);
				std::map<std::string, std::vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(reverse_center);//The gain value of its member is checked and the one with higher gain is kept.
				if(it1 == alignments_in_a_cluster.end()){
				//	if(it->second.size() != 0){
				//		al_of_a_ccs.insert(*it);
				//	}
					alignments_in_a_cluster.insert(*it);// Globally saves all the als.	
					num_cluster_members_al += 1 + it->second.size();
					num_cluster_seq+= 1 + result_it->second.size();
					global_results.insert(*result_it);
				}
				else if(it1->second.size() == 0){
					alignments_in_a_cluster.insert(*it);
					alignments_in_a_cluster.erase(it1);
					num_cluster_members_al += it->second.size();
					num_cluster_seq+= result_it->second.size();
					global_results.insert(*result_it);
					global_results.erase(rev_result_it);
				}
				else{
				//	std::map<std::string, std::vector<pw_alignment> >::iterator it2 = al_of_a_ccs.find(reverse_center);
				//	assert(it2 != al_of_a_ccs.end());
					double sum = 0;
					double reverse_sum = 0;
					for(size_t i =0; i < it->second.size();i++){//forward
						pw_alignment p = it->second.at(i);
						std::string center = it->first;
						sum += m.get_the_gain(p,center);
					}
					for(size_t i =0; i < it1->second.size();i++){//reverse
						pw_alignment p = it1->second.at(i);
						std::string center = it1->first;
						reverse_sum +=m.get_the_gain(p,center);
					}
					if(reverse_sum > sum){
						for(size_t i =0; i < it->second.size();i++){
							pw_alignment & al = it->second.at(i);
							pw_alignment p;
							al.get_reverse(p);
							it1->second.push_back(p);
						//	it2->second.push_back(p);
						}
						num_cluster_members_al += it->second.size();
						for(size_t i =0; i < result_it->second.size();i++){
							std::string mem = result_it->second.at(i);
							std::string rev_mem = get_reverse(mem);
							rev_result_it->second.push_back(rev_mem);
						}
					}else if(reverse_sum< sum){
						alignments_in_a_cluster.insert(*it);
						global_results.insert(*result_it);
					//	if(it->second.size() != 0){
					//		al_of_a_ccs.insert(*it);
					//	}
						std::map<std::string , std::vector<std::string> >::iterator gl_result = global_results.find(it->first);
						assert(gl_result != global_results.end());
					//	std::map<std::string, std::vector<pw_alignment> >::iterator it4 = al_of_a_ccs.find(it->first);
					//	assert(it4 != al_of_a_ccs.end());
						std::map<std::string, std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(it->first);
						assert(it3 != alignments_in_a_cluster.end());
						for(size_t i =0; i < it1->second.size();i++){
							pw_alignment & al = it1->second.at(i);
							pw_alignment p;
							al.get_reverse(p);
							it3->second.push_back(p);
						//	it4->second.push_back(p);
						}
						for(size_t i =0 ; i < rev_result_it->second.size();i++){
							std::string mem = rev_result_it->second.at(i);
							std::string rev_mem = get_reverse(mem);
							gl_result->second.push_back(rev_mem);
						}
						global_results.erase(rev_result_it);
					//	al_of_a_ccs.erase(it2);
						alignments_in_a_cluster.erase(it1);
						num_cluster_members_al += it->second.size();	
						num_cluster_seq+= result_it->second.size();			
					}else{//if their gain valuea are the same we just simply count the number of als and keep the one with the higher number of als
						size_t al_number = it->second.size();
						size_t reverse_al_number = it1->second.size();
						if(al_number >= reverse_al_number){
							alignments_in_a_cluster.insert(*it);
							global_results.insert(*result_it);
					//		al_of_a_ccs.insert(*it);
					//		std::map<std::string, std::vector<pw_alignment> >::iterator it4 = al_of_a_ccs.find(it->first);
					//		assert(it4 != al_of_a_ccs.end());
							std::map<std::string, std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(it->first);
							assert(it3 != alignments_in_a_cluster.end());
							for(size_t i =0; i < it1->second.size();i++){
								pw_alignment & al = it1->second.at(i);
								pw_alignment p;
								al.get_reverse(p);
								it3->second.push_back(p);
							//	it4->second.push_back(p);
							}
							std::map<std::string ,  std::vector<std::string> >::iterator gl_result = global_results.find(it->first);
							for(size_t i =0; i < rev_result_it->second.size();i++){
								std::string mem = rev_result_it->second.at(i);
								std::string rev_mem = get_reverse(mem);
								gl_result->second.push_back(rev_mem);
							}
							global_results.erase(rev_result_it);
							alignments_in_a_cluster.erase(it1);
						//	al_of_a_ccs.erase(it2);
							num_cluster_members_al += it->second.size();	
							num_cluster_seq+= result_it->second.size();			
						}else{
							for(size_t i =0; i < it->second.size();i++){
								pw_alignment & al = it->second.at(i);
								pw_alignment p;
								al.get_reverse(p);
								it1->second.push_back(p);
							//	it2->second.push_back(p);
							}
							num_cluster_members_al += it->second.size();
							for(size_t i = 0; i < result_it->second.size();i++){
								std::string mem = result_it->second.at(i);
								std::string rev_mem = get_reverse(mem);
								rev_result_it->second.push_back(rev_mem);
							}
							num_cluster_seq+= result_it->second.size();
						}
					}
				}
			}
			std::cout << alignments_in_a_cluster.size()<<std::endl;
			std::cout << "num_cluster " << num_cluster_members_al <<std::endl;
			std::cout << global_results.size()<<std::endl;
			std::cout << "num_cluster_seq "<<num_cluster_seq<<std::endl;			
}
//UP to here we already removed centers that might be different with each other only by their direction. Thus from here on there is no need to check for the direction. Notice that there is also no single member with two different directions!


//#pragma omp critical//I can not do them in parallel anymore.
//{
/*			for(std::map<std::string,std::vector<string> >::iterator it=cluster_result.begin();it != cluster_result.end();it++){
				std::string center = it->first;
				num_cluster_seq+= 1 + it->second.size();
				std::map<std::string,std::vector<string> >::iterator it1 = global_results.find(center);
				if(it1 == global_results.end()){
					global_results.insert(make_pair(center, std::vector<std::string>()));
					it1 = global_results.find(center);
				}
				for(size_t i =0; i < it->second.size();i++){
					it1->second.push_back(it->second.at(i));
				}
			}
			std::cout << "num_cluster_seq "<<num_cluster_seq<<std::endl;*/
//}

		} // for cluster_in set  
		std::cout<< "end of cluster in! "<<std::endl;


		//Long centers are created:
		finding_centers centers(data);
		size_t num_seq = data.numSequences();
	//	suffixTree tree(num_seq, centers);
	//	std::vector<size_t> center_with_highest_gain;
	//	size_t highest_index = 0;
	//	tree.make_tree(center_with_highest_gain , highest_index);

//		suffix_tree tree(data,centers);
//		merging_centers merg(data, centers, tree);
	//From this point on orientation matters! It shoud be taken into account in counting the number of each center happening when we are looking for the long centers!
		centers.center_frequency(al_of_a_ccs,centerOnSequence);//It detects all the centers on each sequence and fill in centerOnSequence which contains all the cneters on each sequence.
	//	merg.adding_new_centers(local_long_centers,centersPositionOnASeq);//After merging a group of centers we iteratively create new suffixes and a new tree based on them, then recalculate the gains again. At the end those with gain>0 go to the long_centers and centersPositionOnASeq.
		// New centers (longer ones) are added to a separate data structure.
//		std::cout << "long centers: " << std::endl;
		for(size_t j  =0; j < local_long_centers.size();j++){
			long_centers.push_back(local_long_centers.at(j));
//			for(size_t k =0; k < local_long_centers.at(j).size();k++){
//				std::cout << local_long_centers.at(j).at(k)<< " ";
//			}
//			std::cout << " " <<std::endl;
		}
		std::cout << "al of a ccs size "<< al_of_a_ccs.size() <<std::endl;//The different orientations are kept in the 'al_of_a_ccs'
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = al_of_a_ccs.begin(); it != al_of_a_ccs.end();it++){
			std::cout << it->first << " " << it->second.size()<<std::endl;
		}
//		merg.create_alignment(local_long_centers,new_centers,al_of_a_ccs,centersPositionOnASeq,centerOnSequence);//Here we fill in the 'new_centers' with longer alignments in the way that center is always the second ref
		std::cout << "long alignments: " <<std::endl;
		for(size_t i =0; i < local_long_centers.size();i++){
			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.find(local_long_centers.at(i));
			assert(it != new_centers.end());
		//	for(size_t i =0; i < it->first.size();i++){
		//		std::cout << it->first.at(i) << " ";
		//	}
		//	std::cout << " " << std::endl;
		//	std::cout << "al size " << it->second.size() << std::endl;
			for(size_t k =0; k < it->second.size();k++){
				pw_alignment al = it->second.at(k);
		//		al.print();
			}
		}
	/*	for(std::set< const pw_alignment *, compare_pw_alignment>::iterator it=cc.begin();it!=cc.end();it++ ){
			const pw_alignment *al = *it;
			size_t index = al ->getreference1();

		}*/
//#pragma omp critical(print) 
//{
	//	std::cout << " initial CC " << i << " done " << std::endl << flush;
//}
	} // for connected components
	std::cout << "end of connected components"<<std::endl;
	std::cout<< alignments_in_a_cluster.size() << " " << long_centers.size() << std::endl;
	//Short centers are added to the long centers container. //TODO maybe i can move it pw_al XXX Note that second ref could be backward
	std::map<std::string, std::vector<pw_alignment> > intermediate;
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it != new_centers.end(); it++){
		std::vector<std::string> connected_center = it->first;
		for(size_t i =0; i < connected_center.size();i++){
			std::string center = connected_center.at(i);
			std::cout << "current center "<< center << std::endl;
			std::vector<std::string> cparts;
			strsep(center, ":", cparts);
			unsigned int center_dir = atoi(cparts.at(0).c_str());
			unsigned int center_ref = atoi(cparts.at(1).c_str());
			unsigned int center_left = atoi(cparts.at(2).c_str());
			std::map<std::string, std::vector<pw_alignment> >::iterator cent = alignments_in_a_cluster.find(center);// XXX Attention! Here only one direction is saved! (It is hard to see if there is something wrong here since there is no such case in ecoli dataset)
			if(cent == alignments_in_a_cluster.end()){
				if(center_dir == 0){
					center_dir = 1;
				}else{
					center_dir = 0;
				}
				stringstream temp;
				temp << center_dir<< ":"<< center_ref<<":"<<center_left;
				center = temp.str();
				std::cout << "CENT "<<center <<std::endl;
				cent = alignments_in_a_cluster.find(center);
				assert(cent != alignments_in_a_cluster.end());
			}
			std::cout << "number of als " << cent->second.size() <<std::endl;
			for(size_t j =0;j < cent->second.size();j++){//Go through the short alignments and check to see if they are part a long one
				pw_alignment p = cent->second.at(j); //Check if it happened on a long alignment!
				std::cout << "p is "<<std::endl;
			//	p.print();
				size_t l1,l2,r1,r2;
				size_t ref1 = p.getreference1();
				size_t ref2 = p.getreference2();
				std::cout << "ref1 "<< ref1 << " ref2 "<<ref2 <<std::endl;
				p.get_lr1(l1, r1);
				p.get_lr2(l2, r2);
				if(l1 == center_left && ref1 == center_ref){//Compares its ref2 with the ref1 of the it->seconds member.
					for(size_t k =0; k < it->second.size(); k++){//Long als
						pw_alignment al = it->second.at(k);
				//		al.print();
						size_t left1,left2,right1,right2;
						size_t ref_1 = al.getreference1();
						size_t ref_2 = al.getreference2();
						al.get_lr1(left1, right1);
						al.get_lr2(left2, right2);
						size_t pos_on_ref = data.getSequence(ref2).length();
						size_t rev_center_dir;
						stringstream rev_center;
						std::cout << center_dir << std::endl;
						if(center_dir == 0){
							rev_center_dir = 1;
						}else{
							rev_center_dir = 0;
						}
						std::cout << rev_center_dir << std::endl;
						rev_center << rev_center_dir<< ":"<< center_ref<<":"<<center_left;
						for(std::map<size_t, std::string >::iterator pos=centerOnSequence.at(ref2).begin() ; pos != centerOnSequence.at(ref2).end(); pos++){
							std::cout << pos->first << " " << pos->second << " " << center << " " << l2 << " " << rev_center.str() << std::endl;
							if((center== pos->second && pos->first == l2)||(rev_center.str()==pos->second && pos->first ==l2)){
								std::cout << "pos on ref" << pos_on_ref <<std::endl;
								pos_on_ref = pos->first;
								break;
							}
						}
						if(ref2 == ref_1 && r2 <= right1 && l2 == pos_on_ref && l2 >= left1){//add it to an intermediate map and then add it later to new_center 
							std::cout << "it is part of a long alignment-cneter on ref1 "<<std::endl;
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 != intermediate.end()){//Here we check to get sure we don't have it already in the 'intermediate'. If it is there, should be taken out of it.
								std::vector<pw_alignment> temp;
								for(size_t i =0; i < it1->second.size();i++){
									pw_alignment p1 = it1->second.at(i);
									size_t p1_left1,p1_left2,p1_right1,p1_right2;
									size_t p1_ref1 = p1.getreference1();
									size_t p1_ref2 = p1.getreference2();
									p1.get_lr1(p1_left1, p1_right1);
									p1.get_lr2(p1_left2, p1_right2);
									if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
											std::cout <<" p = p1 " <<std::endl;
									}else{
										temp.push_back(p1);
									}
								}
								it1->second = temp;
							}
							break;
						}else{
							std::cout << "it is not on any long al yet"<<std::endl;
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 == intermediate.end()){
								intermediate.insert(make_pair(center, std::vector<pw_alignment>()));
								it1 = intermediate.find(center);
							}
							bool AlreadyThere = false;
							std::cout << "it1->second size "<<it1->second.size()<<std::endl;
							for(size_t k =0; k < it1->second.size();k++){
								pw_alignment p1 = it1->second.at(k);
								size_t p1_left1,p1_left2,p1_right1,p1_right2;
								size_t p1_ref1 = p1.getreference1();
								size_t p1_ref2 = p1.getreference2();
								p1.get_lr1(p1_left1, p1_right1);
								p1.get_lr2(p1_left2, p1_right2);
								if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
									std::cout <<" p = p1 " <<std::endl;
									AlreadyThere = true;
									break;
								}
							}
							if(AlreadyThere == false){
								it1->second.push_back(p);
							}
						}
					}
				}else{//Compare its ref1 with ref of the long al
					for(size_t k =0; k < it->second.size(); k++){
						pw_alignment al = it->second.at(k);
						size_t left1,left2,right1,right2;
						size_t ref_1 = al.getreference1();
						size_t ref_2 = al.getreference2();
						al.get_lr1(left1, right1);
						al.get_lr2(left2, right2);//Again add the position on seq
						size_t pos_on_ref = data.getSequence(ref1).length();
						size_t rev_center_dir;
						stringstream rev_center;
						if(center_dir == 0){
							rev_center_dir = 1;
						}else{
							rev_center_dir = 0;
						}
						rev_center << rev_center_dir<< ":"<< center_ref<<":"<<center_left;
						for(std::map<size_t, std::string >::iterator pos=centerOnSequence.at(ref1).begin() ; pos != centerOnSequence.at(ref1).end(); pos++){
								if((center== pos->second && pos->first == l1)||(rev_center.str()== pos->second && pos->first == l1)){
								pos_on_ref = pos->first;
								break;
							}
						}
						if(ref1 == ref_1&& l1 == pos_on_ref && l1 >= left1 && r1 <= right1){//add it to an intermediate map and then add it later to new_center							
							std::cout << "it is part of a long alignment "<<std::endl; //Find it in intermediate and delete it
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 != intermediate.end()){
								std::vector<pw_alignment> temp;
								for(size_t i =0; i < it1->second.size();i++){
									pw_alignment p1 = it1->second.at(i);
									size_t p1_left1,p1_left2,p1_right1,p1_right2;
									size_t p1_ref1 = p1.getreference1();
									size_t p1_ref2 = p1.getreference2();
									p1.get_lr1(p1_left1, p1_right1);
									p1.get_lr2(p1_left2, p1_right2);
									if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
											std::cout <<" p = p1 " <<std::endl;
									}else{
										temp.push_back(p1);
									}
								}
								it1->second = temp;
							}
							break;
						}else{
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 == intermediate.end()){
								intermediate.insert(make_pair(center, std::vector<pw_alignment>()));
								it1 = intermediate.find(center);
							}
							bool AlreadyThere = false;
							std::cout << "it1->second size "<<it1->second.size()<<std::endl;
							for(size_t k =0; k < it1->second.size();k++){
								pw_alignment p1 = it1->second.at(k);
						//		p1.print();
								size_t p1_left1,p1_left2,p1_right1,p1_right2;
								size_t p1_ref1 = p1.getreference1();
								size_t p1_ref2 = p1.getreference2();
								p1.get_lr1(p1_left1, p1_right1);
								p1.get_lr2(p1_left2, p1_right2);
								if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
									std::cout <<" p = p1 " <<std::endl;
									AlreadyThere = true;
									break;
								}
							}
							if(AlreadyThere == false){
								it1->second.push_back(p);
							}
						}
					}
				}
			}
		}
	}//Add rest of the centers to the new_center
	std::set<std::string> current_centers_in_long_center;
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin() ; it != new_centers.end();it++){
		for(size_t i =0; i < it->first.size();i++){
			current_centers_in_long_center.insert(it->first.at(i));
		}
	}
	for(std::map<std::string , std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it != alignments_in_a_cluster.end(); it++){
		std::string center = it->first;
		std::set<std::string>::iterator it1 = current_centers_in_long_center.find(center);
		if(it1 == current_centers_in_long_center.end()){
		std::vector<std::string> cparts;
		strsep(center, ":", cparts);
		unsigned int center_dir = atoi(cparts.at(0).c_str());
		unsigned int center_ref = atoi(cparts.at(1).c_str());
		unsigned int center_left = atoi(cparts.at(2).c_str());
		stringstream rev_center;
		if(center_dir == 0){
			center_dir = 1;
		}else{
			center_dir = 0;
		}
		rev_center << center_dir<< ":"<< center_ref<<":"<<center_left;
		std::map<std::string, std::vector<pw_alignment> >::iterator cent = intermediate.find(center);
		std::map<std::string, std::vector<pw_alignment> >::iterator cent1 = intermediate.find(rev_center.str());
		if(cent == intermediate.end()&& cent1 == intermediate.end()){
			std::vector<std::string> new_c;
			new_c.push_back(center);
			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it1 = new_centers.find(new_c);
			if(it1 == new_centers.end()){
				new_centers.insert(make_pair(new_c,std::vector<pw_alignment>()));
				it1 = new_centers.find(new_c);
			}
			for(size_t i = 0; i < it->second.size();i++){
				it1->second.push_back(it->second.at(i));
			}			
		}
}
	}
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it=intermediate.begin(); it!= intermediate.end();it++){
		if(it->second.size() != 0){//Remove centers with 0 alignments
			std::vector<std::string> new_c;
			new_c.push_back(it->first);
			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it1 = new_centers.find(new_c);
			if(it1 == new_centers.end()){
				new_centers.insert(make_pair(new_c,std::vector<pw_alignment>()));
				it1 = new_centers.find(new_c);
			}
			for(size_t i = 0; i < it->second.size();i++){
				it1->second.push_back(it->second.at(i));
			}
		}
	}
	//At the end all the alignemnts are swapped in the way that center be on the second ref. It is used for the graph maf & arith encoding
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it!=new_centers.end();it++){
		if(it->first.size()==1){
			std::string center = it->first.at(0);
			std::vector<std::string> center_parts;
			strsep(center, ":" , center_parts);
			unsigned int center_dir = atoi(center_parts.at(0).c_str());
			unsigned int center_ref = atoi(center_parts.at(1).c_str());
			unsigned int center_left = atoi(center_parts.at(2).c_str());
			for(size_t i =0; i < it->second.size();i++){
				size_t ref1 = it->second.at(i).getreference1();
				size_t ref2 = it->second.at(i).getreference2();
				size_t left1, right1, left2, right2;
				it->second.at(i).get_lr1(left1, right1);
				it->second.at(i).get_lr2(left2, right2);
				unsigned int dir;
				if(it->second.at(i).getbegin1() < it->second.at(i).getend1()){
					dir = 0;
				}else{ 
					dir = 1;
				}
				// swapped alignment is written
				std::stringstream centerref_al;
				std::stringstream otherref_al;
				if(ref1 == center_ref && left1 == center_left && center_dir == dir){	
					for(size_t j=0; j<it->second.at(i).alignment_length(); ++j) {
						char c1, c2;
						it->second.at(i).alignment_col(j, c1, c2);
						centerref_al << c1;
						otherref_al << c2;
					}
					pw_alignment newal(otherref_al.str(),centerref_al.str(),it->second.at(i).getbegin2(),it->second.at(i).getbegin1(),it->second.at(i).getend2(),it->second.at(i).getend1(),
					it->second.at(i).getreference2(),it->second.at(i).getreference1());
					it->second.at(i) = newal;

				} else {
					assert(ref2 == center_ref && left2 == center_left);
				}
			}				
		}
	}
	std::cout<< alignments_in_a_cluster.size() << " " << new_centers.size() << std::endl;
	std::cout<< "global long centers result "<<std::endl;
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it != new_centers.end(); it++){
		std::cout<< "center size: "<< it->first.size() << std::endl;
		for(size_t j =0; j < it->first.size();j++){
			std::cout<<it->first.at(j)<<std::endl;
		}
	//	std::cout << "als size: "<< it->second.size() <<std::endl;
	//	for(size_t i =0;i< it->second.size();i++){
	//		it->second.at(i).print();
	//		std::cout<< " " <<std::endl;		
	//	}
	}

//	for(size_t i =0; i < data.numSequences();i++){
		std::cout<< "centers on sequence "<< 0 << std::endl;
		for(std::map<size_t, std::string>::iterator centPos = centerOnSequence.at(0).begin(); centPos != centerOnSequence.at(0).end();centPos++){
			std::cout << centPos->first << " " << centPos->second <<std::endl;
		}
//	}
	clock_t graph_ma_time = clock();
	// TODO better separation of the different applications of our program: create/read model, compress/decompress sequences, create graph
	// write graph result in maf format
	write_graph_maf(graphout, alignments_in_a_cluster, data);//includes clusters with no associated member.
//	write_graph_maf_with_long_centers(graphout,new_centers,data);//TODO Future work: It can be changed in a way that includes long centers themselves.
	graph_ma_time = clock() - graph_ma_time;
	clock_t arithmetic_encoding_time = clock();
//Remove centers with no other member rather than themselves:
	//(On original centers:)
	std::set<std::string> intermediate_center;
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
		if(it->second.size()==0){
			std::cout << "with 0 al "<< it->first <<std::endl;
			intermediate_center.insert(it->first);
		}
	}
	for(std::set<std::string>::iterator it = intermediate_center.begin();it !=intermediate_center.end();it++){
		std::string cent = *it;
		std::map<std::string, std::vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(cent);
		alignments_in_a_cluster.erase(it1);				
	}

//Defining weights of global clustering results! 
	//On original centers:
	std::map<std::vector<std::string>, size_t> numberOfACenter;
	std::map<std::string, unsigned int>weight;
	size_t max_bit = 8;	
	size_t max_members = 0;//Returns the largest cluster
	std::cout << "size of the original centers map "<< global_results.size() <<std::endl;
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){	
	       	if(it->second.size()> max_members){		
			max_members = it->second.size();
		}
		
	}
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){
		std::string cluster = it ->first;
//		std::cout<< "global result size: "<< it->second.size() << std::endl;
		std::map<std::string, std::vector<pw_alignment> >::iterator cent = alignments_in_a_cluster.find(cluster);
//		if(it->second.size() !=1){
		if(cent != alignments_in_a_cluster.end()){//Removing the centers with no associated member
			std::vector<std::string> temp ;
			temp.push_back(it->first);
			numberOfACenter.insert(make_pair(temp,it->second.size()));
			std::map<std::string, unsigned int >::iterator it1 = weight.find(cluster);
			if(it1 == weight.end()){
				weight.insert(make_pair(cluster,0));
				it1 = weight.find(cluster);
			}
			it1->second = (unsigned int)(it->second.size()*((m.get_powerOfTwo().at(max_bit)-1)/(double)max_members));
			if(it1 ->second == 0){
				it1->second = 1;
			}
		}else continue;
	}
	//On concatenated ones:
	std::map<std::vector<std::string> , unsigned int> long_center_weight;
	size_t maximum_mem = 0;//At the end should be equal to the number of center occurs more than all the others
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it != new_centers.end();it++){
		if(it->first.size() != 1){
			std::vector<std::string> longCenter = it->first;
			size_t number =0;
			number = it->second.size();
			std::cout << "n "<< number <<std::endl;
			numberOfACenter.insert(make_pair(longCenter,number));
		}
	}
	for(std::map<std::vector<std::string>, size_t >::iterator it = numberOfACenter.begin(); it!=numberOfACenter.end(); it++){
		if(it->second > maximum_mem){
			maximum_mem = it->second;
		}
	}
	std::cout<< "maximum_mem "<< maximum_mem << std::endl;	
	for(std::map<std::vector<std::string>, size_t >::iterator it = numberOfACenter.begin(); it!=numberOfACenter.end(); it++){
		std::map<std::vector<std::string>, unsigned int >::iterator it1 = long_center_weight.find(it->first);
		if(it1 == long_center_weight.end()){
			long_center_weight.insert(make_pair(it->first,0));

			it1 = long_center_weight.find(it->first);
		}
		it1->second = (unsigned int)(it->second*((m.get_powerOfTwo().at(max_bit)-1)/(double)maximum_mem));
		if(it1 ->second == 0){
			it1->second = 1;
		}
	}
	std::cout << "size of weight "<<long_center_weight.size() <<std::endl;
/*	for(std::map<std::vector<std::string>, unsigned int >::iterator it = long_center_weight.begin() ; it != long_center_weight.end() ; it++){
		std::cout<< it->second << " ";
	}*/
	size_t no_al = 0;
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin();it !=global_results.end();it++){
			no_al += it->second.size();
	}
//	Filling in the maps includes all memebers of each original cluster center
	std::map<std::string, std::vector<string> >membersOfCluster;//first string represents center and vector of strings are associated members
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
	//	std::cout << "number of als is "<< it->second.size() << std::endl;
		std::map<std::string, std::vector<string> >::iterator it1 = membersOfCluster.find(it->first);
		for(size_t j =0; j < it->second.size();j++){
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			pw_alignment p = it->second.at(j);
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			if(it1 == membersOfCluster.end()){
				membersOfCluster.insert(make_pair(it->first,std::vector<std::string>()));
				it1 = membersOfCluster.find(it->first);
			}
			unsigned int al_dir1;
			unsigned int al_dir2;
			if(p.getbegin1()< p.getbegin2()){
				al_dir1 = 0;
			}else	al_dir1 = 1;
			if(p.getbegin2()<p.getend2()){
				al_dir2 = 0;
			}else	al_dir2 = 1;
			std::stringstream sample1;
			std::stringstream sample2;
			sample1<<al_dir1 << ":" << ref1 << ":" << left_1;
			sample2<<al_dir2 << ":" << ref2 << ":" << left_2;
			it1->second.push_back(sample1.str());
			it1->second.push_back(sample2.str());
		}
	}
	size_t al_size = 0;
	std::map<std::string, std::string>member_of_cluster;// first string is a associated one and the second one is its center/ removed the direction from assocciated memebrs before arith. encoding. directions are not used in arith.encoding
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
		al_size += it->second.size() ;
		for(size_t j =0; j < it->second.size();j++){
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			pw_alignment p = it->second.at(j);
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			unsigned int al_dir1;
			unsigned int al_dir2;
			if(p.getbegin1()< p.getend1()){
				al_dir1 = 0;
			}else	al_dir1 = 1;
			if(p.getbegin2()<p.getend2()){
				al_dir2 = 0;
			}else	al_dir2 = 1;
			std::stringstream sample1;
			std::stringstream sample2;
			sample1 << al_dir1 << ":" <<ref1 << ":" << left_1;
			sample2 << al_dir2 << ":" <<ref2 << ":" << left_2;
			std::stringstream mem1;
			std::stringstream mem2;
			mem1 <<ref1 << ":" << left_1;
			mem2 <<ref2 << ":" << left_2;
			if(sample1.str() != it->first){
				member_of_cluster.insert(make_pair(mem1.str(),it->first));//smaple one after removing dir
			}
			if(sample2.str() != it->first){
				member_of_cluster.insert(make_pair(mem2.str(), it->first));//sample2 after removing dir
			}
			std::cout << "s1 "<< sample1.str() << " s2 "<<sample2.str() << " center " << it->first << std::endl;
			assert(sample1.str()==it->first || sample2.str() == it->first);
		}
		std::string center = it->first;
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		std::stringstream cent;
		cent<< center_ref << ":" << center_left;
		member_of_cluster.insert(make_pair(cent.str(), it->first));//center after removing dir		
	}
//	for(std::map<std::string, std::string>::iterator it = member_of_cluster.begin();it != member_of_cluster.end();it++){
//		std::cout << it->second << std::endl;
//	}
	std::cout<< "size of memeber_of_cluster is : "<< member_of_cluster.size()<<std::endl;
	std::cout << "al size "<< al_size << std::endl;


//	c.calculate_similarity();
//	c.update_values();
//	c.update_clusters();

//	initial_alignment_std::set<model> ias(data, m);
//	ias.compute(o);
//	std::cout << "There are " << o.size() << " alignment parts without partial overlap" << std::endl;
//	std::ifstream in("encode",std::ifstream::binary);
//	m.std::set_patterns(in);
//	for(std::map<std::string, std::vector<unsigned int> >::const_iterator it = m.get_high(0).begin(); it != m.get_high(0).end(); it++){
//		for(size_t k =0; k <5; k++){
//							std::cout<< it->second.at(k)<< " ; ";
//						}
//						std::cout<< " "<< std::endl;

//	}
/*	for(size_t i = 0; i < data.numSequences(); i++){
		std::cout<< "long centers on sequence " << i << " are "<<std::endl;
		for(std::map<size_t, std::vector<std::string> >::iterator it = centersPositionOnASeq.at(i).begin(); it !=  centersPositionOnASeq.at(i).end();it++){
			for(size_t j =0; j < it->second.size();j++){
				std::cout << it->second.at(j)<< " ";
			}
			std::cout<< " "<< std::endl;
		}
	}*/
/*	const dnastring & sequence = data.getSequence(0);
	for(size_t i =336681; i < 336903;i++){
		std::cout << sequence.at(i)<< " ";
	}
	std::cout << " " << std::endl;*/
//TODO set the accession of the new alignments to the one occurs the most (for the begining use a fixed one eg. 0)
//Data compression:
	use_encode en(data,m,wrap);
	if(0!=encoding_out.compare("noencode")) {
	std::cout<< "weight size: "<< weight.size()<<std::endl;
//	en.arithmetic_encoding_alignment(weight,member_of_cluster,alignments_in_a_cluster,outs);
//	en.write_to_stream(alignments_in_a_cluster,outs);
//	std::ofstream al_encode("align_encode",std::ofstream::binary);
	dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
//	en.add_flag_to_stream(outs);
//	en.arithmetic_encoding_seq(outs);
//	en.calculate_high_in_partition(weight,alignments_in_a_cluster);
//	en.add_center_to_stream(outs);
//	en.add_partition_high_to_stream(outs);
//	en.arithmetic_encoding_centers(alignments_in_a_cluster,outs,*enc);
//	en.arithmetic_encoding_alignment(weight,member_of_cluster,alignments_in_a_cluster,outs,*enc);
//	en.write_to_stream(alignments_in_a_cluster,outs);
//	en.al_encode_with_long_center(centerOnSequence,long_center_weight,alignments_in_a_cluster, centersPositionOnASeq, member_of_cluster ,outs, *enc,new_centers);
//	en.weight_flags(new_centers,alignments_in_a_cluster,outs);
	en.al_encoding(weight,member_of_cluster,alignments_in_a_cluster,outs,*enc);
//	en.al_encode_long_center_optimized_flags(centerOnSequence,long_center_weight,alignments_in_a_cluster, centersPositionOnASeq, member_of_cluster,long_centers,new_centers ,outs, *enc);
	delete enc;
	outs.close();
	}

	arithmetic_encoding_time = clock() - arithmetic_encoding_time;
	
	std::cout << "Initial alignments sets summary:" << std::endl;
	std::cout << "Of " << data.numAlignments() << " pairwise alignments with total gain of " << all_gain << " we could use " << all_used_input_alignments << std::endl;
	std::cout << "to create " << all_npo_alignments << " pieces without pairwise overlap, with a total gain of " << all_gain_in_ias << std::endl; 
	std::cout << "inital alignment sets efficiency is " << all_gain_in_ias/all_gain<<std::endl;
	std::cout << "Clustering summary: " << std::endl;
	std::cout << "Input: " << num_cluster_inputs_al << " pw alignments on " << num_cluster_seq <<" sequence pieces " <<std::endl;
	std::cout << "Output: " << num_clusters << " clusters containing " << num_cluster_members_al << " sequence pieces " << std::endl;
#if TIMING
	std::cout << "Time overview: " << std::endl;
	std::cout << "Read data " << (double)read_data_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Train models " << (double)train_models_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute initial connected components " << (double)initial_cc_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute initial alignments set + remove partial overlap " << (double)ias_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Test functions " << (double)test_function_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute second connected components " << (double)second_cc_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Affinity propagation clustering " << (double)ap_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write MSA-maf graph " << (double)graph_ma_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write compressed file " << (double)arithmetic_encoding_time/CLOCKS_PER_SEC << std::endl;
#endif



	return 0;
}
int do_dynamic_mc_model(int argc, char * argv[]) {
//	typedef overlap overlap_type;
	typedef overlap_interval_tree overlap_type;

	typedef compute_cc_avl<overlap_type> compute_cc_type;
//	typedef compute_cc compute_cc_type;
//	typedef compute_cc_with_interval_tree compute_cc_type;

	typedef dynamic_mc_model use_model;
	typedef initial_alignment_set<use_model,overlap_type> use_ias;
	typedef affpro_clusters<use_model,overlap_type> use_affpro;
	typedef dynamic_encoder<use_model> use_encode;
	if(argc < 6) {
		usage();
		cerr << "Program: dy_model" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* output maf file for the graph" << std::endl;
		cerr << "* output binary compressed file (use 'noencode' to skip encoding step)" << std::endl;
		cerr << "* number of threads to use (optional, default 1)" << std::endl;
	}

	clock_t read_data_time = clock();
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string graphout(argv[4]);
	size_t num_threads = 1;
	if(argc == 7) {
		num_threads = atoi(argv[6]);
	}
	std::string encoding_out(argv[5]);
	std::ofstream outs(encoding_out.c_str(),std::ofstream::binary);

// Read all data
	all_data data;
	// TODO also allow for reading sam
	data.read_fasta_maf(fastafile, maffile);
	overlap_type ol(data);
	wrapper wrap;
	test_encoder test;
	read_data_time = clock() - read_data_time;
// Train the model on all data
	clock_t train_models_time = clock();
	use_model m(data,wrap, num_threads);
	m.train(outs);
	train_models_time = clock() - train_models_time;
	size_t total = 0;
	std::cout << "model is trained! "<<std::endl;
// base cost to use an alignment (information need for its adress)
	double cluster_base_cost = 0; // XXX base cost is now estimated in the model and added to modify/gain costs, we can still add something here to select fewer alignments

// Find connected components of alignments with some overlap // Moved it here after training to be able to remove small als first and then connect them to each other.
	std::set<const pw_alignment*, compare_pointer_pw_alignment> al_with_pos_gain; 
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for(size_t i =0; i < data.numAlignments();i++){
		const pw_alignment & al = data.getAlignment(i);
		double g1 ,g2;
		m.gain_function(al,g1,g2);
		double av_gain = (g1+g2)/2 - cluster_base_cost;
		if(av_gain > 0){
#pragma omp critical(gain)
{
			al_with_pos_gain.insert(&al);
}
		}
	}
	//Filtering the alignments(Optional) If one doesn't use filtering 'filtered_als' should be replaced by 'al_with_pos_gain'. In fact one can even mix it with previous step and do a single loop for both of them. Although I keep them seperately for now just to make easier to do one or both.


	size_t num_seq = data.numSequences();
	filter_als filter(m,al_with_pos_gain, num_seq, num_threads);
	filter.find_overlapped_references();
	std::set<const pw_alignment*,  compare_pointer_pw_alignment> filtered_als;
	filtered_als = filter.get_filtered_als();



	//Makes components of partially overlapped alignments
	std::cout << "al with postive gains are kept "<< al_with_pos_gain.size() << " filtered als size "<< filtered_als.size() <<std::endl;
	clock_t initial_cc_time = clock();
//	compute_cc cccs(al_with_pos_gain, data.numSequences(),num_threads);
//	compute_cc_with_interval_tree cccs(al_with_pos_gain, data.numSequences());
	compute_cc_type cccs(al_with_pos_gain, data.numSequences(),num_threads);//In case of not usinf filtering replace 'filtered_als' with 'al_with_pos_gain'
	std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> > ccs; 
	cccs.compute(ccs); //fill in ccs, ordered by size(Notice that they are just connected to each other if they have overlap!)
	initial_cc_time = clock() - initial_cc_time;
/*	for(size_t i=0; i<ccs.size()-1; i++) {//A test function that chckes overlap between components. There shouldn't be any overlap between them.
		std::cout << "Connected component "<< i << " contains " << ccs.at(i).size() << " alignments" << std::endl;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.at(i).begin();it != ccs.at(i).end();it++){
			const pw_alignment * p = *it;
		//	p.print();
			for(size_t j = i+1; j < ccs.size();j++){
				std::set< const pw_alignment* , compare_pointer_pw_alignment> ccs_j = ccs.at(j);
				ol.test_no_overlap_between_ccs(*p, ccs_j);//XXX this is just a test function, comment it the first run!
				data.checkAlignmentRange(*p);//XXX this is just a test function, comment it after the first run!
			}
		}
	}*/

	std::cout<< "ccs are made!"<<std::endl;
// Select an initial alignment set for each connected component (in parallel)
	std::vector<overlap_type> cc_overlap(ccs.size(), overlap_type(data));// An overlap class object is needed for each connected component
	std::map<std::string, std::vector<std::string> > global_results;
	std::map<std::string, std::vector<pw_alignment> > alignments_in_a_cluster;//string ---> center of a cluster, vector ---> alignments with that center
	std::map<vector<std::string>, std::vector<pw_alignment> > new_centers;//Equivalent to alignments_in_a_cluster for long centers.
	std::vector<std::map<size_t, std::vector<std::string> > >centersPositionOnASeq(data.numSequences());//all the long centers of each sequence and their position
	std::vector<std::multimap<size_t, std::string> > centerOnSequence(data.numSequences());//all the centers on a seq and their position(It is used when we create long centers)
	vector<vector<std::string> > long_centers;

	size_t num_clusters = 0;
	size_t num_cluster_members_al = 0;
	size_t num_cluster_seq = 0;
	size_t num_cluster_inputs_al = 0;
#if TIMING
	clock_t ias_time = 0;
	clock_t test_function_time = 0;
	clock_t second_cc_time = 0;
	clock_t ap_time = 0;
#endif
	double all_gain = 0;
	double all_gain_in_ias = 0;
	size_t all_used_input_alignments = 0;
	size_t all_npo_alignments = 0;
	std::cout <<"size of ccs " << ccs.size() <<std::endl;//each member of it has alignments that are partially overlapped.
	for(size_t i=0; i<ccs.size(); ++i) {
		std::cout << " on initial CC " << i << " size " << ccs.at(i).size() << std::endl;
	}	
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for(size_t i=0; i<ccs.size(); ++i) {
#pragma omp critical(print)
{
		std::cout << " on initial CC " << i << " size " << ccs.at(i).size() << std::endl;
}
		clock_t ias_time_local = clock();
		std::set< const pw_alignment*, compare_pointer_pw_alignment> & cc = ccs.at(i);
	//	for(std::set< const pw_alignment*, compare_pointer_pw_alignment>::iterator setit= cc.begin();setit !=cc.end();setit++){
	//		const pw_alignment * al = *setit;
	//		double g1 ,g2;
	//		m.gain_function(*al,g1,g2);
	//		double av_gain = (g1+g2)/2 - cluster_base_cost;
	//		std::cout<< av_gain << std::endl;
	//		assert(av_gain >=0);
	//	}
	//	std::cout<<"all of them had positive gain!"<<std::endl;
		use_ias ias(data,cc, m, cluster_base_cost,num_threads);
		ias.compute(cc_overlap.at(i));//Cutting partial overlaps! XXX Slow!
		ias_time_local = clock() - ias_time_local;
		std::cout<< "DONE"<<std::endl;
		clock_t test_time_local = clock();
#ifndef NDEBUG
		cc_overlap.at(i).test_partial_overlap(); 
#endif
		test_time_local = clock() - test_time_local;
#pragma omp critical(gain_statistics)
{
		all_gain+=ias.get_max_gain();
		all_gain_in_ias+=ias.get_result_gain();
		all_used_input_alignments+=ias.get_used_alignments();
		all_npo_alignments += cc_overlap.at(i).size();
}


		clock_t second_cc_time_local = clock();
		std::vector< std::set<const pw_alignment* , compare_pointer_pw_alignment> > cc_cluster_in; //have shorter length than original alignemnts
		std::cout<< "using the cc again!"<<std::endl;
		size_t num_seq = data.numSequences();
		compute_cc_type occ(cc_overlap.at(i), num_seq ,1);
		occ.compute(cc_cluster_in);//Here we find fully overlapped pieces and save them in a single component, They are oredered by their gain value //TODO!!!!  returns several times the same al
		second_cc_time_local = clock() - second_cc_time_local;
	/*	for(size_t i=0; i<cc_cluster_in.size()-1; i++) {//A test function that chckes overlap between components. There shouldn't be any overlap between them.
			std::cout << "cc cluster in  "<< i << " contains " << cc_cluster_in.at(i).size() << " alignments" << std::endl;
			for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc_cluster_in.at(i).begin();it != cc_cluster_in.at(i).end();it++){
				const pw_alignment * p = *it;
				p->print();
				for(size_t j = i+1; j < cc_cluster_in.size();j++){
					std::set< const pw_alignment* , compare_pointer_pw_alignment> ccs_j = cc_cluster_in.at(j);
					ol.test_no_overlap_between_ccs(*p, ccs_j);//XXX this is just a test function, comment it the first run!
					data.checkAlignmentRange(*p);//XXX this is just a test function, comment it after the first run!
				}
			}
		}
		std::cout << "cc cluster in  "<< cc_cluster_in.size()-1 << " contains " << cc_cluster_in.at(cc_cluster_in.size()-1).size() << " alignments" << std::endl;*/
		for(size_t i =0; i < cc_cluster_in.size();i++){
			std::cout << "cc cluster in  "<< i << " contains " << cc_cluster_in.at(i).size() << " alignments" << std::endl;
			for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = cc_cluster_in.at(i).begin(); it != cc_cluster_in.at(i).end();it++){
				const pw_alignment * p = *it;
			//	p->print();
//				double g1,g2;
//				m.gain_function(*p,g1,g2);
//				std::cout << "g1 " << g1 << " g2 " << g2 << std::endl;
			}
		}
#if TIMING
#pragma omp critical(time)
{
		ias_time += ias_time_local;
		test_function_time += test_time_local;
		second_cc_time += second_cc_time_local;
}
#endif
		map<string, vector<pw_alignment> >al_of_a_ccs;//It is filled in with the original centers and is used for creating long centers & their als.(string ---> center, vector<pw_al> ---> alignments of a cluster)
		vector<vector<std::string> > local_long_centers;
		for(size_t j=0; j<cc_cluster_in.size(); ++j) {
			clock_t ap_time_local = clock();
			std::cout << "cc cluster in at " << j << std::endl;
			use_affpro uaf(cc_cluster_in.at(j), m, cluster_base_cost);
			std::map<std::string, std::vector<string> >cluster_result;
			//Fill in local al in a cluster
			std::map<std::string, std::vector<pw_alignment> > local_al_in_a_cluster;

			std::cout << "afp clustering "<<std::endl;
			uaf.run(cluster_result,local_al_in_a_cluster);
			ap_time_local = clock() - ap_time_local;
#if TIMING
#pragma omp critical(aptime) 
{
			ap_time += ap_time_local;
}

#endif
			size_t counter =0;
			for(std::map<std::string,std::vector<string> >::iterator it=cluster_result.begin();it !=cluster_result.end();it++){
				counter++;
				std::string center = it->first;
				std::cout << "center is "<<center << std::endl;
				std::map<std::string, std::vector<pw_alignment> >::iterator it1 = local_al_in_a_cluster.find(center);
				if(it1 == local_al_in_a_cluster.end()){
					local_al_in_a_cluster.insert(make_pair(center, std::vector<pw_alignment>()));
					it1 = local_al_in_a_cluster.find(center);
				}
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int dir = atoi(center_parts.at(0).c_str());
				unsigned int ref = atoi(center_parts.at(1).c_str());
				unsigned int left = atoi(center_parts.at(2).c_str());
				for(size_t i =0; i < it->second.size();i++){
				//	std::cout<< "member is "<<it->second.at(i)<<std::endl;
					std::vector<std::string> member_parts;
					strsep(it->second.at(i), ":" , member_parts);
					unsigned int mem_dir = atoi(member_parts.at(0).c_str());
					unsigned int mem_ref = atoi(member_parts.at(1).c_str());
					unsigned int mem_left = atoi(member_parts.at(2).c_str());
					// look at all input alignments and determine which cluster it belongs to
					for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = cc_cluster_in.at(j).begin(); it2!=cc_cluster_in.at(j).end(); ++it2){
						const pw_alignment * al = *it2;
						unsigned int al_dir1;
						if(al->getbegin1() < al->getend1()){
							al_dir1 =0;
						}else{
							al_dir1 = 1;
						}
						unsigned int al_dir2;
						if(al->getbegin2() < al->getend2()){
							al_dir2 =0;
						}else{
							al_dir2 = 1;
						}
						size_t ref1 = al->getreference1();
						size_t ref2 = al->getreference2();
						size_t left1;
						size_t left2;
						size_t right1;
						size_t right2;
						al->get_lr1(left1,right1);
						al->get_lr2(left2,right2);
						if(al_dir1 == dir && ref1 == ref && left1 == left && al_dir2 == mem_dir && ref2 == mem_ref && left2 == mem_left){
							it1->second.push_back(*al);
							double gain = m.get_the_gain(*al, center);
							assert(gain > 0);
						//	std::cout << "gain "<<gain << std::endl;
						}
						if(al_dir2 == dir && ref2 == ref && left2 == left && al_dir1==mem_dir&& ref1 == mem_ref && left1 == mem_left ){
							it1->second.push_back(*al);
							double gain = m.get_the_gain(*al, center);
							assert(gain > 0);
						//	std::cout << "gain "<<gain << std::endl;
						}
					} // for cluster_in set
				}
			}
#pragma omp critical 
{
			num_clusters += cluster_result.size();
			num_cluster_inputs_al += cc_cluster_in.at(j).size();
			//Check for those that exist in both directions and keep only one direction
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it = local_al_in_a_cluster.begin(); it!=local_al_in_a_cluster.end(); ++it) {
				if(it->second.size() != 0){
					al_of_a_ccs.insert(*it);//Here we keep two directions separately because we need them later in creating long alignments
				}
				std::string reverse_center = get_reverse(it->first);
				std::cout << "center " << it->first <<" and its reverse "<< reverse_center <<std::endl;
				std::map<std::string,std::vector<string> >::iterator result_it = cluster_result.find(it->first);//local map
				assert(result_it != cluster_result.end());
				std::map<std::string,std::vector<string> >::iterator rev_result_it = global_results.find(reverse_center);
				std::map<std::string, std::vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(reverse_center);//The gain value of its member is checked and the one with higher gain is kept.
				if(it1 == alignments_in_a_cluster.end()){
					alignments_in_a_cluster.insert(*it);// Globally saves all the als.	
					num_cluster_members_al += 1 + it->second.size();
					num_cluster_seq+= 1 + result_it->second.size();
					global_results.insert(*result_it);
				}
				else if(it1->second.size() == 0){
					alignments_in_a_cluster.insert(*it);
					alignments_in_a_cluster.erase(it1);
					num_cluster_members_al += it->second.size();
					num_cluster_seq+= result_it->second.size();
					global_results.insert(*result_it);
					assert(rev_result_it != global_results.end());
					global_results.erase(rev_result_it);
				} else{ //Reverse of the center exists and has couple of members
					double sum = 0;
					double reverse_sum = 0;
					for(size_t i =0; i < it->second.size();i++){//Memeber of the center itself
						pw_alignment p = it->second.at(i);
						std::string center = it->first;
						sum += m.get_the_gain(p,center);
					}
					for(size_t i =0; i < it1->second.size();i++){//Members of its reverse
						pw_alignment p = it1->second.at(i);
						std::string center = it1->first;
						reverse_sum +=m.get_the_gain(p,center);
					}
					if(reverse_sum > sum){
						std::cout << "rev sum > sum "<<std::endl;
						for(size_t i =0; i < it->second.size();i++){
							pw_alignment & al = it->second.at(i);
							pw_alignment p;
							al.get_reverse(p);
							it1->second.push_back(p);
						}
						num_cluster_members_al += it->second.size();
						for(size_t i =0; i < result_it->second.size();i++){
							std::string mem = result_it->second.at(i);
							std::string rev_mem = get_reverse(mem);
							rev_result_it->second.push_back(rev_mem);
						}
					}else if(reverse_sum< sum){
						std::cout << "rev sum < sum " <<std::endl;
						alignments_in_a_cluster.insert(*it);
						global_results.insert(*result_it);
						std::map<std::string , std::vector<std::string> >::iterator gl_result = global_results.find(it->first);
						assert(gl_result != global_results.end());
						std::map<std::string, std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(it->first);
						assert(it3 != alignments_in_a_cluster.end());
						for(size_t i =0; i < it1->second.size();i++){
							pw_alignment & al = it1->second.at(i);
							pw_alignment p;
							al.get_reverse(p);
							it3->second.push_back(p);
						}
						for(size_t i =0 ; i < rev_result_it->second.size();i++){
							std::string mem = rev_result_it->second.at(i);
							std::string rev_mem = get_reverse(mem);
							gl_result->second.push_back(rev_mem);
						}
						global_results.erase(rev_result_it);
						alignments_in_a_cluster.erase(it1);
						num_cluster_members_al += it->second.size();	
						num_cluster_seq+= result_it->second.size();			
					}else{//if their gain valuea are the same we just simply count the number of als and keep the one with the higher number of als
						std::cout<<"equal gain "<<std::endl;
						assert(reverse_sum == sum);
						size_t al_number = it->second.size();
						size_t reverse_al_number = it1->second.size();
						if(al_number >= reverse_al_number){
							alignments_in_a_cluster.insert(*it);
							global_results.insert(*result_it);
							std::map<std::string, std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(it->first);
							assert(it3 != alignments_in_a_cluster.end());
							for(size_t i =0; i < it1->second.size();i++){
								pw_alignment & al = it1->second.at(i);
								pw_alignment p;
								al.get_reverse(p);
								it3->second.push_back(p);
							}
							std::map<std::string ,  std::vector<std::string> >::iterator gl_result = global_results.find(it->first);
							for(size_t i =0; i < rev_result_it->second.size();i++){
								std::string mem = rev_result_it->second.at(i);
								std::string rev_mem = get_reverse(mem);
								gl_result->second.push_back(rev_mem);
							}
							global_results.erase(rev_result_it);
							alignments_in_a_cluster.erase(it1);
							num_cluster_members_al += it->second.size();	
							num_cluster_seq+= result_it->second.size();			
						}else{
							for(size_t i =0; i < it->second.size();i++){
								pw_alignment & al = it->second.at(i);
								pw_alignment p;
								al.get_reverse(p);
								it1->second.push_back(p);
							}
							num_cluster_members_al += it->second.size();
							for(size_t i = 0; i < result_it->second.size();i++){
								std::string mem = result_it->second.at(i);
								std::string rev_mem = get_reverse(mem);
								rev_result_it->second.push_back(rev_mem);
							}
							num_cluster_seq+= result_it->second.size();
						}
					}
				}		
			}
			
			std::cout << " alignments_in_a_cluster "<<alignments_in_a_cluster.size()<<std::endl;
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin() ; it != alignments_in_a_cluster.end() ; it++){
				std::cout << it->first <<std::endl;
			}
			std::cout << "num_cluster " << num_cluster_members_al <<std::endl;
			std::cout << global_results.size()<<std::endl;
			std::cout << "num_cluster_seq "<<num_cluster_seq<<std::endl;	
}
		} // for cluster_in set  
		std::cout<< "end of cluster in! "<<std::endl;
		//XXX Long centers are created here!
		finding_centers centers(data);
		//From this point on orientation matters! It shoud be taken into account in counting the number of each center happening when we are looking for the long centers!
		centers.center_frequency(al_of_a_ccs,centerOnSequence);//All the centers of each sequence are detected and centerOnSequence is filled in. It contains all the cneters on each sequence. 
	//	std::cout<< "centers on sequence 35 "<<std::endl;
	//	for(std::multimap<size_t, std::string >::iterator test_it = centerOnSequence.at(35).begin(); test_it != centerOnSequence.at(35).end(); test_it++){
	//			std::cout<<test_it->first<< " "<< test_it->second<< " "<<std::endl;
	//	}
	/*	std::cout <<"Building suffix tree: "<<std::endl;
		std::vector<std::vector<std::vector<size_t> > > all_current_centers(num_seq);
		for(size_t seq_id = 0; seq_id < num_seq; seq_id++){
			all_current_centers.at(seq_id) = centers.get_centers(seq_id);
		}
		suffixTree tree(num_seq, centers, all_current_centers);
		std::vector<size_t> center_with_highest_gain;
		size_t highest_index = 0;
		tree.make_tree(center_with_highest_gain , highest_index);*/
		//Trees are created iteratively and centers are merged:
		merging_centers merg(data, centers);

		merg.adding_new_centers(local_long_centers,centersPositionOnASeq);//After merging a group of centers we iteratively create new suffixes and a new tree based on them, then recalculate the gains again. At the end those with gain>0 go to the long_centers and centersPositionOnASeq.
		std::cout<< "new centers are made! " << local_long_centers.size() <<std::endl;
		for(size_t j  =0; j < local_long_centers.size();j++){
			long_centers.push_back(local_long_centers.at(j));
		}
		std::cout << "long center size is "<< long_centers.size() << std::endl;
		std::cout<< "long center on 3 "<<std::endl;
		for(std::map<size_t, std::vector<std::string> >::iterator test_it = centersPositionOnASeq.at(3).begin(); test_it != centersPositionOnASeq.at(3).end(); test_it++){
			for(size_t k =0; k < test_it->second.size();k++){
				std::cout<< test_it->second.at(k)<< " ";
			}
			std::cout << " " <<std::endl;
		}
		merg.create_alignment(local_long_centers,new_centers,al_of_a_ccs,centersPositionOnASeq,centerOnSequence);//Here we fill in the 'new_centers' with longer alignments in the way that center is always the second ref

	} // for connected components
	//Short centers are added to the long centers container. //TODO Think about simplifying it! XXX Note that second both ref can be backwards
	std::cout << "number of long centers "<< long_centers.size() << " " <<new_centers.size() <<std::endl;
	std::map<std::string, std::vector<pw_alignment> > intermediate;
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it != new_centers.end(); it++){
		std::vector<std::string> connected_center = it->first;
		for(size_t i =0; i < connected_center.size();i++){
			std::string center = connected_center.at(i);
			std::cout << "current center "<< center << std::endl;
			std::vector<std::string> cparts;
			strsep(center, ":", cparts);
			unsigned int center_dir = atoi(cparts.at(0).c_str());
			unsigned int center_ref = atoi(cparts.at(1).c_str());
			unsigned int center_left = atoi(cparts.at(2).c_str());
			std::map<std::string, std::vector<pw_alignment> >::iterator cent = alignments_in_a_cluster.find(center);// XXX Attention! Here only one direction was saved!
			if(cent == alignments_in_a_cluster.end()){
				if(center_dir == 0){
					center_dir = 1;
				}else{
					center_dir = 0;
				}
				stringstream temp;
				temp << center_dir<< ":"<< center_ref<<":"<<center_left;
				center = temp.str();
				std::cout << "CENT "<<center <<std::endl;
				cent = alignments_in_a_cluster.find(center);
				assert(cent != alignments_in_a_cluster.end());
			}
			std::cout << "number of als " << cent->second.size() <<std::endl;
			for(size_t j =0;j < cent->second.size();j++){//Go through the short alignments and check to see if they are part a long one
				pw_alignment p = cent->second.at(j); //Check if it happened on a long alignment!
				std::cout << "p is "<<std::endl;
			//	p.print();
				size_t l1,l2,r1,r2;
				size_t ref1 = p.getreference1();
				size_t ref2 = p.getreference2();
				std::cout << "ref1 "<< ref1 << " ref2 "<<ref2 <<std::endl;
				p.get_lr1(l1, r1);
				p.get_lr2(l2, r2);
				if(l1 == center_left && ref1 == center_ref){//Compares its ref2 with the ref1 of the it->seconds member.
					for(size_t k =0; k < it->second.size(); k++){//Long als
						pw_alignment al = it->second.at(k);
				//		al.print();
						size_t left1,left2,right1,right2;
						size_t ref_1 = al.getreference1();
						size_t ref_2 = al.getreference2();
						al.get_lr1(left1, right1);
						al.get_lr2(left2, right2);
						size_t pos_on_ref = data.getSequence(ref2).length();
						size_t rev_center_dir;
						stringstream rev_center;
						std::cout << center_dir << std::endl;
						if(center_dir == 0){
							rev_center_dir = 1;
						}else{
							rev_center_dir = 0;
						}
						std::cout << rev_center_dir << std::endl;
						rev_center << rev_center_dir<< ":"<< center_ref<<":"<<center_left;
						for(std::multimap<size_t, std::string >::iterator pos=centerOnSequence.at(ref2).begin() ; pos != centerOnSequence.at(ref2).end(); pos++){
							std::cout << pos->first << " " << pos->second << " " << center << " " << l2 << " " << rev_center.str() << std::endl;
							if((center== pos->second && pos->first == l2)||(rev_center.str()==pos->second && pos->first ==l2)){
								std::cout << "pos on ref" << pos_on_ref <<std::endl;
								pos_on_ref = pos->first;
								break;
							}
						}
						if(ref2 == ref_1 && r2 <= right1 && l2 == pos_on_ref && l2 >= left1){//add it to an intermediate map and then add it later to new_center 
							std::cout << "it is part of a long alignment-cneter on ref1 "<<std::endl;
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 != intermediate.end()){//Here we check to get sure we don't have it already in the 'intermediate'. If it is there, should be taken out of it.
								std::vector<pw_alignment> temp;
								for(size_t i =0; i < it1->second.size();i++){
									pw_alignment p1 = it1->second.at(i);
									size_t p1_left1,p1_left2,p1_right1,p1_right2;
									size_t p1_ref1 = p1.getreference1();
									size_t p1_ref2 = p1.getreference2();
									p1.get_lr1(p1_left1, p1_right1);
									p1.get_lr2(p1_left2, p1_right2);
									if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
											std::cout <<" p = p1 " <<std::endl;
									}else{
										temp.push_back(p1);
									}
								}
								it1->second = temp;
							}
							break;
						}else{
							std::cout << "it is not on any long al yet"<<std::endl;
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 == intermediate.end()){
								intermediate.insert(make_pair(center, std::vector<pw_alignment>()));
								it1 = intermediate.find(center);
							}
							bool AlreadyThere = false;
							std::cout << "it1->second size "<<it1->second.size()<<std::endl;
							for(size_t k =0; k < it1->second.size();k++){
								pw_alignment p1 = it1->second.at(k);
								size_t p1_left1,p1_left2,p1_right1,p1_right2;
								size_t p1_ref1 = p1.getreference1();
								size_t p1_ref2 = p1.getreference2();
								p1.get_lr1(p1_left1, p1_right1);
								p1.get_lr2(p1_left2, p1_right2);
								if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
									std::cout <<" p = p1 " <<std::endl;
									AlreadyThere = true;
									break;
								}
							}
							if(AlreadyThere == false){
								it1->second.push_back(p);
							}
						}
					}
				}else{//Compare its ref1 with ref of the long al
					for(size_t k =0; k < it->second.size(); k++){
						pw_alignment al = it->second.at(k);
						size_t left1,left2,right1,right2;
						size_t ref_1 = al.getreference1();
						size_t ref_2 = al.getreference2();
						al.get_lr1(left1, right1);
						al.get_lr2(left2, right2);//Again add the position on seq
						size_t pos_on_ref = data.getSequence(ref1).length();
						size_t rev_center_dir;
						stringstream rev_center;
						if(center_dir == 0){
							rev_center_dir = 1;
						}else{
							rev_center_dir = 0;
						}
						rev_center << rev_center_dir<< ":"<< center_ref<<":"<<center_left;
						for(std::multimap<size_t, std::string >::iterator pos=centerOnSequence.at(ref1).begin() ; pos != centerOnSequence.at(ref1).end(); pos++){
								if((center== pos->second && pos->first == l1)||(rev_center.str()== pos->second && pos->first == l1)){
								pos_on_ref = pos->first;
								break;
							}
						}
						if(ref1 == ref_1&& l1 == pos_on_ref && l1 >= left1 && r1 <= right1){//add it to an intermediate map and then add it later to new_center							
							std::cout << "it is part of a long alignment "<<std::endl; //Find it in intermediate and delete it
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 != intermediate.end()){
								std::vector<pw_alignment> temp;
								for(size_t i =0; i < it1->second.size();i++){
									pw_alignment p1 = it1->second.at(i);
									size_t p1_left1,p1_left2,p1_right1,p1_right2;
									size_t p1_ref1 = p1.getreference1();
									size_t p1_ref2 = p1.getreference2();
									p1.get_lr1(p1_left1, p1_right1);
									p1.get_lr2(p1_left2, p1_right2);
									if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
											std::cout <<" p = p1 " <<std::endl;
									}else{
										temp.push_back(p1);
									}
								}
								it1->second = temp;
							}
							break;
						}else{
							std::map<std::string, std::vector<pw_alignment> >::iterator it1 = intermediate.find(center);
							if(it1 == intermediate.end()){
								intermediate.insert(make_pair(center, std::vector<pw_alignment>()));
								it1 = intermediate.find(center);
							}
							bool AlreadyThere = false;
							std::cout << "it1->second size "<<it1->second.size()<<std::endl;
							for(size_t k =0; k < it1->second.size();k++){
								pw_alignment p1 = it1->second.at(k);
						//		p1.print();
								size_t p1_left1,p1_left2,p1_right1,p1_right2;
								size_t p1_ref1 = p1.getreference1();
								size_t p1_ref2 = p1.getreference2();
								p1.get_lr1(p1_left1, p1_right1);
								p1.get_lr2(p1_left2, p1_right2);
								if(p1_ref1== ref1 && p1_ref2 == ref2 && p1_left1 == l1&& p1_right1==r1 &&p1_left2==l2 && p1_right2==r2){
									std::cout <<" p = p1 " <<std::endl;
									AlreadyThere = true;
									break;
								}
							}
							if(AlreadyThere == false){
								it1->second.push_back(p);
							}
						}
					}
				}
			}
		}
	}//Add rest of the centers to the new_center
	std::set<std::string> current_centers_in_long_center;
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin() ; it != new_centers.end();it++){
		for(size_t i =0; i < it->first.size();i++){
			current_centers_in_long_center.insert(it->first.at(i));
		}
	}
	std::cout<< "rest of centers are added : "<<std::endl;
	for(std::map<std::string , std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it != alignments_in_a_cluster.end(); it++){
		std::string center = it->first;
		std::set<std::string>::iterator it1 = current_centers_in_long_center.find(center);
		if(it1 == current_centers_in_long_center.end()){
			std::cout << center << std::endl;
			std::vector<std::string> cparts;
			strsep(center, ":", cparts);
			unsigned int center_dir = atoi(cparts.at(0).c_str());
			unsigned int center_ref = atoi(cparts.at(1).c_str());
			unsigned int center_left = atoi(cparts.at(2).c_str());
			stringstream rev_center;
			if(center_dir == 0){
				center_dir = 1;
			}else{
				center_dir = 0;
			}
			rev_center << center_dir<< ":"<< center_ref<<":"<<center_left;
			std::map<std::string, std::vector<pw_alignment> >::iterator cent = intermediate.find(center);
			std::map<std::string, std::vector<pw_alignment> >::iterator cent1 = intermediate.find(rev_center.str());
			if(cent == intermediate.end()&& cent1 == intermediate.end()){
				std::vector<std::string> new_c;
				new_c.push_back(center);
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it1 = new_centers.find(new_c);
				if(it1 == new_centers.end()){
					new_centers.insert(make_pair(new_c,std::vector<pw_alignment>()));
					it1 = new_centers.find(new_c);
				}
				for(size_t i = 0; i < it->second.size();i++){
					it1->second.push_back(it->second.at(i));
				}			
			}
		}
	}
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it=intermediate.begin(); it!= intermediate.end();it++){
		if(it->second.size() != 0){//Remove centers with 0 alignments
			std::vector<std::string> new_c;
			new_c.push_back(it->first);
			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it1 = new_centers.find(new_c);
			if(it1 == new_centers.end()){
				new_centers.insert(make_pair(new_c,std::vector<pw_alignment>()));
				it1 = new_centers.find(new_c);
			}
			for(size_t i = 0; i < it->second.size();i++){
				it1->second.push_back(it->second.at(i));
			}
		}
	}
	//At the end all the alignemnts are swapped in the way that center be on the second ref. It is used for the graph maf & arith encoding
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it!=new_centers.end();it++){
		if(it->first.size()==1){
			std::string center = it->first.at(0);
			std::vector<std::string> center_parts;
			strsep(center, ":" , center_parts);
			unsigned int center_dir = atoi(center_parts.at(0).c_str());
			unsigned int center_ref = atoi(center_parts.at(1).c_str());
			unsigned int center_left = atoi(center_parts.at(2).c_str());
			for(size_t i =0; i < it->second.size();i++){
				size_t ref1 = it->second.at(i).getreference1();
				size_t ref2 = it->second.at(i).getreference2();
				size_t left1, right1, left2, right2;
				it->second.at(i).get_lr1(left1, right1);
				it->second.at(i).get_lr2(left2, right2);
				unsigned int dir;
				if(it->second.at(i).getbegin1() < it->second.at(i).getend1()){
					dir = 0;
				}else{ 
					dir = 1;
				}
				// swapped alignment is written
				std::stringstream centerref_al;
				std::stringstream otherref_al;
				if(ref1 == center_ref && left1 == center_left && center_dir == dir){	
					for(size_t j=0; j<it->second.at(i).alignment_length(); ++j) {
						char c1, c2;
						it->second.at(i).alignment_col(j, c1, c2);
						centerref_al << c1;
						otherref_al << c2;
					}
					pw_alignment newal(otherref_al.str(),centerref_al.str(),it->second.at(i).getbegin2(),it->second.at(i).getbegin1(),it->second.at(i).getend2(),it->second.at(i).getend1(),
					it->second.at(i).getreference2(),it->second.at(i).getreference1());
					it->second.at(i) = newal;

				} else {
					assert(ref2 == center_ref && left2 == center_left);
				}
			}				
		}
	}
	clock_t graph_ma_time = clock();

	// write graph result in maf format
	write_graph_maf(graphout, alignments_in_a_cluster, data);//includes clusters with no associated member.
//	write_graph_maf_with_long_centers(graphout,new_centers,data);//TODO It can be changed in a way that includes long centers themselves as well.

	graph_ma_time = clock() - graph_ma_time;


	clock_t arithmetic_encoding_time = clock();
	//Remove centers with no associated member.
	std::set<std::string> intermediate_center;
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
		if(it->second.size()==0){
			intermediate_center.insert(it->first);
		}
	}
	for(std::set<std::string>::iterator it = intermediate_center.begin();it !=intermediate_center.end();it++){
		std::string cent = *it;
		std::map<std::string, std::vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(cent);
		alignments_in_a_cluster.erase(it1);				
	}
	//Defining weights of global clustering results! 
	//On original centers:
	std::map<std::vector<std::string>, size_t> numberOfACenter;
	std::map<std::string, unsigned int>weight;
	size_t max_bit = 8;	
	size_t max_members = 0;//Returns the largest cluster
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){	
        	if(it->second.size()> max_members){		
			max_members = it->second.size();
		}
		
	}	
	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){
		std::string cluster = it ->first;
		std::map<std::string, std::vector<pw_alignment> >::iterator cent = alignments_in_a_cluster.find(cluster);
		if(cent != alignments_in_a_cluster.end()){//Removing the centers with no associated member
			std::map<std::string, unsigned int >::iterator it1 = weight.find(cluster);
			std::vector<std::string> temp ;
			temp.push_back(it->first);
			numberOfACenter.insert(make_pair(temp,it->second.size()));

			if(it1 == weight.end()){
				weight.insert(make_pair(cluster,0));
				it1 = weight.find(cluster);
			}
			size_t power_of_two = 1<<(max_bit);
			it1->second = (unsigned int)(it->second.size()*(power_of_two-1)/(double)max_members);
			if(it1 ->second == 0){
				it1->second = 1;
			}
		}else continue;
	}
	//On concatenated ones:
	std::map<std::vector<std::string> , unsigned int> long_center_weight;
	size_t maximum_mem = 0;//At the end should be equal to the number of center occurs more than all the others
	for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it != new_centers.end();it++){
		if(it->first.size() != 1){
			std::vector<std::string> longCenter = it->first;
			size_t number =0;
			number = it->second.size();
			std::cout << "n "<< number <<std::endl;
			numberOfACenter.insert(make_pair(longCenter,number));
		}
	}
	for(std::map<std::vector<std::string>, size_t >::iterator it = numberOfACenter.begin(); it!=numberOfACenter.end(); it++){
		if(it->second > maximum_mem){
			maximum_mem = it->second;
		}
	}
	std::cout<< "maximum_mem "<< maximum_mem << std::endl;	
	for(std::map<std::vector<std::string>, size_t >::iterator it = numberOfACenter.begin(); it!=numberOfACenter.end(); it++){
		std::map<std::vector<std::string>, unsigned int >::iterator it1 = long_center_weight.find(it->first);
		if(it1 == long_center_weight.end()){
			long_center_weight.insert(make_pair(it->first,0));

			it1 = long_center_weight.find(it->first);
		}
		it1->second = (unsigned int)(it->second*((1<<(max_bit))-1)/(double)maximum_mem);
		if(it1 ->second == 0){
			it1->second = 1;
		}
	}
	std::cout << "size of weight "<<long_center_weight.size() <<std::endl;

//	size_t no_al = 0;
//	for(std::map<std::string, std::vector<string> >::iterator it=global_results.begin();it !=global_results.end();it++){
//			no_al += it->second.size();
//	}
	std::map<std::string, std::vector<string> >membersOfCluster;//first string represents center and vector of strings are associated members
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
	//	std::cout << "number of als is "<< it->second.size() << std::endl;
		std::map<std::string, std::vector<string> >::iterator it1 = membersOfCluster.find(it->first);
		for(size_t j =0; j < it->second.size();j++){
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			pw_alignment p = it->second.at(j);
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			if(it1 == membersOfCluster.end()){
				membersOfCluster.insert(make_pair(it->first,std::vector<std::string>()));
				it1 = membersOfCluster.find(it->first);
			}
			unsigned int al_dir1;
			unsigned int al_dir2;
			if(p.getbegin1()< p.getbegin2()){
				al_dir1 = 0;
			}else	al_dir1 = 1;
			if(p.getbegin2()<p.getend2()){
				al_dir2 = 0;
			}else	al_dir2 = 1;
			std::stringstream sample1;
			std::stringstream sample2;
			sample1<<al_dir1 << ":" << ref1 << ":" << left_1;
			sample2<<al_dir2 << ":" << ref2 << ":" << left_2;
			it1->second.push_back(sample1.str());
			it1->second.push_back(sample2.str());
		}
	}
	std::cout << "number of centers "<< membersOfCluster.size() << std::endl;

	size_t al_size = 0;
	std::map<std::string,std::string>member_of_cluster;//first string is a associated one and the second one is its center.The direction from assocciated memebrs are removed before arith.encoding and are not used in arith.encoding
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
		al_size += it->second.size() ;
		for(size_t j =0; j < it->second.size();j++){
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			pw_alignment p = it->second.at(j);
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			unsigned int al_dir1;
			unsigned int al_dir2;
			if(p.getbegin1()< p.getend1()){
				al_dir1 = 0;
			}else	al_dir1 = 1;
			if(p.getbegin2()<p.getend2()){
				al_dir2 = 0;
			}else	al_dir2 = 1;
			std::stringstream sample1;
			std::stringstream sample2;
			sample1 << al_dir1 << ":" <<ref1 << ":" << left_1;
			sample2 << al_dir2 << ":" <<ref2 << ":" << left_2;
			std::stringstream mem1;
			std::stringstream mem2;
			mem1 <<ref1 << ":" << left_1;
			mem2 <<ref2 << ":" << left_2;
			if(sample1.str() != it->first){
				member_of_cluster.insert(make_pair(mem1.str(),it->first));//smaple one after removing dir
			}
			if(sample2.str() != it->first){
				member_of_cluster.insert(make_pair(mem2.str(), it->first));//sample2 after removing dir
			}
			std::cout << "s1 "<< sample1.str() << " s2 "<<sample2.str() << " center " << it->first << std::endl;
			assert(sample1.str()==it->first || sample2.str() == it->first);
		}
		std::string center = it->first;
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		std::stringstream cent;
		cent<< center_ref << ":" << center_left;
		member_of_cluster.insert(make_pair(cent.str(), it->first));//center after removing dir		
	}
	for(std::map<std::string, std::string>::iterator it = member_of_cluster.begin();it != member_of_cluster.end();it++){
		std::cout << it->second << std::endl;
	}
	std::cout<< "size of memeber_of_cluster is : "<< member_of_cluster.size()<<std::endl;
	std::cout << "al size "<< al_size << std::endl;

/*	std::ofstream al_high_encode("al_high_encode.txt"); //It is just for testing the high values
	m.write_al_high_onstream(al_high_encode);
	al_high_encode.close();*/

//Data compression:
	use_encode en(data,m,wrap);
	if(0!=encoding_out.compare("noencode")){
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
	//	en.al_encode(weight,member_of_cluster,alignments_in_a_cluster,outs,*enc);
		en.al_encode_with_long_center(centerOnSequence,long_center_weight,alignments_in_a_cluster, centersPositionOnASeq, member_of_cluster ,outs, *enc,new_centers);
		delete enc;
		outs.close();
	}

	arithmetic_encoding_time = clock() - arithmetic_encoding_time;

	std::cout<< "size of memeber_of_cluster is : "<< member_of_cluster.size()<<std::endl;
	std::cout << "Initial alignments sets summary:" << std::endl;
	std::cout << "Of " << data.numAlignments() << " pairwise alignments with total gain of " << all_gain << " we could use " << all_used_input_alignments << std::endl;
	std::cout << "to create " << all_npo_alignments << " pieces without pairwise overlap, with a total gain of " << all_gain_in_ias << std::endl; 
	std::cout << "inital alignment sets efficiency is " << all_gain_in_ias/all_gain<<std::endl;
	std::cout << "Clustering summary: " << std::endl;
	std::cout << "Input: " << num_cluster_inputs_al << " pw alignments on " << num_cluster_seq <<" sequence pieces " <<std::endl;
	std::cout << "Output: " << num_clusters << " clusters containing " << num_cluster_members_al << " alignments " << std::endl;
#if TIMING
	std::cout << "Time overview: " << std::endl;
	std::cout << "Read data " << (double)read_data_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Train models " << (double)train_models_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute initial connected components " << (double)initial_cc_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute initial alignments set + remove partial overlap " << (double)ias_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Test functions " << (double)test_function_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Compute second connected components " << (double)second_cc_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Affinity propagation clustering " << (double)ap_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write MSA-maf graph " << (double)graph_ma_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write compressed file " << (double)arithmetic_encoding_time/CLOCKS_PER_SEC << std::endl;
#endif 
	return 0;

}

// current function
int do_dynamic_mc_model_with_two_edge(int argc, char * argv[]) {
	std::srand(time(0));
	typedef overlap_interval_tree overlap_type;
	typedef dynamic_mc_model use_model;
	typedef compute_cc_avl_fraction cc_type;
	typedef initial_alignment_set<use_model,overlap_type> use_ias;
	if(argc < 8) {
		usage();
		cerr << "Program: new_dy_model" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* output maf file for the graph" << std::endl;
		cerr << "* output binary compressed file (use 'noencode' to skip encoding step)" << std::endl;
		cerr << "* output fasta file for the graph" <<std::endl; //They are neede for mapping reads against the graph
		cerr << "* output dot file represents adjacencies "<<std::endl;
		cerr << "* number of threads to use (optional, default 1)" << std::endl;
		return 1;
	}
	clock_t read_data_time = clock();
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string graphout(argv[4]);
	std::string encoding_out(argv[5]);
	std::ofstream outs(encoding_out.c_str(),std::ofstream::binary);
	std::string graphfasta(argv[6]);
	std::string dotfile(argv[7]);
	std::ofstream dotout(dotfile.c_str());
	size_t num_threads = 1;
	if(argc == 9) {
		num_threads = atoi(argv[8]);
	}
// Read all data
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	overlap_type ol(data);
	wrapper wrap;
	test_encoder test;
	read_data_time = clock() - read_data_time;
// Train the model on all data
	clock_t train_models_time = clock();
	use_model m(data,wrap, num_threads);
	m.train(outs);
	train_models_time = clock() - train_models_time;
	std::cout << "model is trained! "<<std::endl;
// base cost to use an alignment (information need for its adress)
	double cluster_base_cost = 0; // XXX base cost is now estimated in the model and added to modify/gain costs, we can still add something here to select fewer alignments
// Find connected components of alignments with some overlap // Moved it here after training to be able to remove small als first and then connect them to each other.
	std::set<const pw_alignment*, compare_pointer_pw_alignment> al_with_pos_gain; 
	std::set<const located_alignment*, compare_pointer_pw_alignment> lal_with_pos_gain; 
	std::set<const pw_alignment*, compare_pointer_pw_alignment> remainders; 
	double sum_of_input_gain = 0;
	double sum_of_input_gain_loc = 0;
#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i < data.numAlignments();i++){
		const pw_alignment & al = data.getAlignment(i);
		double g1 ,g2;
		m.gain_function(al,g1,g2);
		double av_gain = (g1+g2)/2 ;


		const located_alignment & lal = data.getLAlignment(i);
		m.gain_function(lal, g1, g2);
		double lav_gain = (g1+g2) /2;

	//	std::cout << " al " << lal.alignment_length() << " old gain " << av_gain << " new gain " << lav_gain << std::endl;
	//	al.print();
		if(av_gain > 0){
#pragma omp critical(al_insert)
{
			al_with_pos_gain.insert(&al);
			sum_of_input_gain += av_gain;
}
		}
		if(lav_gain > 0) {
#pragma omp critical(lal_insert)
{
			lal_with_pos_gain.insert(&lal);
			sum_of_input_gain_loc += lav_gain;
}
		}

	}
	std::cout << "al with positive gains are kept " << al_with_pos_gain.size() << " total gain " << sum_of_input_gain << std::endl;
	std::cout << "located al with positive gains are kept " << al_with_pos_gain.size() << " total gain " << sum_of_input_gain_loc << std::endl;


//Filtering the alignments(Optional)
	size_t num_seq = data.numSequences();
//	filter_als filter(m,al_with_pos_gain, num_seq, num_threads);
//	filter.find_overlapped_references();


	size_t gap_in_long_centers = m.estimate_gap_in_long_centers(); // TODO use this to compute long centers
//	alignment_filter alfilter(m, data, al_with_pos_gain, gap_in_long_centers, num_threads);//  TODO XXX 

	clock_t new_ias_time = clock();

	std::vector<std::vector<pw_alignment> > new_overlaps;


	double gain_after_ias = 0;
/*
	divide_and_conquer_alignments dca(data, al_with_pos_gain, m, num_seq, num_threads);
	dca.run(new_overlaps);
	for(size_t i=0; i<new_overlaps.size(); ++i) {
		for(size_t j=0; j<new_overlaps.at(i).size(); ++j) {
			double g1, g2;
			m.gain_function(new_overlaps.at(i).at(j), g1, g2);
			double gain = (g1+g2)/2.0;
			gain_after_ias += gain;
		}

	}
*/
	new_ias_time = clock() - new_ias_time;

	std::cout << " Gain before IAS " << sum_of_input_gain << " after IAS " << gain_after_ias << std::endl;
	std::cout << " Times new " << (double)new_ias_time/CLOCKS_PER_SEC << std::endl;
	

	clock_t old_ias_time = clock();

//	std::set<const pw_alignment*,  compare_pointer_pw_alignment> filtered_als; //XXX OPTIONAL in case of using filtered als all al_with_pos_gain s should be replaces by filtered_als.
//	filtered_als = filter.get_filtered_als();
////Makes components of partially overlapped alignments
//	std::cout << "filtered als size "<< filtered_als.size() << std::endl;

	std::set< pw_alignment, compare_pw_alignment> mixed_als; 
	std::cout << " now we compute connected components " << std::endl;
	clock_t cc_init_time = clock();
	cc_type cccs(al_with_pos_gain, data.numSequences(), num_threads);//The CC one with the fraction
	cc_init_time = clock() - cc_init_time;
	std::cout << "CC Init time " << (double)cc_init_time/CLOCKS_PER_SEC << std::endl;
	remainders = al_with_pos_gain;
	clock_t cc_comp_time = clock();
	double fraction = 0.95;
	while(fraction > 0.60){
		std::set<const pw_alignment* , compare_pointer_pw_alignment> ccs; 
		std::cout << "fraction is "<< fraction << " run on " << remainders.size()<<std::endl;
		cccs.non_recursive_compute(ccs,remainders,fraction); //A graph of als with 95% overlap is made.
		remainders = cccs.get_remaining_als();//XXX is faster if fill it in at the same time with running "non_recursive_compute"
		cc_comp_time = clock() - cc_comp_time;

		std::cout << "CC Init time " << (double)cc_init_time/CLOCKS_PER_SEC << std::endl;
		std::cout << "CC Comp time " << (double)cc_comp_time/CLOCKS_PER_SEC << std::endl;
		std::cout<< "ccs are made!" << ccs.size() <<std::endl; //Initial ccs with 95% overlap
		std::cout << "mixed als size1 "<< mixed_als.size() << std::endl;
		std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment > > stacks;
		std::set< const pw_alignment* , compare_pointer_pw_alignment> & temp = ccs;
		std::cout << "ccs size "<< ccs.size() <<std::endl;
		two_edge_cc twoEdgeCC(cccs,ccs);
		twoEdgeCC.find_bridges(stacks,mixed_als);
		std::cout << "mixed als size2 "<< mixed_als.size() << std::endl;
		size_t Counter =0;
		std::cout << "stack size in main is " << stacks.size()<< std::endl;
		for(size_t i=0; i<stacks.size(); ++i) {
			std::set< const pw_alignment*, compare_pointer_pw_alignment> & stack = stacks.at(i);
			for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = stack.begin();it != stack.end();it++){
				Counter ++;
			}
		}
		std::cout << "total number of stacks' member is "<< Counter <<std::endl;		
//Cutting partial overlaps in 'stacks' 
		std::vector<overlap_type> cc_overlap(stacks.size(), overlap_type(data));
		size_t COUNT = 0;
		for(size_t i=0; i<stacks.size(); ++i) {//TODO Order them by size before cutting, It can happen in parallel!
			std::set< const pw_alignment*, compare_pointer_pw_alignment> & stack = stacks.at(i);
			std::cout << "als number "<< stack.size() <<std::endl;
			for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = stack.begin();it != stack.end();it++){ //TODO it is not thoughtful! Fill in the remainder while fill in the stacks!
				std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it1 = temp.find(*it);
				assert(it1 != temp.end());
				temp.erase(*it1);
			}

			use_ias ias(data,stack, m, cluster_base_cost,num_threads);
			ias.compute(cc_overlap.at(i));
			const std::set< pw_alignment, compare_pw_alignment> & als_in_overlap = cc_overlap.at(i).get_all();
			COUNT += als_in_overlap.size();
			for(std::set<pw_alignment,compare_pw_alignment>::const_iterator it = als_in_overlap.begin(); it != als_in_overlap.end();it++){
				mixed_als.insert(*it);
			}
		}
		std::cout << "total returned numbers after cutting is " << COUNT <<std::endl;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = temp.begin(); it != temp.end();it++){
			mixed_als.insert(**it);
		}
		std::cout << "mixed als size3 "<< mixed_als.size() << std::endl;
		fraction -= 0.05;
	}
	for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = remainders.begin(); it != remainders.end();it++){
		mixed_als.insert(**it);
	}
	std::cout<< "mixed als size is "<< mixed_als.size() << std::endl; //Now we need to run cc without fraction and cut them afterwards
	for(std::set< pw_alignment , compare_pw_alignment>::iterator it = mixed_als.begin(); it != mixed_als.end();it++){
		data.checkAlignmentRange(*it);
	}
//Creating connected components with 0.0 fraction rate over all of them and cut them afterwards.
	typedef compute_cc_avl<overlap_type> ccc_type;
	ccc_type no_fraction_cccs(mixed_als, data.numSequences(), num_threads);
	std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> > no_fraction_ccs;
	no_fraction_cccs.compute(no_fraction_ccs); //TODO check that there is no overlap between the components
	std::cout << "number of components is " << no_fraction_ccs.size() <<std::endl;
	size_t NUM = 0;
	for(size_t i=0; i<no_fraction_ccs.size(); ++i) {
		std::cout << " on CC " << i << " size " << no_fraction_ccs.at(i).size() << std::endl;
			NUM+=no_fraction_ccs.at(i).size();
	}
	std::cout << "NUM " << NUM <<std::endl;

//Cutting partial overlaps: (Fully)
	std::vector<overlap_type> ccc_overlap(no_fraction_ccs.size(), overlap_type(data));

	std::map<std::string, std::vector<std::string> > global_results;
	std::map<std::string, std::vector<pw_alignment> > alignments_in_a_cluster;//string ---> center of a cluster, vector ---> alignments with that center
	std::map<vector<std::string>, std::vector<pw_alignment> > new_centers;//Equivalent to alignments_in_a_cluster for long centers.
	std::vector<std::map<size_t, std::vector<std::string> > >centersPositionOnASeq(data.numSequences());//all the long centers of each sequence and their position
	std::vector<std::multimap<size_t, std::string> > centerOnSequence(data.numSequences());//all the centers on a seq and their position(It is used when we create long centers and adjacencies)
	vector<vector<std::string> > long_centers;

	double sum_of_gain_before_clust = 0;


	size_t num_clusters = 0;
	size_t num_cluster_members_al = 0;
	size_t num_cluster_seq = 0;
	size_t num_cluster_inputs_al = 0;
	clock_t ap_time = 0;
	typedef affpro_clusters<use_model,overlap_type> use_affpro;
	std::map<std::string, std::vector<pw_alignment> > merged_als;
	for(size_t i=0; i<no_fraction_ccs.size(); ++i) {
		std::set< const pw_alignment*, compare_pointer_pw_alignment> & cc = no_fraction_ccs.at(i);
		use_ias ias(data,cc, m, cluster_base_cost,num_threads);
		ias.compute(ccc_overlap.at(i));
		const std::set<pw_alignment, compare_pw_alignment> check_gain = ccc_overlap.at(i).get_all();
		for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = check_gain.begin(); it != check_gain.end(); it++){
			const pw_alignment pal = *it;
			double g1 ,g2;
			m.gain_function(pal,g1,g2);
			double av_gain = (g1+g2)/2 ;
			assert(av_gain > 0);
			if(g1 < 0){
				std::cout << "g1 is less than zero! "<<std::endl;
				pal.print();
			}
			if(g2 < 0){
				std::cout << "g2 is less than zero! "<<std::endl;
				pal.print();
			}
		}
		std::vector< std::set<const pw_alignment* , compare_pointer_pw_alignment> > cc_cluster_in; //have shorter length than original alignemnts
		std::cout<< "using the cc again!"<<std::endl;
		ccc_type occ(ccc_overlap.at(i), num_seq ,1);
		occ.compute(cc_cluster_in);
	/*	for(size_t i=0; i<cc_cluster_in.size()-1; i++) {//A test function that chckes overlap between components. There shouldn't be any overlap between them.
			std::cout << "cc cluster in  "<< i << " contains " << cc_cluster_in.at(i).size() << " alignments" << std::endl;
			for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc_cluster_in.at(i).begin();it != cc_cluster_in.at(i).end();it++){
				const pw_alignment * p = *it;
			//	p->print();
				for(size_t j = i+1; j < cc_cluster_in.size();j++){
					std::set< const pw_alignment* , compare_pointer_pw_alignment> ccs_j = cc_cluster_in.at(j);
					ol.test_no_overlap_between_ccs(*p, ccs_j);//XXX this is just a test function, comment it the first run!
					data.checkAlignmentRange(*p);//XXX this is just a test function, comment it after the first run!
				}
			}
		}*/

		map<string, vector<pw_alignment> >al_of_a_ccs;//It is filled in with the original centers and is used for creating long centers & their als.(string ---> center, vector<pw_al> ---> alignments of a cluster)
		vector<vector<std::string> > local_long_centers;



		// compute sum of gains that go into clustering:
		for(size_t j=0; j<cc_cluster_in.size(); ++j) {
			for(std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc_cluster_in.at(j).begin(); it!=cc_cluster_in.at(j).end(); ++it) {
				const pw_alignment * al = *it;
				double g1,g2;
				m.gain_function(*al, g1, g2);
				double gain = (g1+g2)/2.0;
				sum_of_gain_before_clust += gain;
			}
		}


//AFP clustering
		for(size_t j=0; j<cc_cluster_in.size(); ++j) {
			clock_t ap_time_local = clock();
			std::cout << "cc cluster in at " << j  << "als number "<< cc_cluster_in.at(j).size()<< std::endl;
			use_affpro uaf(cc_cluster_in.at(j), m, cluster_base_cost);
			std::map<std::string, std::vector<string> >cluster_result;
			//Fill in a local container with centers and their alignments
			std::map<std::string, std::vector<pw_alignment> > local_al_in_a_cluster;
			std::cout << "afp clustering "<<std::endl;
			uaf.run(cluster_result,local_al_in_a_cluster);
			ap_time_local = clock() - ap_time_local;

			ap_time += ap_time_local;

			size_t counter =0;
			//Find alignments of each cluster:
		/*	for(std::map<std::string,std::vector<string> >::iterator it=cluster_result.begin();it !=cluster_result.end();it++){
				counter++;
				std::string center = it->first;
				std::cout << "center is "<<center << std::endl;
				std::map<std::string, std::vector<pw_alignment> >::iterator it1 = local_al_in_a_cluster.find(center);
				if(it1 == local_al_in_a_cluster.end()){
					local_al_in_a_cluster.insert(make_pair(center, std::vector<pw_alignment>()));
					it1 = local_al_in_a_cluster.find(center);
10637				}
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int dir = atoi(center_parts.at(0).c_str());
				unsigned int ref = atoi(center_parts.at(1).c_str());
				unsigned int left = atoi(center_parts.at(2).c_str());
				for(size_t i =0; i < it->second.size();i++){
					std::vector<std::string> member_parts;
					strsep(it->second.at(i), ":" , member_parts);
					unsigned int mem_dir = atoi(member_parts.at(0).c_str());
					unsigned int mem_ref = atoi(member_parts.at(1).c_str());
					unsigned int mem_left = atoi(member_parts.at(2).c_str());
					// look at all input alignments and determine which cluster it belongs to
					for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = cc_cluster_in.at(j).begin(); it2!=cc_cluster_in.at(j).end(); ++it2){
						const pw_alignment * al = *it2;
						unsigned int al_dir1;
						if(al->getbegin1() < al->getend1()){
							al_dir1 =0;
						}else{
							al_dir1 = 1;
						}
						unsigned int al_dir2;
						if(al->getbegin2() < al->getend2()){
							al_dir2 =0;
						}else{
							al_dir2 = 1;
						}
						size_t ref1 = al->getreference1();
						size_t ref2 = al->getreference2();
						size_t left1;
						size_t left2;
						size_t right1;
						size_t right2;
						al->get_lr1(left1,right1);
						al->get_lr2(left2,right2);
						if(al_dir1 == dir && ref1 == ref && left1 == left && al_dir2 == mem_dir && ref2 == mem_ref && left2 == mem_left){
							it1->second.push_back(*al);
							double gain = m.get_the_gain(*al, center);
							assert(gain > 0);
							break;
						}
						if(al_dir2 == dir && ref2 == ref && left2 == left && al_dir1==mem_dir&& ref1 == mem_ref && left1 == mem_left ){
							it1->second.push_back(*al);
							double gain = m.get_the_gain(*al, center);
							assert(gain > 0);
							break;
						}
					} // for cc_cluster_in_at set
				}
			}*/
			num_clusters += cluster_result.size();
			num_cluster_inputs_al += cc_cluster_in.at(j).size();
						
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it = local_al_in_a_cluster.begin(); it!=local_al_in_a_cluster.end(); ++it) {
				assert(it->second.size() != 0);
				if(it->second.size() != 0){
					al_of_a_ccs.insert(*it);//Here we keep two directions separately because we need them later in creating long alignments
					alignments_in_a_cluster.insert(*it);// Globally saves all the als.
					std::map<std::string,std::vector<string> >::iterator result_it = cluster_result.find(it->first);//local map
					num_cluster_members_al += 1 + it->second.size();
					num_cluster_seq+= 1 + result_it->second.size();
					global_results.insert(*result_it);
				}
	
			}
		} // for cluster_in set  




		std::cout<< "end of cluster in! " << al_of_a_ccs.size() <<std::endl;
// Long centers are created here!
		finding_centers centers(data);
		//From this point on orientation matters! It shoud be taken into account in counting the number of each center happening when we are looking for the long centers!
		centers.center_frequency(al_of_a_ccs,centerOnSequence);//All the centers of each sequence are detected and centerOnSequence is filled in. It contains all the cneters on each sequence. 
		std::cout << "centers are detected!"<<std::endl;
//Trees are created iteratively and centers are merged:
		merging_centers merg(data, centers);

	//	std::cout << "centers on seq 35 "<<std::endl;
	//	for(std::map<size_t, std::string>::iterator test= centerOnSequence.at(35).begin(); test != centerOnSequence.at(35).end(); test++){
	//		std::cout<< test->first << " "<< test->second << std::endl;
	//	}
	//XXX	merg.adding_new_centers(local_long_centers,centersPositionOnASeq);//After merging a group of centers we iteratively create new suffixes and a new tree based on them, then recalculate the gains again. At the end those with gain>0 go to the long_centers and centersPositionOnASeq.
		std::cout<< "new centers are made! " << local_long_centers.size() <<std::endl;
		for(size_t j  =0; j < local_long_centers.size();j++){
			long_centers.push_back(local_long_centers.at(j));
		}
		std::cout << "long center size is "<< long_centers.size() << std::endl;
	//XXX	merg.create_alignment(local_long_centers,new_centers,al_of_a_ccs,centersPositionOnASeq,centerOnSequence);//Here we fill in the 'new_centers' with longer alignments in the way that center is always the second ref
	//XXX	merg.get_merged_als(merged_als);
	}

	old_ias_time = clock() - old_ias_time;

	std::cout << " gain using old method " << sum_of_gain_before_clust << std::endl;
	std::cout << " Gain before IAS " << sum_of_input_gain << " after IAS new method " << gain_after_ias << std::endl;
	std::cout << " Times old " << (double)old_ias_time/CLOCKS_PER_SEC << " new " << (double)new_ias_time/CLOCKS_PER_SEC << std::endl;
// TODO  XXX 
//	exit(0);

	std::cout << "end of no fraction ccs"<<std::endl;
//Short centers are added to the long centers container.
	mixing_centers mix(data);
	mix.add_original_centers_to_long_centers(new_centers, alignments_in_a_cluster, merged_als);
	std::cout<<"original als were added! "<<std::endl;
	mix.add_nonaligned_regions(centerOnSequence,alignments_in_a_cluster);
	std::cout << "find_adjacencies_with_direction "<<std::endl;
	mix.find_adjacencies_with_direction(adjacencies_with_dir);
//	mix.test(alignments_in_a_cluster);
	std::cout<< "adjacencies: "<<std::endl;
	for(std::map<int, std::set<int> >::iterator it = adjacencies_with_dir.begin(); it != adjacencies_with_dir.end(); it++){
		std::set<int> nodes = it->second;
		std::cout<< it->first << " to: ";
		for(std::set<int>::iterator it1 = nodes.begin(); it1 != nodes.end() ; it1++){
			int node = *it1;
			std::cout << node << " ";
		}
		std::cout << " "<<std::endl;
	}
//Save the graph
	clock_t graph_ma_time = clock();
	// write graph result in dot, fasta and maf format
	write_graph wrgraph(data,mix);
	wrgraph.write_graph_dot(adjacencies_with_dir, dotout);
	wrgraph.write_graph_fasta(graphfasta,alignments_in_a_cluster);//Includes non-aligned regions
	write_graph_maf(graphout, alignments_in_a_cluster, data);//Includes non-aligned regions.
//	write_graph_maf_with_long_centers(graphout,new_centers,data);//TODO It can be changed in a way that includes long centers themselves as well.

	graph_ma_time = clock() - graph_ma_time;


	clock_t arithmetic_encoding_time = clock();
	//Remove centers with no associated member.
//	std::set<std::string> intermediate_center;
	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin(); it!=alignments_in_a_cluster.end();it++){
		assert(it->second.size() != 0);
	//	if(it->second.size()==0){
	//		intermediate_center.insert(it->first);
	//	}
	}
//	for(std::set<std::string>::iterator it = intermediate_center.begin();it !=intermediate_center.end();it++){
//		std::string cent = *it;
//		std::map<std::string, std::vector<pw_alignment> >::iterator it1 = alignments_in_a_cluster.find(cent);
//		alignments_in_a_cluster.erase(it1);				
//	}
//Defining weights of global clustering results! 
	std::map<std::string, unsigned int> weight;
	std::map<std::vector<std::string>, size_t> numberOfACenter;
	mix.calculate_centers_weight(alignments_in_a_cluster, global_results, weight,numberOfACenter);//Only short centers
	std::map<std::vector<std::string> , unsigned int> long_center_weight;
	mix.calculate_long_centers_weight(new_centers,long_center_weight,numberOfACenter);//On mixture of long and short centers //I need to keep all the short centers weight even if they are not in new_centers. 'cus theyare the ones which are encoded.

	std::map<std::string,std::string>member_of_cluster;//first string is a associated one and the second one is its center.The direction from assocciated memebrs are removed before arith.encoding and are not used in arith.encoding
	mix.find_members_of_clusters(member_of_cluster,alignments_in_a_cluster);//Only on short centers


//Data compression:
	typedef dynamic_encoder<use_model> use_encode;
	use_encode en(data,m,wrap);
	if(0!=encoding_out.compare("noencode")){
		dlib::entropy_encoder_kernel_1 * enc = new dlib::entropy_encoder_kernel_1();
		en.al_encode(weight,member_of_cluster,alignments_in_a_cluster,outs,*enc);//TODO fixed up to here.
	//	en.al_encode_with_long_center(centerOnSequence,long_center_weight,alignments_in_a_cluster, centersPositionOnASeq, member_of_cluster ,outs, *enc,new_centers);
		delete enc;
		outs.close();
	}

	arithmetic_encoding_time = clock() - arithmetic_encoding_time;

	std::cout<< "size of memeber_of_cluster is : "<< member_of_cluster.size()<<std::endl;
	std::cout << "Initial alignments sets summary:" << std::endl;
	std::cout << "Clustering summary: " << std::endl;
	std::cout << "Input: " << num_cluster_inputs_al << " pw alignments on " << num_cluster_seq <<" sequence pieces " <<std::endl;
	std::cout << "Output: " << num_clusters << " clusters containing " << num_cluster_members_al << " alignments " << std::endl;
	std::cout << "Time overview: " << std::endl;
	std::cout << "Read data " << (double)read_data_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Train models " << (double)train_models_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Affinity propagation clustering " << (double)ap_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write MSA-maf graph " << (double)graph_ma_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "Write compressed file " << (double)arithmetic_encoding_time/CLOCKS_PER_SEC << std::endl; 

	std::cout << "CC Init time " << (double)cc_init_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "CC Comp time " << (double)cc_comp_time/CLOCKS_PER_SEC << std::endl;
	return 0;
}
int do_decoding(int argc, char * argv[]){
	typedef mc_model use_model;
	typedef decoder<use_model> use_decode;
	if(argc < 4){
		usage();
		cerr << "Program: decoding" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* input binary compressed file from model" << std::endl;
		cerr << " * output txt file for saving retrievd sequences" << std::endl;		
	}
	std::string encoding_out(argv[2]);
	std::string decoding_out(argv[3]);
	std::ifstream in(encoding_out.c_str(),std::ifstream::binary);
	std::ofstream out(decoding_out.c_str());
	all_data data;
	use_model m(data);
//	test_encoder test;
	decoding_wrapper wrap;
	use_decode dec(data,m,wrap);
//Decompression
//	m.set_patterns(in);
//	std::cout << "patterns are decoded ";
//	m.set_alignment_pattern(in);
//	std::cout << "patterns are decoded! ";
	dlib::entropy_decoder_kernel_1  decode;
//	dec.set_flag_from_stream(in);
//	arithmetic_decoding_centers(in,dec);
//	dec.al_decode_with_long_center(in,decode,out);
//	dec.al_decode_long_center_optimized_flag(in,decode,out);
	dec.al_decoding(in,decode,out);
	cout<< "decoding is done!"<<endl;
	return 0;	
}
int do_dynamic_decoding(int argc, char * argv[]){
	typedef dynamic_mc_model use_model;
	typedef dynamic_decoder<use_model> use_decode;
	if(argc < 4){
		usage();
		cerr << "Program: dy_decoding" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* input binary compressed file from model" << std::endl;
		cerr << " * output file for saving retrievd sequences" << std::endl;		
	}
	std::string encoding_out(argv[2]);
	std::string decoding_out(argv[3]);
	std::ifstream in(encoding_out.c_str(),std::ifstream::binary);
	std::ofstream out(decoding_out.c_str());
	all_data data;
	size_t num_threads = 1;
	wrapper wrapp;//TODO Think a bout a solution for capturing the context without calling wrapp in the model class. It makes the encode text file empty and one has to run encoding one more time before comparing. It is possible to add the context in to encoding and remove it from the model. But isnt it any better way?
	use_model m(data,wrapp,num_threads);
//	test_encoder test;
	decoding_wrapper wrap;
	use_decode dec(m,wrap);
	//Decompression
//	m.set_patterns(in);
//	std::cout << "patterns are decoded ";
//	m.set_alignment_pattern(in);
//	std::cout << "patterns are decoded! ";
	dlib::entropy_decoder_kernel_1  decode;
//	dec.set_flag_from_stream(in);
//	dec.arithmetic_decoding_centers(in,decode);
//	dec.al_decode_with_long_center(in,decode,out);
//	dec.al_decode_long_center_optimized_flag(in,decode,out);
//	dec.al_decoding(in,decode,out);
	dec.al_decode_with_long_center(in,decode,out);
	cout<< "decoding is done!"<<endl;
	out.close();

	return 0;	
}

int do_compare_sequences(int argc, char* argv[]){
	if(argc < 5){
		usage();
		cerr << "Program: compare" <<std::endl;
		cerr << "Parameters: "<<std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* file contains decoded sequences" <<std::endl;
	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string decoding_out(argv[4]);
	std::ifstream in(decoding_out.c_str(),std::ifstream::binary);
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	data.compare_seq_with_decoding(in);
	std::cout<< "Comparing is done!"<<std::endl;
	return 0;
}
int do_compare_high(int argc, char* argv[]){//XXX with the dynamic model, encoding needs to be run again after decoding :(
	test_encoder t;
	if(argc < 4){
		usage();
		cerr << "Program: Highcompare" <<std::endl;
		cerr << "Parameters: "<<std::endl;
		cerr << "* file with low and high from encoder" << std::endl;
		cerr << "* file with low and high from decoder" << std::endl;
	}
	std::string encodefile(argv[2]);
	std::string decodefile(argv[3]);
	std::ifstream in(encodefile.c_str());
	std::ifstream in1(decodefile.c_str());
	t.compare(in,in1);//checking the high values in encoding and decoding
//	t.context_compare(in,in1);//checking encoded and decoded bases and modification patterns.
//	t.al_high_compare(in,in1);
	return 0;
}

int do_read_data(int argc, char * argv[]) {
	if(argc < 4) {
		usage();
		cerr << "Program: read_data" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;

	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);

//	std::string samfile(argv[3]);

// Read all data
	all_data data;
	typedef overlap overlap_type;
	data.read_fasta_maf(fastafile, maffile);

//	data.read_accknown_fasta_sam(fastafile,samfile);
//	std::vector<pw_alignment> alignments = data.getAlignments();
	size_t seq_length = 0;
	for(size_t i =0; i < data.numSequences();i++){
		seq_length += data.getSequence(i).length();
		if(i == 4){
			for(size_t j = 518; j < 584; j++){
				std::cout << data.getSequence(i).at(j)<< " " ;
			}
			std::cout<< "" <<std::endl;
		}
	}
	std::cout << "total base number is " << seq_length <<std::endl;
	return 0;
}
int do_simple_test_on_new_cc(int argc, char* argv[]){
	if(argc < 2) {
		cerr << "Program: simple_test_cc" << std::endl;
	}
	std::vector<std::pair<size_t,size_t> > common_int;
	common_int.push_back(make_pair(0,3)); //form left to right +1
	common_int.push_back(make_pair(2,5));
//	common_int.push_back(make_pair(1,4));
	common_int.push_back(make_pair(5,8));
	common_int.push_back(make_pair(4,7));
	common_int.push_back(make_pair(8,11));
	std::vector<std::set<size_t> > ccs;
	compute_cc_with_icl<overlap> cccs(common_int);
	cccs.compute_test(ccs);
	size_t counter = 0;
	for(size_t i=0; i<ccs.size(); ++i) {
		std::cout << " on initial CC " << i << " size " << ccs.at(i).size() << std::endl;
		for(std::set<size_t>::iterator it1 = ccs.at(i).begin(); it1 != ccs.at(i).end(); it1++){
				std::cout << *it1 <<std::endl;
		}
		counter += ccs.at(i).size();
	}	
	std::cout<< "returned number is "<< counter <<std::endl;

	return 0;

}
int do_test_new_cc(int argc, char * argv[]) {
	typedef overlap_interval_tree overlap_type;
	typedef dynamic_mc_model use_model;
//	typedef compute_cc_with_interval_tree<overlap_type> cc_type; //Eric Garisson implementation
//	typedef compute_cc_with_icl cc_type;
//	typedef compute_cc_avl<overlap_type> cc_type;
	typedef compute_cc_avl_fraction cc_type;
	typedef initial_alignment_set<use_model,overlap_type> use_ias;
	if(argc < 5) {
		usage();
		cerr << "Program: test_cc" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* output binary compressed file " << std::endl;
		cerr << "* number of threads to use (optional, default 1)" << std::endl;
		return 1;

	}
	size_t num_threads = 1;
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string encoding_out(argv[4]);
	std::ofstream outs(encoding_out.c_str(),std::ofstream::binary);
	if(argc>5) {
		num_threads = atoi(argv[5]);
	}
// Read all data
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	overlap_type ol(data);
	wrapper wrap;
// Train the model on all data
	use_model m(data,wrap, num_threads);
	m.train(outs);
	std::cout << "model is trained! "<<std::endl;
// Find connected components of alignments with some overlap // Moved it here after training to be able to remove small als first and then connect them to each other.
	std::set<const pw_alignment*, compare_pointer_pw_alignment> al_with_pos_gain; 
	std::set<const pw_alignment*, compare_pointer_pw_alignment> remainders; 
#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i < data.numAlignments();i++){
//	for(size_t i =0; i < 50;i++){
//	if(i%1000==0) std::cout << ".";
//	for(size_t i =0; i < 16;i++){
		const pw_alignment & al = data.getAlignment(i);
		double g1 ,g2;
		m.gain_function(al,g1,g2);
		double av_gain = (g1+g2)/2 ;
		if(av_gain > 0){
#pragma omp critical(al_insert)
{
		al_with_pos_gain.insert(&al);
	//	std::cout<< &al<<std::endl;
	//	al.print();

}
		}
	}
	std::cout << "al with positive gains are kept " << al_with_pos_gain.size() <<std::endl;
	std::cout << std::endl;
	size_t num_seq = data.numSequences();
	filter_als filter(m,al_with_pos_gain, num_seq, num_threads);
	filter.find_overlapped_references();//XXX slow!
	std::set<const pw_alignment*,  compare_pointer_pw_alignment> filtered_als;
	filtered_als = filter.get_filtered_als();
	//Makes components of partially overlapped alignments
	std::cout << "filtered als size "<< filtered_als.size() << std::endl;
//	cc_type component(al_with_pos_gain,data.numSequences());
//	compute_cc_with_icl cccs(al_with_pos_gain, data.numSequences(),num_threads);
	std::set< pw_alignment, compare_pw_alignment> mixed_als; 
	std::cout << " now we compute connected components " << std::endl;
	clock_t cc_init_time = clock();
	cc_type cccs(filtered_als, data.numSequences(), num_threads);//The new one with fraction
	cc_init_time = clock() - cc_init_time;
	std::cout << "CC Init time " << (double)cc_init_time/CLOCKS_PER_SEC << std::endl;
//	std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> > ccs; 
	size_t cluster_base_cost =0;
	
	remainders = filtered_als;
	clock_t cc_comp_time = clock();
	double fraction = 0.95;
//	std::set< pw_alignment, compare_pw_alignment> al_in_overlap; //TODO

	while(fraction > 0.60){
//	for(std::set<pw_alignment,compare_pw_alignment>::iterator it = al_in_overlap.begin(); it != al_in_overlap.end();it++){
//		remainders.insert(*it); // changed the 'remainders' to reference!
//	}
//	for(std::set<pw_alignment,compare_pw_alignment>::iterator it = remainders.begin(); it != remainders.end();it++){
//		const pw_alignment & al = *it;
//		double g1 ,g2;
//		m.gain_function(al,g1,g2);
//		double av_gain = (g1+g2)/2 ;
//		assert(av_gain>0);
//	}
//	al_in_overlap.clear();
//	cc_type cccs(remainders, data.numSequences(), num_threads);
	std::set<const pw_alignment* , compare_pointer_pw_alignment> ccs; 
//	cccs.compute(ccs); //fill in ccs, ordered by size(Notice that they are only connected to each other if they have overlap more than 95% overlap!)
	std::cout << "fraction is "<< fraction << " run on " << remainders.size()<<std::endl;
	cccs.non_recursive_compute(ccs,remainders,fraction); //A graph of als with 95% overlap is made.
	remainders = cccs.get_remaining_als();//XXX is faster if fill it in at the same time with running "non_recursive_compute"

//	mixed_als = cccs.get_remaining_als();//TODO add mixed als to the previous function and fill it in at the same time
	cc_comp_time = clock() - cc_comp_time;

	std::cout << "CC Init time " << (double)cc_init_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "CC Comp time " << (double)cc_comp_time/CLOCKS_PER_SEC << std::endl;
//	component.compute(ccs);
	std::cout<< "ccs are made!" << ccs.size() <<std::endl; //Initial ccs with 95% overlap
	std::cout << "mixed als size1 "<< mixed_als.size() << std::endl;
	size_t counter = 0;
//	for(size_t i=0; i<ccs.size(); ++i) { // Only those have more than 95% overlap
//	for(size_t i=0; i<ccs.size()-1; ++i) {
//		std::cout << " on initial CC " << i << " size " << ccs.at(i).size() << std::endl;
//		counter += ccs.at(i).size();
	/*	for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.at(i).begin();it != ccs.at(i).end();it++){
			std::cout << *it << std::endl;

			const pw_alignment * p = *it;
//			p->print();
			for(size_t j = i+1; j < ccs.size();j++){
				std::set< const pw_alignment* , compare_pointer_pw_alignment> ccs_j = ccs.at(j);
				cccs.test_redundancy(p,ccs_j);
//				ol.test_no_overlap_between_ccs(*p, ccs_j);//XXX this is just a test function, comment it the first run!
				data.checkAlignmentRange(*p);//XXX this is just a test function, comment it after the first run!
			}
		}*/
//	}

//	std::cout<< "returned number is "<< counter <<std::endl;
	//Finding biconnected components
/*	biconnected_component bicomp(cccs,ccs);
	for(size_t i=0; i<ccs.size(); ++i) { //
		std::cout << "ccs at " << i << std::endl;
		bicomp.creat_biconnected_component(i);
		std::cout << "D"<<std::endl;
	}*/
	std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment > > stacks;
	std::set< const pw_alignment* , compare_pointer_pw_alignment> & temp = ccs;
//	for(size_t i=0; i<ccs.size(); ++i) {
		std::cout << "ccs size "<< ccs.size() <<std::endl;
		two_edge_cc twoEdgeCC(cccs,ccs);
		twoEdgeCC.find_bridges(stacks,mixed_als);
		std::cout << "mixed als size2 "<< mixed_als.size() << std::endl;


//	}
	size_t Counter =0;
	std::cout << "stack size in main is " << stacks.size()<< std::endl;
	for(size_t i=0; i<stacks.size(); ++i) {
	//	std::cout << "at "<< i <<std::endl;
		std::set< const pw_alignment*, compare_pointer_pw_alignment> & stack = stacks.at(i);
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = stack.begin();it != stack.end();it++){
		//	std::cout << *it << std::endl;
	//		const pw_alignment * al = *it;
	//		double g1 ,g2;
	//		m.gain_function(*al,g1,g2);
	//		double av_gain = (g1+g2)/2 ;
	//		assert(av_gain > 0);
			Counter ++;
		}
		//	std::cout << "at "<< i << " "<< Counter <<std::endl;		
	}
	std::cout << "total number of stacks' member is "<< Counter <<std::endl;

			
	//Cutting partial overlaps in 'stacks' 
	std::vector<overlap_type> cc_overlap(stacks.size(), overlap_type(data));
//	std::set<const pw_alignment*, compare_pointer_pw_alignment> al_pointers_in_overlap;
	size_t COUNT = 0;
	for(size_t i=0; i<stacks.size(); ++i){//TODO Order them by size before cutting
		std::cout << "at stck " << i << " with "<< stacks.at(i).size() << " als "<<std::endl;
		std::set< const pw_alignment*, compare_pointer_pw_alignment> & stack = stacks.at(i);
	//	std::cout << "temp size before "<< temp.size()<<std::endl;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = stack.begin();it != stack.end();it++){ //TODO it is not thoughtful! Fillin the remainder while fill in the stacks!
			std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it1 = temp.find(*it);
			assert(it1 != temp.end());
			temp.erase(*it1);
		}
	//	std::cout << "temp size after "<< temp.size()<<std::endl;

	//	std::cout << "stack size at " << i << " is " << stack.size()<<std::endl;
		use_ias ias(data,stack, m, cluster_base_cost,num_threads);
		ias.compute(cc_overlap.at(i));
		const std::set< pw_alignment, compare_pw_alignment> & als_in_overlap = cc_overlap.at(i).get_all();
		COUNT += als_in_overlap.size();
	//	std::cout<< "number of returned als after cutting is "<< als_in_overlap.size() << std::endl;
		for(std::set<pw_alignment,compare_pw_alignment>::const_iterator it = als_in_overlap.begin(); it != als_in_overlap.end();it++){
		//	std::cout<< &*it<<std::endl;
		//	al_pointers_in_overlap.insert(&*it);
			mixed_als.insert(*it);
		//	remainders.insert(&*it);
		//	al_in_overlap.insert(*it);
		}
	}
	std::cout << "total returned numbers after cutting is " << COUNT <<std::endl;
	for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = temp.begin(); it != temp.end();it++){
		mixed_als.insert(**it);
	}
	std::cout << "mixed als size3 "<< mixed_als.size() << std::endl;
//	remainders = cccs.get_remainders();
	fraction -= 0.05;
}//Add the remaining als to mixed_als(can use remainder container in main)
	for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = remainders.begin(); it != remainders.end();it++){
		mixed_als.insert(**it);
	}
//	for(std::set<pw_alignment,compare_pw_alignment>::iterator it = al_in_overlap.begin(); it != al_in_overlap.end();it++){
//		mixed_als.insert(*it);
//	}

	std::cout<< "mixed als size is "<< mixed_als.size() << std::endl; //Now we need to run cc without fraction and cut them afterwards
	for(std::set< pw_alignment , compare_pw_alignment>::iterator it = mixed_als.begin(); it != mixed_als.end();it++){
		data.checkAlignmentRange(*it);
	}

	//Creating connected components with no fraction rate over all of them and cut them afterwards.
	typedef compute_cc_avl<overlap_type> ccc_type;
	ccc_type no_fraction_cccs(mixed_als, data.numSequences(), num_threads);
	std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> > no_fraction_ccs;
	no_fraction_cccs.compute(no_fraction_ccs); 
	std::cout << "number of components is " << no_fraction_ccs.size() <<std::endl;
	size_t NUM = 0;
	for(size_t i=0; i<no_fraction_ccs.size(); ++i) {
		std::cout << " on CC " << i << " size " << no_fraction_ccs.at(i).size() << std::endl;
	//	for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = no_fraction_ccs.at(i).begin();it != no_fraction_ccs.at(i).end();it++){
		//	std::cout << *it << std::endl;
			NUM+=no_fraction_ccs.at(i).size();
//			const pw_alignment * p = *it;
		//	for(size_t j = i+1; j < no_fraction_ccs.size();j++){
			//	std::set< const pw_alignment* , compare_pointer_pw_alignment> ccs_j = no_fraction_ccs.at(j);
			//	ol.test_no_overlap_between_ccs(*p, ccs_j);//XXX this is just a test function, comment it the first run!
			//	data.checkAlignmentRange(*p);//XXX this is just a test function, comment it after the first run!
		//	}
	//	}
	}
	std::cout << "NUM " << NUM <<std::endl;

	//Cutting partial overlaps:
	std::vector<overlap_type> ccc_overlap(no_fraction_ccs.size(), overlap_type(data));
	for(size_t i=0; i<no_fraction_ccs.size(); ++i) {
		std::set< const pw_alignment*, compare_pointer_pw_alignment> & cc = no_fraction_ccs.at(i);
		use_ias ias(data,cc, m, cluster_base_cost,num_threads);
		ias.compute(ccc_overlap.at(i));
		std::vector< std::set<const pw_alignment* , compare_pointer_pw_alignment> > cc_cluster_in; //have shorter length than original alignemnts
		std::cout<< "using the cc again!"<<std::endl;
	//	size_t num_seq = data.numSequences();
		ccc_type occ(ccc_overlap.at(i), num_seq ,1);
		occ.compute(cc_cluster_in);
	}
	//Cutting partial overlaps:
/*	std::vector<overlap_type> cc_overlap(no_fraction_ccs.size(), overlap_type(data));
	size_t cluster_base_cost =0;
	for(size_t i=0; i<no_fraction_ccs.size(); ++i) {
		std::set< const pw_alignment*, compare_pointer_pw_alignment> & cc = no_fraction_ccs.at(i);
		std::cout << "this cc size is "<< cc.size() << std::endl;
	//	for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc.begin();it != cc.end();it++){
	//		const pw_alignment * p = *it;
	//		p->print();
	//	}
		detect_overlap_rate ov_rate(data.numSequences(),cc);
	//	std::vector<std::multimap<size_t , size_t> > adjucents; // each nodes and its adjucents
		std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > strongly_connected;
		ov_rate.compute(strongly_connected);
		std::cout<< "strongly connected size is " << strongly_connected.size() <<std::endl;
		for(size_t i =0;i< strongly_connected.size();i++){
			std::cout << "at i " << i <<std::endl;
			for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = strongly_connected.at(i).begin();it != strongly_connected.at(i).end();it++){
				const pw_alignment * p = *it;
				p->print();
				size_t id = ov_rate.get_id(p);
				std::cout << id << std::endl;
			}
		}

		exit(0);
//		assert(adjucents.size() == strongly_connected.size());
		for(size_t j =0; j < strongly_connected.size();j++){
			biconnected_component bicomp(ov_rate,strongly_connected.at(j));
		}
//		use_ias ias(data,cc, m, cluster_base_cost,num_threads);
		std::set< const pw_alignment* , compare_pointer_pw_alignment> less_redundent;
	//	std::cout <<"repeats size is "<< repeats.size()<<std::endl;
	//	for(size_t j =0; j< repeats.size();j++){
	//		std::cout<< repeats.at(j).size()<<std::endl;
	//		for(std::set<const pw_alignment*>::iterator it = repeats.at(j).begin();it != repeats.at(j).end();it++){
	//			const pw_alignment* al = *it;
	//			al->print();
	//		}
	//		use_ias ias(data,repeats.at(j), m, cluster_base_cost,num_threads);
	//		ias.compute(cc_overlap.at(i));
			std::cout<< "DONE on CC "<< i << " with " << cc_overlap.at(i).size() << " fragments "<<std::endl;
		//	std::set<pw_alignment, compare_pw_alignment> al_in_overlap = cc_overlap.at(i).get_all();
		//	for(std::set<pw_alignment,compare_pw_alignment>::iterator it = al_in_overlap.begin(); it != al_in_overlap.end();it++){
		//		const pw_alignment & p1 = *it;
		//		p1.print();
		//	}
	//	}
	//	ias.compute(cc_overlap.at(i));
	//	std::cout<< "DONE on CC "<< i << " with " << cc_overlap.at(i).size() << " fragments "<<std::endl;
	//	cc_overlap.at(i).test_partial_overlap();
		std::cout << "no partial overlap! "<<std::endl;
		std::vector< std::set<const pw_alignment* , compare_pointer_pw_alignment> > cc_cluster_in; 
	//	cc_type occ(cc_overlap.at(i), data.numSequences(),1);
	//	occ.compute(cc_cluster_in);*/
	//	std::cout << "cluster in size: " << cc_cluster_in.size() <<std::endl;
	/*	for(size_t j =0; j < cc_cluster_in.size();j++){
			for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = cc_cluster_in.at(j).begin(); it != cc_cluster_in.at(j).end();it++){
				const pw_alignment * p = *it;
				p->print();
				double g1,g2;
				m.gain_function(*p,g1,g2);
				assert((g1+g2)/2 >= 0);
				std::cout << "g1 " << g1 << " g2 " << g2 << std::endl;
			}
		}*/

//	}
	std::cout << "cutting is over!"<<std::endl;
	std::cout << "CC Init time " << (double)cc_init_time/CLOCKS_PER_SEC << std::endl;
	std::cout << "CC Comp time " << (double)cc_comp_time/CLOCKS_PER_SEC << std::endl;
//	std::cout<< "returned number is "<< counter <<std::endl;
	return 0;
}


int do_model_simple(int argc, char * argv[]){//simple model
/*	typedef model use_model;
	typedef initial_alignment_set<use_model> use_ias;
	typedef affpro_clusters<use_model> use_affpro;
	if(argc < 5) {
		usage();
		cerr << "Program: simple_model" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;
		cerr << "* output maf file for the graph" << std::endl;
		cerr << "* number of threads to use (optional, default 10)" << std::endl;
	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);
	std::string graphout(argv[4]);
	size_t num_threads = 1;
	if(argc == 6) {
		num_threads = atoi(argv[5]);
	}
//Reading data:
	all_data data;
	data.read_fasta_maf(fastafile, maffile);
	overlap ol(data);
	std::ofstream outs;
// Find connected components of alignments with some overlap
	compute_cc cccs(data);
	std::vector<std::set< pw_alignment*, compare_pointer_pw_alignment> > ccs;
	cccs.compute(ccs);
//Train all the sequences:
	use_model m(data);
	m.train();
	std::vector<overlap> cc_overlap(ccs.size(), overlap(data));
	// base cost to use an alignment
//	double cluster_base_cost = log2(data.numAlignments());
	double cluster_base_cost = 5.0;
	std::cout << " base cost " << cluster_base_cost << std::endl;

// Select an initial alignment std::set for each connected component
	for(size_t i=0; i<ccs.size(); ++i) {
		std::set< pw_alignment*, compare_pointer_pw_alignment> & cc = ccs.at(i);
		use_ias ias(data,cc, m, cluster_base_cost);
		ias.compute(cc_overlap.at(i));
	std::vector< std::set<pw_alignment , compare_pw_alignment> > cc_cluster_in; //has shorter length than original sequence pieces
	compute_cc occ(cc_overlap.at(i), data.numSequences());
	occ.compute(cc_cluster_in);
	std::cout<< "cc_cluster_in size: "<< cc_cluster_in.size()<<std::endl;
	std::cout << "cc " << i << ": from " << cc.size() << " original als with total gain " << ias.get_max_gain() << " we made " << cc_overlap.at(i).size() << " pieces with total gain " << ias.get_result_gain() <<  std::endl; 
	std::vector<std::string> all_centers;
	map<string, vector<pw_alignment> >al_of_a_ccs;// it is used for creating long centers (string ---> center, vector<pw_al> ---> alignments of a cluster)
	for(size_t j=0; j<cc_cluster_in.size(); ++j) {
		std::cout << " run affpro at " << j << " on " << cc_cluster_in.at(j).size()<<std::endl;	
		std::cout << " in al std::set size " << cc_cluster_in.at(j).size() << std::endl;
		size_t nums=0;
		for(std::set<pw_alignment, compare_pw_alignment>::iterator ait=cc_cluster_in.at(j).begin(); ait!=cc_cluster_in.at(j).end(); ++ait) {
		//	std::cout << " nums " << nums << std::endl;
			const pw_alignment & cc_al = *ait;
			nums++;
		//	cc_al.print();
			std::cout << std::endl;
			double c1;
			double c2;
			double m1;
			double m2;
			m.cost_function(cc_al, c1, c2, m1,m2);
			std::cout << "c1 " << c1 << " c2 " << c2 << " m1 " << m1 << " m2 "<< m2 <<std::endl;
		}
		
		use_affpro uaf(cc_cluster_in.at(j), m, cluster_base_cost,outs);
		std::map<std::string, std::vector<string> >cluster_result;
		uaf.run(cluster_result);
		for(std::map<std::string,std::vector<string> >::iterator it=cluster_result.begin();it != cluster_result.end();it++){
			std::cout << "center " << it->first << std::endl;
		}
	}
}




	*/
	return 0;
}
int do_daligning(int argc, char * argv[]){//TODO need to be reviewed! for some reads more than one las file is produced. have to figure it out why
	/* 1. Read the read file, includes all the reads.
	   2. Align one read at the time over the given graph(Daligner)*/
	typedef dynamic_mc_model use_model;
	if(argc < 6){
		usage();//TODO need to change it!!
		cerr << "Program: run_daligner" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* trimmed read fasta file" << std::endl;
		cerr << "* fasta file contains the reference graph" << std::endl;
		cerr << "* output sam file includes all the alignment" << std::endl;
		cerr << "* read number "<<std::endl;
		cerr << "num_threads" <<std::endl;
		return 1;
	}
	size_t num_threads = 1;
	if(argc ==7){
		num_threads = atoi(argv[6]);
	}
	std::string readfastafile(argv[2]);
	std::string graphfasta(argv[3]);
	std::string samfile(argv[4]);
	size_t read_number = atoi(argv[5]);
	all_data data;
	daligner dalign;
	data.read_fasta(readfastafile);//Read all the trimmed reads.
	std::vector<std::string> reads;
	reads = data.get_reads();
	std::cout << "reads size is "<< reads.size() << std::endl;

//Create refdb here:
	std::string cmline = "fasta2DAM";
	cmline += " ";
	cmline +="./dbref";
	cmline +=" ";
	cmline += graphfasta.c_str();
	system(cmline.c_str());
	//TODO make it works over all the reads!!
	//make a read_db for a single read
	dalign.make_dbread(reads.at(read_number),read_number);
	std::string cmline_daligner= "runDaligner.sh";
	cmline_daligner +=" ";
	cmline_daligner +=samfile.c_str();
	system(cmline_daligner.c_str());
	////XXX just a test to see if the sam file is not empty
	//std::string checkSamFile =  "grep -cvq '^@' ";
	//checkSamFile += samfile.c_str();
	//if(system(checkSamFile.c_str())){
	//	std::cout << "Sam file empty "<<std::endl;
	//	continue;
	//}
	return 0;
}

int do_read_map(int argc, char * argv[]){
	/*Gets read fasta, graph fasta, alignment maf file and read+graph fasta as an input
	  For each read makes alignment graph and detect the best path on the graph.
	  Note that all the fasta files are output of fasta_prepare*/
	typedef dynamic_mc_model use_model;
	typedef overlap_interval_tree overlap_type;
	if(argc < 8){
		usage();//TODO need to change it!!
		cerr << "Program: map_read" << std::endl;
		cerr << "Parameters:" << std::endl;
		cerr << "* Reads fasta file" << std::endl; //after running fasta_prepare //Note that reads names should be unique
		cerr << "* Ref(graph) fasta file" << std::endl; //after running fasta_prepare
		cerr << "* Ref(graph) dot file" << std::endl;
		cerr << "* Read+Ref fasta file" << std::endl; //after running fasta_prepare //Attention! it has to be made from original ref and read file before fasta_prepare
		cerr << "* alignment maf file "<<std::endl; //I used LAST
		cerr << "* output a fasta file " << std::endl; 
		cerr << "* num_threads(optional, defult is 1)" <<std::endl;
		return 1;
	}
	size_t num_threads = 1;
	if(argc ==9){
		num_threads = atoi(argv[8]);
	}
	std::string readfasta(argv[2]);
	std::string reffasta(argv[3]);
	std::string refdotfile(argv[4]);
	std::string allfasta(argv[5]);
	std::string alignmentmaf(argv[6]);
	std::string outputfasta(argv[7]);
	std::ofstream output(outputfasta.c_str());
	all_data data = all_data();
	data.read_fasta_maf(allfasta,alignmentmaf);
	wrapper wrapp;
	use_model m(data,wrapp,num_threads);
	std::ofstream out("noencode", std::ofstream::binary);
	m.train(out);
	std::cout << data.numSequences()<<std::endl;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> al_with_pos_gain; 
	std::string ref_accession;
#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i < data.numAlignments();i++){
		const pw_alignment & al = data.getAlignment(i);
		if(i == 0){
			size_t acc = data.accNumber(al.getreference1());
			std::cout << "ref accession is " << acc <<std::endl;
			ref_accession = data.get_acc(acc);
		}
		double g1 ,g2;
		m.gain_function(al,g1,g2);
		double av_gain = (g1+g2)/2 ;
		if(av_gain > 0){
#pragma omp critical(al_insert)
{
			al_with_pos_gain.insert(&al);
}
		}
	}
//	needleman<use_model> needle(data, m);
	std::map< std::string, size_t> longname2seqidx;
	longname2seqidx = data.getLongname2seqidx();
	std::cout<<longname2seqidx.size()<<std::endl;
	std::cout<< "read dot file "<<std::endl;
	Graph graph;
	ref_graph rgraph(graph,data);
	rgraph.read_dot_file(refdotfile,ref_accession);//Fill in the vetix and edge container with the reference graph nodes and adjacencies
	graph.print();
	std::cout << "read read.fasta "<<std::endl;
	data.read_fasta(readfasta);//Read all the reads.
	std::vector<std::string> reads;
	reads = data.get_reads();
	std::cout << "reads size is "<< reads.size() << std::endl;
	std::vector<std::vector<pw_alignment> > all_als(reads.size());
	for(size_t i =0; i< data.numAlignments();i++){//Read all the als and fill in all_als vector with als of each read
		const pw_alignment & p = data.getAlignment(i);
		size_t ref2 = p.getreference2();
		size_t acc = data.accNumber(p.getreference1());
		if(acc == 1){
			all_als.at(ref2).push_back(p);	
		}else{
			assert(acc == 0);
			size_t num = data.get_seq_number_of_acc(acc);
			all_als.at(ref2-num).push_back(p);
		}
	}
	for(size_t i =0; i < all_als.size();i++){
		std::cout <<"read i  "<< i <<std::endl;
		std::vector<pw_alignment> als = all_als.at(i);
		std::cout<<als.size()<<std::endl;
	}
//Goes through all the reads and for each of them find the best path using Dijkestra graph, the best path should be save in the output
	for(size_t i =0; i < all_als.size();i++){
		std::cout <<"read i  "<< i <<std::endl;
		std::vector<pw_alignment> als = all_als.at(i);
		if(als.size()>100){//If the read aligned on more than 100 places on the ref , we reduce the number of als by choosing the best ones only
			std::multimap<double, pw_alignment> sorting_als;
			std::vector<pw_alignment> chosen_als;
			for(size_t j =0; j < als.size();j++){
				const pw_alignment & p = als.at(j);
				sorting_als.insert(std::make_pair(p.get_modify1()/data.get_seq_size(i), p));
			}
			size_t limit = 70;
			size_t count =0;
			for(std::multimap<double, pw_alignment>::reverse_iterator it = sorting_als.rbegin(); it!= sorting_als.rend();it++){
				while(count < limit){
					chosen_als.push_back(it->second);
					count++;
				}
			}
			als = chosen_als;
		}
		//TODO from here on! 
		std::map<int, std::set<int> > adjacencies = rgraph.get_adjacencies();
		deep_first df(data, adjacencies);
		dnastring current_read = data.getSequence(i);
		std::cout << "als size "<< als.size() <<std::endl;
	/*	for(size_t j =0; j < als.size(); j++){
			pw_alignment p = als.at(j);
			size_t ref1 = p.getreference1();
			std::string seqname = data.get_seq_name(ref1);
			int name = std::stoi(seqname);
			if(p.getbegin1()>p.getend1()){
				name = -1*name;
			}
			std::cout << "name is " << name <<std::endl;
			df.deep_first_search(name, ref_accession);
			std::vector<std::vector<int> > paths = df.get_paths();
			std::cout<< "path size is "<< paths.size()<<std::endl;
			for(size_t i =0; i < paths.size(); i++){
				for(size_t j =0; j < paths.at(i).size();j++){
					std::cout << paths.at(i).at(j)<< " ";
				}
				std::cout<< " "<<std::endl;
			}

		}*/
		std::cout << "here!!!!!"<<std::endl;
		typedef compute_cc_avl<overlap_type> ccc_type;
		ccc_type no_fraction_cccs(als, data.numSequences(), num_threads);
		std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> > no_fraction_ccs;
		no_fraction_cccs.compute(no_fraction_ccs); 
		std::cout << "number of components is " << no_fraction_ccs.size() <<std::endl;

		als_components ccs(data, df, no_fraction_ccs);//check the distance between the last right and the first left of the next component
		ccs.merg_components();
		ccs.find_als_on_paths();

		exit(1);
		Graph newgraph;
		model m1(data);//XXX temp, because Marion used simple model, later one i need to change everything in her code make it work with the dynamic model(just for getting her code works)
		graph.prepare_read(newgraph,als, current_read,m1,data,output);
	}
	std::cout << "done! "<<std::endl;




	return 0;
}
int needleman_wunsch(int argc, char * argv[]){
	if(argc < 4){
		usage();
		cerr << "Program: needleman_test" << std::endl;
		cerr << "* fasta file from fasta_prepare" << std::endl;
		cerr << "* maf file containing alignments of sequences contained in the fasta file" << std::endl;

		
		return 1;
	}
	std::string fastafile(argv[2]);
	std::string maffile(argv[3]);

	size_t num_threads = 1;
	std::string ref = "ACTGCTGATTC";
	std::string read = "GACTGGATTC";
	all_data data;
	data.read_fasta_maf(fastafile, maffile);

	wrapper wrap;
	std::ofstream outs("testoutput");
	dynamic_mc_model m(data,wrap, num_threads);
	m.train(outs);
	needleman<dynamic_mc_model> n(data, m, read, ref);
	size_t readid = 1;
	size_t refid = 0;
	n.compute_matrix(readid, refid);
	n.print_path();
	size_t type=1;//TODO have to fix it for the other types as soon as i mix it with the other part of the mapping code and running it on real data.
	std::string readstr;
	std::string refstr;
	n.find_the_best_path(type, readstr, refstr);



	std::map<int, std::set<int> > adjacencies;
//	adjacencies.insert(make_pair(1,2));
//	adjacencies.insert(make_pair(1,3));
//	adjacencies.insert(make_pair(3,4));
//	adjacencies.insert(make_pair(3,5));
//	adjacencies.insert(make_pair(4,5));
//	adjacencies.insert(make_pair(4,6));
//	adjacencies.insert(make_pair(4,7));




	deep_first df(data, adjacencies);

	return 0;
}

#if !TEST
	
int main(int argc, char * argv[]) {

	if(argc < 2) {
		usage();
		exit(1);
	}
	std::string program(argv[1]);

	if(0==program.compare("fasta_prepare")) {
		return do_fasta_prepare(argc, argv);
	}
	else if(0==program.compare("model")){
		return do_mc_model(argc, argv);
	}
	else if(0==program.compare("decoding")){
		return do_decoding(argc,argv);
	}
	else if(0==program.compare("dy_decoding")){
		return do_dynamic_decoding(argc,argv);
	}
	else if(0==program.compare("simple_model")){
		return do_model_simple(argc,argv);
	}
	else if(0==program.compare("dy_model")){
		return do_dynamic_mc_model(argc,argv);
	}
	else if(0==program.compare("new_dy_model")){
		return do_dynamic_mc_model_with_two_edge(argc,argv);
	}
	else if(0 ==program.compare("compare")){
		return do_compare_sequences(argc,argv);
	}
	else if(0 ==program.compare("Highcompare")){
		return do_compare_high(argc,argv);
	}
	else if (0== program.compare("read_data")){
		return do_read_data(argc,argv);
	}
	else if(0 ==program.compare("test_cc")){
		return do_test_new_cc(argc,argv);
	}
	else if(0 == program.compare("simple_test_cc")){
		return do_simple_test_on_new_cc(argc,argv);
	}
	else if(0==program.compare("run_daligner")){
		return do_daligning(argc,argv);
	}
	else if(0==program.compare("map_read")){
		return do_read_map(argc,argv);
	}
	else if(0==program.compare("needleman_test")){
		return needleman_wunsch(argc,argv);
	}
	else {
		usage();
	}

	
	return 1;
}



#endif


#include "model.cpp"
#include "overlap.cpp"
#include "encoder.cpp"
#include "connected_component.cpp"
#include "intervals.cpp"
#include "dynamic_encoder.cpp"
#include "dynamic_decoder.cpp"
#include "needleman_wunsch.cpp"
#include "filter_alignments.cpp"
#include "alignment_index.cpp"





















