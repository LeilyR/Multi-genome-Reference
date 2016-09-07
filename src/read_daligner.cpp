#include "read_daligner.hpp"
#ifndef DALIGNER_CPP
#define DALIGNER_CPP



	void daligner::make_dbread(std::string & single_read, size_t & read_id){

		std::ofstream single_read_fasta;
		single_read_fasta.open("single_read.fasta");
		single_read_fasta << "> "<<read_id <<std::endl;
		size_t rest = 0;
		if(single_read.size()>70){
			for(size_t i =0; i < single_read.size()-70;i++){
				for(size_t j =i; j< i+70;j++){
					single_read_fasta<<single_read.at(j);
				}
				i +=69;
				rest = i+1;
				single_read_fasta<<std::endl;
			}
			for(size_t i =rest; i < single_read.size(); i++){
				single_read_fasta<<single_read.at(i);
			}
			single_read_fasta<<std::endl;			
		}else{
			single_read_fasta<<single_read<<std::endl;
		}
		single_read_fasta.close();

		std::string cmline = "fasta2DAM";
		cmline += " ";
		cmline +="./dbread";
		cmline +=" ";
		cmline += "single_read.fasta";
		system(cmline.c_str());

	}

#endif
