#include "needleman_wunsch.hpp"

#ifndef NEEDLEMAN_WUNSCH.CPP
#define NEEDLEMAN_WUNSCH.CPP

//	template<typename T>	
//	needleman<T>::needleman(all_data & d, T & m):data(d),model(m){
//
//	}

	template<typename T>
	needleman<T>::needleman(const all_data & d, T & m , std::string & read , std::string & ref):data(d),model(m),mod_matrix(read.length()+1 , std::vector<double>(ref.length()+1)),score_matrix(read.length()+1 , std::vector<double>(ref.length()+1)), context_matrix(read.length()+1 , std::vector<std::string>(ref.length()+1)), keep_number(read.length()+1 , std::vector<size_t>(ref.length()+1)), delete_number(read.length()+1 , std::vector<size_t>(ref.length()+1)), path(read.length()+1 , std::vector<size_t>(ref.length()+1)){
		this->read = read;
		this->ref = ref;
		std::cout << "this ref "<< this->ref << std::endl;
		std::cout<< "this read" << this->read <<std::endl;
	}

	template<typename T>
	needleman<T>::~needleman(){}

	template<typename T>
	void needleman<T>::run_needleman(size_t & readacc , size_t & refacc, size_t & type, std::string & out_read, std::string & out_ref){
		compute_matrix(readacc , refacc);
		find_the_best_path(type,out_read,out_ref);
	}
	template<typename T>
	void needleman<T>::compute_matrix(size_t & readacc, size_t & refacc){
	//	size_t readacc = data.accNumber(read_id);
	//	size_t refacc = data.accNumber(ref_id);
		const std::map< std::string, std::vector<double> > model_cost = model.get_al_cost(refacc, readacc);
		//Fill in the first cell:
		std::string context = get_first_context();
		mod_matrix.at(0).at(0) = 0.0;
		score_matrix.at(0).at(0) = 0.0;
		keep_number.at(0).at(0) = 0;
		delete_number.at(0).at(0) = 0;
		context_matrix.at(0).at(0) = context;
		path.at(0).at(0) = 0;
		//Filling in the first row(Adding deletes at each position)
		std::string first_row;//updated pattern of the first row
		std::cout << "ref length "<< ref.length()<<std::endl;
		for(size_t j = 1; j <= ref.length() ; j++){
		//	std::cout<< "at " << j <<std::endl;
			//Find the context + the last base:
			context+=(char) dnastring::base_to_index(ref.at(j-1));
			std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(context);
			assert(it != model_cost.end());
			//Find delete of length j+1
			char numdelete;
			size_t delete_length = j;
			count_first_row_delete(delete_length, numdelete,first_row); //Update the delete length return the last one as numdelete //TODO make it more efficient!
			keep_number.at(0).at(j)= 0;
			delete_number.at(0).at(j)= j;
			mod_matrix.at(0).at(j)= it->second.at((size_t)numdelete);
			double delcost = 0.0;
			std::string last_context_before_delete = get_first_context();
			last_context_before_delete += (char) dnastring::base_to_index(ref.at(0));
		//	std::cout<<"first pattern "<< last_context_before_delete <<std::endl;
			delete_cost(model_cost,last_context_before_delete ,first_row,delcost);
			score_matrix.at(0).at(j) = delcost;
		//	for(size_t k =0; k < first_row.size();k++){
		//		std::cout<< (size_t)first_row.at(k)<< " ";
		//	}
		//	std::cout<<""<<std::endl;

			context = first_row.at(first_row.size()-MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
			for(size_t i = 1; i < MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL; i++){
				context += first_row.at(first_row.size()-MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL+i);
			}
			assert(context.size()== MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
			context_matrix.at(0).at(j) = context;
			path.at(0).at(j) = 2;
		//	std::cout << "j is "<< j << std::endl;
		}
	//	std::cout<< "first row was filled"<<std::endl;
	//	print_context_matrix();	
	//	print_score_matrix();
		context = get_first_context();
		//Filling in the first column (Adding insertion at position i)
		std::string first_column;//Pattern on first column
 		first_column =get_first_context();
		for(size_t i = 1; i <= read.length();i++){
			context += (char) dnastring::base_to_index(ref.at(0));//There is no actual base on the ref before this insertion. For the moment i keep the first base on ref.
			std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(context);
			assert(it != model_cost.end());
			char insert = model.modification_character(-1, -1, dnastring::base_to_index(read.at(i-1)),-1);	
			mod_matrix.at(i).at(0) = it->second.at((size_t)insert);
			score_matrix.at(i).at(0) = mod_matrix.at(i).at(0) + score_matrix.at(i-1).at(0);
			first_column += insert;
			//calculate the next context 
			context = first_column.at(first_column.size()-MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
			for(size_t k = 1; k < MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL; k++){
				context += first_column.at(first_column.size()-MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL+k);
			}
			context_matrix.at(i).at(0) = context;
			keep_number.at(i).at(0) = 0;
			delete_number.at(i).at(0) = 0;
			path.at(i).at(0) = 3;
		}
	//	std::cout<<"first column was filled "<<std::endl;
		//Filling in the rest of matrices:
		for(size_t i=1; i <= read.size() ; i++){
		//	std::cout << " at i = "<<i<<std::endl;
			for(size_t j =1; j<=ref.size(); ++j){
				std::vector<double> mod;
				std::vector<double> score;
				std::vector<std::string> temp_context;
				size_t deletenumber =0;
				size_t keepnumber = 0;
			//	std::cout<< "at j = "<< j<<std::endl;
				bool ifkeep = false;
				if(read.at(i-1)==ref.at(j-1)){//If match
				//	std::cout<<"match!"<<std::endl;
					double previous_value =mod_matrix.at(i-1).at(j-1);	
					std::string previous_context = context_matrix.at(i-1).at(j-1);
				//	std::cout << (size_t)previous_context.at(previous_context.length()-1)<<std::endl;
					if(previous_context.at(previous_context.length()-1) >= NUM_KEEP_DYN + NUM_DELETE_DYN + 5 || previous_context.at(previous_context.length()-1)< NUM_DELETE_DYN + 5){//If the last one is not a keep
						ifkeep = true;
						keepnumber++;
						context = previous_context;
						context +=(char) dnastring::base_to_index(ref.at(j-1));
						std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(context);
						assert(it != model_cost.end());
						double value = it->second.at(NUM_DELETE_DYN + 5);//keep of length 1
						mod.push_back(value);//First save them localy, then the maximum is saved globaly
						score.push_back(value+score_matrix.at(i-1).at(j-1));
					//	std::cout<< "score "<< value << " + " << score_matrix.at(i-1).at(j-1)<<std::endl;
						previous_context.erase(previous_context.begin());
						std::string current =previous_context;
						current +=model.modification_character(-1,-1, -1, 0);
						temp_context.push_back(current);
					//	std::cout << "match context"<<std::endl;
					//	for(size_t k =0; k < current.size();k++){
					//		std::cout<< (size_t)current.at(k)<<" ";
					//	}
					//	std::cout << " " <<std::endl;
					}else{//if the last one is a keep
						ifkeep = true;
						keepnumber = keep_number.at(i-1).at(j-1);
						keepnumber++;//Don't add it to keep_number because not sure if this is the best context yet.
						std::string last_nonkeep_context = context_matrix.at(i-keepnumber).at(j-keepnumber);//It shouldn't be ended by keep. Check it:
					//	std::cout << "keep number " << keepnumber << " last char " << (size_t)last_nonkeep_context.at(last_nonkeep_context.size()-1) << std::endl;
						assert(last_nonkeep_context.at(last_nonkeep_context.size()-1)< NUM_DELETE_DYN + 5 || last_nonkeep_context.at(last_nonkeep_context.size()-1) >= NUM_DELETE_DYN + NUM_KEEP_DYN + 5);
						char current_mod;
						std::string current_context;
						std::string all_keeps;
						size_t keepnum = keepnumber;
						update_keep(keepnum, last_nonkeep_context, current_context , current_mod, all_keeps);//update the keep length return the current context.
						context = previous_context;
						context += (char) dnastring::base_to_index(ref.at(j-1));
						std::map< std::string, std::vector<double> >::const_iterator it= model_cost.find(context);
						assert(it != model_cost.end());
						double value = it->second.at((size_t)current_mod);
						double keepcost = 0.0;
						last_nonkeep_context += (char) dnastring::base_to_index(ref.at(j-keepnumber));
						keep_cost(model_cost,last_nonkeep_context,all_keeps,keepcost);
						score.push_back(keepcost + score_matrix.at(i-keepnumber).at(j-keepnumber));
						temp_context.push_back(current_context);
						mod.push_back(value);
					}
				}
				else{ //If mismatch
				//	std::cout<<"mismatch!"<<std::endl;
					assert(read.at(i-1)!=ref.at(j-1));
					double previous_value =mod_matrix.at(i-1).at(j-1);
					std::string previous_context = context_matrix.at(i-1).at(j-1);	
					std::string context;
					context+=previous_context;
					context += (char)dnastring::base_to_index(ref.at(j-1));
					std::map< std::string, std::vector<double> >::const_iterator it= model_cost.find(context);
					assert(it != model_cost.end());
					double value = it->second.at(dnastring::base_to_index(read.at(i-1)));//Mismatches are the first four rows of the it->second;
					score.push_back(value + score_matrix.at(i-1).at(j-1));
					previous_context.erase(previous_context.begin());
					std::string current_context = previous_context;
					current_context += dnastring::base_to_index(read.at(i-1));
					temp_context.push_back(current_context);
					mod.push_back(value);
				}
				//if deletion(horizontal, simillar to the match case)
				double previous_value =mod_matrix.at(i).at(j-1);	
				std::string previous_context = context_matrix.at(i).at(j-1);
			//	std::cout << (size_t)previous_context.at(previous_context.length()-1)<<std::endl;
				if(previous_context.at(previous_context.length()-1) >= NUM_DELETE_DYN + 5 || previous_context.at(previous_context.length()-1)< 5){//If the last one is not a delete
				//	std::cout<< "no del before it"<<std::endl;
					deletenumber++;
					context = previous_context;
					context +=(char) dnastring::base_to_index(ref.at(j-1));
					std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(context);
					assert(it != model_cost.end());
					double value = it->second.at(5);//Delete of length 1
					mod.push_back(value);//First save them localy, then the maximum is saved globaly
					score.push_back(value+score_matrix.at(i).at(j-1));
				//	std::cout<< "score "<< value << " + " << score_matrix.at(i).at(j-1)<<std::endl;
					previous_context.erase(previous_context.begin());
					std::string current =previous_context;
					current +=model.modification_character(-1,0, -1, -1);
					temp_context.push_back(current);
				//	std::cout <<"del context "<<std::endl;
				//	for(size_t k =0; k < current.size();k++){
				//		std::cout<< (size_t)current.at(k)<<" ";
				//	}
				//	std::cout << " " <<std::endl;

				}else{
					assert(previous_context.at(previous_context.length()-1) < NUM_DELETE_DYN + 5 && previous_context.at(previous_context.length()-1) >= 5); //previous context ended with delete
					deletenumber = delete_number.at(i).at(j-1);
				//	std::cout<<"del num "<<deletenumber << " " << delete_number.at(i).at(j-1)<<std::endl;
					deletenumber++;//Don't add it to delete_number because not sure if this is the best context yet.
					std::string last_nondelete_context = context_matrix.at(i).at(j-deletenumber);//It shouldn't be ended by delete. Check it:
				//	std::cout << "i " << i <<" j "<< j<< " del num " << deletenumber<< " last context "<< (size_t)last_nondelete_context.at(last_nondelete_context.size()-1) <<std::endl;
					assert(last_nondelete_context.at(last_nondelete_context.size()-1) >= NUM_DELETE_DYN + 5 || last_nondelete_context.at(last_nondelete_context.size()-1) < 5);
					char current_mod;
					std::string current_context;
					std::string all_deletes;
					size_t delnum = deletenumber;
				//	std::cout<<"delete of length "<<delnum<<std::endl;
					update_delete(delnum, last_nondelete_context, current_context , current_mod, all_deletes);//updates the delete length return the current context.
					context = previous_context;
					context += (char)dnastring::base_to_index(ref.at(j-1));
					std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(context);
					assert(it != model_cost.end());
					double value = it->second.at((size_t)current_mod);
					last_nondelete_context += (char)dnastring::base_to_index(ref.at(j-deletenumber));
					double long_delete_cost = 0.0;
					delete_cost(model_cost,last_nondelete_context,all_deletes,long_delete_cost);
					score.push_back(long_delete_cost + score_matrix.at(i).at(j-deletenumber));
					temp_context.push_back(current_context);
				//	std::cout <<"del context "<<std::endl;
				//	for(size_t k =0; k < current_context.size();k++){
				//		std::cout<< (size_t)current_context.at(k)<<" ";
				//	}
				//	std::cout << " " <<std::endl;

					mod.push_back(value);
				}
				//if insertion(vertical)
				previous_value =mod_matrix.at(i-1).at(j);	
				previous_context = context_matrix.at(i-1).at(j);
				context = previous_context;
				context +=(char) dnastring::base_to_index(ref.at(j-1));
				std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(context);
				assert(it != model_cost.end());
				double value = it->second.at(5+NUM_DELETE_DYN+NUM_KEEP_DYN+dnastring::base_to_index(read.at(i-1)));
				score.push_back(value+score_matrix.at(i-1).at(j));
			//	std::cout<< "score "<< value << " + " << score_matrix.at(i-1).at(j)<<std::endl;
				previous_context.erase(previous_context.begin());
				std::string current_context = previous_context;
				current_context += dnastring::base_to_index(read.at(i-1));
				temp_context.push_back(current_context);
				mod.push_back(value);
			//	std::cout<<"insert context "<<std::endl;
			//	for(size_t k =0; k < current_context.size();k++){
			//		std::cout<< (size_t)current_context.at(k)<<" ";
			//	}
			//	std::cout << " " <<std::endl;

				//Keep the max of these three values
				add_minimums(i,j,score,mod,temp_context,keepnumber,deletenumber, ifkeep);
			//	print_score_matrix();
			//	std::cout << " "<<std::endl;
			//	print_context_matrix();
			}
		}
		if(ref.length()==1){
			print_score_matrix();
			std::cout << " "<<std::endl;
			print_context_matrix();
		}
	}
	template<typename T>
	void needleman<T>::add_minimums(size_t & row , size_t & column , std::vector<double> & score , std::vector<double> & mod , std::vector<std::string> & context, size_t & keep , size_t & del, bool& ifkeep){
		double score_min = score.at(0);
		size_t min =0;
	//	std::cout<<"scores: "<<std::endl;
		for(size_t i = 0; i < score.size();i++){
		//	std::cout<< score.at(i)<<std::endl;
			if(score.at(i)<score_min){
				score_min = score.at(i);
				min = i;
			}
		}
	//	std::cout << "min is "<<min <<std::endl;
		score_matrix.at(row).at(column)=score_min;
		mod_matrix.at(row).at(column) = mod.at(min);
		context_matrix.at(row).at(column) = context.at(min);
		if(ifkeep == true && min == 0){
		//	std::cout << "keep at "<< row << " at " << column << " : " <<keep <<std::endl;

			keep_number.at(row).at(column) = keep;
		}
		if(min == 0 && ifkeep == false){
			keep_number.at(row).at(column) = 0;
		}
		if(min == 1){
		//	std::cout << "del at "<< row << " at " << column << " : " <<del <<std::endl;
			delete_number.at(row).at(column)=del;
		}
		if(min!=1){
			delete_number.at(row).at(column)=0;
		}
		add_path(row,column,min);
	}
	template<typename T>
	void needleman<T>::add_path(size_t & row, size_t & column, size_t& type){
		if(type == 0){
			path.at(row).at(column) = 1; //keep/modification(diagonal)
		}
		if(type == 1){
			path.at(row).at(column) = 2;//delete(horizontal)
		}
		if(type == 2){
			path.at(row).at(column) = 3;//insert(vertical)
		}
	}
	template<typename T>
	std::string needleman<T>::get_first_context(){//Changed from keeps as we had in the model class to modifications, because Keeps can cause an issue in the length of the pattern by merging short keeps into one longer one.
		std::string first_context;
		for(size_t i= MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL; i>0; --i) {
			char modc = model.modification_character(0, -1, -1, -1);
			first_context.append(1, modc);
		}	
		return first_context;
	}
	template<typename T>
	void needleman<T>::count_first_row_delete(size_t & delete_length , char & numdelete,std::string & first_row){//TODO make it better by adding the largest possible deletes, keeping them and not counting them again.
	//	std::cout<< "length of del is "<< delete_length<<std::endl;
		first_row=get_first_context();
		while(delete_length > (1<<(NUM_DELETE_DYN-1))){
			first_row += model.modification_character(-1, NUM_DELETE_DYN-1 , -1, -1);
			delete_length -= (1<<(NUM_DELETE_DYN-1));
		}
		for(size_t i = NUM_DELETE_DYN; i>0;i--){
			if(delete_length == 0) break;
			if((delete_length & 1<<(i-1)) != 0){
				numdelete = model.modification_character(-1, i-1 , -1, -1);
				delete_length-=1<<(i-1);
				first_row += numdelete;
			}
		}
	}
	
	template<typename T>
	void needleman<T>::update_keep(size_t & numkeep , std::string & nonkeep_context, std::string & context, char & last_keep, std::string & all_keeps){
		//Find the length of previous keep, add one to it and update the value
		all_keeps += nonkeep_context;
		count_keep(numkeep, all_keeps);
		last_keep = all_keeps.at(all_keeps.size()-1);
		assert(all_keeps.size() > MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
		for(size_t i = MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL; i > 0 ; i--){
			context += all_keeps.at(all_keeps.size()-i);
		}	
	}
	template<typename T>
	void needleman<T>::count_keep(size_t & keep,std::string & previous){
		//encode the length
		while(keep > (1<<(NUM_KEEP_DYN-1))){
			previous += model.modification_character(-1, -1, -1, NUM_KEEP_DYN-1);
			keep -= 1<<(NUM_KEEP_DYN-1);
		}
		for(size_t i = NUM_KEEP_DYN; i >0;i--){
			if((keep & (1<<(i-1)))!=0){
				previous += model.modification_character(-1, -1, -1, i-1);
				keep-=1<<(i-1);
				if(keep == 0) break;
			}
		}
		
	}
	template<typename T>
	void needleman<T>::update_delete(size_t & deletenumber, std::string & last_nondelete_context, std::string & current_context , char & current_mod, std::string & all_deletes){
		//Find the length of previous delete, add one to it and update the value
		all_deletes += last_nondelete_context;
		count_delete(deletenumber, all_deletes);
	//	std::cout<<"all deletes"<<std::endl;
	//	for(size_t i = 0; i < all_deletes.size();i++){
		//	std::cout<< (size_t)all_deletes.at(i)<< " ";
	//	}
	//	std::cout<<""<<std::endl;
		current_mod = all_deletes.at(all_deletes.size()-1);
		assert(all_deletes.size() > MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL);
		for(size_t i = MAX_ALIGNMENT_MARKOV_CHAIN_LEVEL; i > 0 ; i--){
			current_context += all_deletes.at(all_deletes.size()-i);
		}	

	}
	template<typename T>
	void needleman<T>::count_delete(size_t & deletenum ,std::string & context){
		//encode the length
		while(deletenum > (1<<(NUM_DELETE_DYN-1))){
			context += model.modification_character(-1, NUM_DELETE_DYN-1, -1, -1);
			deletenum -= 1<<(NUM_DELETE_DYN-1);
		}
		for(size_t i = NUM_DELETE_DYN; i >0;i--){
			if((deletenum & (1<<(i-1)))!=0){
			//	std::cout <<"inja "<< i-1<<std::endl;
				context += model.modification_character(-1, i-1, -1, -1);
				deletenum-=1<<(i-1);
				if(deletenum == 0) break;
			}
		}	
	}
	
	template<typename T>
	void needleman<T>::delete_cost(const std::map< std::string, std::vector<double>  >& model_cost, std::string & last_context, std::string & pattern, double & cost){
		cost = 0;
	//	std::cout << "last context: "<<last_context << " its size "<< last_context.size() <<std::endl;
		std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(last_context);
		assert(it != model_cost.end());
	//	std::cout<< pattern.size()<<std::endl;
		for(size_t i = pattern.size(); i > 0; i--){
		//	std::cout << (size_t)pattern.at(i-1)<<std::endl; 
			if((size_t)pattern.at(i-1) >=5 && (size_t)pattern.at(i-1)<5+NUM_DELETE_DYN){
				cost += it->second.at((size_t)pattern.at(i-1));
			}
			else break;

		}


	}
	template<typename T>
	void needleman<T>::keep_cost(const std::map< std::string, std::vector<double>  >& model_cost, std::string & last_nonkeep_context, std::string & all_keeps, double & keepcost){
		std::map< std::string, std::vector<double>  >::const_iterator it= model_cost.find(last_nonkeep_context);
		assert(it != model_cost.end());
		for(size_t i = all_keeps.size(); i > 0; i--){
			if(all_keeps.at(i-1) >=5+NUM_DELETE_DYN && all_keeps.at(i-1)<5+NUM_DELETE_DYN+NUM_KEEP_DYN){
				keepcost += it->second.at((size_t)all_keeps.at(i-1));
			}
			else break;
		}
	}

	template<typename T>
	void needleman<T>::find_the_best_path(size_t & type, std::string & from_read , std::string & from_ref){
		std::vector<std::vector<double> > copy_score = score_matrix;
		std::cout << "read length "<< read.length() << " ref length "<< ref.length() << std::endl;
		size_t row,column;
		if(type == 1 || type == 2){//when both start and end or only the end is fixed. We need to start from the last cell in the matrix.
			row = read.length();
			column = ref.length();
		}else{
			assert(type ==3 || type == 4); //when the start on read is fixed but the end is not.
			//choose the smallest score on the last row
			double min = score_matrix.at(read.length()).at(0);
			column = 0;
			for(size_t i =1; i < ref.length()+1; i++){
				std::cout << score_matrix.at(read.length()).at(i) <<std::endl;
				if(score_matrix.at(read.length()).at(i)<min){
					min = score_matrix.at(read.length()).at(i);
					column = i;
				} 
			}
			row = read.length();
			std::cout << "start here! "<<column<<std::endl;
		}
	/*	if(type == 3 && column > 0){//XXX Get sure if it is necessary to keep the gap on read !
			for(size_t i =ref.length()-1; i > column-1 ; i--){
				from_read += '-';
				from_ref += ref.at(i);	
			}
		}*/
		if(type == 2){
			//set the first row to zero
			for(size_t i =0; i < ref.length()+1; i++){
				copy_score.at(0).at(0) = 0.0;
			}
		}
		size_t last_part;
		while(row != 0){
		//	std::cout << path.size() << " " << path.at(0).size() <<" row "<< row << " column "<< column << std::endl;
			last_part= path.at(row).at(column);
			if(last_part == 1){
			//	std::cout<< "match/mismatch"<<std::endl;
				from_read += read.at(row-1);
				from_ref += ref.at(column-1);
				row--;
				column--;
			}
			if(last_part == 2){
			//	std::cout << "del"<<std::endl;
				from_read += '-';
				from_ref += ref.at(column-1);
				column--;
			}
			if(last_part == 3){
			//	std::cout<< "insertion"<<std::endl;
				from_read += read.at(row-1);
				from_ref += '-';
				row--;
			}
		}
		if(row == 0 && type ==2){
			std::cout << "at type two! "<<std::endl;
		}
		if(row == 0 && (type ==1 || type ==3 || type == 4)){
			std::cout<< "HERE AT TYPE 3 "<<  "column "<<column  << " row " << row <<std::endl;
			while(column != 0){
				from_read += '-';
				from_ref += ref.at(column-1);
				column--;
			}
		}
		std::reverse(from_ref.begin(), from_ref.end());
		std::reverse(from_read.begin(), from_read.end());

		std::cout<<from_ref <<std::endl;
		std::cout << from_read<<std::endl;
	}
	template<typename T>
	void needleman<T>::print_score_matrix(){
		std::cout<< "\t"<<"\t";
		for(size_t i =0; i < ref.length();i++){
			std::cout<<ref.at(i)<< "\t";
		}
		std::cout << "" <<std::endl;
		for(size_t i =0; i < score_matrix.size();i++){
			if(i != 0){
				std::cout<<read.at(i-1)<< "\t";
			}else{
				std::cout<< "\t";
			}
			for(size_t j =0; j < score_matrix.at(i).size(); j++){
				std::cout<< score_matrix.at(i).at(j) << "\t";
			}
			std::cout << "" <<std::endl;
		}

	}
	template<typename T>
	void needleman<T>::print_context_matrix(){
		std::cout<< "\t"<<"\t";
		for(size_t i =0; i < ref.length();i++){
			std::cout<<ref.at(i)<< "\t";
		}
		std::cout << "" <<std::endl;
		for(size_t i =0; i < context_matrix.size();i++){
			if(i != 0){
				std::cout<<read.at(i-1)<< "\t";
			}else{
				std::cout<< "\t";
			}
			for(size_t j =0; j < context_matrix.at(i).size(); j++){
				std::string str = context_matrix.at(i).at(j);
				for(size_t k =0; k<str.length();k++){
					std::cout<< (size_t)str.at(k) << " ";
				}
				std::cout<<"\t";
			}
			std::cout << "" <<std::endl;
		}

	}
	template<typename T>
	void needleman<T>::print_path(){
		std::cout<< "\t"<<"\t";
		for(size_t i =0; i < ref.length();i++){
			std::cout<<ref.at(i)<< "\t";
		}
		std::cout << "" <<std::endl;
		for(size_t i =0; i < path.size();i++){
			if(i != 0){
				std::cout<<read.at(i-1)<< "\t";
			}else{
				std::cout<< "\t";
			}
			for(size_t j =0; j < path.at(i).size(); j++){
				std::cout<< path.at(i).at(j) << "\t";
			}
			std::cout << "" <<std::endl;
		}

	}





#endif
