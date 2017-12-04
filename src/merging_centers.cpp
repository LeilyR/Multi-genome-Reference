#include "merging_centers.hpp"

#ifndef MERGING_CENTERS_CPP
#define MERGING_CENTERS_CPP


	void merging_centers::updating_centers(std::vector<int> & center_string, int & index){//replace 'center_string' with its new 'index' where ever 'center_sting' occurs on the tree and updates the number of happening. 
		//First we make a new tree ! 
		size_t num_seq = data.numSequences();
		suffixTree tree(num_seq, centers, all_current_centers);
		tree.make_tree(center_string,index);
		for(size_t seq_id = 0; seq_id < num_seq; seq_id++){
			std::cout << seq_id <<std::endl;
			all_current_centers.at(seq_id) = tree.get_current_centers(seq_id);
		}
		std::map<size_t,std::vector<int> > edges = tree.get_edges();
		std::map<std::vector<size_t> , size_t> counts = tree.get_branches();
		std::map<std::vector<int>, int> intermediate;
		std::map<std::vector<int>, size_t> centers_frequency;
		for(std::map<std::vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			vector<size_t> br = it->first;//list of nodes of a path
			size_t number = it->second;//number of the path happening
			std::vector<int> seriesOfCenters;
			std::vector<int> reverse_centers;
		//	std::cout<< "centers on these nodes are "<<std::endl;
		//	std::cout<< "series of centers: " <<std::endl;
			for(size_t j = 0 ; j < br.size(); j++){
				std::map<size_t, std::vector<int> >::iterator nodes= edges.find(br.at(j));
				assert(nodes != edges.end());
				for(size_t k =0; k < nodes->second.size();k++){
					if(nodes->second.at(k) != 1 <<31){
					//	std::cout<< nodes->second.at(k) <<std::endl;
						seriesOfCenters.push_back(nodes->second.at(k));//push back all the centers of all the nodes on the path br
						reverse_centers.push_back(-1 * nodes->second.at(k));
					}
				}
			}
			std::reverse(reverse_centers.begin(),reverse_centers.end());
		/*	if(seriesOfCenters.at(seriesOfCenters.size()-1) == 1<<31){//The extra ending character is removed 
				std::cout << "romve the last char!"<<std::endl;
				seriesOfCenters.pop_back();	
			}*/

			std::map<std::vector<int>, int>::iterator it1 = gains.find(seriesOfCenters);
			std::map<std::vector<int>, int>::iterator rev_it1 = gains.find(reverse_centers);
			if(it1 != gains.end()){//It is kept as it is
				if(rev_it1 == gains.end()){
					std::cout << "rev: "<< reverse_centers << std::endl;
					std::cout << "cent: " << seriesOfCenters << std::endl;
				}
				assert(rev_it1 == gains.end());
				intermediate.insert(make_pair(it1->first,it1->second));
			}else if(rev_it1 != gains.end()){
				assert(it1 == gains.end());
				intermediate.insert(make_pair(rev_it1->first,rev_it1->second));
			}else{
				assert(it1 == gains.end() && rev_it1 == gains.end());
			//	size_t number_of_new_center = 0;
			//	for(size_t i = 0; i < seriesOfCenters.size(); i++){
			//		if(seriesOfCenters.at(i)== index){//??
			//			number_of_new_center = number_of_new_center + 1;
			//		}
			//	}
		//		size_t old_length = seriesOfCenters.size()+ (number_of_new_center*center_string.size());
		//		std::cout << "old_length " << old_length << std::endl;
//				size_t gain = number*old_length - (number + seriesOfCenters.size());
			/*	size_t gain = number*seriesOfCenters.size() - (number + seriesOfCenters.size());
				intermediate.insert(make_pair(seriesOfCenters,gain)); */
			//	gains.insert(make_pair(seriesOfCenters,gain));
				if(seriesOfCenters.size() >1){
					std::map<std::vector<int> , size_t>::iterator freq=centers_frequency.find(seriesOfCenters);
					std::map<std::vector<int> , size_t>::iterator freq1=centers_frequency.find(reverse_centers);
					if(freq == centers_frequency.end() && freq1 == centers_frequency.end()){ //non of the 'centers' and its reverse have been saved. 
						centers_frequency.insert(std::make_pair(seriesOfCenters,number));
					}else if(freq == centers_frequency.end() && freq1 != centers_frequency.end()){
						freq1->second += number;
					}else{
						assert(freq != centers_frequency.end() && freq1 == centers_frequency.end());
						freq->second += number;
					}
				}	
			}
		}
		gains.clear();
		for(std::map<std::vector<int>,int>::iterator it = intermediate.begin(); it != intermediate.end(); it++){
			gains.insert(make_pair(it->first,it->second));
		}
		for(std::map<std::vector<int> , size_t>::iterator it = centers_frequency.begin(); it != centers_frequency.end(); it++){
			int gain = (it->second * it->first.size()) - (it->second + it->first.size());
			gains.insert(std::make_pair(it->first , gain));
		}
	}

	void merging_centers::merg_gain_value(suffixTree & tree){//calculates the gain value and builds second tree and so on iteratively
		std::map<size_t,std::vector<int> > edges = tree.get_edges();
		std::cout<< "size of the tree " << edges.size()<<std::endl;
		size_t original_center_numbers = centers.get_number_of_centers();
		std::cout << "original center number: "<< original_center_numbers << std::endl; 
		std::map<std::vector<size_t> , size_t> counts = tree.get_branches();//vector<size_t> shows a path(size_t s are node indices) and size_t is its number of happening. 
		//calculating the initial gain values:
		std::map<std::vector<int> , size_t> centers_frequency;
		for(std::map<std::vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			std::vector<size_t> br = it->first;//list of nodes of a path
			size_t number = it->second;
			std::vector<int> centers;
			std::vector<int> rev_centers;
			int gain = 0;
			for(size_t j = 0 ; j < br.size(); j++){
				std::map<size_t, std::vector<int> >::iterator nodes= edges.find(br.at(j));
				assert(nodes != edges.end());
				for(size_t k =0; k < nodes->second.size();k++){
					if(nodes->second.at(k) !=  1<<31){
						centers.push_back(nodes->second.at(k));//push back all the centers of all the nodes on the path br
						rev_centers.push_back( -1*(nodes->second.at(k)));
					}
				}
			}
			std::reverse(rev_centers.begin(),rev_centers.end());
		/*	if(centers.at(centers.size()-1) == 1<<31){//The extra ending character is removed 
				centers.pop_back();	
			}*/
			std::cout<< "centers " ;
			for(size_t j =0; j < centers.size();j++){
				std::cout<< centers.at(j) << " " ;
			}
			std::cout << " " <<std::endl;
		//	std::cout<< "center size: " << centers.size() << "  number of happening: "<< number << std::endl;
		/*	gain = number*(centers.size())-(number+centers.size());
			if(centers.size()!= 0){// We may get zero when the original path only had the extra endign char.
				gains.insert(make_pair(centers,gain));//centers-->index of centers, gain of the long center
				std::cout<< "gain "<<gain <<std::endl;
			} */
			if(centers.size() > 1){
				std::map<std::vector<int> , size_t>::iterator freq=centers_frequency.find(centers);
				std::map<std::vector<int> , size_t>::iterator freq1=centers_frequency.find(rev_centers);
				if(freq == centers_frequency.end() && freq1 == centers_frequency.end()){ //non of the 'centers' and its reverse have been saved. 
					centers_frequency.insert(std::make_pair(centers,number));
					std::cout << "here_freq!"<<std::endl;
				}else if(freq == centers_frequency.end() && freq1 != centers_frequency.end()){
					freq1->second += number;
					std::cout << "here_freq1!"<<std::endl;
				}else{
					assert(freq != centers_frequency.end() && freq1 == centers_frequency.end());
					freq->second += number;
					std::cout << "here_freq2!"<<std::endl;
				}
			}
		}
		for(std::map<std::vector<int> , size_t>::iterator it = centers_frequency.begin(); it != centers_frequency.end(); it++){
			int gain = (it->second * it->first.size()) - (it->second + it->first.size());
			gains.insert(std::make_pair(it->first , gain));
		}
		int highest_gain=0;
		std::vector<int> highest_path;
		std::vector<int> seen1;
		find_highest_gain(highest_path, highest_gain,seen1);
		if(highest_gain > 0){
			merged_centers.insert(make_pair(highest_path, original_center_numbers+1));//Making the first new center with a new index which is 'original_center_numbers+1'
		}
		int center_numbers;
		center_numbers = original_center_numbers + 1;//it will be used when we are going to insert the next megerd center to the merged_centers map.
		while(highest_gain > 0){
			std::cout << "highest gain: "<< highest_gain << std::endl; 
			std::cout <<" highest gain path " ;
			for(size_t i = 0 ; i < highest_path.size(); i++){
				std::cout << highest_path.at(i)<< " ";
			}
			std::cout << "" << std::endl;
			updating_centers(highest_path, center_numbers);//update all the strings, their number of happpening and gains! Notice that new trees are made in this function!
			highest_gain=0;
			std::map<std::vector<int>, size_t>::iterator hi = merged_centers.find(highest_path);
			assert(hi != merged_centers.end());
			std::vector<int> seen;
			std::cout << "hi second "<<hi->second << std::endl;
			seen.push_back(hi->second);//Because we dont want to use the same center again as the one with the highest gain //XXX Shall i look for all the previous ones?!!
			find_highest_gain(highest_path, highest_gain,seen);
			if(gains.size() != 0){
			/*	std::cout << "check highest path: " <<std::endl;
				for(size_t j = 0; j < highest_path.size(); j++){
					std::cout << highest_path.at(j)<< " ";
				}*/
			//	std::cout << " " <<std::endl;
			}
			center_numbers = center_numbers +1;
			merged_centers.insert(make_pair(highest_path, center_numbers));
		}
		std::cout << "final result " << std::endl;
//		for(map<std::vector<size_t>, int>::iterator it = gains.begin(); it != gains.end(); it++){
//			for(size_t j =0; j < it->first.size(); j++){
//				std::cout<< it->first.at(j)<< " ";
//			}
//			std::cout<< " gain is " << it->second << std::endl;
//		}
//		std::cout << "new centers: "<<std::endl;
//		for(std::map<std::vector<size_t>, size_t>::iterator it = merged_centers.begin(); it != merged_centers.end(); it++){
//			if(it != merged_centers.end()){
//				for(size_t i = 0; i < it->first.size(); i ++){
//					std::cout << it->first.at(i)<< " ";
//				}
//				std::cout << " its index is " << it->second <<std::endl;
//			}else {std::cout << "there is no merged center! " <<std::endl;}
//			
//		}		
	}
	void merging_centers::find_highest_gain(std::vector<int> & highest_path, int & highest_gain, std::vector<int> & seen){
		for(std::map<std::vector<int>, int>::iterator it = gains.begin(); it != gains.end(); it++){
			if(it->second > highest_gain && it->first != seen && it->first.size() != 1){
				highest_path = it->first;
				highest_gain = it->second;
			}else continue;
		}
	}
	void merging_centers::adding_new_centers(std::vector<std::vector<std::string> > & long_centers, std::vector<std::map<size_t, std::vector<std::string> > > & centersPositionOnASeq , std::vector<std::map<size_t, std::string> > & centerOnSequence){//Filling in the long centers vector and centersPositionOnASeq which contain all the long centers 
		//First the initial tree is built:
		std::vector<int> h_gain;//it is only used for the first time of making tree. for the next times the center with the highest gain is used.
		int index = 0; // The index of the center with the highest gain is used frome the next rounds.
		size_t num_seq = data.numSequences();
		for(size_t seq_id = 0; seq_id < num_seq; seq_id++){
			all_current_centers.at(seq_id) = centers.get_centers(seq_id);
		}
		suffixTree tree(num_seq, centers, all_current_centers);
		tree.make_tree(h_gain,index);
		merg_gain_value(tree);
		size_t biggest_index = 0;
		std::vector<std::string> sequence_of_centers;
	//	std::cout << "merged centers size: "<< merged_centers.size()<<std::endl;
		for(std::map<std::vector<int>,size_t>::iterator it = merged_centers.begin(); it != merged_centers.end(); it++){//The new indices are converted to their name
		/*	std::cout << " merged_center: " << it->second << std::endl;
			for(size_t j =0;j< it->first.size();j++){ std::cout<< it->first.at(j)<< " ";}
			std::cout<< " "<<std::endl;*/
		//	std::vector<size_t> updated_center = it->first;
			std::vector<int> list;
			list = it->first;
			int id = list.at(0);	
			bool ThereIsStillABigID = true;
			while (ThereIsStillABigID == true){
			//	std::cout << "here!" << std::endl;
			//	std::cout << "number of original centers: "<< centers.get_number_of_centers()<<std::endl;
				for(size_t j =0; j < list.size(); j++){
					id = list.at(j);
				//	std::cout << "id " << id <<std::endl;
					int id1 = id;
					if(id < 0){
						id1 = (-1)*id;
					}
					if(id1 > centers.get_number_of_centers()){
						for(std::map<std::vector<int>,size_t>::iterator it1 = merged_centers.begin(); it1 != merged_centers.end(); it1++){
							if(it1->second == id1){
								std::vector<int> temp;
							//	std::cout << "j "<< j << std::endl;
								if(j != 0){
									for(size_t i = 0; i < j;i++){
										temp.push_back(list.at(i));	
									}
								//	std::cout<< "temp0: "<<std::endl;
								//	for(size_t i = 0; i < temp.size(); i ++){
								//		std::cout << temp.at(i)<< " ";
								//	}
								//	std::cout << " " <<std::endl;
								}
								for(size_t i = 0; i < it1->first.size();i++){
									temp.push_back(it1->first.at(i));
								}
							//	std::cout<< "temp1: "<<std::endl;
							//	for(size_t i = 0; i < temp.size(); i ++){
							//		std::cout << temp.at(i)<< " ";
							//	}
							//	std::cout << " " <<std::endl;
								for(size_t i = j+1; i < list.size();i++){
									temp.push_back(list.at(i));
								}
							//	std::cout<< "temp: "<<std::endl;
							//	for(size_t i = 0; i < temp.size(); i ++){
							//		std::cout << temp.at(i)<< " ";
							//	}
								list = temp;
								j = j + it1->first.size()-1;
							//	std::cout << "j1 "<< j <<std::endl;
								break;
							} else continue;
						}
					}
				}
				ThereIsStillABigID = false;
				for(size_t j =0; j < list.size(); j++){
					id = list.at(j);//TODO Check how you got negetive id!!
				//	std::cout << "j2 " << j <<std::endl;
				//	std::cout << "id "<< id <<" num of centers "<< centers.get_number_of_centers() <<std::endl;
					int id1 = id;
					if(id < 0){
						id1 = (-1)*id;
					}
					if(id1 > centers.get_number_of_centers()){
						ThereIsStillABigID = true;
						std::cout << "Still a new center! "<<std::endl;
						break;
					}else continue;
				}
			}
		//	std::cout << "list size is "<< list.size() <<std::endl;
			for(size_t j =0; j < list.size(); j++){
				int id = list.at(j);
				std::string center = centers.find_center_name(id);
				sequence_of_centers.push_back(center);
			//	std::cout << center << std::endl;
			}
			long_centers.push_back(sequence_of_centers);
			for(size_t i = 0; i < data.numSequences(); i++){
				find_new_centers(it->second,sequence_of_centers,i, centersPositionOnASeq);
			}
			sequence_of_centers.clear();
		}
	}
	void merging_centers::index_centers(std::map<std::string , std::vector<pw_alignment> > & al_of_a_ccs){
		int id = 1;
		center_index.clear();
		for(std::map<std::string , std::vector<pw_alignment> >::iterator it = al_of_a_ccs.begin();it != al_of_a_ccs.end();it++){
			std::string center = it->first;
			std::vector<std::string> cent_parts;
			strsep(center, ":" , cent_parts);
			unsigned int cent_dir = atoi(cent_parts.at(0).c_str());
			unsigned int cent_ref = atoi(cent_parts.at(1).c_str());
			unsigned int cent_left = atoi(cent_parts.at(2).c_str());
			std::stringstream rev_cent;
			if(cent_dir==0){
				rev_cent <<1<<":" << cent_ref <<":"<< cent_left;
			}else{
				rev_cent <<0<<":" << cent_ref <<":"<< cent_left;
			}
			std::string reverse = rev_cent.str();
			std::map<std::string,int>::iterator it1 = center_index.find(reverse);
			if(it1 == center_index.end()){
				std::map<std::string, int>::iterator it2 = center_index.find(center);
				assert(it2 == center_index.end());
				if(cent_dir == 0){
					center_index.insert(make_pair(center,id));
				}else{
					center_index.insert(make_pair(center, (-1*(id))));
				}
				id ++;
			}else{
				center_index.insert(make_pair(center,-1*(it1->second)));
			}			
		}
	/*	std::cout << "indices"<<std::endl;
		for(std::map<std::string, int>::iterator it = center_index.begin(); it!=center_index.end();it++){
			std::cout << it->first << " " << it->second <<std::endl;
		}	*/	
	}
	void merging_centers::add_long_centers_to_map(std::vector<std::vector<std::string> > & long_centers, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers){
		for(size_t i =0; i < long_centers.size(); i ++){
			vector<std::string> long_center = long_centers.at(i);
			std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(long_center);
			if(new_cent == new_centers.end()){
				new_centers.insert(make_pair(long_center,std::vector<pw_alignment>()));
			}
		}
	}

	void merging_centers::create_alignment(std::vector<std::vector<std::string> > & long_centers, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string , std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::vector<std::string> > > & centersPositionOnASeq, std::vector<std::map<size_t , std::string> > & centerOnSequence,size_t & gap_in_long_centers){
		std::cout << "create long als "<< std::endl;
//		index_centers(alignments_in_a_cluster);
		center_index = centers.get_index();
		size_t artificial_ref = data.numSequences()+new_centers.size();
		add_long_centers_to_map(long_centers,new_centers);
		for(size_t i =0; i < long_centers.size();i++){//Includes local long centers
			std::cout << "current long center " << std::endl;
			for(size_t k =0; k < long_centers.at(i).size();k++){
				std::cout << long_centers.at(i).at(k)<< " ";
			}
			std::cout << " " <<std::endl;
			for(size_t j = 0 ; j < data.numSequences(); j++){
				std::cout << "num seq "<< j << " number of long centers on it " << centersPositionOnASeq.at(j).size() <<std::endl;
				for(std::map<size_t , std::vector<std::string> >::iterator it = centersPositionOnASeq.at(j).begin(); it != centersPositionOnASeq.at(j).end(); it++){//It includes all the long centers on a sequence //TODO make a better container that one doesn't need to loop over!
					size_t reverse = 0;
					if(j == 3){
						std::cout << "this pos "<< it->first <<std::endl;
					}
					if(it->second == long_centers.at(i)){//If the long center is on that sequence
						//first ref is considered on the forward strand of a sequence 
						for(size_t k =0; k < it->second.size(); k++){
							std::cout<< it->second.at(k)<<std::endl;
						}
						std::vector<bool> sample1;
						std::vector<bool> sample2;
						pw_alignment al;
						pw_alignment self_al;
						al.setreference1(j);
						al.setbegin1(it->first);
						size_t left_one = it->first; 
					//	size_t right_of_last_piece=0;
					//	size_t length = 0;
						size_t position = it->first;
						std::cout << "position " << position << std::endl;
					//	find_long_center_length(it->second, alignments_in_a_cluster,j,position,length,right_of_last_piece, gap_in_long_centers);//It returns length of the long center and its right position on the sequence //TODO
					//	std::cout << "end of last piece "<< right_of_last_piece << " length "<< length << std::endl;
					//	al.setend1(right_of_last_piece); 
					//	size_t right_one = right_of_last_piece; 
						//Up to now first reference of an pw-alignment is set. Center is always considered as the second reference.
						std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(it->second);
						assert(new_cent != new_centers.end());
						if(new_cent->second.size() == 0){//when we are at the first sequence that this long center is located	
							al.setreference2(artificial_ref);
							artificial_ref = artificial_ref+1;
							al.setbegin2(0);//Long center is considered as an artificial sequence
						//	al.setend2(length-1);
						}else{
							pw_alignment p1 = new_cent->second.at(0);
							al.setreference2(p1.getreference2());
							assert(p1.getreference2() >=data.numSequences());
						//	al.setbegin2(p1.getbegin2());//Long center is considered as an artificial sequence
						//	al.setend2(p1.getend2());
							al.setbegin2(0);//Long center is considered as an artificial sequence
						//	al.setend2(length-1); 

						}
						std::cout<< "begin2 "<< al.getbegin2() << " end2 " << al.getend2() << std::endl;
						//Second reference is already set.
						//We are setting samples here:
						size_t right = 0; //end of the previous center on a long center
						std::cout<< "it second size " << it->second.size() <<std::endl;
					//	set_smaples(it->second,j);
						size_t long_center_length = 0;
						size_t left_of_a_sample = it->first;
						for(size_t i =0; i < it->second.size();i++){
							std::cout<< " center "<< it->second.at(i)<<std::endl;
							std::string center = it->second.at(i);
							std::vector<std::string> cent_parts;
							strsep(center, ":" , cent_parts);
							unsigned int cent_dir = atoi(cent_parts.at(0).c_str());
							unsigned int cent_ref = atoi(cent_parts.at(1).c_str());
							unsigned int cent_left = atoi(cent_parts.at(2).c_str());
							std::stringstream cent;
							cent<< cent_ref<<":"<<cent_left;
							std::map<std::string , std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(cent.str());
							assert(it3 != alignments_in_a_cluster.end());
							std::cout << "left of this center is "<< left_of_a_sample <<std::endl;
							std::multimap<size_t , std::string>::iterator this_left = centerOnSequence.at(j).find(left_of_a_sample);
							assert(this_left != centerOnSequence.at(j).end());
							size_t nextleft = 0;
							if(i != it->second.size()-1){
							//	std::multimap<size_t , std::string>::iterator next_left = this_left++;
								std::advance(this_left , 1);
								assert(this_left != centerOnSequence.at(j).end());
								nextleft = this_left->first;
								std::cout << "next left is "<< nextleft <<std::endl;
							}
							assert(nextleft != 0 || i == it->second.size()-1);
							std::cout << "number of als " << it3->second.size()<<std::endl;
							//Setting samples:
							size_t cent_length = 0;
							set_samples(it3->second,reverse, sample1,sample2,left_of_a_sample,j,cent_ref,cent_left,right, cent_length);
							std::cout << "right is " << right << " cent_length " <<cent_length<<std::endl;
							long_center_length +=cent_length;
							if(i == it->second.size() -1){
								al.setend1(right);
								al.setend2(long_center_length-1); 
							}
							if(i != it->second.size()-1 && right != nextleft -1){
								assert(right < nextleft-1);
								vector<bool> middle_part_of_sample;
								vector<bool> gap_sample;
								std::cout << "right+1 "<< right+1 << " nextleft "<< nextleft <<std::endl;
								for(size_t m = right+1 ; m < nextleft ; m++){
									char base = data.getSequence(j).at(m);
									char gap = '-';
									vector<bool> bits;
									vector<bool> gap_bit;
									pw_alignment::get_bits(base,bits);
									pw_alignment::get_bits(gap,gap_bit);
									for(size_t n = 0; n < 3; n++){
										middle_part_of_sample.push_back(bits.at(n));
										gap_sample.push_back(gap_bit.at(n));
									}
								}
								std::cout<< "size of middle part is " << middle_part_of_sample.size() << std::endl;
								for(size_t m = 0 ; m < middle_part_of_sample.size();m++){
									sample1.push_back(middle_part_of_sample.at(m));
									sample2.push_back(gap_sample.at(m));
								}

							}
							left_of_a_sample = nextleft;
						}
						std::cout << "sam1 "<< sample1.size() << " sam2 " << sample2.size() << std::endl;
						al.set_alignment_bits(sample1,sample2);
						std::cout << "alignment is: " <<std::endl;
						al.print();
						std::cout << "reverse is " << reverse << " and long center size is " << it->second.size()<<std::endl;
						if(reverse==it->second.size()){//XXX it is not the most thoughtful way, the better way is making it in a right direction from the scratch, it is just to see if it works this way
						//	pw_alignment p;
						//	al.get_reverse(p);
							std::cout<< "reverse is added "<<std::endl;
							al.print();
						//	assert((p.getbegin2()== length -1) && (p.getend2()==0));
						//	p.setbegin2(0);
						//	p.setend2(length-1);
						//	new_cent->second.push_back(p);
						//	al.setbegin2(length-1);
							al.setbegin2(long_center_length -1);
							al.setend2(0);
							new_cent->second.push_back(al);	
						}else if(reverse > 0 && reverse != it->second.size()){
								std::cout<< "centers are happened in different direction so the al is not added."<<std::endl;
						}else{
							assert(reverse == 0);
							new_cent->second.push_back(al);
							std::cout << "al is added! "<<std::endl;
							al.print();
						//	if(al.getreference2() == 9){//XXX
						//		std::cout << al.get_al_ref1() << std::endl;
						//		std::cout << al.get_al_ref2() << std::endl;
						//	}
						}
						if(it->second.size()==3 && it->second.at(0)=="0:1:1941043"){
							std::cout << al.get_al_ref1() << std::endl;
							std::cout << al.get_al_ref2() << std::endl;
						}
					}// no break is needed at the end of this if loop becasue a long center can occur more than one time on a sequence
				}
			}
		}
		remove_fully_reverse_refs(long_centers, new_centers);//Edits new_centers
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it=merging_als.begin(); it!=merging_als.end();it++){
			std::cout << "center is "<< it->first<<std::endl;
			std::cout<<"members are "<<std::endl;
			for(size_t i =0; i < it->second.size();i++){
				it->second.at(i).print();
			}
		}
	}
	void merging_centers::remove_fully_reverse_refs(vector<vector<std::string> > & long_centers, std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers){
	//Translate centers to their indices and check for negetive reverse of each center. If it exists one of them is chosen and the other one is considered as the reverse of that reference.	
		std::cout << "remove reverse centers "<< std::endl;
		std::map<std::vector<int> , std::vector<std::string> >long_center_index;
		for(size_t j =0;j< long_centers.size();j++){
			std::vector<std::string> center = long_centers.at(j);
			std::vector<int> indices;
			for(size_t i =0; i < center.size();i++){
				std::map<std::string,int>::iterator it1 = center_index.find(center.at(i));
			//	std::cout << "cent " << it1->first << " index " << it1->second << std::endl;
				assert(it1!=center_index.end());
				indices.push_back(it1->second);
			}
			long_center_index.insert(std::make_pair(indices, center));
		}
		for(std::map<std::vector<int> , std::vector<std::string> >::iterator it = long_center_index.begin(); it != long_center_index.end();it++){
			for(size_t i =0 ; i < it->first.size();i++){
				std::cout << it->first.at(i) << " ";
			}
			std::cout << " " << std::endl;

		}
		for(std::map<std::vector<int> , std::vector<std::string> >::iterator it = long_center_index.begin(); it != long_center_index.end();it++){
			//Pick a vector make its negetive reverse and find it in "long_center_index" if it exists check for the one who has more als since it helps us use less reverse_flag in encoding for the others use the same virtual ref in the opposite direction
			std::vector<int> rev_indices;
			for(size_t i = it->first.size() ; i >0 ;i--){
				int rev = -1 * it->first.at(i-1);
				std::cout << rev << std::endl;
				rev_indices.push_back(rev);
			}
			for(size_t i =0; i < rev_indices.size();i++){
				std::cout << rev_indices.at(i) << " ";
			}
			std::cout << " " << std::endl;
			std::map<std::vector<int>, std::vector<std::string> >::iterator it1 = long_center_index.find(rev_indices);
			if(it1 != long_center_index.end()){//if the negative reverse exists
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it2 = new_centers.find(it1 ->second);
				assert(it2 != new_centers.end());
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it3 = new_centers.find(it->second);
				assert(it3 != new_centers.end());
				if(it2->second.size() > it3->second.size()){
					pw_alignment al = it2->second.at(0);
					for(size_t i =0; i < it3->second.size(); i++){
						pw_alignment p = it3->second.at(i);
						p.setreference2(al.getreference2());
						size_t begin = p.getbegin2();
						size_t end = p.getend2();
						p.setbegin2(end);
						p.setend2(begin);
						std::cout<< "negative reverse exists: "<<std::endl;
						p.print();
					}
				}else{
					pw_alignment al = it3->second.at(0);
					for(size_t i =0; i < it2->second.size(); i++){
						pw_alignment p = it2->second.at(i);
						p.setreference2(al.getreference2());
						size_t begin = p.getbegin2();
						size_t end = p.getend2();
						p.setbegin2(end);
						p.setend2(begin);
						std::cout<< "negative reverse exists: "<<std::endl;
						p.print();
					}
				}
			}
		}
	}
//The following function finds long centers on dna sequences
	void merging_centers::find_new_centers(size_t & center_indices, std::vector<std::string > & current_long_center, size_t & seq_id , std::vector<std::map<size_t, std::vector<std::string> > > & centersPositionOnASeq){
	//	std::map<size_t,std::vector<size_t> > all_centers = tree.get_center_on_a_sequence(seq_id);
		std::vector<size_t> positions = centers.get_long_center_position(seq_id , current_long_center);// All the positions a long center occur on a sequence, even if was part of a longer suffix it will be found
	//	std::vector<std::vector<std::vector<size_t> > > all_suffix = tree.get_suffix_of_sequence(seq_id);
	//	std::cout << "all the suffixes on "<< seq_id << std::endl;
	//	for(size_t i =0;i <all_suffix.size();i++){
	//		for(size_t j = 0; j < all_suffix.at(i).size();j++){
	//			for(size_t k =0; k < all_suffix.at(i).at(j).size();k++){
	//				std::cout <<all_suffix.at(i).at(j).at(k)<< " ";
	//			}
	//			std::cout << " "<<std::endl;
	//		}
	//		std::cout << " "<<std::endl;
	//	}
		std::cout<<"sequence is "<< seq_id << std::endl;
		std::cout << "current index "<< center_indices<< std::endl;
		std::cout << "long center is "<<std::endl;
		for(size_t i =0; i < current_long_center.size();i++){
			std::cout << current_long_center.at(i) << " ";
		}
		std::cout << " " <<std::endl;
		for(size_t  i =0; i < positions.size(); i++){
			std::cout << "at position "<< positions.at(i) << std::endl;
			centersPositionOnASeq.at(seq_id).insert(make_pair(positions.at(i),current_long_center));			
		}
	//	for(std::map<size_t,std::vector<size_t> >::iterator it1 = all_centers.begin(); it1!= all_centers.end();it1++){
	//			std::cout << " " <<std::endl;
	//			std::cout << "at position "<< it1->first << std::endl;
	//			for(size_t i =0; i < it1->second.size(); i++){
	//				std::cout << it1->second.at(i)<< " ";
	//			}
	//			std::cout << " " <<std::endl;
	//			if(it1->second.at(0) == center_indices){//XXX Why at(0)?
	//				centersPositionOnASeq.at(seq_id).insert(make_pair(it1->first,current_long_center));
	//				std::cout << "seq "<<seq_id << "position "<<it1->first << std::endl;
	//			}
	//			for(size_t i =0; i < it1->second.size();i++){
	//				if(it1->second.at(i)==center_indices){
//
//					}
//				}
//		}
	}
	void merging_centers::find_long_center_length(std::vector<std::string> & centers, std::map<std::string,std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & seq_id,size_t & position ,size_t & center_length, size_t & end_of_last_piece, size_t & gap_in_long_centers){
//		end_of_last_piece = 0;
		size_t current_position  = position;
		center_length = 0;
		std::map<size_t , size_t > first_left_and_right;
		for(size_t i =0; i < centers.size();i++){
			std::string center = centers.at(i);
			std::vector<std::string> center_parts;
			strsep(center, ":" , center_parts);
			unsigned int center_dir = atoi(center_parts.at(0).c_str());
			unsigned int center_ref = atoi(center_parts.at(1).c_str());
			unsigned int center_left = atoi(center_parts.at(2).c_str());
			std::stringstream cent;
			cent<<center_ref<<":"<<center_left;
			std::map<std::string , std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.find(cent.str());//removed the dir
			std::cout << "center is "<< cent.str() <<std::endl;
			assert(it != alignments_in_a_cluster.end());
		//	std::cout<< "al size "<<it->second.size() << std::endl;
		//	bool selfAligned = true;
			for(size_t j =0; j < it->second.size();j++){//TODO there could be more than one al fits the condition the closest one should be chosen
				pw_alignment p = it->second.at(j);
				size_t r1,r2,l1,l2;
				p.get_lr1(l1,r1);
				p.get_lr2(l2,r2);
				size_t ref1 = p.getreference1();
				size_t ref2 = p.getreference2();
		//		size_t centerPosition = data.getSequence(seq_id).length();
				std::cout << "l1 " << l1 << " r1 "<< r1 << " l2 "<< l2 << " r2 " << r2 << " ref1 "<<ref1 << " ref2 "<< ref2<<std::endl;
				std::cout << "seq id " << seq_id << " cent ref "<< center_ref << " current pos " << current_position << " cent lef "<< center_left << " allowed gap is " << gap_in_long_centers<< std::endl;
				//If center_ref is not the seq_id 
		/*		for ( std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(seq_id).begin() ; it1 != centerOnSequence.at(seq_id).end();it1++){
					if(it1->second == centers.at(i)){//It is wrong cus a single center may happen more than once on a seq
						centerPosition = it1->first;
						std::cout << "center position " << centerPosition << std::endl;
						break;
					}
				}
				assert(centerPosition != data.getSequence(seq_id).length());*/
				//If we need to align the same piece against itself
				if(ref1 == seq_id && seq_id == center_ref && center_left == l1 && center_left >= current_position && center_left <= current_position + gap_in_long_centers){
					center_length +=r1-l1+1;
					current_position = r1;
					std::cout << "center length here! "<< center_length  <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r1;
						std::cout << "end point2 "<< end_of_last_piece << " r1 "<< r1 <<std::endl;
					}
					break;
				}
				else if(ref2 == seq_id && seq_id == center_ref && center_left == l2 && center_left >= current_position && center_left <= current_position + gap_in_long_centers){
					center_length +=r2-l2+1;
					current_position = r2;
					std::cout << "center length there! " << center_length <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r2;
						std::cout << "end point3 "<< end_of_last_piece << " r2 "<< r2 <<std::endl;

					}
					break;

				}

				else if((ref1== seq_id && ref2== center_ref && l1 >= current_position && l1 <= current_position + gap_in_long_centers && l2 ==center_left)){//when center is on ref2
					std::cout<< "center on ref 2"<<std::endl;
					center_length += r2-l2+1;
					current_position = r1;
					std::cout << "cur pos " << current_position << " length "<< center_length << std::endl;
					std::map<std::string,std::vector<pw_alignment> >::iterator find_al = merging_als.find(cent.str());
					if(find_al==merging_als.end()){
						merging_als.insert(std::make_pair(cent.str(),std::vector<pw_alignment>()));
						find_al=merging_als.find(cent.str());
					}
					find_al->second.push_back(p);
					if(i == centers.size()-1){
						end_of_last_piece = r1;
						std::cout << "end point0 "<< end_of_last_piece << " r1 "<< r1 <<std::endl;
					}
					break;
				}
				else if(ref2== seq_id && ref1== center_ref && l2 >= current_position && l2<= current_position + gap_in_long_centers && l1 ==center_left){//when center is on ref1
					std::cout<< "center on ref 1"<<std::endl;
					center_length +=r1-l1+1;
					current_position = r2;
					std::cout << "cur pos " << current_position << " length "<< center_length << std::endl;
					std::map<std::string,std::vector<pw_alignment> >::iterator find_al = merging_als.find(cent.str());
					if(find_al==merging_als.end()){
						merging_als.insert(std::make_pair(cent.str(),std::vector<pw_alignment>()));
						find_al=merging_als.find(cent.str());
					}
					find_al->second.push_back(p);
					if(i == centers.size()-1){
						end_of_last_piece = r2;
						std::cout << "end point1 "<< end_of_last_piece << " r2 "<< r2 <<std::endl;

					}
					break;
				}
			

			}
		}
	}
	void merging_centers::set_samples(std::vector<pw_alignment>& als,size_t & reverse, std::vector<bool> & sample1, std::vector<bool> & sample2, size_t & left_of_a_sample, size_t & seq_id, unsigned int & cent_ref, unsigned int & cent_left, size_t & right, size_t & cent_length){
		for(size_t i = 0; i < als.size();i++){
			pw_alignment p = als.at(i);
			size_t r1,r2,l1,l2;
			p.get_lr1(l1,r1);
			p.get_lr2(l2,r2);
			std::vector<bool> sample1_p = p.getsample1();
			std::vector<bool> sample2_p = p.getsample2();
			std::cout << " sample1_p.size() " << sample1_p.size()<< " " << sample2_p.size()<< " " << p.alignment_length() <<std::endl;
			std::cout<< "ref 1 "<< p.getreference1() <<" ref2 "<< p.getreference2()<<" seq_id " << seq_id << " l1 " << l1 << " l2 " << l2 << " left_of_a_sample "<< left_of_a_sample << std::endl;
			std::vector< std::vector<bool> >revSample;
			p.get_reverse_complement_sample(revSample);
			std::cout << "r1 " << r1 << " r2 " << r2 <<std::endl;
			if(p.getreference1()== seq_id && l1 == left_of_a_sample && p.getreference2()== cent_ref && l2 == cent_left){
				right = r1;
				std::cout << "center on the second ref " << right<<std::endl;
				if(right == 893735 || right == 896025){
					std::cout << p.get_al_ref1() <<std::endl;
					std::cout << p.get_al_ref2() <<std::endl;
				}
				cent_length = r2-l2+1;
				if(p.getbegin1() < p.getend1()){
					for(size_t m =0; m < sample1_p.size(); m++){
						sample1.push_back(sample1_p.at(m));
						sample2.push_back(sample2_p.at(m));
					}
				}else{//I should always think about turning the corresponding center!
					for(size_t m = 0; m < revSample.at(0).size();m++){
						sample1.push_back(revSample.at(0).at(m));
						sample2.push_back(revSample.at(1).at(m));
					}
					reverse++;
				}
				std::stringstream cent;//XXX JUST ADDED
				cent<<cent_ref<<":"<<cent_left;
				std::map<std::string,std::vector<pw_alignment> >::iterator find_al = merging_als.find(cent.str());
				if(find_al==merging_als.end()){
					merging_als.insert(std::make_pair(cent.str(),std::vector<pw_alignment>()));
					find_al=merging_als.find(cent.str());
				}
				find_al->second.push_back(p);
				std::cout << "add to merg" << std::endl;
				p.print();

				break;
			}else if(p.getreference2() == seq_id && l2 == left_of_a_sample && p.getreference1() == cent_ref && l1 == cent_left){
				std::cout << "center on the first ref" << std::endl;
				right = r2;
				cent_length = r1-l1+1;
				if(p.getbegin2()<p.getend2()){
					std::cout << "forward "<<std::endl;
					for(size_t m =0; m < sample2_p.size();m++){
						sample1.push_back(sample2_p.at(m));
						sample2.push_back(sample1_p.at(m));
					}
				}else{
					std::cout << "reverse "<<std::endl;
					for(size_t m = 0; m < revSample.at(1).size();m++){
						sample1.push_back(revSample.at(1).at(m));
						sample2.push_back(revSample.at(0).at(m));
					}
					reverse ++;
				}
				std::stringstream cent;//XXX JUST ADDED
				cent<<cent_ref<<":"<<cent_left;
				std::map<std::string,std::vector<pw_alignment> >::iterator find_al = merging_als.find(cent.str());
				if(find_al==merging_als.end()){
					merging_als.insert(std::make_pair(cent.str(),std::vector<pw_alignment>()));
					find_al=merging_als.find(cent.str());
				}
				find_al->second.push_back(p);
				std::cout << "add to merg" << std::endl;
				p.print();
				break;
			}else if(p.getreference2() == seq_id && cent_ref == seq_id && cent_left == l2 && cent_left == left_of_a_sample){//Instead of using samples sequence bases should be used! Samples are including gaps! 
				std::cout << "center on the ref - second sample "<<std::endl;
				right = r2;
				cent_length = r2-l2 + 1;
				std::vector<bool> intermediate;
				for(size_t m =l2; m <= r2; m++){
					char base = data.getSequence(seq_id).at(m);
					std::vector<bool> bits;
					pw_alignment::get_bits(base,bits);
					for(size_t n = 0; n < 3; n++){
						intermediate.push_back(bits.at(n));
					}
				}
				std::cout << intermediate.size() << std::endl;
				for(size_t m =0 ; m < intermediate.size();m++){
					sample1.push_back(intermediate.at(m));
					sample2.push_back(intermediate.at(m));
				}
			/*	std::stringstream cent;//XXX JUST ADDED
				cent<<cent_ref<<":"<<cent_left;
				std::map<std::string,std::vector<pw_alignment> >::iterator find_al = merging_als.find(cent.str());
				if(find_al==merging_als.end()){
					merging_als.insert(std::make_pair(cent.str(),std::vector<pw_alignment>()));
					find_al=merging_als.find(cent.str());
				}
				find_al->second.push_back(p);
				std::cout << "add to merg" << std::endl;
				p.print();*/
				break;
			}else if(p.getreference1()== seq_id && cent_ref == seq_id && cent_left == l1 && cent_left == left_of_a_sample ){
				std::cout << "center on the ref - first sample "<<std::endl;
				right = r1;
				cent_length = r1 - l1 + 1;
				std::vector<bool> intermediate;
				for(size_t m =l1; m <= r1; m++){
					char base = data.getSequence(seq_id).at(m);
					vector<bool> bits;
					pw_alignment::get_bits(base,bits);
					for(size_t n = 0; n < 3; n++){
						intermediate.push_back(bits.at(n));
					}
				}
				for(size_t m =0 ; m < intermediate.size();m++){
					sample1.push_back(intermediate.at(m));
					sample2.push_back(intermediate.at(m));
				}

			/*	std::stringstream cent;//XXX JUST ADDED
				cent<<cent_ref<<":"<<cent_left;
				std::map<std::string,std::vector<pw_alignment> >::iterator find_al = merging_als.find(cent.str());
				if(find_al==merging_als.end()){
					merging_als.insert(std::make_pair(cent.str(),std::vector<pw_alignment>()));
					find_al=merging_als.find(cent.str());
				}
				find_al->second.push_back(p);
				std::cout << "add to merg" << std::endl;
				p.print();*/
				break;
			}
		//	std::cout << sample1.size() << " " << sample2.size() << std::endl;
		}
	}
	//Short centers are added to the long centers container.
/*	void mixing_centers::add_original_centers_to_long_centers(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster,std::vector<std::multimap<size_t, std::string > > & centerOnSequence){//TODO !!
	// XXX Note that second both ref can be backwards
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
			std::stringstream str;
			str << center_ref<<":"<<center_left;
			std::map<std::string, std::vector<pw_alignment> >::iterator cent = alignments_in_a_cluster.find(str.str());
			assert(cent != alignments_in_a_cluster.end());
		//	if(cent == alignments_in_a_cluster.end()){
		//		if(center_dir == 0){
		//			center_dir = 1;
		//		}else{
		//			center_dir = 0;
		//		}
		//		std::stringstream temp;
		//		temp << center_dir<< ":"<< center_ref<<":"<<center_left;
		//		center = temp.str();
		//		std::cout << "CENT "<<center <<std::endl;
		//		cent = alignments_in_a_cluster.find(center);
		//		assert(cent != alignments_in_a_cluster.end());
		//	}
			std::cout << "number of als " << cent->second.size() <<std::endl;
			for(size_t j =0;j < cent->second.size();j++){//Go through the short alignments and check to see if they are part a long one
				pw_alignment p = cent->second.at(j); //Check if it happened on a long alignment!
				std::cout << "p is "<<std::endl;
				size_t l1,l2,r1,r2;
				size_t ref1 = p.getreference1();
				size_t ref2 = p.getreference2();
				std::cout << "ref1 "<< ref1 << " ref2 "<<ref2 <<std::endl;
				p.get_lr1(l1, r1);
				p.get_lr2(l2, r2);
				if(l1 == center_left && ref1 == center_ref){//Compares its ref2 with the ref1 of the it->seconds member.
					for(size_t k =0; k < it->second.size(); k++){//Long als
						pw_alignment al = it->second.at(k);
						size_t left1,left2,right1,right2;
						size_t ref_1 = al.getreference1();
						size_t ref_2 = al.getreference2();
						al.get_lr1(left1, right1);
						al.get_lr2(left2, right2);
						size_t pos_on_ref = data.getSequence(ref2).length();
						size_t rev_center_dir;
						std::stringstream rev_center;
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
						std::stringstream rev_center;
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
	}
	//Add rest of the centers to the new_center
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
			std::stringstream rev_center;
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
	swap_references(new_centers);
}*/
	void mixing_centers::add_original_centers_to_long_centers(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster,std::map<std::string, std::vector<pw_alignment> > & merged_als, std::vector<std::map<size_t, std::vector<std::string> > > & long_centers_on_a_seq, std::vector<std::map<size_t, std::string> > & centerOnSequence){
		all_centers_on_a_seq = long_centers_on_a_seq; //has both long and short centers
		std::cout<<"new centers before: "<< new_centers.size() <<std::endl;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it= new_centers.begin() ; it != new_centers.end(); it++){
			std::cout << "cent: "<< it->first <<std::endl;
			for(size_t i = 0; i < it->second.size();i++){
				it->second.at(i).print();
			}
		}
		std::cout<< "original cluster in: "<< alignments_in_a_cluster.size() <<std::endl;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it= alignments_in_a_cluster.begin() ; it != alignments_in_a_cluster.end() ;it++){
			assert(it->second.size() != 0);
			std::cout << "cent: "<< it->first <<std::endl;
			for(size_t i = 0; i < it->second.size();i++){
				it->second.at(i).print();
			}
			std::stringstream str;
			str<<0<<":"<<it->first;//Forward direction is added to all the centers
			std::vector<std::string> temp;
			temp.push_back(str.str());
			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator new_cent=new_centers.find(temp);
			assert(new_cent==new_centers.end());
			new_centers.insert(std::make_pair(temp, std::vector<pw_alignment>()));
			new_cent=new_centers.find(temp);
			
			std::map<std::string,std::vector<pw_alignment> >::iterator find_al= merged_als.find(it->first);
			if(find_al==merged_als.end()){
				//Add it to new_centers
				std::cout << "not part of a long center " << it->second.size() <<std::endl;
			//	new_cent->second = it->second;
				//Add this centers to long centers that happen on a seq
				for(size_t i = 0; i < it->second.size();i++){
				/*	bool direction = is_forward(it->first, it->second.at(i));
					if(direction == true){//Check coordinate and see if forward or backward
						new_cent->second.push_back(it->second.at(i));
					}else{
						assert(direction == false);
						std::stringstream str1;
						str1<<1<<":"<<it->first;//Reverse direction is added to all the centers
						std::vector<std::string> temp1;
						temp1.push_back(str1.str());
						std::cout<< str1.str()<<std::endl;
						std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator rnew_cent=new_centers.find(temp1);
						if(rnew_cent==new_centers.end()){
							new_centers.insert(std::make_pair(temp1, std::vector<pw_alignment>()));
							rnew_cent=new_centers.find(temp1);
						}
						rnew_cent->second.push_back(it->second.at(i));
					}*/
					add_to_long_center_map(it->first,it->second.at(i), centerOnSequence , new_centers);
				}
			}else{
				//Add the rest of als but not the existing ones
				for(size_t k =0; k < it->second.size();k++){
					bool equal = false;
					for(size_t j = 0; j < find_al->second.size();j++){
						equal =it->second.at(k).equals(find_al->second.at(j));
						if(equal == true){
							break;
						}
					}
					if(equal == false){
					//	std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator new_cent=new_centers.find(temp);
					//	if(new_cent==new_centers.end()){
					//		new_centers.insert(std::make_pair(temp, std::vector<pw_alignment>()));
					//		new_cent=new_centers.find(temp);
					//	}
						std::cout << "original al is saved "<<std::endl;
						it->second.at(k).print();
					/*	bool direction = is_forward(it->first, it->second.at(k));
						if(direction == true){//Check coordinate and see if forward or backward
							new_cent->second.push_back(it->second.at(k));
						}else{
							assert(direction == false);
							std::stringstream str1;
							str1<<1<<":"<<it->first;//Reverse direction is added to all the centers
							std::vector<std::string> temp1;
							temp1.push_back(str1.str());
							std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator rnew_cent=new_centers.find(temp1);
							if(rnew_cent==new_centers.end()){
								new_centers.insert(std::make_pair(temp1, std::vector<pw_alignment>()));
								rnew_cent=new_centers.find(temp1);
							}
							rnew_cent->second.push_back(it->second.at(k));
						}*/
						add_to_long_center_map(it->first,it->second.at(k), centerOnSequence, new_centers);
					}
				}
				if(new_cent->second.size()==0){
					//Cent left and cent ref are added to the all_centers_on_a_seq.at(cent_ref)
					std::string center = it->first;
					std::vector<std::string> center_parts;
					strsep(center, ":" , center_parts);
					unsigned int center_ref = atoi(center_parts.at(0).c_str());
					unsigned int center_left = atoi(center_parts.at(1).c_str());
					std::stringstream str;
					str << 0<<":"<<center_ref<<":"<<center_left;
					std::vector<std::string> temp;
					temp.push_back(str.str());
					//Find it and assert it reaches the end
					std::map<size_t, std::vector<std::string> >::iterator it1 = all_centers_on_a_seq.at(center_ref).find(center_left);
					if(it1 == all_centers_on_a_seq.at(center_ref).end()){
						all_centers_on_a_seq.at(center_ref).insert(std::make_pair(center_left,temp));
						std::cout<< "center is always part of a long center but not on its own ref "<< center <<std::endl;

					}else{
						std::cout<< "center is always part of a long center "<< center <<std::endl;
					}
				}
			}
		}
		swap_references(new_centers);
		std::cout<<"new centers after: "<< new_centers.size()<<std::endl;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it= new_centers.begin() ; it != new_centers.end(); it++){
			std::cout << "cent: "<< it->first <<std::endl;
			for(size_t i = 0; i < it->second.size();i++){
				it->second.at(i).print();
			}
		}
		std::cout << "all the centers on seq 0: "<<std::endl;
		for(std::map<size_t, std::vector<std::string> >::iterator it = all_centers_on_a_seq.at(0).begin() ; it != all_centers_on_a_seq.at(0).end() ; it++){
			std::cout << "at " << it->first << " is "<<std::endl;
			std::cout << it->second << std::endl;
		}
	}
	bool mixing_centers::is_forward(const std::string & center , pw_alignment & p){
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_ref = atoi(center_parts.at(0).c_str());
		unsigned int center_left = atoi(center_parts.at(1).c_str());
		unsigned int ref1 = p.getreference1();
		unsigned int ref2 = p.getreference2();
		size_t left1, right1, left2, right2;
		p.get_lr1(left1, right1);
		p.get_lr2(left2, right2);
		std::cout << p.get_al_ref1()  << " " << p.get_al_ref2() << std::endl;
		if(ref1 == center_ref && left1 == center_left){	
			std::cout << p.get_al_ref1()  << " " << p.get_al_ref2() << " "<< data.extract_reverse_seq_part(ref1,left1,right1)<< " " << data.extract_seq_part(ref1,left1,right1) <<std::endl;
			std::cout << data.extract_reverse_seq_part(ref2,left2,right2)<< " " << data.extract_seq_part(ref2,left2,right2) <<std::endl;

		//	if(p.getbegin1() > p.getend1() && p.getbegin2() < p.getend2()){ //XXX OR (begin==end && thisbase == reverse of seq at this position)
			if((p.getbegin1()>p.getend1()) || (p.getbegin1() == p.getend1() && p.get_al_ref1() == data.extract_reverse_seq_part(ref1,left1,right1))){
				std::cout << "rev1"<<std::endl;
				return false;
			}
		}
		if(ref2 == center_ref && left2 == center_left){	
			std::cout << "rev2"<<std::endl;
		//	if(p.getbegin2() > p.getend2() && p.getbegin1() < p.getend1()){
			if((p.getbegin2()>p.getend2()) || (p.getbegin2() == p.getend2() && p.get_al_ref2() == data.extract_reverse_seq_part(ref2,left2,right2))){

				return false;
			}
		}
		return true;
	}
	void mixing_centers::swap_references(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers){
		//At the end all the alignemnts are swapped in the way that center be on the second ref. It is used for the graph maf & arith encoding
		std::cout<< "new cent in swap ref: "<< std::endl;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it!=new_centers.end();it++){
			std::cout<< it->first <<std::endl;
			if(it->first.size()==1 && it->second.size()!=0){
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
				//	unsigned int dir;
				//	if(it->second.at(i).getbegin1() < it->second.at(i).getend1()){
				//		dir = 0;
				//	}else{ 
				//		dir = 1;
				//	}
					// swapped alignment is written
					std::stringstream centerref_al;
					std::stringstream otherref_al;
					if(ref1 == center_ref && left1 == center_left){	
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
	}
	

	void mixing_centers::add_to_long_center_map(const std::string & center,const pw_alignment & al, std::vector<std::map<size_t, std::string> > & centerOnSequence, std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers){
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_ref = atoi(center_parts.at(0).c_str());
		unsigned int center_left = atoi(center_parts.at(1).c_str());

		unsigned int ref1 = al.getreference1();
		unsigned int ref2 = al.getreference2();
		size_t left1, right1, left2, right2;
		al.get_lr1(left1, right1);
		al.get_lr2(left2, right2);
		if(al.alignment_length()==1){ //XXX error point in ecoli16!
			std::cout << "center of length 1 "<< center_ref << " " << center_left << std::endl;
			std::cout << al.get_al_ref1()  << " " << al.get_al_ref2() << " "<< data.extract_reverse_seq_part(ref1,left1,right1)<< " " << data.extract_seq_part(ref1,left1,right1) <<std::endl;
			std::cout << data.extract_reverse_seq_part(ref2,left2,right2)<< " " << data.extract_seq_part(ref2,left2,right2) <<std::endl;
		}
		if(ref1 ==center_ref && left1 == center_left){
			std::map<size_t, std::vector<std::string> >::iterator it = all_centers_on_a_seq.at(ref1).find(left1);
			if(it==all_centers_on_a_seq.at(ref1).end()){
				std::stringstream temp;
				temp<<0<<":"<<center;
				std::vector<std::string> this_cent;
				this_cent.push_back(temp.str());
				all_centers_on_a_seq.at(ref1).insert(std::make_pair(left1,this_cent));	
			}
			std::map<size_t, std::vector<std::string> >::iterator it1 = all_centers_on_a_seq.at(ref2).find(left2);
			if(it1 == all_centers_on_a_seq.at(ref2).end()){
				std::cout<< "left "<< left2 << std::endl;
			}
			assert(it1 == all_centers_on_a_seq.at(ref2).end());
			std::vector<std::string> this_cent;
			std::multimap<size_t, std::string>::iterator short_cent=centerOnSequence.at(ref2).find(left2);
			assert(short_cent != centerOnSequence.at(ref2).end());
			std::cout<< "short cent: " <<short_cent->second <<std::endl;
			this_cent.push_back(short_cent->second);
			all_centers_on_a_seq.at(ref2).insert(std::make_pair(left2,this_cent));
			//ADD IT TO NEW_CENTERS:
		//	std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator long_cent = new_centers.find(this_cent);
			std::stringstream temp;
			temp<<0<<":"<<center;
			std::vector<std::string>Center;
			Center.push_back(temp.str());
			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator long_cent = new_centers.find(Center);

			if(long_cent == new_centers.end()){
				new_centers.insert(std::make_pair(Center,std::vector<pw_alignment>()));
				long_cent = new_centers.find(Center);
			}
			long_cent->second.push_back(al);
		}
		else{
			assert(ref2 == center_ref && left2 == center_left);
			std::map<size_t, std::vector<std::string> >::iterator it = all_centers_on_a_seq.at(ref2).find(left2);
			if(it == all_centers_on_a_seq.at(ref2).end()){
				std::stringstream temp;
				temp<<0<<":"<<center;
				std::vector<std::string> this_cent;
				this_cent.push_back(temp.str());
				all_centers_on_a_seq.at(ref2).insert(std::make_pair(left2,this_cent));	
			}
			std::map<size_t, std::vector<std::string> >::iterator it1 = all_centers_on_a_seq.at(ref1).find(left1);
			if(it1 != all_centers_on_a_seq.at(ref1).end()){
				std::cout << it1->first << " " << it1->second << std::endl;
			}
			assert(it1==all_centers_on_a_seq.at(ref1).end());
			std::vector<std::string> this_cent;
			std::multimap<size_t, std::string>::iterator short_cent=centerOnSequence.at(ref1).find(left1);
			assert(short_cent != centerOnSequence.at(ref1).end());
			std::cout<< "short cent: " <<short_cent->second <<std::endl;
			this_cent.push_back(short_cent->second);
			all_centers_on_a_seq.at(ref1).insert(std::make_pair(left1,this_cent));
			//ADD IT TO NEW_CENTERS:
		//	std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator long_cent = new_centers.find(this_cent);
			std::stringstream temp;
			temp<<0<<":"<<center;
			std::vector<std::string>Center;
			Center.push_back(temp.str());

			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator long_cent = new_centers.find(Center);
			if(long_cent == new_centers.end()){
				new_centers.insert(std::make_pair(Center,std::vector<pw_alignment>()));
				long_cent = new_centers.find(Center);
			}
			long_cent->second.push_back(al);

		}
	}
	void mixing_centers::get_direction(const pw_alignment & p, size_t & dir){
	/*	if(p.getbegin1()< p.getend1() && p.getbegin2() < p.getend2()){
			dir = 0;
		}
		if(p.getbegin1()< p.getend1() && p.getbegin2() > p.getend2()){
			dir = 1;
		}
		if(p.getbegin1()> p.getend1() && p.getbegin2() < p.getend2()){
			dir = 1;
		}
		if(p.getbegin1()> p.getend1() && p.getbegin2() > p.getend2()){
			dir = 0;
		}*/
	}
	//computing the weight of all the cluster centers
	void mixing_centers::calculate_centers_weight(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, unsigned int> & weight){//On original centers:
		size_t max_bit = 8;	
		size_t max_members = 0;//Returns the size of the largest cluster
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.begin(); it !=alignments_in_a_cluster.end(); it++){	
	        	if(it->second.size()+1> max_members){		
				max_members = it->second.size()+1;
			}
		
		}	
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.begin(); it !=alignments_in_a_cluster.end(); it++){
			std::string center = it ->first;
			std::map<std::string, unsigned int >::iterator it1 = weight.find(center);
			if(it1 == weight.end()){
				weight.insert(std::make_pair(center,0));
				it1 = weight.find(center);
			}
			size_t power_of_two = 1<<(max_bit);
			it1->second = (unsigned int)((it->second.size()+1)*(power_of_two-1)/(double)max_members);
			if(it1 ->second == 0){
				it1->second = 1;
			}
		}
	}
	void mixing_centers::calculate_long_centers_weight(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::vector<std::string> , unsigned int> & long_center_weight){//On concatenated ones:
		size_t max_bit = 8;	
		size_t max_members = 0;//At the end should be equal to the number of center occurs more than all the others
		std::cout << "here new center size is "<< new_centers.size()<<std::endl;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it != new_centers.end();it++){
			if(it->second.size()+1> max_members){		
				max_members = it->second.size()+1;
			}	
		}
		std::cout<< "maximum_mem "<< max_members << std::endl;	
		std::cout<<"number of a center"<<std::endl;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.begin(); it!=new_centers.end(); it++){
			std::cout<< it->first << std::endl;
			std::vector<std::string> this_center = it->first;
			std::map<std::vector<std::string>, unsigned int >::iterator it1 = long_center_weight.find(this_center);
			if(it1 == long_center_weight.end()){
				long_center_weight.insert(make_pair(this_center,0));

				it1 = long_center_weight.find(this_center);
			}
			it1->second = (unsigned int)((it->second.size()+1)*((1<<(max_bit))-1)/(double)max_members);
			if(it1 ->second == 0){
				it1->second = 1;
			}
		}
		std::cout << "size of weight "<<long_center_weight.size() <<std::endl;
	}

	void mixing_centers::find_members_of_clusters(std::map<std::string,std::string> & member_of_cluster, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster){
		size_t al_size = 0;
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
				std::stringstream sample1;
				std::stringstream sample2;
				sample1 <<ref1 << ":" << left_1;
				sample2 <<ref2 << ":" << left_2;
				if(sample1.str() != it->first){
					member_of_cluster.insert(make_pair(sample1.str(),it->first));
				}
				if(sample2.str() != it->first){
					member_of_cluster.insert(make_pair(sample2.str(), it->first));
				}
				std::cout << "s1 "<< sample1.str() << " s2 "<<sample2.str() << " center " << it->first << std::endl;
				assert(sample1.str()==it->first || sample2.str() == it->first);
			}
			std::string center = it->first;
			member_of_cluster.insert(make_pair(it->first, it->first));		
		}
		for(std::map<std::string, std::string>::iterator it = member_of_cluster.begin();it != member_of_cluster.end();it++){
			std::cout << it->second << std::endl;
		}
		std::cout<< "size of memeber_of_cluster is : "<< member_of_cluster.size()<<std::endl;
		std::cout << "al size "<< al_size << std::endl;
	}
	void  mixing_centers::find_adjacencies(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::multimap<size_t, std::string > > & centerOnSequence, std::map<size_t, std::set<size_t> > &adjacencies){
		std::vector<std::string> centers_with_members;
		std::vector<std::set<size_t> >related_references;
		std::map<std::string, size_t> center_id;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
			if(it->second.size()>0){
				std::set<size_t> related_ref;
				centers_with_members.push_back(it->first);
				center_id.insert(std::make_pair(it->first,related_references.size()));
				for(size_t i = 0; i < it->second.size(); i++){
					size_t direction;
					size_t l1,r1,l2,r2;
					pw_alignment al = it->second.at(i);
					size_t ref1 = al.getreference1();
					size_t ref2 = al.getreference2();
					al.get_lr1(l1,r1);
					al.get_lr2(l2,r2);
					std::string center = it->first;
					std::vector<std::string> center_parts;
					strsep(center, ":" , center_parts);
					unsigned int center_dir = atoi(center_parts.at(0).c_str());
					unsigned int center_ref = atoi(center_parts.at(1).c_str());
					unsigned int center_left = atoi(center_parts.at(2).c_str());

					related_ref.insert(ref1);
					related_ref.insert(ref2);
				}
				related_references.push_back(related_ref);
			}
		}
		assert(centers_with_members.size()== related_references.size());
		std::cout<< "index size is "<< centers_with_members.size()<<std::endl;
		for(size_t i = 0; i < related_references.size();i++){
			adjacencies.insert(std::make_pair(i, std::set<size_t>()));
			std::string center = centers_with_members.at(i);
			for(std::set<size_t>::iterator it = related_references.at(i).begin(); it!= related_references.at(i).end(); it++){
				std::string previous;
				size_t ref = *it;
				for(std::multimap<size_t, std::string >::iterator it1 = centerOnSequence.at(ref).begin(); it1 != centerOnSequence.at(ref).end();it1++){
					std::cout << "pos "<< it1->first <<std::endl;
					if(previous == center){
						std::cout << "adj pos " << it1->first<<std::endl;
						std::map<std::string, size_t>::iterator id = center_id.find(it1->second);
						if(id != center_id.end()){
							std::map<size_t, std::set<size_t> >::iterator adj = adjacencies.find(i);
							assert(adj != adjacencies.end());
							adj->second.insert(id->second);			
							break;
						}else{std::cout << "need ot look for the reverse!"<<std::endl;}
					}else{
						previous = it1->second;
					}
				}
			}
		}
		
	}
	void mixing_centers::compute_centers_length(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string , size_t> & center_length){
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.begin() ; it!= alignments_in_a_cluster.end(); it++){
			pw_alignment p= it->second.at(0);
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			std::stringstream sample1;
			std::stringstream sample2;
			sample1 <<ref1 << ":" << left_1;
			sample2 <<ref2 << ":" << left_2;
			if(sample1.str() != it->first){
				center_length.insert(std::make_pair(it->first,(right_2-left_2+1)));
			}
			if(sample2.str() != it->first){
				center_length.insert(std::make_pair(it->first,(right_1-left_1+1)));
			}
		}
	}
	const size_t mixing_centers::find_right_on_seq(std::string & center,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & left , size_t & ref, std::map<std::string, size_t> & non_aligned_right)const{
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		std::stringstream str;
		str<<center_ref<<":"<<center_left;

		std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.find(str.str());
		if(it!= alignments_in_a_cluster.end()){
			assert(it->second.size() > 0);
			std::map<std::string, size_t>::iterator it1 =non_aligned_right.find(center);
			assert(it1 == non_aligned_right.end());
			for(size_t i = 0 ; i < it->second.size();i++){
				pw_alignment p= it->second.at(i);
				size_t left_1; 
				size_t left_2;
				size_t right_1;
				size_t right_2;
				size_t ref1;
				size_t ref2;
				p.get_lr1(left_1,right_1);
				p.get_lr2(left_2,right_2);
				ref1 = p.getreference1();
				ref2 = p.getreference2();
				if((ref1 == 1 && ref2 == 4) || (ref1 == 4 && ref2 == 1)) p.print();
				std::stringstream sample1;
				std::stringstream sample2;
				sample1 <<ref1 << ":" << left_1;
				sample2 <<ref2 << ":" << left_2;
				if(sample1.str()==it->first && left_1 == left && ref1 == ref){
					return right_1;
				}
				if(sample2.str() == it->first && left_2 == left && ref2 == ref){
					return right_2;
				}
				if(sample2.str()==it->first && left_1 == left && ref1 == ref){
					return right_1;
				}
				if(sample1.str() == it->first && left_2 == left && ref2 == ref){
					return right_2;
				}
			}
		}else{
			std::cout << "is on a non aligned region! " <<std::endl;
			std::map<std::string , size_t>::const_iterator it1 = non_aligned_right.find(center);
			assert(it1 != non_aligned_right.end());
			size_t right = it1->second;
			return right;
		}
	}
const size_t mixing_centers::find_right_on_seq(std::string & center,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & left , size_t & ref, std::map<std::string, size_t> & non_aligned_right, pw_alignment & al, size_t & gap_on_long_centers)const{
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		std::stringstream str;
		str<<center_ref<<":"<<center_left;

		std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.find(str.str());
		assert(it!= alignments_in_a_cluster.end());
		assert(it->second.size() > 0);
		std::map<std::string, size_t>::iterator it1 =non_aligned_right.find(center);
		assert(it1 == non_aligned_right.end());
		std::map<size_t , size_t> left_and_right;
		for(size_t i = 0 ; i < it->second.size();i++){
			pw_alignment p= it->second.at(i);
			p.print();
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			std::stringstream sample1;
			std::stringstream sample2;
			sample1 <<ref1 << ":" << left_1;
			sample2 <<ref2 << ":" << left_2;
			if(sample1.str()==it->first && left_1 - left >=0 && left_1 - left < gap_on_long_centers && ref1 == ref){//center on the seq (on ref1)
				al = p;
				std::cout << "r1" << std::endl;
				if(left_1 == left){
					return right_1;
				}else{
					left_and_right.insert(std::make_pair(left_1, right_1));
				}
			}
			if(sample2.str() == it->first && left_2 - left >=0 && left_2 - left < gap_on_long_centers && ref2 == ref){//center on the seq (on ref2)
				al = p;
				std::cout << "r2" << std::endl;
				if(left_2 == left){
					return right_2;
				}else{
					left_and_right.insert(std::make_pair(left_2,right_2));
				}
			}
			if(sample2.str()==it->first && left_1 - left >= 0 && left_1 - left < gap_on_long_centers && ref1 == ref){
				al = p;
				std::cout << "r3" << std::endl;
				if(left_1 == left){
					return right_1;
				}else{
					left_and_right.insert(std::make_pair(left_1, right_1));
				}
			}
			if(sample1.str() == it->first && left_2 - left >=0 && left_2 - left < gap_on_long_centers && ref2 == ref){
				al = p;
				std::cout << "r4" << std::endl;
				if(left_2 == left){
					return right_2;
				}else{
					left_and_right.insert(std::make_pair(left_2,right_2));
				}
			}
		}
		std::cout << "r5" << std::endl;	
		std::map<size_t , size_t>::iterator it2=left_and_right.begin();
		assert(it2 != left_and_right.end());
		left = it2->first;
		return it2->second;
	}

	const size_t mixing_centers::find_right_of_long_center_on_seq(std::vector<std::string>& center ,  std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers, size_t & position , size_t & reference, const std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster)const{
		std::vector<std::string> this_cent = center;
		std::map<vector<std::string>, std::vector<pw_alignment> >::iterator it = new_centers.find(this_cent);
		if(it==new_centers.end()){//Now I keep the right drection of centers then this condition can be removed
			this_cent = get_reverse(this_cent);
			std::cout << "got reversed " << this_cent <<std::endl;
			it = new_centers.find(this_cent);
			assert(this_cent.size()==1);
		}
		assert(it != new_centers.end());
		std::vector<pw_alignment> als = it->second;
		std::cout << "als size "<< als.size() <<std::endl;
		//if ref1 is the currecnt seq and left1 is the postition the return the right1;
		if(als.size()>0){
			for(size_t i =0; i < als.size(); i++){
				als.at(i).print();
				size_t ref1 = als.at(i).getreference1();
				size_t ref2 = als.at(i).getreference2();
				size_t left1, right1, left2, right2;
				als.at(i).get_lr1(left1, right1);
				als.at(i).get_lr2(left2, right2);

				if(left1 == position && ref1 == reference){
				//	std::cout<< "right1 "	<< right1<<std::endl;
				/*	if(center.at(0)== "0:0:729890" || center.at(0) == "0:0:2066158" || center.at(0) == "0:0:2288918"){
						std::cout << "al_length "<<als.at(i).alignment_length() << std::endl;
						for(size_t j = 0 ; j < als.at(i).alignment_length(); j++){//XXX very slow test
							std::cout << als.at(i).get_al_ref1().at(j) << " "<< als.at(i).get_al_ref2().at(j)<<std::endl;
							assert(als.at(i).get_al_ref1().at(j) == als.at(i).get_al_ref2().at(j));
						}
					}*/
					return right1;
				}
				if(left2 == position && ref2 == reference){//long cneter itself is on the sequence
				//	std::cout<< "right2 "	<< right2<<std::endl;
					return right2;
				}
			
			}
		}else{//If als.size()==0 check in als_in_cluster
			assert(this_cent.size()==1);
			std::cout << "no als left! "<<std::endl;
			std::vector<std::string> center_parts;
			strsep(this_cent.at(0), ":" , center_parts);
			unsigned int center_dir = atoi(center_parts.at(0).c_str());
			unsigned int center_ref = atoi(center_parts.at(1).c_str());
			unsigned int center_left = atoi(center_parts.at(2).c_str());
			std::stringstream str;
			str<<center_ref<<":"<<center_left;
			std::map<std::string, std::vector<pw_alignment> >::const_iterator it1=alignments_in_a_cluster.find(str.str());
			assert(it1 != alignments_in_a_cluster.end());
			assert(it1->second.size() != 0);
			std::cout<< center_ref << " "<< reference << " "<< center_left << " "<< position <<std::endl;
			assert(center_ref == reference && center_left == position);

			const pw_alignment al = it1->second.at(0);
			size_t ref1 = al.getreference1();
			size_t ref2 = al.getreference2();
			size_t left1, right1, left2, right2;
			al.get_lr1(left1, right1);
			al.get_lr2(left2, right2);
			if(left1 == position && ref1 == reference){
			//	std::cout<< "right1 "	<< right1<<std::endl;
				return right1;
			}
			if(left2 == position && ref2 == reference){
			//	std::cout<< "right2 "	<< right2<<std::endl;
				return right2;
			}

		}
	}

/*	void mixing_centers::add_nonaligned_regions(std::vector<std::multimap<size_t, std::string > > & centerOnSequence, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster){
		for(size_t i =0; i < data.numSequences(); i++){
			std::cout << "centers on sequence " << i << " are :"<<std::endl;
			size_t pre_pos = 0;//later on will be updated to right+1
			for(std::multimap<size_t, std::string >::iterator it = centerOnSequence.at(i).begin(); it != centerOnSequence.at(i).end();it++){
				std::cout<< it->second << " ";
				size_t pos = it->first;
				std::cout<< "prepos "<< pre_pos << " pos "<<pos <<std::endl;
				assert(pre_pos <=pos);
				if(pre_pos != pos){//There is a non aligned region in between
 					std::stringstream str;
					str<<0<<":"<<i<<":"<<pre_pos;
					//Get the context from pre_pos to pos-1
					size_t to = pos-1;
					unsigned int ref = i;
					std::string this_part = data.extract_seq_part(ref,pre_pos,to);
					std::map<std::string, std::string>::iterator context=non_aligned_context.find(this_part);//TODO infact i can also check for the reverse - complement
					if(context==non_aligned_context.end()){
						std::cout<<"if context doesn't exist"<< str.str()<<std::endl;
						non_aligned_context.insert(std::make_pair(this_part, str.str()));//TODO cen be created inside this function!
						non_aligned_clusters.insert(std::make_pair(str.str(),std::vector<std::string>()));//TODO it is not necessary , can be removed!
						all_pieces.at(i).insert(std::make_pair(pre_pos,str.str()));//TODO Create it in main and used it while makin suffix tree
						non_aligned_right.insert(std::make_pair(str.str(),pos-1));
					}else{
						std::cout<< "if context exists "<< context->second << std::endl;
						std::map<std::string, std::vector<std::string> >::iterator pseudo_center=non_aligned_clusters.find(context->second);
						assert(pseudo_center != non_aligned_clusters.end());
						pseudo_center->second.push_back(str.str());
						all_pieces.at(i).insert(std::make_pair(pre_pos,pseudo_center->first));
						non_aligned_right.insert(std::make_pair(pseudo_center->first,pos-1));

					}
				}
				pre_pos = find_right_on_seq(it->second,alignments_in_a_cluster, pos , i) + 1;
				all_pieces.at(i).insert(std::make_pair(it->first,it->second));//Add the aligned region to the all_pieces;

			}
			if(pre_pos < data.get_seq_size(i)-1){
				std::cout <<"there is a last piece! "<<std::endl;
				//add the last piece!
				std::stringstream str;
				str<<0<<":"<<i<<":"<<pre_pos;
				all_pieces.at(i).insert(std::make_pair(pre_pos,str.str()));
				non_aligned_right.insert(std::make_pair(str.str(),data.get_seq_size(i)-1));
			}
			std::cout << std::endl;
		}
	}*/
	void mixing_centers::remove_inclusive_long_centers(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster , std::vector<std::map<size_t, std::vector<std::string> > > & long_centers_on_seq, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string, size_t> & non_aligned_right,size_t & gap_in_long_centers){
		//Go over the sequences, at each position on a seq check if there any short or long als
		std::vector<std::set<size_t> > intermediate(data.numSequences());
		std::map<std::string, size_t> seen;
		for(size_t i =0; i < data.numSequences(); i++){
			std::cout << "centers on sequence " << i <<" with length "<< data.getSequence(i).length()<< " are :"<<std::endl;
			size_t pre_pos = 0;//later on will be updated to right+1
			size_t right = 0;
			size_t previous_right = 0;
			for(std::map<size_t, std::vector<std::string> >::iterator it = all_centers_on_a_seq.at(i).begin(); it!= all_centers_on_a_seq.at(i).end(); it++){
				std::cout<< it->second << " ";
				size_t pos = it->first;
				std::cout<< "prepos "<< pre_pos << " pos "<<pos <<std::endl;
				right = find_right_of_long_center_on_seq(it->second,new_centers, pos , i, alignments_in_a_cluster);//It takes into account the case that center itself is on the sequence
				assert(right != 0);
				std::cout<< "pos "<< pos << " pre right "<< previous_right << " right "<< right <<std::endl;
			//	if(pos != 0){
			//		assert(previous_right < pos);
			//	}
				if(pos> previous_right || pos==0){
					assert(pre_pos <=pos);
					if(pre_pos != pos){//There is a non aligned region in between
 						std::stringstream str;
						str<<0<<":"<<i<<":"<<pre_pos;
						std::cout<< "there is a non aligned region " << str.str()<<std::endl; //It checks if this piece is already in the non aligned container 
						//Get the context from pre_pos to pos-1
					//	size_t to = pos-1;
					//	unsigned int ref = i;
					//	std::string this_part = data.extract_seq_part(ref,pre_pos,to);
						std::map<std::string, size_t>::iterator it1= non_aligned_right.find(str.str());
						assert(it1 != non_aligned_right.end());
						seen.insert(std::make_pair(it1->first,it1->second));
					/*	if(it1==non_aligned_right.end()){
							std::stringstream str1;
							str1<<i<<":"<<pre_pos;
							std::map<std::string, std::vector<pw_alignment> >::iterator it2 = alignments_in_a_cluster.find(str1.str());
							assert(it2 != alignments_in_a_cluster.end());//Add it to the non_aligned_right!
							pw_alignment p = it2->second.at(0);
							size_t l1,l2,r1,r2;
							p.get_lr1(l1,r1);
							p.get_lr2(l2,r2);
							size_t this_right;
							if(p.getreference1()==i && l1 == pre_pos){
								this_right = r1;
							}else{
								assert(p.getreference2()==i && l2 == pre_pos);
								this_right = r2;
							}
							non_aligned_right.insert(std::make_pair(str.str(),this_right));
							seen.insert(std::make_pair(str.str(),this_right));
						}else{
							seen.insert(std::make_pair(it1->first,it1->second));
						}*/
						std::vector<std::string> temp;
						temp.push_back(str.str());
						all_pieces_long_centers.at(i).insert(std::make_pair(pre_pos,temp));

						//std::cout<<"this part of seq "<< this_part <<std::endl;
					/*	std::map<std::string, std::string>::iterator context=non_aligned_context.find(this_part);
						if(context==non_aligned_context.end()){
							std::cout<<"if context doesn't exist"<< str.str()<<std::endl;
							non_aligned_context.insert(std::make_pair(this_part, str.str()));
							non_aligned_clusters.insert(std::make_pair(str.str(),std::vector<std::string>()));
							std::vector<std::string> temp;
							temp.push_back(str.str());
							all_pieces_long_centers.at(i).insert(std::make_pair(pre_pos,temp));
							non_aligned_right.insert(std::make_pair(str.str(),pos-1));
						}else{
							std::cout<< "if context exists "<< context->second << std::endl;
							std::map<std::string, std::vector<std::string> >::iterator pseudo_center=non_aligned_clusters.find(context->second);
							assert(pseudo_center != non_aligned_clusters.end());
							pseudo_center->second.push_back(str.str());
							std::vector<std::string> temp;
							temp.push_back(pseudo_center->first);
							all_pieces_long_centers.at(i).insert(std::make_pair(pre_pos,temp));
							non_aligned_right.insert(std::make_pair(pseudo_center->first,pos-1));


						}*/
					}
					pre_pos = right +1;
					previous_right = right;
					std::cout << "this pre pos is "<< pre_pos <<std::endl;
					assert(pre_pos>0 && pre_pos<= data.getSequence(i).length());
					all_pieces_long_centers.at(i).insert(std::make_pair(it->first,it->second));//Add the aligned region to the all_pieces;
				}
				else if(pos< previous_right && right > previous_right){//XXX need to get checked
					std::cout << "overlap case "<<std::endl;
					//Add the remaining ones one by one. Maybe you can check for the longer one and keep that one and break the others. 
					for(size_t j = 0; j< it->second.size(); j++){
						pw_alignment p;
						pre_pos = pos;
						right = find_right_on_seq(it->second.at(j), alignments_in_a_cluster, pos, i, non_aligned_right, p, gap_in_long_centers);
						std::cout<< "this center right "<< right << std::endl;
						if(right <= previous_right){
							pos = right +1;
						}
						else{
							if(pre_pos == pos){
							}else{//Add the non aligned ! 
								assert(pos > pre_pos);
								//from pre_pos to pos-1 TODO complecated cus need to check the content!! :(((( 
								unsigned int this_ref  = i;
								size_t to = pos-1;
								std::string seq = data.extract_seq_part(this_ref, pre_pos, to);
								std::stringstream str;
								str << 0<<":"<<i<<":"<<pre_pos;
								std::vector<std::string> temp;
								temp.push_back(str.str());
								std::cout << "gap is added to non aligned "<< str.str() << std::endl;
								seen.insert(std::make_pair(str.str(), pos-1));
								
								all_pieces_long_centers.at(i).insert(std::make_pair(pre_pos,temp));
							}
							std::vector<std::string> temp;
							std::vector<std::string> temp1;

							std::vector<std::string> center_parts;
							strsep(it->second.at(j),":", center_parts);
							unsigned int center_dir = atoi(center_parts.at(0).c_str());
							unsigned int center_ref = atoi(center_parts.at(1).c_str());
							unsigned int center_left = atoi(center_parts.at(2).c_str());
							std::stringstream center;
							center<< 0<<":"<<center_ref<<":"<<center_left;
							temp1.push_back(it->second.at(j));
							temp.push_back(center.str());
							std::cout << "added: "<< pos << " to "<< right << std::endl;
							all_pieces_long_centers.at(i).insert(std::make_pair(pos,temp1));
							std::map<vector<std::string>, std::vector<pw_alignment> >::iterator add_al = new_centers.find(temp);// Only the forward direction is kept
							if(add_al == new_centers.end()){
								new_centers.insert(std::make_pair(temp,std::vector<pw_alignment>()));
								add_al = new_centers.find(temp);
							}
							assert(add_al != new_centers.end());
							//Swap ref if necessary!
							swap_refs(p, it->second.at(j));
							add_al->second.push_back(p);
							pos = right + 1;
						}
					}
					pre_pos = right +1;
					previous_right = right;
				}
				else{ //Add to an intermediate map and delete from long centers(There are centers that are included fully in a longer ones and seems I entered both to the long center container. The shorter one has to be removed.)
					assert(right <= previous_right && pos < previous_right);
					std::cout << "add it to the intermediate! "<<std::endl;
					std::set<size_t>::iterator it1= intermediate.at(i).find(it->first);
					assert(it1 == intermediate.at(i).end());
					intermediate.at(i).insert(it->first);
				}
			}//Add the last piece:
			if(pre_pos < data.get_seq_size(i)){
				std::cout <<"there is a last piece! "<<std::endl;
				std::stringstream str;
				str<<0<<":"<<i<<":"<<pre_pos;
				std::cout<< str.str() << std::endl;
				std::map<std::string, size_t>::iterator it1= non_aligned_right.find(str.str());
				assert(it1 != non_aligned_right.end());
				seen.insert(std::make_pair(it1->first,it1->second));
				std::vector<std::string> temp;
				temp.push_back(str.str());
				all_pieces_long_centers.at(i).insert(std::make_pair(pre_pos,temp));
			//	non_aligned_right.insert(std::make_pair(str.str(),data.get_seq_size(i)-1));
			}
		}
 		for(size_t i =0; i < data.numSequences(); i++){
			std::cout<<"size of intermediate at "<< i << " is "<< intermediate.at(i).size()<<std::endl;
			for(std::set<size_t>::iterator it = intermediate.at(i).begin() ; it != intermediate.at(i).end() ; it++){
				std::map<size_t, std::vector<std::string> >::iterator it1 = all_centers_on_a_seq.at(i).find(*it);
				assert(it1 != all_centers_on_a_seq.at(i).end());
				std::cout << "this inclusive "<< it1->second <<std::endl;
				std::map<size_t, std::vector<std::string> >::iterator longcent = long_centers_on_seq.at(i).find(*it);
				if(it1->second.size() > 1){
					assert(longcent != long_centers_on_seq.at(i).end());
					long_centers_on_seq.at(i).erase(longcent);//Remove from both all_ and long_
				}
				std::vector<std::string> center = it1->second;
				std::map<vector<std::string>, std::vector<pw_alignment> >::iterator cent = new_centers.find(center);
				if(cent==new_centers.end()){
					center = get_reverse(center);
					cent = new_centers.find(center);
				}
				assert(cent != new_centers.end());
				std::vector<pw_alignment> als = cent->second;//Delete the corresponding alignment
				if(als.size()!=0){//Shouldnt delete it if als.size == 0
				//if ref1 is the currecnt seq and left1 is the postition:
					for(size_t j =0; j < als.size(); j++){
						als.at(j).print();
						size_t ref1 = als.at(j).getreference1();
						size_t ref2 = als.at(j).getreference2();
						size_t left1, right1, left2, right2;
						als.at(j).get_lr1(left1, right1);
						als.at(j).get_lr2(left2, right2);

						if(left1 == *it && ref1 == i){//We know that centers are on the second reference
							std::cout<<"at "<< j << "right: "<< right1<<std::endl;
							als.erase (als.begin()+j);
							break;
						}
					}
				}
				all_centers_on_a_seq.at(i).erase(it1);
			}
			std::cout << "size of the long center container for seq " << i << "  is " << long_centers_on_seq.at(i).size() << std::endl;
			std::cout << "size of all centers on seq "<< i << " is "<< all_centers_on_a_seq.at(i).size() << std::endl;
		}
		non_aligned_right = seen;
	}
	void mixing_centers::swap_refs(pw_alignment & al , const std::string & center){
		std::vector<std::string> center_parts;
		strsep(center , ":", center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());

		size_t left1, right1, left2, right2;
		al.get_lr1(left1,right1);
		al.get_lr2(left2,right2);
		if(al.getreference1()== center_ref && left1 == center_left){
			//Swap!
			pw_alignment new_al(al.get_al_ref2(), al.get_al_ref1(), al.getbegin2() , al.getbegin1() , al.getend2() , al.getend1() , al.getreference2(), al.getreference1());
			al = new_al;
		}else{
			//assert!
			assert(al.getreference2()== center_ref && left2 == center_left);
		}
	}
	void mixing_centers::make_index(std::vector<std::map<size_t , std::string> > & all_pieces){//Ids are made regardless of the direction!
		size_t num = 1;
		for(size_t i =0; i < data.numSequences();i++){
			std::cout<< "all pieces on seq "<< i << std::endl;
			for(std::map<size_t , std::string>::iterator it = all_pieces.at(i).begin() ; it != all_pieces.at(i).end() ; it++){
				std::string center = it->second;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int center_dir = atoi(center_parts.at(0).c_str());
				unsigned int center_ref = atoi(center_parts.at(1).c_str());
				unsigned int center_left = atoi(center_parts.at(2).c_str());
				std::stringstream str;
				str<<center_ref<<":"<<center_left;
				std::cout<< str.str() <<std::endl;
				std::map<std::string, size_t>::iterator id = center_id.find(str.str());
				if(id == center_id.end()){
					center_id.insert(std::make_pair(str.str(),num));
					num++;
				}
			}
		}
	}
	void mixing_centers::make_long_cent_index(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::vector<std::string>, std::vector<pw_alignment> > & long_centers, std::map<std::string, size_t> & non_aligned_right){
	//	std::cout << "indices of long centers are between 1 and ";//Walking over the new_center  + non aligned container and set an index for each of them.
		size_t num = 1;
	/*	for(size_t i =0; i < data.numSequences();i++){
			for(std::map<size_t , std::vector<std::string> >::iterator it = all_pieces_long_centers.at(i).begin() ; it != all_pieces_long_centers.at(i).end() ; it++){
				std::map<std::vector<std::string>, size_t>::iterator it1 = long_center_id.find(it->second);
				if(it1==long_center_id.end()){
					long_center_id.insert(std::make_pair(it->second,num));
					num++;
				}
			}
		}*/
		for(size_t i =0; i < data.numSequences();i++){ //XXX TEST
			std::cout << "on sequence "<< i << ":" <<std::endl;
			for(std::map<size_t , std::vector<std::string> >::iterator it = all_pieces_long_centers.at(i).begin() ; it != all_pieces_long_centers.at(i).end() ; it++){
				std::cout << "at "<< it->first << std::endl;
				size_t at = it->first;
				std::vector<std::string> this_center = it->second;
				std::vector<std::string> reverse = get_reverse(this_center);
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it1= long_centers.find(it->second);
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator rev_cent= long_centers.find(reverse);//Since only forward direction of short centers are kept in new_center container
				std::cout << "this center "<< this_center.size() <<std::endl;
				std::cout << this_center << std::endl;
				if(it1 != long_centers.end()){
					size_t right = find_right_of_long_center_on_seq(it->second,long_centers, at , i, alignments_in_a_cluster);
					std::cout<< "to" << right << std::endl;
					std::cout <<"here!!"<<std::endl;
				}else if(it->second.size()==1 && rev_cent != long_centers.end()){
					size_t right = find_right_of_long_center_on_seq(it->second,long_centers, at , i, alignments_in_a_cluster);
					std::cout << "to "<<right << std::endl;
					std::cout << "was a reverse center! "<<std::endl;
				}
				else{//NON ALIGNED
					assert(it->second.size()==1);
					std::map<std::string, size_t>::iterator it2= non_aligned_right.find(it->second.at(0));
					if(it2 == non_aligned_right.end()){
						std::cout << it->second.at(0)<<std::endl;
					}
					assert(it2 != non_aligned_right.end());
					std::cout << "to "<< it2->second << std::endl;
				}
			}
		}
		std::cout<< "make indices "<<std::endl;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it = long_centers.begin() ; it != long_centers.end() ; it++){
			std::cout << it->first<<std::endl;
			std::vector<std::string> this_center = it->first;
			std::vector<std::string> reverse_center = get_reverse(this_center);
			std::map<std::vector<std::string>, size_t>::iterator it1 = long_center_id.find(it->first);
			std::map<std::vector<std::string>, size_t>::iterator rev = long_center_id.find(reverse_center);
			assert(it1==long_center_id.end());
			if(it->first.size() == 1 && rev == long_center_id.end()){//TODO check for the reverse!! //TODO update! If it works with one direction then this is unnecessay
				long_center_id.insert(std::make_pair(it->first,num));
				num++;
			}
			if(it->first.size() > 1){
				assert(rev==long_center_id.end());
				long_center_id.insert(std::make_pair(it->first,num));
				num++;
			}

		}
		for(std::map<std::string, size_t>::iterator it= non_aligned_right.begin(); it != non_aligned_right.end() ; it++){
			std::vector<std::string> temp;
			temp.push_back(it->first);
			std::map<std::vector<std::string>, size_t>::iterator it1 = long_center_id.find(temp);
			assert(it1==long_center_id.end());
			long_center_id.insert(std::make_pair(temp,num));
			num++;
		}
		std::cout << "indices of long centers are between 1 and ";
		std::cout<< num -1 << std::endl;
	}
	void mixing_centers::find_adjacencies_with_direction(std::map<int, std::set<int> > &adjacencies, std::vector<std::map<size_t , std::string> > & all_pieces){
		make_index(all_pieces);
		for(size_t i =0; i < data.numSequences(); i++){
			std::cout << "on seq " << i << std::endl;
			std::string previous;
			int pre_id = 0;
			// go over the all pieces 
			for(std::map<size_t, std::string >::iterator it1 = all_pieces.at(i).begin(); it1 != all_pieces.at(i).end();it1++){
				size_t pos = it1->first;
				std::string center = it1->second;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int center_dir = atoi(center_parts.at(0).c_str());
				unsigned int center_ref = atoi(center_parts.at(1).c_str());
				unsigned int center_left = atoi(center_parts.at(2).c_str());
				std::stringstream str;
				str<<center_ref<<":"<<center_left;
				std::cout << "pos "<< it1->first << " "<< it1->second <<std::endl;
				if(it1 == all_pieces.at(i).begin()){
					std::cout << "first center on the seq"<<std::endl;
					assert(previous.length()==0);
					assert(it1->first == 0);
					std::map<std::string, size_t>::iterator id = center_id.find(str.str());
					assert(id != center_id.end());
					if(center_dir == 0){
						pre_id = id->second;
					}else{
						pre_id = (-1)*id->second;
					}
					std::cout << "first id "<<pre_id<<std::endl;

					previous = it1->second;
				//	std::map< int , std::set<int> >::iterator adj = adjacencies.find(pre_id);
				//	if(adj == adjacencies.end()){
				//		adjacencies.insert(std::make_pair(pre_id,std::set<int>()));
				//	}
				}
				else{
					assert(previous.length() !=0);
					assert(it1 !=  all_pieces.at(i).begin());
					std::cout << "adj pos " << it1->first<<std::endl;
					std::map<std::string, size_t>::iterator id = center_id.find(str.str());
					assert(id != center_id.end());//If center itself is in center_id
					std::cout << "pre id "<<pre_id<<std::endl;
					std::map< int , std::set<int> >::iterator adj = adjacencies.find(pre_id);
					if(adj == adjacencies.end()){
						adjacencies.insert(std::make_pair(pre_id,std::set<int>()));
						adj = adjacencies.find(pre_id);
					}
					if(center_dir == 0){
						pre_id = id->second;
					}else{
						pre_id = (-1)*id->second;
					}
					std::cout<<" this id is "<< pre_id<<std::endl;
					adj->second.insert(pre_id);
					path_on_a_seq.at(i).push_back(pre_id);
					previous = id->second;
				}
			}
			if(pre_id != 0){
				path_on_a_seq.at(i).push_back(pre_id);

				std::map< int , std::set<int> >::iterator adj = adjacencies.find(pre_id);
				if(adj == adjacencies.end()){
					adjacencies.insert(std::make_pair(pre_id,std::set<int>()));
				}
			}
		}
	}
	void mixing_centers::find_adjacencies_with_direction_long_centers(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<int, std::set<int> > &adjacencies, std::map<std::vector<std::string>, std::vector<pw_alignment> > & long_centers, std::map<std::string, size_t> & non_aligned_right){//To figure out if a long center happened forward or reverse on a sequence, I check the long_centers container alignments.
		make_long_cent_index(alignments_in_a_cluster, long_centers, non_aligned_right);
		for(size_t i =0; i < data.numSequences(); i++){
			std::cout << "on seq " << i << std::endl;
			int pre_id = 0;
			for(std::multimap<size_t, std::vector<std::string> >::iterator it = all_pieces_long_centers.at(i).begin(); it != all_pieces_long_centers.at(i).end();it++){//Attention! all_piece_long_centers have the exact same direction that exists on that specific sequence, it could be different from what actually is in new_centers and long_center_id
				std::cout << "at "<< it->first << " " << it->second << std::endl;
				int this_id =0;
				size_t pos = it->first;
				std::vector<std::string> this_center = it->second;
				std::map<std::vector<std::string>, size_t>::iterator id = long_center_id.find(this_center);
				size_t direction = 0;
				this_id = id->second;
				if(id == long_center_id.end()){
					this_center = get_reverse(this_center);
					id = long_center_id.find(this_center);
					assert(id != long_center_id.end());
					direction = 1;
					this_id = (-1)*id->second;
				}
				path_on_a_seq.at(i).push_back(this_id);
			//	assert(id != long_center_id.end());
				if(it == all_pieces_long_centers.at(i).begin()){
					assert(pos == 0);
				}
				//Find the direction from long_centers map
			//	size_t direction = 2;
			//	find_direction_of_center(i,pos,direction,this_center,long_centers);				
			//	if(direction == 0){
			//		this_id = id->second;
			//	}else{
			//		assert(direction == 1);
			//		this_id = (-1)*id->second;
			//	}
				std::cout << "previous id "<<pre_id<<std::endl;
				if(it == all_pieces_long_centers.at(i).begin()){
					std::cout << "first long center on the seq"<<std::endl;
					assert(pos ==0);
					pre_id = this_id;
				}else{
					std::map< int , std::set<int> >::iterator adj = adjacencies.find(pre_id);
					if(adj == adjacencies.end()){
						adjacencies.insert(std::make_pair(pre_id,std::set<int>()));
						adj = adjacencies.find(pre_id);
					}
					adj->second.insert(this_id);
					pre_id = this_id;
				}

			}//The last piece on a seq:
			if(pre_id != 0){//pre_id == 0 could be so unlikely but might happen if the seuqnece has only one segment.
				std::map< int , std::set<int> >::iterator adj = adjacencies.find(pre_id);
				if(adj == adjacencies.end()){
					adjacencies.insert(std::make_pair(pre_id,std::set<int>()));
				}
			}

		}
	}
	void mixing_centers::find_direction_of_center(size_t & this_ref, size_t & this_left ,size_t & direction , std::vector<std::string> & center, std::map<std::vector<std::string>, std::vector<pw_alignment> > & long_centers){
	/*	std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it=long_centers.find(center);
			if(it != long_centers.end()){
				std::vector<pw_alignment> als = it->second;
				bool thisRefExists = false;
				for(size_t i =0; i < als.size(); i++){
					pw_alignment p = als.at(i);
					size_t left1,right1;
					p.get_lr1(left1,right1);
					if(p.getreference1()== this_ref && left1 == this_left){
						if(p.getbegin2()< p.getend2()){
							direction = 0;
						}else direction = 1;
						thisRefExists = true;
						break;
					}
		
				}
				assert(thisRefExists==true);
			}else{
				direction = 0;
			}*/

	}
	bool mixing_centers::find_coordinate(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & position , size_t & ref, std::string & center){
		std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.find(center);
		assert(it != alignments_in_a_cluster.end());
		assert(it->second.size() > 0);
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		std::cout << "center is "<< center <<std::endl;
		for(size_t i =0; i < it->second.size();i++){
			pw_alignment p = it->second.at(i);
			p.print();
			size_t l1,l2,r1,r2;
			p.get_lr1(l1,r1);
			p.get_lr2(l2,r2);
			if(p.getreference1()== ref && l1 == position){
				if(ref == center_ref && l1 ==center_left){
					if(p.getbegin1() < p.getend1()){
						return false;
					}else{
						return true;
					}
				}else{
					assert(p.getreference2()== center_ref && l2 == center_left);
					if(p.getbegin2() < p.getend2()){
						return false;
					}else{
						return true;
					}
				}
			}
			else if(p.getreference2() == ref && l2 == position){
				if(ref == center_ref && l2 == center_left){
					if(p.getbegin2() < p.getend2()){
						return false;
					}else{
						return true;
					}
				}else{
					assert(p.getreference1()== center_ref && l1 == center_left);
					if(p.getbegin1() < p.getend1()){
						return false;
					}else{
						return true;
					}
				}
			}
		}

	}

	void mixing_centers::test(std::map<int, std::set<int> > & adjacencies, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::string> > & all_pieces, std::map<std::string, size_t> & non_aligned_right){
		std::set<int> adj_nodes;
		for(size_t i = 0; i < data.numSequences(); i++){
			std::cout << "on seq "<< i << " :" <<std::endl;
			for(std::multimap<size_t , std::string>::iterator it = all_pieces.at(i).begin() ; it != all_pieces.at(i).end(); it++){
				size_t left = it->first;
				
				std::cout<<"piece "<< it->second << " at "<<it->first<<std::endl;
				size_t right = find_right_on_seq(it->second, alignments_in_a_cluster, left, i, non_aligned_right);
				std::vector<std::string> center_parts;
				strsep(it->second, ":" , center_parts);
				unsigned int center_dir = atoi(center_parts.at(0).c_str());
				unsigned int center_ref = atoi(center_parts.at(1).c_str());
				unsigned int center_left = atoi(center_parts.at(2).c_str());
				std::stringstream str;
				str<<center_ref<<":"<<center_left;
				std::map<std::string , size_t>::iterator it1 = center_id.find(str.str());
				assert(it1 != center_id.end());
				std::cout<<" from "<< it->first << " to "<< right << " is " << it->second << " with id "<< it1->second <<std::endl;
				std::map<int, std::set<int> >::iterator adj;
				if(center_dir == 0){
					adj = adjacencies.find(it1->second);
					assert(adj != adjacencies.end());
					if( it != all_pieces.at(i).begin()){
						std::set<int>::iterator it2 = adj_nodes.find(it1->second);
						assert(it2 != adj_nodes.end());
					}
					adj_nodes = adj->second;
				}else{
					assert(center_dir == 1);
					adj = adjacencies.find(-1*it1->second);	
					if( it != all_pieces.at(i).begin()){
						std::set<int>::iterator it2 = adj_nodes.find(-1*it1->second);
						assert(it2 != adj_nodes.end());
					}
					adj_nodes = adj->second;				
				}
			}
			std::cout<<std::endl;
		}


	}
	const std::map<std::string, size_t> mixing_centers::get_center_id()const{
		return center_id;
	}
	
	const std::map<std::vector<std::string>, size_t> mixing_centers::get_long_center_id()const{
		return long_center_id;
	}
	const size_t mixing_centers::get_right(std::string & center_name,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right)const{
		//if it is in alignments_in_a_cluster
		std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.find(center_name);
		if(it != alignments_in_a_cluster.end()){
			assert(it->second.size() != 0);
			pw_alignment p = it->second.at(0);
			size_t left_1; 
			size_t left_2;
			size_t right_1;
			size_t right_2;
			size_t ref1;
			size_t ref2;
			p.get_lr1(left_1,right_1);
			p.get_lr2(left_2,right_2);
			ref1 = p.getreference1();
			ref2 = p.getreference2();
			std::stringstream sample1;
			std::stringstream sample2;
			sample1 <<ref1 << ":" << left_1;
			sample2 <<ref2 << ":" << left_2;
			if(sample1.str()==center_name){
				std::cout<< right_1 <<std::endl;
				return right_1;
			}else{
				assert(sample2.str()==center_name);
				std::cout<< right_2 <<std::endl;

				return right_2;
			}
		}else{
			//else look for it in non_aligned right
		//	std::cout<< "non-aligned " <<std::endl;
			std::stringstream str;
			str<<0<<":"<<center_name;
		//	std::cout<< "non-aligned "<< str.str() <<std::endl;			
			std::map<std::string, size_t>::const_iterator it1 = non_aligned_right.find(str.str());
			assert(it1 != non_aligned_right.end());
		//	std::cout<< it1->second<<std::endl;
			return it1->second;
		}

	}

	std::vector<std::string> mixing_centers::get_reverse(std::vector<std::string> & center)const{
		std::vector<std::string> this_center;
		for(size_t i =0; i < center.size();i++){
			std::vector<std::string> center_parts;
			strsep(center.at(i), ":" , center_parts);
			unsigned int center_dir = atoi(center_parts.at(0).c_str());
			unsigned int center_ref = atoi(center_parts.at(1).c_str());
			unsigned int center_left = atoi(center_parts.at(2).c_str());
			std::stringstream temp;
			if(center_dir == 0){
				temp<< 1<<":"<<center_ref<<":"<<center_left;
			}else{
				assert(center_dir = 1);
				temp<< 0<<":"<<center_ref<<":"<<center_left;
			}
			this_center.push_back(temp.str());
		}
		std::reverse(this_center.begin(),this_center.end());
		return this_center;
	}
	void write_graph::write_graph_dot(std::map<int, std::set<int> > & adjacencies , std::ostream & dotfile){
		dotfile << "/Graph referenceGraph {"<<std::endl;
		for(std::map<int, std::set<int> >::iterator it = adjacencies.begin(); it != adjacencies.end(); it++){
			std::set<int> nodes = it->second;
			if(it->second.size()>0){
				for(std::set<int>::iterator it1 = nodes.begin(); it1 != nodes.end() ; it1++){
					if(it->first > 0){
						dotfile << it->first <<"+ -> ";
					}else{
						assert(it->first<0);
						dotfile << -1*it->first <<"- -> ";
					}
					if(*it1 > 0){
						dotfile << *it1 << "+ ;"<<std::endl;
					}else{
						assert(*it1 < 0);
						dotfile << -1*(*it1) << "- ;"<<std::endl;					
					}
				}
			}else{
				if(it->first > 0){
					dotfile << it->first<<"+ ;"<<std::endl;
				}else{		
					dotfile << -1*it->first<<"- ;"<<std::endl;
				}
			}
		}
		dotfile << "/};"<<std::endl;
	}
	void write_graph::write_graph_fasta(const std::string & graphout, std::ofstream & txtout, std::map<std::string ,std::vector<pw_alignment> > & alignments_in_a_cluster,std::map<std::string, size_t> & non_aligned_right){//All the cneters are saved from their forward strand
		std::cout<< "write fasta file "<<std::endl;
	//	std::ofstream longid;
	//	longid.open("centersid_name.txt");//This text file includes all the information about each center. In fact it shows where each center came from.
		std::map<std::string, size_t> center_id = mixcenter.get_center_id();//All pieces on sequences(centers and non aligned regions) names and id with no direction.
		std::ofstream gout(graphout.c_str());
		std::map< std::string , std::string> centers;
		for(std::map<std::string, size_t>::const_iterator it = center_id.begin(); it!=center_id.end(); ++it){
			std::vector<std::string> center_parts;
			strsep(it->first, ":" , center_parts);
			unsigned int center_ref = atoi(center_parts.at(0).c_str());
			unsigned int center_left = atoi(center_parts.at(1).c_str());
			std::string seqname = data.get_seq_name(center_ref);
			//find its right, save it form its left to its right in centers map.
			std::string center = it->first;
			std::cout << center << " with id " << it->second << " at "<< center_left << std::endl;
			size_t right = mixcenter.get_right(center,alignments_in_a_cluster, non_aligned_right);
			//assert(right != 0); Right can be 0 if the non aligned region had length one and includes only the first base on the sequence
			//extract the seq from left to the included right!
			std::string sub;
			std::stringstream id;
			id<<it->second;
			size_t length = right - center_left + 1;
			txtout << it->second << ":"<<seqname<<":"<<center_left<<":"<<length<<std::endl;
			dnastring seq = data.getSequence(center_ref);
			for(size_t i = center_left ; i <=right ; i++){
				sub+=seq.at(i);
			}
			centers.insert(std::make_pair(id.str(), sub));
		}
		make_fasta(gout,centers);
		txtout.close();
	}
	void write_graph::write_graph_fasta_with_long_centers(const std::string & graphout, std::ofstream & txtout,std::map<vector<std::string>, std::vector<pw_alignment> > & long_centers, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right){
	//Go over long_centers if the length is one, its id will be the id from center_index if not it will get an id which center_id size + the number of already seen long centers
		std::cout<< "write fasta file "<<std::endl;
	//	std::ofstream longid;
	//	longid.open("centersid_name.txt");//This text file includes all the information about each center. In fact it shows where each center came from.
		std::map< std::string , std::string> centers; //first string is the id, second string is the context		
		std::map<std::vector<std::string>, size_t> center_id = mixcenter.get_long_center_id();//All long centers and non aligned regions in between
		std::cout<< "center id size is "<< center_id.size() << std::endl;
		std::ofstream gout(graphout.c_str());
		for(std::map<vector<std::string>, size_t >::iterator it = center_id.begin() ; it != center_id.end(); it++){//Go over all the pieces, using long_center_id. 
			std::vector<std::string> this_center = it->first;
			txtout <<it->second<<":";//Id of a long center is added to the file
			std::string sequence;
			size_t length = 0;
			for(size_t i = 0; i < this_center.size(); i++){
			//	std::cout<< "this center is "<< this_center.at(i) << std::endl;
				std::vector<std::string> center_parts;
				strsep(this_center.at(i), ":" , center_parts);
				unsigned int center_dir = atoi(center_parts.at(0).c_str());
				unsigned int center_ref = atoi(center_parts.at(1).c_str());
				unsigned int center_left = atoi(center_parts.at(2).c_str());
				std::string seqname = data.get_seq_name(center_ref);
			//	longid <<seqname<<":"<<center_left<<":";
				std::stringstream temp;
				temp<<center_ref<<":"<<center_left;
				std::string cent = temp.str();
				size_t right = mixcenter.get_right(cent,alignments_in_a_cluster, non_aligned_right);
				size_t left = center_left;
				length +=(right-left)+1;
			//	longid<<"r"<<right;
				if(this_center.size()==1) assert(center_dir == 0);
				if(center_dir==0){
				 	sequence.append(data.extract_seq_part(center_ref,left,right));
					txtout <<seqname<<":"<<center_left<<"r"<<right<<":";

				}else{
					assert(center_dir == 1);
					std::cout << temp.str() <<"has dir = 1" << std::endl;
					sequence.append(data.extract_reverse_seq_part(center_ref,left,right));
					txtout <<seqname<<":"<<right<<"r"<<left<<":";

				}
			}
			txtout <<"l"<<length<<std::endl;
			std::stringstream this_id;
			this_id<<it->second;
			centers.insert(std::make_pair(this_id.str(), sequence));
		}
		make_fasta(gout,centers);
		txtout.close();
	}
	void write_graph::make_fasta(std::ostream & graphfasta, std::map<std::string, std::string> & centers){
		for(std::map<std::string, std::string>::iterator it = centers.begin() ; it != centers.end() ; it++){
			size_t rest = 0;
			graphfasta<<">"<<it->first <<std::endl;
			if(it->second.size() > 70){
				for(size_t i = 0; i < it->second.size()-70; i++){
					for(size_t j = i; j < i+70;j++){
						graphfasta<<it->second.at(j);
					}
					i +=69;
					rest = i+1;
					graphfasta<<std::endl;
				}
				size_t length = 0;
				for(size_t i = rest ; i < it->second.size();i++){
					graphfasta<<it->second.at(i);
					length++;
				}
				graphfasta<<std::endl;
				assert(length <= 70);
			}else{
				graphfasta<<it->second<<std::endl;
			}
		}

	}
/*
    src -- The name of one of the source sequences for the alignment. For sequences that are resident in a browser assembly, the form 'database.chromosome' allows automatic creation of links to other assemblies. Non-browser sequences are typically reference by the species name alone.
    start -- The start of the aligning region in the source sequence. This is a zero-based number. If the strand field is '-' then this is the start relative to the reverse-complemented source sequence.
    size -- The size of the aligning region in the source sequence. This number is equal to the number of non-dash characters in the alignment text field below.
    strand -- Either '+' or '-'. If '-', then the alignment is to the reverse-complemented source.
    srcSize -- The size of the entire source sequence, not just the parts involved in the alignment.
    text -- The nucleotides (or amino acids) in the alignment and any insertions (dashes) as well.
*/

void write_graph::write_maf_record(std::ostream & out, const std::string & src, size_t start, size_t size, char strand, size_t srcSize, const std::string & alignment_part) {
	std::cout << "s " << src << " " << start << " " << size << " " << strand << " " << srcSize << std::endl;

	out << "s " << src << " " << start << " " << size << " " << strand << " " << srcSize << " " << alignment_part << std::endl;
}


/*
	write half a pairwise alignment to a maf file 
	(this write the s-line only)

*/
void write_graph::write_maf_record(std::ostream & out, const all_data & data, const pw_alignment & al, size_t reference) {
	assert(reference==0 || reference==1);
	// get the requested half alignment
	size_t print_seq;
	size_t print_start;
	size_t print_end;
	std::string print_al;
	//print size not end
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
	std::string accname("longcenter");
	std::string seqname;
	size_t size;
	if(reference == 1 && al.getreference2()>=data.numSequences()){
		std::stringstream temp;
		temp<<print_seq;
		seqname = temp.str();
		size = print_end - print_start + 1;
		if(print_end < print_start) {
			size = print_start - print_end + 1;
		}
	}else{
		al.print();
		std::cout<< data.numSequences()<<std::endl;
		size_t acc = data.accNumber(print_seq);
		accname = data.get_acc(acc);
		seqname = data.get_seq_name(print_seq);
		size = data.get_seq_size(print_seq);
	}	
	std::stringstream write_longname;
	assert(accname.length()>0);
	write_longname << accname << ':' << seqname;
	char strand = '+';
	size_t al_on_ref_size = print_end - print_start + 1;
	if(print_end < print_start) {
		al_on_ref_size = print_start - print_end + 1;
		strand = '-';
		print_start = size - print_start -1;//It reads the backward sequence


		
	} 
	write_maf_record(out, write_longname.str(), print_start, al_on_ref_size, strand, size, print_al);
	

}





/*
	compute star phylogeny multiple sequence alignment
	modifies input alignments by inserting gaps
*/
	void write_graph::msa_star_alignment(const std::string & center, std::vector<pw_alignment> & alignments) {
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

void write_graph::msa_star_al_with_long_centers(const std::vector<std::string> & longCenter, std::vector<pw_alignment> & alignment){//center is always on the second reference! 
	size_t center_length = 0; // length of the cluster center on its reference
	std::vector<size_t> gaps_before;
	std::cout<< "al size "<< alignment.size() <<std::endl;
	for(size_t i=0; i<alignment.size(); ++i) {
		size_t ref1 = alignment.at(i).getreference1();
		size_t ref2 = alignment.at(i).getreference2();
		size_t left1, right1, left2, right2;
		alignment.at(i).get_lr1(left1, right1);
		alignment.at(i).get_lr2(left2, right2);
		std::cout<< "ref2: "<< ref2<<std::endl;
		pw_alignment p1 = alignment.at(i);
		p1.print();
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
	void write_graph::write_graph_maf(std::ofstream & gout, const std::map<std::string, std::vector<pw_alignment> > & cluster_result_al, std::map<std::string, size_t> & center_id) {	
	//	std::ofstream gout(graphout.c_str());
	//	gout << "##maf version=1 scoring=probability" << std::endl;
	//	gout << "##maf version=1 cluster=number of cluster" << std::endl;

	//	gout << "# graph version " << VERSION << std::endl;
		// TODO add more metadata
		gout << "# Coordinates are 0-based.  For - strand matches, coordinates" << std::endl;
		gout << "# in the reverse complement of the 2nd sequence are used." << std::endl;
		gout << "#" << std::endl;
		gout << "# name start alnSize strand seqSize alignment" << std::endl;
		gout << "#" << std::endl;
		for(std::map<std::string, std::vector<pw_alignment> >::const_iterator it = cluster_result_al.begin(); it!=cluster_result_al.end(); ++it) {
			std::string center = it->first;
			std::vector<pw_alignment> als = it->second;
			msa_star_alignment(center, als);
			if(als.size()>0) {
				std::map<std::string, size_t>::iterator centid = center_id.find(it->first);
				assert(centid != center_id.end());
			//	gout << "# cluster " << centid->second<<std::endl;			
			//	gout << ")" << std::endl;
			//	gout << "a score=1" << std::endl; 

				gout << "a cluster="<<centid->second << std::endl; //Since score is optional and there is no requirement for the name nor the value i changed it to a cluster
				std::cout<< "for center "<< centid->second << " als are "<<std::endl;
				write_maf_record(gout, data, als.at(0), 0);
				for(size_t i=0; i<als.size(); ++i) {
					als.at(i).print();
					write_maf_record(gout, data, als.at(i), 1);
				}	
			}
		}
		gout.close();
}

	void write_graph::write_graph_maf_with_long_centers(std::ofstream & gout ,const std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::vector<std::string>, size_t>& long_center_id){
	//	std::ofstream gout(graphout.c_str());
	//	gout << "##maf version=1 scoring=probability" << std::endl;
	//	gout << "# graph version " << VERSION << std::endl;
		gout << "# Coordinates are 0-based.  For - strand matches, coordinates" << std::endl;
		gout << "# in the reverse complement of the 2nd sequence are used." << std::endl;
		gout << "#" << std::endl;
		gout << "# name start alnSize strand seqSize alignment" << std::endl;
		gout << "#" << std::endl;
		for(std::map<std::vector<std::string> , std::vector<pw_alignment> >::const_iterator it = new_centers.begin(); it != new_centers.end(); it++){//It includes all the centers
			std::vector<std::string> longCenter = it->first;
			std::cout<< "this long center: "<<std::endl;
			std::cout<< longCenter << std::endl;
			std::vector<pw_alignment> als = it->second;
		//	assert(als.size()>0);
			if(als.size()==0){
				std::cout<< "center with no memebr after merging "<< it->first << std::endl;
			//	assert(it->first.size() == 1); XXX NOT BE THE CASE FOR THE MOMENT SINCE I EXCLUDED INCLUSIVE LONG CENTERS!
			}
			std::map<std::vector<std::string>, size_t>::iterator it1 = long_center_id.find(longCenter);
			assert(it1 != long_center_id.end());
		//	for(size_t i =0; i < longCenter.size();i++){
		//		std::cout << longCenter.at(i)<<std::endl;
		//	}
			std::cout << "cent size " << longCenter.size() << std::endl;
			msa_star_al_with_long_centers(longCenter,als);
			if(als.size()>0){
				gout << "a cluster="<< it1->second << std::endl;
				//Add the center string itself
				write_maf_record(gout, data, als.at(0), 1);
				//Add one half of all the als
				for(size_t i=0; i<als.size(); ++i) {
					write_maf_record(gout, data, als.at(i), 0); 
				}	
			}
		}
	}
	void write_graph::write_graph_gfa(std::map<int, std::set<int> > & adjacencies , std::ostream & graph_gfa, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right){
		std::cout << "write gfa "<< std::endl;
		graph_gfa << "H\tVN:Z:1.0" << std::endl;
		std::map<size_t , std::string> centers; //size_t is the id of a center and std::string is the sequence content of it
		std::map<std::string, size_t> short_center_id = mixcenter.get_center_id();//All the centers and non aligned regions in between them
		std::cout << "short center id size "<< short_center_id.size()<<std::endl;
		for(std::map<std::string, size_t >::iterator it = short_center_id.begin() ; it != short_center_id.end(); it++){//Go over all the pieces, using long_center_id. 
			std::string this_center = it->first;
		//	std::cout << "this center "<<std::endl;
			assert(it->second> 0);
		//	std::cout << this_center << std::endl;
			graph_gfa << "S\t"<<it->second<<"\t";
			std::string sequence;
			size_t length = 0;
			std::vector<std::string> center_parts;
			strsep(this_center, ":" , center_parts);
			unsigned int center_ref = atoi(center_parts.at(0).c_str());
			unsigned int center_left = atoi(center_parts.at(1).c_str());
			std::string seqname = data.get_seq_name(center_ref);
			std::stringstream temp;
			temp<<center_ref<<":"<<center_left;
			std::string cent = temp.str();
			size_t right = mixcenter.get_right(cent,alignments_in_a_cluster, non_aligned_right);
			size_t left = center_left;
			length +=(right-left)+1;
		 	sequence.append(data.extract_seq_part(center_ref,left,right));
			std::stringstream this_id;
			this_id<<it->second;
			centers.insert(std::make_pair(it->second, sequence));
			graph_gfa<<sequence<<std::endl;
			int this_node_id = it->second;
			std::map<int, std::set<int> >::iterator adj = adjacencies.find(this_node_id);
			if(adj->second.size()>0){
				std::set<int> nodes = adj->second;	
				for(std::set<int>::iterator it1 = nodes.begin(); it1 != nodes.end() ; it1++){
					assert(adj->first > 0);
					graph_gfa << "L\t"<<adj->first<<"\t+\t";
					if(*it1 > 0){
						graph_gfa <<*it1<<"\t+\t0M"<<std::endl;							
					}else{
						graph_gfa <<std::abs(*it1)<<"\t-\t0M"<<std::endl;			
					}
				}
			}
			int reverse_id = -1 * this_node_id;
			adj = adjacencies.find(reverse_id);
			if(adj != adjacencies.end() && adj->second.size()>0){
			//	std::cout << "neg. node" << reverse_id <<std::endl;
				std::set<int> nodes = adj->second;	
				for(std::set<int>::iterator it1 = nodes.begin(); it1 != nodes.end() ; it1++){
					assert(adj->first < 0);
					graph_gfa << "L\t"<<std::abs(adj->first)<<"\t-\t";
					if(*it1 > 0){
						graph_gfa <<*it1<<"\t+\t0M"<<std::endl;							
					}else{
						graph_gfa <<std::abs(*it1)<<"\t-\t0M"<<std::endl;			
					}
				}
			}
		}

		size_t number = 0; //TODO instead have to write the name of the seq
		std::vector<std::vector<int> > path_on_sequences = mixcenter.get_path();
		std::cout<< path_on_sequences.size() << std::endl;
		for(size_t i =0; i < path_on_sequences.size() ; i++){
			std::cout << "on path "<< i<<std::endl;
			std::cout << path_on_sequences.at(i).size()<<std::endl;
			graph_gfa<<"P\t"<<number<<"\t";
			for(size_t j = 0 ; j < path_on_sequences.at(i).size() ; j++){
				if(path_on_sequences.at(i).at(j)> 0){
					graph_gfa<< std::abs(path_on_sequences.at(i).at(j))<<"+";
				}else{
					graph_gfa<< std::abs(path_on_sequences.at(i).at(j))<<"-";
				}
				if(j != path_on_sequences.at(i).size()-1){
					graph_gfa<<",";
				}
			}
			graph_gfa<<"\t";
			for(size_t j = 0 ; j < path_on_sequences.at(i).size()-2 ; j++){
				graph_gfa<<"0M,";
			}
			graph_gfa <<"0M"<<std::endl;
			number++;
		}
	}


	void write_graph::write_graph_gfa_with_long_centers(std::map<int, std::set<int> > & adjacencies , std::ostream & graph_gfa, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right){
		std::cout << "write gfa "<< std::endl;
		graph_gfa << "H\tVN:Z:1.0" << std::endl;
	//	graph_gfa << "H" << std::endl;
		std::map<size_t , std::string> centers; //size_t is the id of a center and std::string is the sequence content of it
		std::map<std::vector<std::string>, size_t> center_id = mixcenter.get_long_center_id();//All long centers and non aligned regions in between
		std::cout << "long center id size "<< center_id.size()<<std::endl;
		for(std::map<vector<std::string>, size_t >::iterator it = center_id.begin() ; it != center_id.end(); it++){//Go over all the pieces, using long_center_id. 
			std::vector<std::string> this_center = it->first;
			std::cout << "this center "<<std::endl;
			std::cout << this_center << std::endl;
			graph_gfa << "S\t"<<std::abs(it->second)<<"\t";
			std::string sequence;
			size_t length = 0;
			for(size_t i = 0; i < this_center.size(); i++){
				std::vector<std::string> center_parts;
				strsep(this_center.at(i), ":" , center_parts);
				unsigned int center_dir = atoi(center_parts.at(0).c_str());
				unsigned int center_ref = atoi(center_parts.at(1).c_str());
				unsigned int center_left = atoi(center_parts.at(2).c_str());
				std::string seqname = data.get_seq_name(center_ref);
				std::stringstream temp;
				temp<<center_ref<<":"<<center_left;
				std::string cent = temp.str();
				size_t right = mixcenter.get_right(cent,alignments_in_a_cluster, non_aligned_right);
				size_t left = center_left;
				length +=(right-left)+1;
				if(center_dir==0){
				 	sequence.append(data.extract_seq_part(center_ref,left,right));
				}else{
					assert(center_dir == 1);
					sequence.append(data.extract_reverse_seq_part(center_ref,left,right));
				}
			}
			std::stringstream this_id;
			this_id<<it->second;
			centers.insert(std::make_pair(it->second, sequence));
			graph_gfa<<sequence<<std::endl;
		}
	
		std::cout << "adj size "<< adjacencies.size()<<std::endl;
		for(std::map<int, std::set<int> >::iterator adj = adjacencies.begin(); adj != adjacencies.end() ; adj++){
			if(adj->second.size()>0){
			//	std::string seq;
			//	std::map<size_t , std::string>::iterator it = centers.find(std::abs(adj->first));
			//	assert(it != centers.end());
			//	char ch = it->second.back();
			//	if(adj->first < 0){
			//		ch = dnastring::complement(it->second.at(0));
			//	}
				std::set<int> nodes = adj->second;	
				for(std::set<int>::iterator it1 = nodes.begin(); it1 != nodes.end() ; it1++){
			//		graph_gfa << "S\t"<< number<< "\t";
			//		seq += ch;
			//		it = centers.find(std::abs(*it1));
			//		assert(it != centers.end());
			//		char ch1 = it->second.at(0);
			//		if(*it1 < 0){
			//			ch1 = dnastring::complement(it->second.back());
			//		}
			//		seq += ch1;
			//		graph_gfa <<seq<<std::endl;
					if(adj->first > 0){
					//	graph_gfa << "L\t"<<adj->first<<"\t+\t"<<number<<"\t+\t1M"<<std::endl;
						graph_gfa << "L\t"<<adj->first<<"\t+\t";
					}else{
					//	graph_gfa << "L\t"<<std::abs(adj->first)<<"\t-\t"<<number<<"\t+\t1M"<<std::endl;
						graph_gfa << "L\t"<<std::abs(adj->first)<<"\t-\t";
					}
					if(*it1 > 0){
					//	graph_gfa <<"L\t"<<number<<"\t+\t"<<*it1<<"\t+\t1M"<<std::endl;	
						graph_gfa <<*it1<<"\t+\t0M"<<std::endl;							
					}else{
					//	graph_gfa <<"L\t"<<number<<"\t+\t"<<(-1*(*it1))<<"\t-\t1M"<<std::endl;						
						graph_gfa <<std::abs(*it1)<<"\t-\t0M"<<std::endl;			
					}
				//	number++;
				}
			}
		}
	//	size_t number = center_id.size();
		size_t number = 0; //TODO Change it the name of sequences
		std::vector<std::vector<int> > path_on_sequences = mixcenter.get_path();
		std::cout<< path_on_sequences.size() << std::endl;
		for(size_t i =0; i < path_on_sequences.size() ; i++){
			std::cout << "on path "<< i<<std::endl;
			std::cout << path_on_sequences.at(i).size()<<std::endl;
			graph_gfa<<"P\t"<<number<<"\t";
			for(size_t j = 0 ; j < path_on_sequences.at(i).size() ; j++){
				if(path_on_sequences.at(i).at(j)> 0){
					graph_gfa<< std::abs(path_on_sequences.at(i).at(j))<<"+";
				}else{
					graph_gfa<< std::abs(path_on_sequences.at(i).at(j))<<"-";
				}
				if(j != path_on_sequences.at(i).size()-1){
					graph_gfa<<",";
				}
			}
			graph_gfa<<"\t";
			for(size_t j = 0 ; j < path_on_sequences.at(i).size()-2 ; j++){
				graph_gfa<<"0M,";
			}
			graph_gfa <<"0M"<<std::endl;
			number++;
		}
	}
	

#endif


