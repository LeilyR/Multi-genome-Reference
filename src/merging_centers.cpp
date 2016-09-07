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
		for(std::map<std::vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			vector<size_t> br = it->first;//list of nodes of a path
			size_t number = it->second;//number of the path happening
			std::vector<int> seriesOfCenters;
		//	std::cout<< "centers on these nodes are "<<std::endl;
			std::cout<< "series of centers: " <<std::endl;
			for(size_t j = 0 ; j < br.size(); j++){
				std::map<size_t, std::vector<int> >::iterator nodes= edges.find(br.at(j));
				assert(nodes != edges.end());
				for(size_t k =0; k < nodes->second.size();k++){
					std::cout<< nodes->second.at(k) <<std::endl;
					seriesOfCenters.push_back(nodes->second.at(k));//push back all the centers of all the nodes on the path br
				}
			}
			if(seriesOfCenters.at(seriesOfCenters.size()-1) == 1<<31){//The extra ending character is removed 
				seriesOfCenters.pop_back();	
			}

			std::map<std::vector<int>, int>::iterator it1 = gains.find(seriesOfCenters);
			if(it1 != gains.end()){//It is kept as it is
				intermediate.insert(make_pair(it1->first,it1->second));
			}else{
			//	size_t number_of_new_center = 0;
			//	for(size_t i = 0; i < seriesOfCenters.size(); i++){
			//		if(seriesOfCenters.at(i)== index){//??
			//			number_of_new_center = number_of_new_center + 1;
			//		}
			//	}
		//		size_t old_length = seriesOfCenters.size()+ (number_of_new_center*center_string.size());
		//		std::cout << "old_length " << old_length << std::endl;
//				size_t gain = number*old_length - (number + seriesOfCenters.size());
				size_t gain = number*seriesOfCenters.size() - (number + seriesOfCenters.size());
				intermediate.insert(make_pair(seriesOfCenters,gain));
			//	gains.insert(make_pair(seriesOfCenters,gain));
			}
		}
		gains.clear();
		for(std::map<std::vector<int>,int>::iterator it = intermediate.begin(); it != intermediate.end(); it++){
			gains.insert(make_pair(it->first,it->second));
		}

	}

	void merging_centers::merg_gain_value(suffixTree & tree){//calculates the gain value and builds second tree and so on iteratively
		std::map<size_t,std::vector<int> > edges = tree.get_edges();
		std::cout<< "size of the tree " << edges.size()<<std::endl;
		size_t original_center_numbers = centers.get_number_of_centers();
		std::cout << "original center number: "<< original_center_numbers << std::endl; 
		std::map<std::vector<size_t> , size_t> counts = tree.get_branches();//vector<size_t> shows a path(size_t s are node indices) and size_t is its number of happening. 
		//calculating the initial gain values:
		for(std::map<std::vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			std::vector<size_t> br = it->first;//list of nodes of a path
			size_t number = it->second;
			std::vector<int> centers;
			int gain = 0;
			for(size_t j = 0 ; j < br.size(); j++){
				std::map<size_t, std::vector<int> >::iterator nodes= edges.find(br.at(j));
				assert(nodes != edges.end());
				for(size_t k =0; k < nodes->second.size();k++){
					centers.push_back(nodes->second.at(k));//push back all the centers of all the nodes on the path br
				}
			}
			if(centers.at(centers.size()-1) == 1<<31){//The extra ending character is removed 
				centers.pop_back();	
			}
		//	std::cout<< "centers " ;
		//	for(size_t j =0; j < centers.size();j++){
		//		cout<< centers.at(j) << " " ;
		//	}
		//	cout << " " <<endl;
			std::cout<< "center size: " << centers.size() << "  number of happening: "<< number << std::endl;
			gain = number*(centers.size())-(number+centers.size());
			if(centers.size()!= 0){// We may get zero when the original path only had the extra endign char.
				gains.insert(make_pair(centers,gain));//centers-->index of centers, gain of the long center
				std::cout<< "gain "<<gain <<std::endl;
			}
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
				std::cout << "check highest path: " <<std::endl;
				for(size_t j = 0; j < highest_path.size(); j++){
					std::cout << highest_path.at(j)<< " ";
				}
			//	std::cout << " " <<std::endl;
			}
			center_numbers = center_numbers +1;
			merged_centers.insert(make_pair(highest_path, center_numbers));
		}
//		std::cout << "final result: " << std::endl;
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
	void merging_centers::adding_new_centers(std::vector<std::vector<std::string> > & long_centers, std::vector<std::map<size_t, std::vector<std::string> > > & centersPositionOnASeq){//Filling in the long centers vector and centersPositionOnASeq which contain all the long centers 
		//First the initial tree is built:
		std::vector<int> h_gain;//it is only used for the first time of making tree. for the next times the center with the highest gain is used.
		int index = 0; // The index of the center with the highest gain is used frome the next rounds.
		size_t num_seq = data.numSequences();
		for(size_t seq_id = 0; seq_id < num_seq; seq_id++){
			all_current_centers.at(seq_id) = centers.get_centers(seq_id);
		}
		suffixTree tree(num_seq, centers,all_current_centers);
		tree.make_tree(h_gain,index);
		merg_gain_value(tree);
		size_t biggest_index = 0;
		std::vector<std::string> sequence_of_centers;
		std::cout << "merged centers size: "<< merged_centers.size()<<std::endl;
		for(std::map<std::vector<int>,size_t>::iterator it = merged_centers.begin(); it != merged_centers.end(); it++){//The new indices are converted to their name
			std::cout << " merged_center: " << it->second << std::endl;
			for(size_t j =0;j< it->first.size();j++){ std::cout<< it->first.at(j)<< " ";}
			std::cout<< " "<<std::endl;
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
					std::cout << "id " << id <<std::endl;
					if(id > centers.get_number_of_centers()){
						for(std::map<std::vector<int>,size_t>::iterator it1 = merged_centers.begin(); it1 != merged_centers.end(); it1++){
							if(it1->second == id){
								std::vector<int> temp;
								std::cout << "j "<< j << std::endl;
								if(j != 0){
									for(size_t i = 0; i < j;i++){
										temp.push_back(list.at(i));	
									}
									std::cout<< "temp0: "<<std::endl;
									for(size_t i = 0; i < temp.size(); i ++){
										std::cout << temp.at(i)<< " ";
									}
									std::cout << " " <<std::endl;
								}
								for(size_t i = 0; i < it1->first.size();i++){
									temp.push_back(it1->first.at(i));
								}
								std::cout<< "temp1: "<<std::endl;
								for(size_t i = 0; i < temp.size(); i ++){
									std::cout << temp.at(i)<< " ";
								}
								std::cout << " " <<std::endl;
								for(size_t i = j+1; i < list.size();i++){
									temp.push_back(list.at(i));
								}
								std::cout<< "temp: "<<std::endl;
								for(size_t i = 0; i < temp.size(); i ++){
									std::cout << temp.at(i)<< " ";
								}
								list = temp;
								j = j + it1->first.size()-1;
								std::cout << "j1 "<< j <<std::endl;
								break;
							} else continue;
						}
					}
				}
				ThereIsStillABigID = false;
				for(size_t j =0; j < list.size(); j++){
					id = list.at(j);
					std::cout << "j2 " << j <<std::endl;
					std::cout << "id "<< id <<" num of centers "<< centers.get_number_of_centers() <<std::endl;
					if(id > centers.get_number_of_centers()){
						ThereIsStillABigID = true;
						std::cout << "Still a new center! "<<std::endl;
						break;
					}else continue;
				}
			}
			std::cout << "list size is "<< list.size() <<std::endl;
			for(size_t j =0; j < list.size(); j++){
				int id = list.at(j);
				std::string center = centers.find_center_name(id);
				sequence_of_centers.push_back(center);
				std::cout << center << std::endl;
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
		std::cout << "indices"<<std::endl;
		for(std::map<std::string, int>::iterator it = center_index.begin(); it!=center_index.end();it++){
			std::cout << it->first << " " << it->second <<std::endl;
		}		
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

	void merging_centers::create_alignment(std::vector<std::vector<std::string> > & long_centers, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string , std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::vector<std::string> > > & centersPositionOnASeq, std::vector<std::multimap<size_t , std::string> > & centerOnSequence){
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
						size_t right_of_last_piece=0;
						size_t length = 0;
						size_t position = it->first;
						std::cout << "position " << position << std::endl;
						find_long_center_length(it->second, alignments_in_a_cluster,j,position,length,right_of_last_piece);//It returns length of the long center and its right position on the sequence
						std::cout << "end of last piece "<< right_of_last_piece << " length "<< length << std::endl;
						al.setend1(right_of_last_piece);
						size_t right_one = right_of_last_piece;
						//Up to now first reference of an pw-alignment is set. Center is always considered as the second reference.
						std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(it->second);
						assert(new_cent != new_centers.end());
						if(new_cent->second.size() == 0){//when we are at the first sequence that this long center is located	
							al.setreference2(artificial_ref);
							artificial_ref = artificial_ref+1;
							al.setbegin2(0);//Long center is considered as an artificial sequence
							al.setend2(length-1);
						}else{
							pw_alignment p1 = new_cent->second.at(0);
							al.setreference2(p1.getreference2());
							assert(p1.getreference2() >=data.numSequences());
						//	al.setbegin2(p1.getbegin2());//Long center is considered as an artificial sequence
						//	al.setend2(p1.getend2());
							al.setbegin2(0);//Long center is considered as an artificial sequence
							al.setend2(length-1);

						}
						std::cout<< "begin2 "<< al.getbegin2() << " end2 " << al.getend2() << std::endl;
						//Second reference is already set.
						//We are setting samples here:
						size_t right = 0; //end of the previous center on a long center
						std::cout<< "it second size " << it->second.size() <<std::endl;
					//	set_smaples(it->second,j);
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
							size_t left_of_a_sample;
							for(std::multimap<size_t , std::string>::iterator leftOfSample = centerOnSequence.at(j).begin(); leftOfSample !=  centerOnSequence.at(j).end(); leftOfSample++){
								if (center == leftOfSample->second && leftOfSample->first >= it->first && leftOfSample->first >= right){
									left_of_a_sample = leftOfSample->first;
									std::cout<< "left is set! " << left_of_a_sample <<std::endl;
									break;
								}
							}
							assert(it3 != alignments_in_a_cluster.end());
							std::cout << "number of als " << it3->second.size()<<std::endl;
							//Setting samples:
							set_samples(it3->second,reverse, sample1,sample2,left_of_a_sample,j,cent_ref,cent_left,right);
						/*	for(size_t k = 0; k < it3->second.size();k++){
								pw_alignment p = it3->second.at(k);
								size_t r1,r2,l1,l2;
								p.get_lr1(l1,r1);
								p.get_lr2(l2,r2);
								std::vector<bool> sample1_p = p.getsample1();
								std::vector<bool> sample2_p = p.getsample2();
								std::cout << " sample1_p.size() " << sample1_p.size()<< " " << sample2_p.size()<< " " << p.alignment_length() <<std::endl;
								std::cout<< "ref 1 "<< p.getreference1() <<" ref2 "<< p.getreference2()<<" j " << j << " l1 " << l1 << " l2 " << l2 << " left_of_a_sample "<< left_of_a_sample << std::endl;
								std::vector< std::vector<bool> >sample;
								p.get_reverse_complement_sample(sample);
								std::cout << "r1 " << r1 << " r2 " << r2 <<std::endl;
								if(p.getreference1()== j && l1 == left_of_a_sample && p.getreference2()== cent_ref && l2 == cent_left){
									std::cout << "center on the second ref "<<std::endl;
									right = r1;
									if(p.getbegin1() < p.getend1()){
										for(size_t m =0; m < sample1_p.size(); m++){
											sample1.push_back(sample1_p.at(m));
											sample2.push_back(sample2_p.at(m));
										}
									}else{//I should always think about turning the corresponding center!
										for(size_t m = 0; m < sample.at(0).size();m++){
											sample1.push_back(sample.at(0).at(m));
											sample2.push_back(sample.at(1).at(m));
										//	reverse++;
										}
											reverse++;

									}
									break;
								}
								else if(p.getreference2() == j && l2 == left_of_a_sample && p.getreference1() == cent_ref && l1 == cent_left){
									std::cout << "center on the first ref" << std::endl;
									right = r2;
									if(p.getbegin2()<p.getend2()){
										std::cout << "forward "<<std::endl;
										for(size_t m =0; m < sample2_p.size();m++){
											sample1.push_back(sample2_p.at(m));
											sample2.push_back(sample1_p.at(m));
										}
									}else{
										std::cout << "reverse "<<std::endl;
										for(size_t m = 0; m < sample.at(1).size();m++){
											sample1.push_back(sample.at(1).at(m));
											sample2.push_back(sample.at(0).at(m));
										//	reverse ++;
										}
											reverse ++;
									}
									break;
								}else if(p.getreference2() == j && cent_ref == j && cent_left == l2 && cent_left == left_of_a_sample){//Instead of using samples sequence bases should be used! Samples are including gaps! 
									std::cout << "center on the ref - second sample "<<std::endl;
									right = r2;
								//	if(cent_dir == 0){ //seq from l2 to r2 for both refs
										std::vector<bool> intermediate;
										for(size_t m =l2; m <= r2; m++){
											char base = data.getSequence(j).at(m);
										//	std::cout << base <<std::endl;
											std::vector<bool> bits;
											pw_alignment::get_bits(base,bits);
										//	std::cout << "bits " << bits.at(0) << " " << bits.at(1) << " " <<bits.at(2) <<std::endl;
											for(size_t n = 0; n < 3; n++){
												intermediate.push_back(bits.at(n));
											}
										}
										std::cout << intermediate.size() << std::endl;
										for(size_t m =0 ; m < intermediate.size();m++){
											sample1.push_back(intermediate.at(m));
											sample2.push_back(intermediate.at(m));
										}
								//	}
								//	else{ //rev_comp of seq from l2 to r2 
								//		for(size_t m = 0; m < sample.at(1).size();m++){
								//			sample1.push_back(sample.at(1).at(m));
								//			sample2.push_back(sample.at(1).at(m));
								//		}
								//			std::cout << "dunno! " <<std::endl;
								//	}
									break;
								}else if(p.getreference1()== j && cent_ref == j && cent_left == l1 && cent_left == left_of_a_sample ){ //XXX Just added cent_left == left_of_a_sample
									std::cout << "center on the ref - first sample "<<std::endl;
									right = r1;
							//		if(cent_dir == 0){
										std::vector<bool> intermediate;
										for(size_t m =l1; m <= r1; m++){
											char base = data.getSequence(j).at(m);
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
								//	}else{
								//		for(size_t m = 0; m < sample.at(0).size();m++){
								//			sample1.push_back(sample.at(0).at(m));
								//			sample2.push_back(sample.at(0).at(m));
								//		}
								//			std::cout << "dunno! " <<std::endl;
								//
								//	}
									break;
								}
								std::cout << sample1.size() << " " << sample2.size() << std::endl;
							}*/
							if(i != it->second.size()-1){//gap between two centers
								size_t left_of_next_center=0;
								for(std::multimap<size_t , std::string>::iterator leftOfNextCenter = centerOnSequence.at(j).begin(); leftOfNextCenter != centerOnSequence.at(j).end(); leftOfNextCenter++){
									size_t left_1,right_1;
									al.get_lr1(left_1,right_1);// Need to check if it is bigger than the previous right !
									if (it->second.at(i+1)== leftOfNextCenter->second && leftOfNextCenter->first >= left_1 && leftOfNextCenter->first <= right_1 && leftOfNextCenter->first > right){//A center can happen couple of times on a sequence, the right one should be chosen!
										left_of_next_center = leftOfNextCenter->first;
										std::cout<< "left of next center is set! " << it->second.at(i+1) <<" "<< left_of_next_center <<std::endl;
										break;
									}
								}
								vector<bool> middle_part_of_sample;
								vector<bool> gap_sample;
								std::cout << "right+1 "<< right+1 << " left "<< left_of_next_center <<std::endl;
							//	size_t counter = 0;
								for(size_t m = right+1 ; m < left_of_next_center ; m++){
							//		counter +=1;
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
							//	std::cout << "counter "<< counter <<std::endl;
								std::cout<< "size of middle part is " << middle_part_of_sample.size() << std::endl;
								for(size_t m = 0 ; m < middle_part_of_sample.size();m++){
									sample1.push_back(middle_part_of_sample.at(m));
									sample2.push_back(gap_sample.at(m));
								}
							}
						}
						std::cout << "sam1 "<< sample1.size() << " sam2 " << sample2.size() << std::endl;
						al.set_alignment_bits(sample1,sample2);
						std::cout << "alignment is: " <<std::endl;
						al.print();
					//	std::vector<char> seq_base;
					//	std::vector<char>al_base;
					//	std::cout << al.alignment_length()<< " " << al.getbegin1() << " " << al.getend1() << std::endl;
					//	for(size_t m = al.getbegin1(); m<=al.getend1();m++){
					//		char base = data.getSequence(j).at(m);
					//		seq_base.push_back(base);
					//	}
					//	for(size_t m =0; m < al.alignment_length();m++){
					//		if(al.get_al_ref1().at(m) != '-'){
					//	//		al_base.push_back(al.get_al_ref1().at(m));
					//		}else{
					//			std::cout<< "there is a gap!"<<std::endl;
					//		}
					//	}
					//	std::cout << "al base size "<< al_base.size() << " seq base size "<< seq_base.size() << std::endl;
					//	assert(seq_base.size() == al_base.size());
					//	for(size_t m =0; m < al_base.size();m++){
					//		assert(al_base.at(m)==seq_base.at(m));
					//	}
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
							al.setbegin2(length-1);
							al.setend2(0);
							new_cent->second.push_back(al);	
						}else if(reverse > 0 && reverse != it->second.size()){
								std::cout<< "centers are happened in different direction so the al is not added."<<std::endl;
						}else{
							assert(reverse == 0);
							new_cent->second.push_back(al);
							std::cout << "al is added! "<<std::endl;
							al.print();
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
	void merging_centers::remove_fully_reverse_refs(vector<vector<std::string> > & local_long_centers, std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers){
	//Translate centers to their indices and check for negetive reverse of each center. If it exists one of them is chosen and the other one is considered as the reverse of that reference.	
		std::cout << "remove reverse centers "<< std::endl;
		std::map<std::vector<int> , std::vector<std::string> >long_center_index;
		for(size_t j =0;j< local_long_centers.size();j++){
			std::vector<std::string> center = local_long_centers.at(j);
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
		std::vector<size_t> positions = centers.get_long_center_position(seq_id , current_long_center);//TODO !!! need to check for its substings as well!!
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
	void merging_centers::find_long_center_length(std::vector<std::string> & centers, std::map<std::string,std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & seq_id,size_t & position ,size_t & center_length, size_t & end_of_last_piece){
//		end_of_last_piece = 0;
		size_t current_position  = position;
		center_length = 0;
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
			for(size_t j =0; j < it->second.size();j++){
				pw_alignment p = it->second.at(j);
				size_t r1,r2,l1,l2;
				p.get_lr1(l1,r1);
				p.get_lr2(l2,r2);
				size_t ref1 = p.getreference1();
				size_t ref2 = p.getreference2();
		//		size_t centerPosition = data.getSequence(seq_id).length();
				std::cout << "l1 " << l1 << " r1 "<< r1 << " l2 "<< l2 << " r2 " << r2 << " ref1 "<<ref1 << " ref2 "<< ref2<<std::endl;
				std::cout << "seq id " << seq_id << " cent ref "<< center_ref << " current pos " << current_position << " cent lef "<< center_left << " allowed gap is " << ALLOWED_GAP<< std::endl;
				//If center_ref is not the seq_id 
		/*		for ( std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(seq_id).begin() ; it1 != centerOnSequence.at(seq_id).end();it1++){
					if(it1->second == centers.at(i)){//It is wrong cus a single center may happen more than once on a seq
						centerPosition = it1->first;
						std::cout << "center position " << centerPosition << std::endl;
						break;
					}
				}
				assert(centerPosition != data.getSequence(seq_id).length());*/
				if((ref1== seq_id && ref2== center_ref && l1 >= current_position && l1 <= current_position + ALLOWED_GAP && l2 ==center_left)){//when center is on ref2
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
				else if(ref2== seq_id && ref1== center_ref && l2 >= current_position && l2<= current_position + ALLOWED_GAP && l1 ==center_left){//when center is on ref1
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
				//If we need to align the same piece against itself
				else if(ref1 == seq_id && seq_id == center_ref && center_left == l1 && center_left >= current_position && center_left <= current_position + ALLOWED_GAP){
					center_length +=r1-l1+1;
					current_position = r1;
					std::cout << "center length here! "<< center_length  <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r1;
						std::cout << "end point2 "<< end_of_last_piece << " r1 "<< r1 <<std::endl;
					}
					break;
				}
				else if(ref2 == seq_id && seq_id == center_ref && center_left == l2 && center_left >= current_position && center_left <= current_position + ALLOWED_GAP){
					center_length +=r2-l2+1;
					current_position = r2;
					std::cout << "center length there! " << center_length <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r2;
						std::cout << "end point3 "<< end_of_last_piece << " r2 "<< r2 <<std::endl;

					}
					break;

				}


			}
		}
	}
	void merging_centers::set_samples(std::vector<pw_alignment>& als,size_t & reverse, std::vector<bool> & sample1, std::vector<bool> & sample2, size_t & left_of_a_sample, size_t & seq_id, unsigned int & cent_ref, unsigned int & cent_left, size_t & right){
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
				std::cout << "center on the second ref "<<std::endl;
				right = r1;
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
				break;
			}else if(p.getreference2() == seq_id && l2 == left_of_a_sample && p.getreference1() == cent_ref && l1 == cent_left){
				std::cout << "center on the first ref" << std::endl;
				right = r2;
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
				break;
			}else if(p.getreference2() == seq_id && cent_ref == seq_id && cent_left == l2 && cent_left == left_of_a_sample){//Instead of using samples sequence bases should be used! Samples are including gaps! 
				std::cout << "center on the ref - second sample "<<std::endl;
				right = r2;
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
				break;
			}else if(p.getreference1()== seq_id && cent_ref == seq_id && cent_left == l1 && cent_left == left_of_a_sample ){ //XXX Just added cent_left == left_of_a_sample
				std::cout << "center on the ref - first sample "<<std::endl;
				right = r1;
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
				break;
			}
			std::cout << sample1.size() << " " << sample2.size() << std::endl;
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
	void mixing_centers::add_original_centers_to_long_centers(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster,std::map<std::string, std::vector<pw_alignment> > & merged_als){
		std::cout<<"new centers before: "<< new_centers.size() <<std::endl;
		for(std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it= new_centers.begin() ; it != new_centers.end(); it++){
			std::cout << "cent: "<< it->first <<std::endl;
			for(size_t i = 0; i < it->second.size();i++){
				it->second.at(i).print();
			}
		}
		std::cout<< "original cluster in: "<< alignments_in_a_cluster.size() <<std::endl;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it= alignments_in_a_cluster.begin() ; it != alignments_in_a_cluster.end() ;it++){
			std::cout << "cent: "<< it->first <<std::endl;
			for(size_t i = 0; i < it->second.size();i++){
				it->second.at(i).print();
			}
			std::stringstream str;
			str<<0<<":"<<it->first;//Forward direction is added to all the centers
			std::vector<std::string> temp;
			temp.push_back(str.str());
			std::map<std::string,std::vector<pw_alignment> >::iterator find_al= merged_als.find(it->first);
			if(find_al==merged_als.end()){
				//Add it to new_centers
				new_centers.insert(std::make_pair(temp,it->second));
			}else{
				//Add the rest of als but not the existing ones
				for(size_t i =0; i < it->second.size();i++){
					bool equal = false;
					for(size_t j = 0; j < find_al->second.size();j++){
						equal =it->second.at(i).equals(find_al->second.at(j));
						if(equal == true){
							break;
						}
					}
					if(equal == false){
						std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator new_cent=new_centers.find(temp);
						if(new_cent==new_centers.end()){
							new_centers.insert(std::make_pair(temp, std::vector<pw_alignment>()));
							new_cent=new_centers.find(temp);
						}
						new_cent->second.push_back(it->second.at(i));
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

	}
	void mixing_centers::swap_references(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers){
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
	}
	
	//computing the weight of all the cluster centers
	void mixing_centers::calculate_centers_weight(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, std::vector<std::string> > & global_results, std::map<std::string, unsigned int> & weight, std::map<std::vector<std::string>, size_t> &  numberOfACenter){//On original centers:
		size_t max_bit = 8;	
		size_t max_members = 0;//Returns the largest cluster
		for(std::map<std::string, std::vector<std::string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){	
	        	if(it->second.size()> max_members){		
				max_members = it->second.size();
			}
		
		}	
		for(std::map<std::string, std::vector<std::string> >::iterator it=global_results.begin(); it !=global_results.end(); it++){
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
	}
	void mixing_centers::calculate_long_centers_weight(std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::vector<std::string> , unsigned int> & long_center_weight, std::map<std::vector<std::string>, size_t> & numberOfACenter){//On concatenated ones:
		size_t max_bit = 8;	
		size_t maximum_mem = 0;//At the end should be equal to the number of center occurs more than all the others
		std::cout << "here new center size is "<< new_centers.size()<<std::endl;
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
	const size_t mixing_centers::find_right_on_seq(std::string & center,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & left , size_t & ref)const{
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		std::stringstream str;
		str<<center_ref<<":"<<center_left;

		std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.find(str.str());
		if(it!= alignments_in_a_cluster.end()){
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
			std::map<std::string , size_t>::const_iterator it1 = non_aligned_right.find(center);
			assert(it1 != non_aligned_right.end());
			size_t right = it1->second;
			return right;
		}
	}
	void mixing_centers::add_nonaligned_regions(std::vector<std::multimap<size_t, std::string > > & centerOnSequence, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster){
		for(size_t i =0; i < data.numSequences(); i++){
			std::cout << "centers on sequence " << i << " are :"<<std::endl;
			size_t pre_pos = 0;//later on will be updated to right+1
			for(std::multimap<size_t, std::string >::iterator it = centerOnSequence.at(i).begin(); it != centerOnSequence.at(i).end();it++){
				std::cout<< it->second << " ";
				size_t pos = it->first;
				if(pre_pos != pos){
 					std::stringstream str;
					str<<0<<":"<<i<<":"<<pre_pos;
					all_pieces.at(i).insert(std::make_pair(pre_pos,str.str()));
					non_aligned_right.insert(std::make_pair(str.str(),pos-1));
				}
				pre_pos = find_right_on_seq(it->second,alignments_in_a_cluster, pos , i) + 1;
				all_pieces.at(i).insert(std::make_pair(it->first,it->second));

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
	}
	void mixing_centers::make_index(){
	//	std::cout << alignments_in_a_cluster.size() << std::endl;
		size_t num = 1;
	/*	for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.begin();it != alignments_in_a_cluster.end();it++){
			assert(it->second.size()>0);
			std::cout << "num "<< num << " center " << it->first <<" number of als "<< it->second.size()<<std::endl;
			center_id.insert(std::make_pair(it->first,num));
			num++;
		}*/
		for(size_t i =0; i < data.numSequences();i++){
			for(std::map<size_t , std::string>::iterator it = all_pieces.at(i).begin() ; it != all_pieces.at(i).end() ; it++){
				std::string center = it->second;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int center_dir = atoi(center_parts.at(0).c_str());
				unsigned int center_ref = atoi(center_parts.at(1).c_str());
				unsigned int center_left = atoi(center_parts.at(2).c_str());
				std::stringstream str;
				str<<center_ref<<":"<<center_left;
				std::map<std::string, size_t>::iterator id = center_id.find(str.str());
				if(id == center_id.end()){
					center_id.insert(std::make_pair(str.str(),num));
					num++;
				}
			}
		}
	}
	void mixing_centers::find_adjacencies_with_direction(std::map<int, std::set<int> > &adjacencies){//TODO think of a center happening itself on a seq!!!
		make_index();
		for(size_t i =0; i < data.numSequences(); i++){
			std::cout << "on seq " << i << std::endl;
			std::string previous;
			int pre_id = 0;
		//	std::cout << "number of its cneter is "<< centerOnSequence.at(i).size()<<std::endl;//TODO go over the all pieces 
		//	for(std::multimap<size_t, std::string >::iterator it1 = centerOnSequence.at(i).begin(); it1 != centerOnSequence.at(i).end();it1++){//In 'centerOnSequence' I kept different directions of centers on a seq
			for(std::multimap<size_t, std::string >::iterator it1 = all_pieces.at(i).begin(); it1 != all_pieces.at(i).end();it1++){
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
					std::map<std::string, size_t>::iterator id = center_id.find(str.str());
					assert(id != center_id.end());
					if(center_dir == 0){
						pre_id = id->second;
					}else{
						pre_id = (-1)*id->second;
					}
					std::cout << "pre id "<<pre_id<<std::endl;

					previous = it1->second;
					std::map< int , std::set<int> >::iterator adj = adjacencies.find(pre_id);
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
					adj->second.insert(pre_id);
					previous = id->second;
				}
			}
			if(pre_id != 0){
				std::map< int , std::set<int> >::iterator adj = adjacencies.find(pre_id);
				if(adj == adjacencies.end()){
					adjacencies.insert(std::make_pair(pre_id,std::set<int>()));
				}
			}
		}
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

	void mixing_centers::test(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster){
		for(size_t i = 0; i < data.numSequences(); i++){
			std::cout << "on seq "<< i << " :" <<std::endl;
			for(std::multimap<size_t , std::string>::iterator it = all_pieces.at(i).begin() ; it != all_pieces.at(i).end(); it++){
				size_t left = it->first;
				size_t right = find_right_on_seq(it->second, alignments_in_a_cluster, left, i);
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
			}
			std::cout<<std::endl;
	
		}


	}
	const std::map<std::string, size_t> mixing_centers::get_center_id()const{
		return center_id;
	}
	
	const size_t mixing_centers::get_right(std::string & center_name,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster)const{
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
				return right_1;
			}else{
				assert(sample2.str()==center_name);
				return right_2;
			}
		}else{
			//else look for it in non_aligned right
			std::stringstream str;
			str<<0<<":"<<center_name;
			std::map<std::string, size_t>::const_iterator it1 = non_aligned_right.find(str.str());
			assert(it1 != non_aligned_right.end());
			return it1->second;
		}

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
				dotfile << it->first<<"+ ;"<<std::endl;
				dotfile << it->first<<"- ;"<<std::endl;
			}
		}
		dotfile << "/};"<<std::endl;
	}
	void write_graph::write_graph_fasta(const std::string & graphout, std::map<std::string ,std::vector<pw_alignment> > & alignments_in_a_cluster){//All the cneters are saved from their forward strand
		std::map<std::string, size_t> center_id = mixcenter.get_center_id();
		std::ofstream gout(graphout.c_str());
		std::map< size_t , std::string> centers;
		for(std::map<std::string, size_t>::const_iterator it = center_id.begin(); it!=center_id.end(); ++it){
			std::vector<std::string> center_parts;
			strsep(it->first, ":" , center_parts);
			unsigned int center_ref = atoi(center_parts.at(0).c_str());
			unsigned int center_left = atoi(center_parts.at(1).c_str());
			//find its right, save it form its left to its right in centers map.
			std::string center = it->first;
			size_t right = mixcenter.get_right(center,alignments_in_a_cluster);
			assert(right != 0);
			//extract the seq from left to right!
			std::string sub;
			dnastring seq = data.getSequence(center_ref);
			for(size_t i = center_left ; i <=right ; i++){
				sub+=seq.at(i);
			}
			centers.insert(std::make_pair(it->second, sub));
		}
		make_fasta(gout,centers);
	}
	void write_graph::make_fasta(std::ostream & graphfasta, std::map<size_t, std::string> & centers){
	for(std::map<size_t, std::string>::iterator it = centers.begin() ; it != centers.end() ; it++){
		size_t rest = 0;
		graphfasta<<">"<<it->first<<std::endl;
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



	

#endif


