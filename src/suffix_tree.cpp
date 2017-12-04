#include "suffix_tree.hpp"


//Finds all the centers the happen on a sequence:
	finding_centers::finding_centers(const all_data & d):data(d),AlignmentsFromClustering(data.numSequences()), centersOfASequence(data.numSequences()),initial_suffixes(data.numSequences()){
		INDEX = 1;
	}
	finding_centers::~finding_centers(){}
	void finding_centers::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//set alignments of each reference
		AlignmentsFromClustering.resize(data.numSequences());
		for (std::map<std::string,std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
			assert(it->second.size() != 0);
			for(size_t j = 0; j < it->second.size();j++){
				pw_alignment p = it->second.at(j);
				size_t left1; 
				size_t left2;
				size_t right1;
				size_t right2;
				p.get_lr1(left1,right1);
				p.get_lr2(left2,right2);
				AlignmentsFromClustering.at(p.getreference1()).insert(std::make_pair(left1,p));
				AlignmentsFromClustering.at(p.getreference2()).insert(std::make_pair(left2,p));
			}
		}
	/*	for(size_t i = 0 ; i < data.numSequences();i++){
			std::cout << "on sequence " << i << ":" << std::endl;
			for(std::map<size_t, pw_alignment>::iterator it1 = AlignmentsFromClustering.at(i).begin(); it1 != AlignmentsFromClustering.at(i).end(); it1++){
				pw_alignment p = it1->second;
				p.print();
				std::cout << "position is  "<< it1->first << std::endl;
			}
		}*/
	}
	void finding_centers::findMemberOfClusters(std::map<std::string,vector<pw_alignment> > & alignmentsOfClusters){//fills in a map with centers and their members and removes the direction from members
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){
			assert(it->second.size() != 0);
		//	std::cout << "no of als " << it->second.size()  << " " << it->first << std::endl;
			for(size_t i = 0; i < it->second.size(); i++){
				size_t ref1;
				size_t ref2;
				size_t left1;
				size_t left2;
				size_t right1;
				size_t right2;
				pw_alignment p =it->second.at(i);
			//	p.print();
				p.get_lr1(left1,right1);
				p.get_lr2(left2,right2);
				ref1 = p.getreference1();
				ref2 = p.getreference2();
				std::stringstream sample1;
				std::stringstream sample2;
				sample1<< ref1 << ":" << left1;
				sample2<< ref2 << ":" << left2;
				std::stringstream member1;
				std::stringstream member2;
				member1 << ref1 << ":"<<left1;
				member2 << ref2 << ":" << left2;
				if(sample1.str()==it->first){
					std::map<std::string,std::string>::iterator mem = memberOfCluster.find(member2.str());
					assert(mem == memberOfCluster.end());
					memberOfCluster.insert(make_pair(member2.str(),it->first));
				}else{
					assert(sample2.str()==it->first);
					std::map<std::string,std::string>::iterator mem = memberOfCluster.find(member1.str());
					assert(mem == memberOfCluster.end());

					memberOfCluster.insert(make_pair(member1.str(),it->first));
				}
			}
			std::map<std::string,std::string>::iterator mem = memberOfCluster.find(it->first);
			assert(mem == memberOfCluster.end());

			memberOfCluster.insert(std::make_pair(it->first,it->first));
		}
	//	std::cout << "memberOfCluster size "<< memberOfCluster.size() << std::endl;
	}
	void finding_centers::center_frequency(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::vector<std::map<size_t, std::string> > & centerOnseq){//it basically returns indices of centers on each sequence
		setOfAlignments(alignmentsOfClusters);//set all the alignments on each sequence
	//	std::cout<<"end of al set"<<std::endl;
		findMemberOfClusters(alignmentsOfClusters);// returns strings that are member of a cluster in memberOfCluster
		//First: fill in the "center_index" vector
	//	size_t index_counter = 1;
	//	std::map<int,bool> if_center;
	/*	for(std::map<std::string, std::vector<pw_alignment> >::iterator it2=alignmentsOfClusters.begin(); it2 != alignmentsOfClusters.end();it2++){ //Create index for both directions of each center by adding direction to their name.
			std::stringstream str1;
			std::stringstream str2;
			std::string center = it2->first;
			str1<<0<<":"<<center;
			str2<<1<<":"<<center;
			center_id.insert(std::make_pair(str1.str(),index_counter));
			center_index.insert(std::make_pair(index_counter, str1.str()));
			if_center.insert(std::make_pair(index_counter,false));
			center_id.insert(std::make_pair(str2.str(),(-1*index_counter)));
			center_index.insert(std::make_pair((-1*index_counter),str2.str()));
			if_center.insert(std::make_pair(-1*index_counter,false));
			index_counter++;	
		}*/
	//	std::cout<< "center index size before " << center_index.size()<<std::endl;
		for(size_t i = 0; i< data.numSequences(); i++){
		//	std::cout << "sequence: " << i << std::endl;
			//Second: Find all the centers on a sequence
			find_seq_centers(i,centerOnseq);
		//	std::cout<< "center size at " << i <<" is "<< centerOnseq.at(i).size()<< std::endl;
		}
	/*	for(std::map<int, bool>::iterator it = if_center.begin(); it != if_center.end(); it++){
			if(it->second == false){
				std::map<int,std::string>::iterator index = center_index.find(it->first);
				assert(index != center_index.end());
				std::string center = index->second;
				std::map<std::string,int>::iterator id=center_id.find(center);
				assert(id!=center_id.end());
				center_index.erase(index);
				center_id.erase(id);
			}
		}*/
	//	std::cout<< "center index size after " << center_index.size() << " " <<center_id.size()<<std::endl;
		//Just an assertion:
		std::map<int, std::string>::reverse_iterator id = center_index.rbegin();
	//	std::cout << id->first <<  " "<< alignmentsOfClusters.size() <<std::endl;
		assert(id->first <= alignmentsOfClusters.size());
	}
	void finding_centers::create_long_centers_candidates(std::vector<std::map<size_t, std::string> > & centerOnseq, size_t & gap_in_long_centers){
		//Find potential candidates for being a long center. 
		//We are going to find those centers that have less than ALLOWED_GAP base pairs  distances. Their indices are saved in a vector.
		for(size_t i = 0; i< data.numSequences(); i++){
			std::cout << "long centers on seq "<< i << std::endl;
			find_long_centers(i, centerOnseq.at(i), gap_in_long_centers);
			std::cout << "check for fully reversed ones" <<std::endl;
			//At this stage i need to check for fully reversed centers and merge two clusters in to one.
			std::map<size_t , std::vector<int> > intermediate;
			for(std::map<size_t,std::vector<int> >::iterator it = centersOfASequence.at(i).begin() ; it != centersOfASequence.at(i).end();it++){
				std::map<size_t,std::vector<int> >::iterator rev_pos = intermediate.find(it->first);
				assert(rev_pos == intermediate.end());
					std::map<int,std::string>::iterator centid = center_index.find(it->second.at(0));
				//	std::cout<< it->second.at(0)<<std::endl;
					assert(centid != center_index.end());
					std::string center = centid->second;
					std::vector<std::string> center_parts;
					strsep(center, ":" , center_parts);
					unsigned int dir = atoi(center_parts.at(0).c_str());
					unsigned int ref = atoi(center_parts.at(1).c_str());
					unsigned int left = atoi(center_parts.at(2).c_str());
					for(std::map<size_t , std::vector<int> >::iterator it1 = it ; it1 != centersOfASequence.at(i).end();it1++){//We can start from it instead of begin
						centid = center_index.find(it1->second.at(it1->second.size()-1));
						std::string cent = centid->second;
						std::vector<std::string> cent_parts;
						strsep(cent, ":" , cent_parts);
						unsigned int dir1 = atoi(center_parts.at(0).c_str());
						unsigned int ref1 = atoi(center_parts.at(1).c_str());
						unsigned int left1 = atoi(center_parts.at(2).c_str());
						bool Fullyreverse = true;
						if(ref == ref1 && left == left1 && dir != dir1 && it1->second.size() == it->second.size()){
							for(size_t i = 0; i < it1->second.size();i++){//We are sure that the length is higher than 1
								// All the next left and ref should be checked!
								centid = center_index.find(it1->second.at(i));
								std::string center1 = centid->second;
								std::vector<std::string> center1_parts;
								strsep(cent, ":" , center1_parts);
								unsigned int dir_c1 = atoi(center_parts.at(0).c_str());
								unsigned int ref_c1 = atoi(center_parts.at(1).c_str());
								unsigned int left_c1 = atoi(center_parts.at(2).c_str());
								centid = center_index.find(it->second.at(it1->second.size()-1+i));
								std::string center2 = centid->second;
								std::vector<std::string> center2_parts;
								strsep(cent, ":" , center2_parts);
								unsigned int dir_c2 = atoi(center_parts.at(0).c_str());
								unsigned int ref_c2 = atoi(center_parts.at(1).c_str());
								unsigned int left_c2 = atoi(center_parts.at(2).c_str());
								if(ref_c1 == ref_c2 && left_c1 == left_c2 && dir_c1 != dir_c2){
								
								}else{
									Fullyreverse = false;
									break;
								}
							}
						}else{
							Fullyreverse = false;
						}
						if(Fullyreverse == true){//The fully reverse exists and we keep one of them
							std::cout << "fully reversed! "<<std::endl;
							intermediate.insert(std::make_pair(it1->first,it1->second));
							break;
						}
					}

			}
			//Remove the intermediate from centersOfASequence
			std::cout <<"intermediate size "<< intermediate.size()<<std::endl;
			for(std::map< size_t , std::vector<int> >::iterator it = intermediate.begin();it != intermediate.end();it++){
				centersOfASequence.at(i).erase(it);
				std::map<size_t , std::vector<std::vector< int> > >::iterator it1 =long_centers_of_a_sequence.find(i);
				for(size_t j = 0; j < it1->second.size();j++){
					if(it1->second.at(j)==it->second){
						it1->second.erase(it1->second.begin()+j);
						break;
					}
				}
			}
			//Print the potential long centers:
		/*	for(std::map< size_t , std::vector<int> >::iterator it = centersOfASequence.at(i).begin();it != centersOfASequence.at(i).end();it++){
				std::cout << "position "<< it->first <<std::endl;
				for(size_t i =0; i < it->second.size();i++){
					std::cout << it->second.at(i)<<  " ";
				}
				std::cout << ""<<std::endl;
			}*/
		}
	}
	void finding_centers::find_seq_centers( size_t & seq_id , std::vector<std::map<size_t, std::string> > & centerOnseq){
		const dnastring & sequence = data.getSequence(seq_id);
		for(size_t n= 0; n < sequence.length(); n++){
			std::map<size_t, pw_alignment>::iterator it=AlignmentsFromClustering.at(seq_id).find(n);
			if(it != AlignmentsFromClustering.at(seq_id).end()){
				std::cout << "at "<< n;
				pw_alignment p = it->second;
				p.print();
				size_t l1,r1,l2,r2;
				p.get_lr1(l1,r1);
				p.get_lr2(l2,r2);
				std::stringstream member;
				member<< seq_id << ":" << n;
				std::cout << "mem" << member.str() <<std::endl;
				std::map<std::string,std::string>::iterator it1 = memberOfCluster.find(member.str());
				assert(it1 != memberOfCluster.end());
				std::string center = it1->second;				
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				bool dir;
				unsigned int cent_ref = atoi(center_parts.at(0).c_str());
				unsigned int cent_left = atoi(center_parts.at(1).c_str());
				if(p.getreference1()==cent_ref && l1==cent_left){
				//	if(p.getbegin1()<p.getend1()){//XXX CHECK!!
					if((p.getbegin1()<p.getend1()&&p.getbegin2()<p.getend2())||(p.getbegin1()>p.getend1()&&p.getbegin2()>p.getend2())){ //Dont check it here, it is checked in encoding if it occurs reverse or forward
						std::stringstream str;
						str<<0<<":"<<center;
						center = str.str();
						dir = false;	
					}else{
						std::stringstream str;
						str<<1<<":"<<center;
						center = str.str();
						dir = true;	
					}
				}else{
					assert(p.getreference2()==cent_ref && l2==cent_left);
				//	if(p.getbegin2()<p.getend2()){
					if((p.getbegin1()<p.getend1()&&p.getbegin2()<p.getend2())||(p.getbegin1()>p.getend1()&&p.getbegin2()>p.getend2())){
						std::stringstream str;
						str<<0<<":"<<center;
						center = str.str();
						dir = false;	
					}else{
						std::stringstream str;
						str<<1<<":"<<center;
						center = str.str();
						dir = true;	
					}
				}
				std::cout<< "center: "<< center << std::endl;//Direction was added to center
			//	centersOnSequence.at(seq_id).insert(std::make_pair(n,center));
			//	std::cout << " is " << center <<std::endl;
				centerOnseq.at(seq_id).insert(std::make_pair(n,center));//This is the global one. It is initialized in main function
				std::map<std::string,int>::iterator id=center_id.find(center);//Check also for the reverese and remove the if_center instead 
				if(id == center_id.end()){
					if(dir == false){
						std::stringstream str;
						str<<1<<":"<<cent_ref<<":"<<cent_left;
						std::map<std::string,int>::iterator revid=center_id.find(str.str());
						if(revid == center_id.end()){
							center_id.insert(std::make_pair(center,INDEX));
							center_index.insert(std::make_pair(INDEX,center));
							INDEX++;
						}else{
							center_id.insert(std::make_pair(center,-1*revid->second));
							center_index.insert(std::make_pair(-1*revid->second,center));
						}
					}else{
						assert(dir==true);
						std::stringstream str;
						str<<0<<":"<<cent_ref<<":"<<cent_left;
						std::map<std::string,int>::iterator revid=center_id.find(str.str());
						if(revid == center_id.end()){
							center_id.insert(std::make_pair(center,-1*INDEX));
							center_index.insert(std::make_pair(-1*INDEX,center));
							INDEX++;
						}else{
							assert(revid->second>0);
							center_id.insert(std::make_pair(center,-1*revid->second));
							center_index.insert(std::make_pair(-1*revid->second,center));
						}
					}
				}
			//	std::cout << id->second <<std::endl;
			//	std::map<int,bool>::iterator it = if_center.find(id->second);
			//	assert(it!= if_center.end());
			//	it->second = true;
			}
		}
	}
	void finding_centers::find_long_centers(size_t & seq_id,  std::map<size_t, std::string> & centerOnseq, size_t & gap_in_long_centers){
		size_t last_position;
		size_t first_position;
		std::map<size_t, std::vector<int> > AllConnectedOnes; //First one is position, second is vector of centers
		vector<int> vectorOfcenters;
		std::vector<size_t> position_of_each_suffix; 
		bool direction = true;
		size_t dir;
		size_t dir1;
		for(std::map<size_t, std::string>::iterator it = centerOnseq.begin(); it != centerOnseq.end();it++){
		//	std::cout<< "at position "<< it->first << " center name is "<< it->second<<std::endl;
			std::map<size_t, pw_alignment>::iterator it1=AlignmentsFromClustering.at(seq_id).find(it->first);
			assert(it1 != AlignmentsFromClustering.at(seq_id).end());
			pw_alignment al = it1->second;
			size_t l1,l2,r1,r2;
			al.get_lr1(l1,r1);
			al.get_lr2(l2,r2);
			assert(l1==it->first || l2 == it->first);
			if(l1 == it->first && al.getreference1()== seq_id){
				if(al.getbegin1()<al.getend1()){
					dir1 = 1;
				}else{
					dir1 = 0;
				}
			}
			else{
				assert(l2 == it->first && al.getreference2() == seq_id);
				if(al.getbegin2()< al.getend2()){
					dir1 = 1;
				}else{
					dir1 = 0;
				}
			}
			if(vectorOfcenters.size()==0){
				last_position = it->first;
				first_position = it->first;
				dir = dir1;
			}
			if(dir ==dir1){
				direction = true;
			}else{
				direction = false;
			}
		//	std::cout << "position " << it->first << " allowed gap "<< gap_in_long_centers<<std::endl;
			if((it->first - last_position) < gap_in_long_centers && direction == true){
				if(vectorOfcenters.size()!=0){
				//	std::cout << "smaller than ALLOWED_GAP! "<< it->second <<std::endl;
				}
				std::map<std::string,int>::iterator id=center_id.find(it->second);
				assert(id != center_id.end());
				int index = id->second;//TODO what if id is negative?
				vectorOfcenters.push_back(index);
				position_of_each_suffix.push_back(it->first);
			//	std::cout << "index "<< index << std::endl;		
			}else{
			/*	std::cout<< "bigger than ALLOWED_GAP: " <<std::endl;
				for(size_t i =0; i < vectorOfcenters.size();i++){
					std::cout << vectorOfcenters.at(i)<< " ";
				}
				std::cout << "-------"<<std::endl*/
				AllConnectedOnes.insert(std::make_pair(first_position, vectorOfcenters));
				//Creates suffix and their related position , then clear the 'position_of_each_suffix' 
				if(vectorOfcenters.size()>1){
					for(size_t j =0; j < vectorOfcenters.size(); j++){
						std::vector<int> suffix;
						for(size_t k =j; k < vectorOfcenters.size();k++){
							suffix.push_back(vectorOfcenters.at(k));
						}
						initial_suffixes.at(seq_id).insert(std::make_pair(suffix,position_of_each_suffix.at(j)));
					}
				}
				position_of_each_suffix.clear();
				vectorOfcenters.clear();
			}
			if(l1 == it->first && al.getreference1()== seq_id){
				last_position = r1;
				if(al.getbegin1()<al.getend1()){
					dir = 1;
				}else{
					dir = 0;
				}
			}
			if(l2 == it->first && al.getreference2() == seq_id){
				last_position = r2;
				if(al.getbegin2()< al.getend2()){
					dir = 1;
				}else{
					dir = 0;
				}
			}
		//	std::cout << "last position "<< last_position << std::endl;
		}
		//If the last center had ALLOWED_GAP with its previous one and there is no more left to go through the else
		AllConnectedOnes.insert(make_pair(first_position, vectorOfcenters));
		if(vectorOfcenters.size()>1){
			for(size_t j =0; j < vectorOfcenters.size(); j++){
				std::vector<int> suffix;
				for(size_t k =j; k < vectorOfcenters.size();k++){
					suffix.push_back(vectorOfcenters.at(k));
				}
				initial_suffixes.at(seq_id).insert(std::make_pair(suffix,position_of_each_suffix.at(j)));
			}
		}
		std::vector<std::vector<int> > temp;
		for(std::map<size_t, std::vector<int> >::iterator it = AllConnectedOnes.begin(); it != AllConnectedOnes.end(); it++){
			if(it->second.size() != 1 && it->second.size() != 0){
			//	std::cout<< "pos "<< it->first <<std::endl;
				centersOfASequence.at(seq_id).insert(std::make_pair(it->first, it->second));//connected centers that can be potential long centers
				temp.push_back(it->second);
			/*	for(size_t i =0; i < it->second.size();i++){
					std::cout << it->second.at(i)<< " ";
				}
				std::cout << "-------"<<std::endl;*/
			}
		}
		long_centers_of_a_sequence.insert(std::make_pair(seq_id, temp));
	}
//	std::multimap< size_t, std::string> finding_centers::get_sequence_centers(size_t& id)const{
//		return centersOnSequence.at(id);
//	}
	std::string finding_centers::find_center_name(int & centerIndex)const{
		std::map<int,std::string>::const_iterator it= center_index.find(centerIndex);
		assert(it != center_index.end());
		return it->second;
	}
	std::vector<size_t> finding_centers::get_long_center_position(size_t & seq_id , std::vector<std::string> & long_center){
		std::vector<size_t> position;
		std::vector<int> centers;
	//	std::cout << "current long cent " <<std::endl;

		for(size_t i =0; i < long_center.size();i++){
			std::map<std::string, int>::iterator it = center_id.find(long_center.at(i));
			assert(it != center_id.end());
		//	std::cout << long_center.at(i)<< " "<< it->second <<std::endl;
			centers.push_back(it->second);
		}
		std::multimap<std::vector<int>,size_t>::iterator suf = initial_suffixes.at(seq_id).find(centers);
		if(suf != initial_suffixes.at(seq_id).end()){
			std::pair<std::multimap<std::vector<int>,size_t>::iterator , std::multimap<std::vector<int>,size_t>::iterator > it = initial_suffixes.at(seq_id).equal_range(centers);
			for(std::multimap<std::vector<int>,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
				position.push_back(it1->second);
			}	
			return position;
		}else{// Need to loop over the initial suffixes multimap which is not so efficient TODO think of a better way
			for(std::multimap<std::vector<int>,size_t>::iterator it = initial_suffixes.at(seq_id).begin() ; it != initial_suffixes.at(seq_id).end();it++){
				if(it->first.at(0)!=centers.at(0) || it->first.size()<centers.size()){
					continue;
				}
				if(it->first.at(0)==centers.at(0) && it->first.size()>centers.size()){
					bool isNot = false;
					for(size_t i =0; i < centers.size();i++){
						if(centers.at(i) != it->first.at(i)){
							isNot = true;
							break;
						}						
					}
					if(isNot == false){
						position.push_back(it->second);
					}
				}
			}
		}
		return position;
	}
	std::map<size_t, std::vector<int> >finding_centers::get_center(size_t & seq_id)const{
		return centersOfASequence.at(seq_id);
	}
	const std::vector<std::vector<int> > & finding_centers::get_centers(size_t & seq_id)const{
		std::vector<std::vector< int> > empty;
		std::map<size_t , std::vector<std::vector<int> > >::const_iterator it =long_centers_of_a_sequence.find(seq_id);
	//	assert(it != long_centers_of_a_sequence.end());
		if(it != long_centers_of_a_sequence.end()){
			return it->second;
		}else return empty;
	}
	const size_t finding_centers::get_number_of_centers(){
		std::map<int, std::string>::reverse_iterator it = center_index.rbegin();
		std::map<int, std::string>::const_iterator it1 = center_index.begin();
		if(std::abs(it->first) > std::abs(it1->first)){
			return std::abs(it->first);
		}else{
			return std::abs(it1->first);
		}
	}
	const std::map<std::string, int> finding_centers::get_index()const{
		return center_id;
	}


void finding_centers::add_nonaligned_regions(std::vector<std::map<size_t, std::string > > & centerOnSequence, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::string> > & all_pieces, std::map<std::string, size_t> & non_aligned_right){//If they gather together and make a group an index should be refered to them!! I better add them to al_in_a_cluster too
		std::map<std::string, std::string> non_aligned_context;//string--> content of seq , string --> name of center (pos:ref)
		std::set<std::string> erased;
		assert(center_index.size()==center_id.size());
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
					std::string reverse_comp_this_string = data.extract_reverse_seq_part(ref,pre_pos,to);
					std::map<std::string, std::string>::iterator context=non_aligned_context.find(this_part);//the reverse complement can also be checked.
					std::map<std::string, std::string>::iterator rev_context=non_aligned_context.find(reverse_comp_this_string);
					if(context==non_aligned_context.end()&& rev_context==non_aligned_context.end()){ //If no other non aligned region has this context or its reverse.
						std::cout<<"if context doesn't exist"<< str.str()<<std::endl;
						non_aligned_context.insert(std::make_pair(this_part, str.str()));
						all_pieces.at(i).insert(std::make_pair(pre_pos,str.str()));
						non_aligned_right.insert(std::make_pair(str.str(),pos-1));
					}else if(context!=non_aligned_context.end()){ //If another non aligned region with the exact same context exist.
						if(rev_context != non_aligned_context.end()){
						//	std::cout << this_part << std::endl;
						//	std::cout << reverse_comp_this_string << std::endl;
						}
						if(this_part != reverse_comp_this_string){//To avoid the Palindromic regions
							assert(rev_context==non_aligned_context.end());
						}
						assert(rev_context == non_aligned_context.end() || this_part == reverse_comp_this_string);
						std::cout<< "if context exists "<< context->second<< std::endl;
						std::vector<std::string> center_parts;
						strsep(context->second, ":" , center_parts);
						unsigned int center_dir = atoi(center_parts.at(0).c_str());
						unsigned int center_ref = atoi(center_parts.at(1).c_str());
						unsigned int center_left = atoi(center_parts.at(2).c_str());
						assert(center_dir == 0);
						std::stringstream str1;
						str1<<1<<":"<<center_ref<<":"<<center_left;

						all_pieces.at(i).insert(std::make_pair(pre_pos,context->second));
						centerOnSequence.at(i).insert(std::make_pair(pre_pos, context->second));
					//	non_aligned_right.insert(std::make_pair(context->second,pos-1));
						std::map<std::string, size_t>::iterator nonaligned = non_aligned_right.find(context->second);
						if(nonaligned == non_aligned_right.end()){
						//	std::cout << "already removed "<<context->second << " "<< context->first <<std::endl;
						}
						if(nonaligned != non_aligned_right.end()){
							non_aligned_right.erase(nonaligned); //It is getting removed from the non_aligned_right container 'cus it is not a non aligned region anymore
							erased.insert(context->second);
						}else{
							std::set<std::string>::iterator test = erased.find(context->second);
							assert(test != erased.end());
						}
						std::map<std::string, int>::iterator id = center_id.find(context->second);//Should check if its reverse exists!
						std::map<std::string , int>::iterator rev_id = center_id.find(str1.str());
						if(id == center_id.end() && rev_id == center_id.end()){
						//	std::cout<< "center index size is "<< center_index.size() << std::endl;
							std::map<int, std::string>::reverse_iterator lastind = center_index.rbegin();
							size_t this_index = lastind->first + 1;
						//	std::cout<< "this index "<< this_index <<std::endl;
							center_index.insert(std::make_pair(this_index, context->second));
							center_id.insert(std::make_pair(context->second, this_index));
						}else if(id == center_id.end()){
							assert(rev_id != center_id.end());
							assert(rev_id->second < 0);
							center_index.insert(std::make_pair(-1*rev_id->second,context->second));
							center_id.insert(std::make_pair(context->second, -1*rev_id->second));
						}else if(rev_id == center_id.end()){
							assert(id != center_id.end());
							assert(id->second > 0);
						}
						std::cout << "when true " << context->second <<std::endl;
						add_to_alignments(centerOnSequence, alignments_in_a_cluster, i, pre_pos, this_part, context->second,true); //this_part --> sequence content of both member and so called center , context->second --> center name
					}else{//If its context doesn't exist itself but its reverse context does exist.
						assert(context==non_aligned_context.end() && rev_context!=non_aligned_context.end());
						std::vector<std::string> center_parts;
						strsep(rev_context->second, ":" , center_parts);
						unsigned int center_dir = atoi(center_parts.at(0).c_str());
						unsigned int center_ref = atoi(center_parts.at(1).c_str());
						unsigned int center_left = atoi(center_parts.at(2).c_str());
						assert(center_dir == 0);
						std::stringstream str1;
						str1<<1<<":"<<center_ref<<":"<<center_left;
						all_pieces.at(i).insert(std::make_pair(pre_pos,str1.str()));
						centerOnSequence.at(i).insert(std::make_pair(pre_pos, str1.str()));
					//	non_aligned_right.insert(std::make_pair(str.str(),pos-1));
						std::map<std::string, size_t>::iterator nonaligned = non_aligned_right.find(rev_context->second);
						if(nonaligned == non_aligned_right.end()){
						//	std::cout<<"if its reverse exists " << rev_context->second <<std::endl;
						//	std::cout << rev_context->first << " "<< this_part <<std::endl;
						}
						if(nonaligned != non_aligned_right.end()){
							non_aligned_right.erase(nonaligned);
							erased.insert(rev_context->second);
						}else{
							std::set<std::string>::iterator test = erased.find(rev_context->second);
							assert(test != erased.end());
						}
						std::map<std::string, int>::iterator id = center_id.find(rev_context->second);
						std::map<std::string, int>::iterator rev_id = center_id.find(str1.str());
						if(id == center_id.end() && rev_id == center_id.end()){
						//	std::cout << "here! "<<std::endl;
							std::map<int, std::string>::reverse_iterator lastind = center_index.rbegin();
							std::map<int,std::string>::iterator firstind = center_index.begin();
							size_t this_index = std::abs(lastind->first) + 1;
							if(std::abs(lastind->first) < std::abs(firstind->first)){
								this_index = std::abs(firstind->first)+1;
							}
							center_index.insert(std::make_pair(-1*this_index, str1.str()));
							center_id.insert(std::make_pair(str1.str(), -1*this_index));
							/////////JUST ADDED
							center_index.insert(std::make_pair(this_index, rev_context->second));
							center_id.insert(std::make_pair(rev_context->second,this_index));

						}else if(id == center_id.end()){ //XXX Is it even happening?? Seems it shouldn't happen 'cus it happend atleast once before.
						//	std::cout << "HERE! "<<std::endl;
							assert(rev_id != center_id.end());
							assert(rev_id->second < 0);

						}else if(rev_id == center_id.end()){
							assert(id != center_id.end());
							assert(id->second > 0);
						//	std::cout << "here1"<<std::endl;
							center_index.insert(std::make_pair(-1*id->second, str1.str()));
							center_id.insert(std::make_pair(str1.str(), -1*id->second));

						}else{
							assert(rev_id != center_id.end());
							assert(rev_id->second < 0);
							assert(id != center_id.end());
							assert(id->second > 0);
						}
						std::string rev_center = str1.str();
						add_to_alignments(centerOnSequence, alignments_in_a_cluster, i, pre_pos, reverse_comp_this_string, rev_center,false);
					}
				}
				pre_pos = find_right_on_seq(it->second,alignments_in_a_cluster, pos , i) + 1;
				all_pieces.at(i).insert(std::make_pair(it->first,it->second));//Add the aligned region to the all_pieces;

			}
			if(pre_pos < data.get_seq_size(i)){
			//	std::cout <<"there is a last piece! "<<std::endl;
				//add the last piece!
				std::stringstream str;
				str<<0<<":"<<i<<":"<<pre_pos;
				all_pieces.at(i).insert(std::make_pair(pre_pos,str.str()));
				non_aligned_right.insert(std::make_pair(str.str(),data.get_seq_size(i)-1));
			}
		//	std::cout << std::endl;
		}
		std::cout<<"center index size after adding non aligned regions "<< center_index.size()<<std::endl;
		//Making potential long centers:

	}
	void finding_centers::add_to_alignments(std::vector<std::map<size_t, std::string > > & centerOnSequence, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & ref, size_t & position, std::string & sequence, std::string & center_name , bool direction){//it adds newly made alignments between identical non aligned regions to alignmnets_in_a_cluster and AlignmentFromClustering
	//	std::cout << "add the al "<<std::endl;
		std::vector<std::string> center_parts;
		strsep(center_name, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		if(direction == false){
			assert(center_dir == 1);
		}else{
			assert(center_dir ==0);
		}
		size_t begin1 = center_left;
		size_t end1 = center_left+ sequence.length()-1;
		std::string sample1 = sequence;
		std::string sample2 = sequence;
		size_t begin2 = position;
		size_t end2 = position+sequence.length()-1;
		if(direction==false){
			begin2 = end2;
			end2 = position;
		//	sample1 = data.extract_reverse_seq_part(center_ref, begin1, end1);
		}
		pw_alignment p(sample1, sample2, begin1, begin2, end1, end2 , center_ref , ref);
		std::map<size_t , pw_alignment >::iterator al = AlignmentsFromClustering.at(ref).find(position);
		assert(al==AlignmentsFromClustering.at(ref).end());
		AlignmentsFromClustering.at(ref).insert(std::make_pair(position, p));
		std::stringstream str;
		str<<center_ref<<":"<<center_left;	
		std::map<std::string, std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.find(str.str());
	//	std::cout << str.str()<<std::endl;
		if(it == alignments_in_a_cluster.end()){
			alignments_in_a_cluster.insert(std::make_pair(str.str(), std::vector<pw_alignment>()));
			it = alignments_in_a_cluster.find(str.str());
			std::map<size_t , pw_alignment >::iterator cent = AlignmentsFromClustering.at(center_ref).find(center_left);
			assert(cent==AlignmentsFromClustering.at(center_ref).end());
			AlignmentsFromClustering.at(center_ref).insert(std::make_pair(center_left, p));
			std::stringstream str1;
			str1<<0<<":"<<center_ref<<":"<<center_left;
		//	std::cout << "center left in add al "<< center_left << " ref " << center_ref << std::endl;
			centerOnSequence.at(center_ref).insert(std::make_pair(center_left, str1.str()));
		}
		it->second.push_back(p);
	//	p.print();
	//	std::cout << p.get_al_ref1() << std::endl;
	//	std::cout << p.get_al_ref2() << std::endl;

	}
	const size_t finding_centers::find_right_on_seq(std::string & center,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & left , size_t & ref)const{
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
		unsigned int center_dir = atoi(center_parts.at(0).c_str());
		unsigned int center_ref = atoi(center_parts.at(1).c_str());
		unsigned int center_left = atoi(center_parts.at(2).c_str());
		std::stringstream str;
		str<<center_ref<<":"<<center_left;

		std::map<std::string, std::vector<pw_alignment> >::iterator it=alignments_in_a_cluster.find(str.str());
		assert(it!= alignments_in_a_cluster.end());
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
				std::cout << "1: "<<std::endl;
				p.print();
				return right_1;
			}
			if(sample2.str() == it->first && left_2 == left && ref2 == ref){
				std::cout << "2: "<<std::endl;
				p.print();

				return right_2;
			}
			if(sample2.str()==it->first && left_1 == left && ref1 == ref){
				std::cout << "3: "<<std::endl;
				p.print();

				return right_1;
			}
			if(sample1.str() == it->first && left_2 == left && ref2 == ref){
				std::cout << "4: "<<std::endl;
				p.print();

				return right_2;
			}
		}
		std::cout<< "shouldnt happen! "<<std::endl;
		return 0;
	}
	suffixTree::suffixTree(size_t & num_seq, finding_centers & cent, std::vector<std::vector<std::vector<int> > > & words):centers(cent){
		this->num_seq = num_seq;
		this->words =words;
		last_node_index=0;
		edges.insert(make_pair(last_node_index,root));
		last_node_index ++;
	}

	suffixTree::~suffixTree(){}

	void suffixTree::add_leaves(size_t & parent , std::vector<int> & new_suffix){
	//	double inf = std::numeric_limits<double>::infinity();
	//	std::cout<< "last node index "<< last_node_index << " new suffix size "<< new_suffix.size() << std::endl;
		edges_relations.insert(std::make_pair(parent, last_node_index));
		edges.insert(std::make_pair(last_node_index,new_suffix));
		last_node_index++;
	//	std::cout << "relation size "<<edges_relations.size()<< std::endl;
	}
	void suffixTree::add_internal_node(size_t & parent_index , std::vector<int> & new_suffix){//equal range for parent
		std::multimap<size_t,size_t>::iterator it2 = edges_relations.find(parent_index);
		assert(it2 != edges_relations.end());
	//	std::cout << last_node_index << " " << parent_index <<std::endl;
		std::vector<size_t> intermediate;
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			intermediate.push_back(it1->second);
		}
		edges_relations.erase(parent_index);
		it2 = edges_relations.find(parent_index);
		assert(it2 == edges_relations.end());
	//	std::cout << "child nodes were deleted"<<std::endl;
		edges_relations.insert(std::make_pair(parent_index, last_node_index));
	//	std::cout<< "intermediate size "<< intermediate.size()<<std::endl;
		for(size_t i =0; i < intermediate.size();i++){
			edges_relations.insert(std::make_pair(last_node_index, intermediate.at(i)));						
		//	std::cout << "test at "<< i << " is " << intermediate.at(i)<<std::endl;
		}
		edges.insert(std::make_pair(last_node_index,new_suffix));
		last_node_index ++;

	}
	void suffixTree::update_node(size_t & node_index, std::vector<int> & new_content){
	//	std::cout<< "node index "<< node_index << " new content "<< new_content << std::endl;
		std::map<size_t, std::vector<int> >::iterator it = edges.find(node_index);
		assert(it != edges.end());
		it->second = new_content;

	}
	void suffixTree::split(size_t & parent_index, std::vector<int> & parent , std::vector<int> & new_suffix, std::vector<int> & common_part, bool & make_a_tree){
		assert(common_part.size() <= parent.size());
		std::vector<int> parent_non_common;
	//	std::cout << "common part size "<< common_part.size() << " parent size "<< parent.size()<<std::endl;
		for(size_t i =common_part.size(); i < parent.size();i++){
			parent_non_common.push_back(parent.at(i));
		}
		if(parent_non_common.size() != 0){//common part shorter than parent's length
		//	std::cout << "if parent_non_common.size() != 0 "<<std::endl;
			std::multimap<size_t, size_t>::iterator it = edges_relations.find(parent_index);
			if(it == edges_relations.end()){
			//	std::cout<< "add leaves "<<std::endl;
				add_leaves(parent_index,parent_non_common); 
			}else{
			//	std::cout<< "add internal "<<std::endl;
				add_internal_node(parent_index,parent_non_common);
			}
		//	std::cout<< "parent index is "<< parent_index<< " " << new_suffix.size() << " >= " << common_part.size() << " " << parent_non_common <<std::endl;
			update_node(parent_index,common_part);//update its content
		//	std::cout<< "1parent index is "<< parent_index<< " " << new_suffix.size() << " >= " << common_part.size() <<std::endl;
			if(new_suffix.size()>= common_part.size()){
				std::vector<int> current_non_common;
				for(size_t i =common_part.size(); i < new_suffix.size();i++){
				//	std::cout << new_suffix.at(i) << " " ;
					current_non_common.push_back(new_suffix.at(i));
				}
			//	std::cout<< "current non common size is "<< current_non_common <<std::endl;
				if(current_non_common.size() != 0){
				//	std::cout<< "add leaves1 "<<std::endl;
					add_leaves(parent_index,current_non_common); 
				}
			}
			make_a_tree = false;
		}else{//common part is equal to parent's length - the new suffix is longer than the parent and fully covers it.
		//	std::cout << "i am here!" <<std::endl;
			assert(parent_non_common.size() == 0);
			if(common_part.size() < new_suffix.size()){ //TRICKY ONE!!!!
				std::vector<int> current_non_common;
				for(size_t i =common_part.size(); i < new_suffix.size();i++){
				//	std::cout<< new_suffix.at(i) << " ";
					current_non_common.push_back(new_suffix.at(i));
				}
			//	std::cout << " "<<std::endl;
				assert(current_non_common.size()!=0);
				std::multimap<size_t, size_t>::iterator it = edges_relations.find(parent_index);
				if(it != edges_relations.end()){//check on the childrens for a common part
				//	std::cout<< "has child node"<<std::endl;
					size_t new_parent_index = edges.size();
					find_next_parent(parent_index, current_non_common, common_part, new_parent_index, parent);
				//	std::cout << new_parent_index <<std::endl;
					if(new_parent_index != edges.size()){
						new_suffix = current_non_common;
						parent_index = new_parent_index;
					//	std::cout << "here "<< std::endl;
					}else{
						add_leaves(parent_index,current_non_common);
						make_a_tree = false;
					}
				}else{//add it as a leave
					assert(it == edges_relations.end());
					add_leaves(parent_index,current_non_common);
					make_a_tree = false;
				}

			}else{
				assert(common_part.size() == new_suffix.size());
				make_a_tree = false;
			}

		}
	}
	void suffixTree::make_suffix(std::vector<int> & word){
		suffixes.clear();
		for(size_t i =0; i < word.size();i++){
		//	std::cout << word.at(i) << " ";
			std::vector<int> suffix;
			for(size_t k =i; k < word.size();k++){
				suffix.push_back(word.at(k));
			}//For convenience add 1<<31 at the end of all the suffixes
			suffix.push_back(1<<31);
			suffixes.push_back(suffix);
		} 
//		std::cout << " "<<std::endl;	
	
	}
	void suffixTree::find_first_edges_after_root(size_t & deepest_parent, std::vector<int> & new_suffix){
		std::map<size_t, size_t >::iterator it = first_edges.find(new_suffix.at(0));
		if(it!= first_edges.end()){
			deepest_parent = it->second;
		}
	
	}
	void suffixTree::find_next_parent(size_t & parent_index, std::vector<int> & current_non_common, std::vector<int> & common_part, size_t & new_parent_index, std::vector<int> & parent){
		common_part.clear();
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
		//	std::cout << "child "<< it1->second << std::endl;
			std::map<size_t, std::vector<int> >::iterator edge=edges.find(it1->second);
			assert(edge != edges.end());
		//	std::cout << edge->second.at(0)<< " "<<current_non_common.at(0)<<std::endl;
			if(edge->second.at(0)==current_non_common.at(0)){
			//	std::cout<< "got it"<< std::endl;
				find_common_part(edge->second, current_non_common,common_part);
				new_parent_index = edge->first;
				parent = edge->second;
			//	std::cout << "new parent index "<< new_parent_index << std::endl;
				break;
			}
		}


	}
	void suffixTree::find_children(std::vector<int> & parent ,std::vector<std::vector<int> >& children){


	}
	void suffixTree::find_common_part(std::vector<int> & parent , std::vector<int> & current, std::vector<int> & common_part){
		size_t length = current.size();
		if(length > parent.size()) length = parent.size();
		for(size_t i =0; i < length; i++){
			if(parent.at(i)==current.at(i)){
				common_part.push_back(parent.at(i));
			}else{
				break;
			}
		}
	}
	void suffixTree::make_tree_for_a_seq(size_t & seq){
		std::cout << "at seq "<< seq <<std::endl;
		for(size_t i = 0; i < words.at(seq).size(); i ++){
			std::vector<int> this_word = words.at(seq).at(i);//current string of centers on a seq
			make_suffix(this_word);
	//		std::cout<< "suffixes size "<<suffixes.size()<<std::endl;
			for(size_t j =0; j< suffixes.size(); j++){
				size_t deepest_parent =edges.size();
				find_first_edges_after_root(deepest_parent, suffixes.at(j));
				if(deepest_parent == edges.size()){
				//	std::cout << "no first parent! "<<std::endl;
					first_edges.insert(std::make_pair(suffixes.at(j).at(0),last_node_index));
					size_t root_index =0;
					add_leaves(root_index,suffixes.at(j));
				}else{
					std::map<size_t, std::vector<int> >::iterator it = edges.find(deepest_parent);
					assert(it != edges.end());
				//	std::cout<< "has a parent"<<std::endl;
					std::vector<int> common_part;
					std::vector<int> parent = it->second;
					std::vector<int> current = suffixes.at(j);
					find_common_part(parent, current,common_part);
				//	for(size_t n =0; n < common_part.size();n++){
				//		std::cout<< common_part.at(n) << " ";
				//	}
				//	std::cout << " "<<std::endl;
					bool make_a_tree = true;
					size_t parent_index = deepest_parent;
					while(make_a_tree == true){
						split(parent_index,parent,current, common_part, make_a_tree);
					}
				}
			}
		}
	//	print_tree();

	}
	void suffixTree::make_tree(std::vector<int> & center_with_highest_gain, int & highest_index){	
	//	std::cout << "words size is "<< words.size() <<std::endl;
		for(size_t seq_id = 0; seq_id < num_seq; seq_id++){
			if(center_with_highest_gain.size() == 0){
				std::cout<< "seq id "<<seq_id <<std::endl;
			//	words.at(seq_id) = centers.get_centers(seq_id);
				std::cout << "its words size is "<< words.at(seq_id).size() <<std::endl;
			}else{
				update_centers(seq_id, center_with_highest_gain, highest_index);
			}
			make_tree_for_a_seq(seq_id);
		}
		count_branches();
	}
	void suffixTree::update_centers(size_t & seq_id, std::vector<int> & edge_with_highest_gain, int & highest_index){ 
		for(size_t i = 0; i < words.at(seq_id).size(); i++){
			std::cout << "before update at " << i <<std::endl;
			std::vector<int> centers = words.at(seq_id).at(i);
			std::cout << centers <<std::endl;
			std::cout << "rev: " << -1 * edge_with_highest_gain.back()<<std::endl;
			for(size_t k =0; k < centers.size();k++){
				if(centers.at(k) == edge_with_highest_gain.at(0)&& ((centers.size()-k)>= edge_with_highest_gain.size())){
					size_t first_common_index;
					std::vector<int> common;
					for(size_t i = 0; i < edge_with_highest_gain.size();i++){
						first_common_index = k;
						if(centers.at(i+k)== edge_with_highest_gain.at(i)){
							common.push_back( edge_with_highest_gain.at(i));
						}else break;
					}
					if(common.size() == edge_with_highest_gain.size()){
						std::vector<int> new_center;
						for(size_t m = 0; m < first_common_index ;m++){
							new_center.push_back(centers.at(m));
						}
						new_center.push_back(highest_index);
						for(size_t m = first_common_index + common.size(); m < centers.size();m++){
							new_center.push_back(centers.at(m));
						}
						centers=new_center;
					}
				}else if(centers.at(k) == -1 * edge_with_highest_gain.back() && ((centers.size()-k)>= edge_with_highest_gain.size())){
					std::cout << "here!"<<std::endl;
					std::vector<int> common;
					common.push_back(-1 * edge_with_highest_gain.back());
					for(size_t i = edge_with_highest_gain.size()-1; i> 0; i--){
						std::cout << " i " << i << std::endl;
						if(centers.at(k+edge_with_highest_gain.size()-i)== -1 * edge_with_highest_gain.at(i-1)){
							std::cout << "cent " << centers.at(k+edge_with_highest_gain.size()-i) << " " << -1 * edge_with_highest_gain.at(i-1) <<std::endl;
							common.push_back(-1 * edge_with_highest_gain.at(i-1));
						}else break;
					}
					if(common.size() == edge_with_highest_gain.size()){
						std::vector<int> new_center;
						for(size_t m = 0; m < k ;m++){
							new_center.push_back(centers.at(m));
						}
						new_center.push_back(-1*highest_index);
						for(size_t m = k + common.size(); m < centers.size();m++){
							new_center.push_back(centers.at(m));
						}
						centers=new_center;
					}
				}
			}
			words.at(seq_id).at(i) = centers;
			std::cout << "after update at "<< i <<std::endl;
			std::cout << centers << std::endl;
		}

	}
	void suffixTree::count_branches(){
		for(size_t seq =0; seq < num_seq; seq++){
			for(size_t w = 0; w< words.at(seq).size(); w++){
				std::vector<int> this_word = words.at(seq).at(w);//current string of centers on a seq
				make_suffix(this_word);
			//	std::cout<< "suffixes size "<<suffixes.size()<<std::endl;
				for(size_t j =0; j< suffixes.size(); j++){
				//	std::cout << "j "<<j <<std::endl;
					std::vector<size_t> branch;
					std::vector<int> current = suffixes.at(j);
					size_t deepest_parent =edges.size();
					find_first_edges_after_root(deepest_parent,current);
					assert(deepest_parent != edges.size());
					std::map<size_t, std::vector<int> >::iterator it=edges.find(deepest_parent);
					assert(it != edges.end());
					if(it->second.size()==current.size()){
						branch.push_back(deepest_parent);
					//	std::cout << "all the suffix is on one edge!"<<std::endl;
					}else{//Note that first parent can not be longer than current.
						std::vector<int> updated_current;
					//	std::cout<<"updated current: "<<std::endl;
						for(size_t k =it->second.size() ; k < current.size(); k++){
						//	std::cout << current.at(k)<< " ";
							updated_current.push_back(current.at(k));
						}
					//	std::cout << " "<<std::endl;
						size_t parent_index = deepest_parent;
						branch.push_back(parent_index);
						current = updated_current;
						while(current.size()>0){
							traverse_the_tree(parent_index,current);
						//	std::cout << "current size "<< current.size() <<std::endl;
							branch.push_back(parent_index);							
						}
					}
					std::vector<size_t> sub_branch;
				//	std::cout<<"branch size "<< branch.size()<<std::endl;
					for(size_t i = 0 ; i < branch.size(); i++){
					//	std::cout << "i "<< i << " " << branch.at(i) <<std::endl;
						sub_branch.push_back(branch.at(i));
						std::map<std::vector<size_t>,size_t>::iterator it1 = branches.find(sub_branch);
						if(it1 == branches.end()){
							branches.insert(make_pair(sub_branch,0));
							it1 = branches.find(sub_branch);
						}
						it1->second = it1->second +1;
					}

				}
			}
		}
	}
	void suffixTree::traverse_the_tree(size_t & parent_index, std::vector<int> & new_suffix){
	//	std::cout<< "parent index is "<< parent_index <<std::endl;
		std::multimap<size_t, size_t>::iterator check = edges_relations.find(parent_index);
		assert(check != edges_relations.end());
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator > it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t, size_t>::iterator it1 = it.first; it1 !=it.second;it1++){
			std::map<size_t,std::vector<int> >::iterator it2 = edges.find(it1->second);
		//	std::cout<< it1->second <<std::endl;
			assert(it2 != edges.end());
			if(it2->second.at(0)== new_suffix.at(0)){
				parent_index = it1->second;
			//	std::cout<< "new parent index " << parent_index<<std::endl;
				std::vector<int> updated_current;
			//	std::cout<< "updated:"<<std::endl;
				for(size_t i = it2->second.size(); i < new_suffix.size(); i++){
				//	std::cout<< new_suffix.at(i)<< " ";
					updated_current.push_back(new_suffix.at(i));
				}		
			//	std::cout << " "<<std::endl;
				new_suffix = updated_current;
				break;
			}
		}
		


	}
	const std::map<size_t, std::vector<int> > suffixTree::get_edges()const{
		return edges;
	}
	const std::map<std::vector<size_t>, size_t > suffixTree::get_branches()const{
		return branches;
	}

	const std::vector<std::vector<int> > suffixTree::get_current_centers(size_t & seqid)const{
		return words.at(seqid);
	}

	void suffixTree::print_tree(){
		for(std::map<size_t, std::vector<int> >::iterator it = edges.begin();it !=edges.end(); it++){
			std::cout<< " node is "<<it->first << " its edge content is ";
			for(size_t i = 0 ; i < it->second.size();i++){
				std::cout << it->second.at(i)<< " ";
			}
			std::cout<< " "<<std::endl;
			std::cout << "its children are "<<std::endl;
			std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it2 = edges_relations.equal_range(it->first);
			for(std::multimap<size_t,size_t>::iterator it1 = it2.first ; it1 != it2.second; it1++){
				std::cout << it1->second << " ";
			}
			std::cout << " " <<std::endl;
		}
	}



