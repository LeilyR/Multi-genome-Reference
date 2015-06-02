#include "graph.hpp"




// ---------------------------
// 		Vertex
// ---------------------------

// -------------------------------------------Functions
//Graph::Vertex::Vertex(){}
Graph::Vertex::Vertex(int id){
	idVertex = id;
	name = "";
	visited = false;
	costScore = std::numeric_limits<double>::infinity(); // initialize at none ? or -1 ?
	distance = std::numeric_limits<double>::infinity();
	startOnRead = -1;
	endOnRead = -1 ;
}
//Graph::Vertex::Vertex(Graph::Vertex const& v){}
//Graph::Vertex::~Vertex(){}
// -------------------------------------------Get
const int& Graph::Vertex::getIdVertex()const{return idVertex;}

//-------------
//Split a string
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
 }


//--------------------------------------------------------- Set and Get
void Graph::Vertex::setStartToStart(const int& v){
	if((std::find(startToStart.begin(),startToStart.end(),v)) == startToStart.end())
		startToStart.push_back(v);
}

void Graph::Vertex::setStartToEnd(const int& v){
	if((std::find(startToEnd.begin(),startToEnd.end(),v) == startToEnd.end()))
		startToEnd.push_back(v);
}

void Graph::Vertex::setEndToStart(const int& v){
	if((std::find(endToStart.begin(),endToStart.end(),v) == endToStart.end()))
		endToStart.push_back(v);
}

void Graph::Vertex::setEndToEnd(const int& v){
	if((std::find(endToEnd.begin(),endToEnd.end(),v) == endToEnd.end()))
		endToEnd.push_back(v);
}
void Graph::Vertex::setCostScore(const double& score){
	costScore = score;
}
void Graph::Vertex::setDistance(const double& d){
	distance = d;
}
void Graph::Vertex::setPrevious(const int& idPrevious){
	if((std::find(previous.begin(),previous.end(),idPrevious) == previous.end()))
		previous.push_back(idPrevious);
}
void Graph::Vertex::setStartOnRead(const int& start){
	startOnRead = start;
}
void Graph::Vertex::setEndOnRead(const int& end){
	endOnRead = end;
}
void Graph::Vertex::setName(const std::string& n){
	name =n;
}
const std::vector<int>& Graph::Vertex::getStartToStart(){return startToStart;}
const std::vector<int>& Graph::Vertex::getStartToEnd(){return startToEnd;}
const std::vector<int>& Graph::Vertex::getEndToStart(){return endToStart;}
const std::vector<int>& Graph::Vertex::getEndToEnd(){return endToEnd;}
const double& Graph::Vertex::getCostScore(){return costScore;}
const double& Graph::Vertex::getDistance(){return distance;}
const std::vector<int>& Graph::Vertex::getPrevious(){return previous;}
const int& Graph::Vertex::getStartOnRead(){	return startOnRead;}
// ---------------------------
// 		Graph
// ---------------------------

// -------------------------------------------Functions
//Graph::Graph(){}
//Graph::Graph(Graph const& g){}
//Graph::~Graph(){}
void Graph::setStart(const int&id, const int&start){
	std::map<int,Graph::Vertex>::iterator it = vertices.find(id);
		it->second.setStartOnRead(start);
}
void Graph::setEnd(const int&id, const int&end){
	std::map<int,Graph::Vertex>::iterator it = vertices.find(id);
		it->second.setEndOnRead(end);
}
const std::map<int,Graph::Vertex>& Graph::getVertices()const{return vertices;}
void Graph::setScore(const int& name, const double& score){
	std::map<int,Graph::Vertex>::iterator it = vertices.find(name);
	it->second.setCostScore(score);
}
const double& Graph::getScore(const int& id){
	std::map<int,Graph::Vertex>::iterator it = vertices.find(id);
	return it->second.getCostScore();
}
const Graph::Vertex& Graph::getVertex(const int& id){
	std::map<int,Graph::Vertex>::iterator it = vertices.find(id);
		return it->second;
}

void Graph::readDotFile(std::string file,std::map<std::string, size_t> longname2seqidx){
	std::cout << "readDotFile "<<std::endl;
	/*
	Achtung : Need spaces in the dot file
	like that : 1 -> 3 [headport=value] [tailport=value] ;
	For the moment node in .dot correspond to the index in the fasta file.
	So node are interger !
	Really important to respect the format :
	write from left to right to know the sens of the arrow
	and precise the sens of the node with + and -
	if + + = EndToStart
	if - + = StartToEnd
	if - - = StartToEnd
	if + - = EndToEnd
	*/
	std::string accession = "allSequencesNodes"; // TODO find a way to get the accession of the file
	std::ifstream dotFileIn(file.c_str());
	if(!dotFileIn){
		std::cerr << "Error : Cannot open " << file.c_str() << std::endl;
		exit(1);
	}
	else{
		std::string line;
		while(getline(dotFileIn,line)){
			if(line[0] == '/'){
				continue;
			}
			else{
				size_t pos = line.find("->");
				if(pos != std::string::npos){
					int origin,destination;
					std::string arrow, tmpOrigin,tmpDestination;
					std::istringstream iss(line);
					iss >> tmpOrigin >> arrow >>tmpDestination;

					std::string tmp1(tmpOrigin.end()-1,tmpOrigin.end());
					std::string tmp2(tmpDestination.end()-1,tmpDestination.end());

					std::string name1(tmpOrigin.begin(),tmpOrigin.end()-1);
					std::string name2(tmpDestination.begin(),tmpDestination.end()-1);
					std::string tmpAccName1 = accession + ":"+name1;
					map<string, size_t>::iterator findseq1 = longname2seqidx.find(tmpAccName1);
					if(findseq1==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence in dot File: " << tmpAccName1 << std::endl;
						exit(1);
							}
					origin = findseq1->second;
					std::string tmpAccName2 = accession + ":"+name2;
					map<string, size_t>::iterator findseq2 = longname2seqidx.find(tmpAccName2);
					if(findseq2==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence in dot File: " << tmpAccName2 << std::endl;
						exit(1);
					}
					destination = findseq2->second;
					//case 1 : startToStart
					if((tmp1 == "-") && (tmp2 == "+")){
						addEdge(1,origin,name1,destination,name2);
						continue;
					}
					//case 2 : startToEnd
					if((tmp1 == "-") && (tmp2 == "-")){
						addEdge(2,origin,name1,destination,name2);
						continue;
					}
					//case 3 : endToStart
					if((tmp1 == "+") && (tmp2 == "+")){
						addEdge(3,origin,name1,destination,name2);
						continue;
					}
					//case 4 : endToEnd
					if((tmp1 == "+") && (tmp2 == "-")){
						addEdge(4,origin,name1,destination,name2);
						continue;
					}
				}
			}
		}
	}
}

void Graph::addVertex(Graph::Vertex v){
	vertices.insert(std::pair<int,Graph::Vertex>(v.getIdVertex(), v));
}

void Graph::addEdge(int typeOfEdge,int id1, std::string name1,int id2,std::string name2){ // create nodes and add edges

	std::map<int,Graph::Vertex>::iterator it = vertices.find(id1);
	if(it == vertices.end()){ // If vertex does not exist
		Graph::Vertex v1= Graph::Vertex(id1);
		switch(typeOfEdge){
			case 1 :
				v1.setStartToStart(id2);
				break;
			case 2 :
				v1.setStartToEnd(id2);
				break;
			case 3 :
				v1.setEndToStart(id2);
				break;
			case 4 :
				v1.setEndToEnd(id2);
				break;
		}
		v1.setName(name1);
		addVertex(v1);

	}
	else{ // if exist, update
		switch(typeOfEdge){
			case 1 :
				it->second.setStartToStart(id2);
				break;
			case 2 :
				it->second.setStartToEnd(id2);
				break;
			case 3 :
				it->second.setEndToStart(id2);
				break;
			case 4 :
				it->second.setEndToEnd(id2);
				break;
		}
	}
	std::map<int,Graph::Vertex>::iterator it2 = vertices.find(id2);
	if(it2 == vertices.end()){// if v2 does not exist, creation
		Graph::Vertex v2 = Graph::Vertex(id2);
		v2.setPrevious(id1);
		v2.setName(name2);
		addVertex(v2);
	}
	else{
		it2->second.setPrevious(id1);
	}

}

std::ostream& operator<<(std::ostream &os, Graph & g)
{
	g.printGraph(os);
	return os;
}

void Graph::printGraph(std::ostream &os){
	os << "There is " << vertices.size() << " nodes" << std::endl;

	for (std::map<int,Graph::Vertex>::iterator it=vertices.begin(); it!=vertices.end(); ++it)
		(it->second).printVertex(os) ;
}

std::ostream& operator<<(std::ostream &os, Graph::Vertex const& v)
{
	v.printVertex(os);
	return os;
}

void Graph::Vertex::printVertex(std::ostream &os)const{
	os << " vertex : " << idVertex << " score "<< costScore << "pos["<<startOnRead<<"-"<<endOnRead<<"] name "<< name<<"\n";
	if(previous.size() != 0)
		os<<"previous node(s) " ;
	for(unsigned int i=0; i<previous.size(); ++i){
		os << previous[i] << " ";
	}
	os <<"\t";
	if(startToStart.size() != 0)
		os <<"startToStart: ";
	for(unsigned int i=0 ; i < startToStart.size(); ++i){
		 os << startToStart[i] <<" ";
	}
	os <<"\t";
	if(startToEnd.size() != 0)
		os <<"startToEnd: ";
	for(unsigned int i=0 ; i < startToEnd.size(); ++i){
		os << startToEnd[i] <<" ";
	}
	os <<"\t";
	if(endToStart.size() != 0)
		os <<"endToStart: ";
	for(unsigned int i=0 ; i < endToStart.size(); ++i){
		os << endToStart[i] <<" ";
	}
	os <<"\t";
	if(endToEnd.size() != 0)
		os <<"endToEnd: ";
	for(unsigned int i=0 ; i < endToEnd.size(); ++i){
		os << endToEnd[i] <<" ";
	}
	os << "\n";
}

bool operator== ( const Graph::Vertex &v1, const Graph::Vertex &v2){
        return v1.getIdVertex() == v2.getIdVertex();
}

void Graph::updateDistance(Graph::Vertex& v1, Graph::Vertex& v2, std::map<int,int>& previous){
	double weight = v2.getCostScore();
	std::cout << "updateDistance" <<  "v1.getDistance() "<<v1.getDistance()  << " weight "<< weight << " v2.getDistance()" << v2.getDistance()<<std::endl;
	if(v1.getDistance()+ weight <= v2.getDistance()){
		v2.setDistance(v1.getDistance()+weight);
		previous[v2.getIdVertex()] = v1.getIdVertex();
	}
}

std::vector<int> Graph::dijkstra(int start,int end){
	//Init
	std::vector<Graph::Vertex> notVisited;
	std::map<int,int> previous;
	std::vector<Graph::Vertex> neighbor;
    for( std::map<int,Graph::Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it ) {
	    	previous.insert(std::pair<int,int>(it->first,-1)); // pair <node, previous>
	}
	//Initialization depending to the start
	notVisited.push_back(vertices.find(start)->second);
	notVisited.begin()->setDistance(notVisited.begin()->getCostScore()); // TODO no need to initialize the start node because it is alone in the Queue at the beginning

	while (!notVisited.empty()){
		std::sort(notVisited.begin(),notVisited.end());
		if(notVisited.begin()->getIdVertex() == end){
			break;
		}
		else {
			neighbor.clear();
			std::vector<int> tmpneighbor = notVisited.begin()->getEndToStart(); // I need to sort neighbor to find the next one with the smaller cost
			for(std::vector<int>::iterator itt2 =tmpneighbor.begin(); itt2 != tmpneighbor.end(); ++itt2 ){
				neighbor.push_back(vertices.find(*itt2)->second);
				if (notVisited.begin()->getIdVertex()==start){ //first run, so we need to update distance of neighbor
					neighbor.back().setDistance(neighbor.back().getCostScore()); //set the distance to the costScore of the node
					previous[neighbor.back().getIdVertex()] = notVisited.begin()->getIdVertex();
				}
			}
			std::sort(neighbor.begin(), neighbor.end()); // sort neighbor by distance
			for(std::vector<Graph::Vertex>::iterator it2 = neighbor.begin(); it2 != neighbor.end(); ++it2 ){ // for each neighbor, updateDistance
				//update distance
				updateDistance(*notVisited.begin(),*it2,previous);
				// add node to the queue
				if ( std::find(notVisited.begin(), notVisited.end(),*it2)== notVisited.end() && it2->getIdVertex() != start) //TODO pb with circle
					notVisited.push_back(*it2);
				std::cout <<notVisited.begin()->getIdVertex() <<" --  "<<it2->getIdVertex()<< " distance " << it2->getDistance()<< std::endl;
			//Before erase neighbor, update vertices
			vertices.find(it2->getIdVertex())->second.setDistance(it2->getDistance());
			}
		}
		notVisited.erase(notVisited.begin());
	}

	if (previous[end] == -1){
		std::cout << "not reachable "<< "start :"<< start << " end "<< end<< std::endl;
		exit(1);

	}

	else{
	// Find Path :
		std::vector<int> path;
	int node = end;
	path.insert(path.begin(),node); // add end at the path
	while(node != start){
		path.insert(path.begin(),previous[node]);
		node = previous[node];
	}
	// print Path
	std::cout << "path ";
	std::copy(path.begin(), path.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout <<std::endl;
	return path;
	}

}
void Graph::initAllCostScore(all_data data,Graph& newGraph){
	if(VERBOSE) std::cout << "initAllCostScore" << std::endl;
	// Train the model on all data
	typedef mc_model use_model;
	use_model m(data);
	ofstream outs("encode",std::ofstream::binary);
	m.train(outs);
	for(size_t i=0 ; i<data.numAlignments(); ++i){
		double c1;
		double c2;
		double m1;
		double m2;
		m.cost_function(data.getAlignment(i),c1,c2,m1,m2,outs);
		if(VERBOSE) std::cout << "cost function with "<<data.getAlignment(i).getreference1()<<" c1 "<<c1<<" m1 "<<m1 <<" c2 "<<c2<<" m2 "<<m2 << std::endl;
		//std::cout << "data.getAlignment(i).getreference1() " << newGraph.getScore(data.getAlignment(i).getreference1()) << " "
		if(newGraph.getScore(data.getAlignment(i).getreference1()) > c1)//look if theformer score was better or replace it
			newGraph.setScore(data.getAlignment(i).getreference1(),c1);//c1 or c2 ?
	}
}


void Graph::parseData(all_data data, Graph& newGraph){
		//init
		size_t maxNumberOfBases = 2000;
		std::vector<pair<std::vector<int>,double>> allPath; // vector of path and distance


		vector<pw_alignment> vectorAl = data.getAlignments(); // TODO copy of vector of pw_alignments, not good !! But I need to sort the alignments
		std::sort(vectorAl.begin(),vectorAl.end());
		dnastring seqOfRead = data.getSequence(data.getAlignment(0).getreference2()); // save read sequence
		std::cout << "length of the read " << seqOfRead.length() << std::endl;

		for(size_t itP =0; itP < data.numAlignments(); ++itP){
			int rest = 0 ;
			//Verify that the 5' end of the read fit to the graph
			if(data.getAlignment(itP).getbegin2() < maxNumberOfBases && data.getAlignment(itP).getbegin2() !=0){
				std::string partOfRead, partOfPreviousNode;
				std::cout << "Look at the start of the read from position "<< data.getAlignment(itP).getbegin2()<< " node "<<data.getAlignment(itP).getreference1() <<std::endl; //call needleman option 2
				partOfRead = extractPartOfSeq(seqOfRead, 0, data.getAlignment(itP).getbegin2());
				Graph::Vertex it = vertices.find(data.getAlignment(itP).getreference1())->second;
				size_t lengthToLook = partOfRead.size()+ maxNumberOfBases; //TODO change that : stats ?

				for(unsigned int i = 0; i < it.getPrevious().size(); ++i){

					Graph::Vertex previousNode = vertices.find(it.getPrevious()[i])->second;
					std::cout << "previous node "<< previousNode.getIdVertex()<<std::endl;
					if(lengthToLook > data.get_seq_size(previousNode.getIdVertex())){
						rest = lengthToLook - data.get_seq_size(previousNode.getIdVertex()); // TODO rerun until rest =0
						lengthToLook = data.get_seq_size(previousNode.getIdVertex());
					}
					while(rest != 0 ){ // !! boucle infinie pour le moment !!
					int startOnPreviousNode =  data.get_seq_size(previousNode.getIdVertex())-lengthToLook;
					partOfPreviousNode = extractPartOfSeq(data.getSequence(previousNode.getIdVertex()), startOnPreviousNode,data.get_seq_size(previousNode.getIdVertex()) );
					std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,partOfPreviousNode,2);
					size_t idx1 = previousNode.getIdVertex();
					size_t idx2 = data.getAlignment(itP).getreference2();
					size_t incl_end1 = data.get_seq_size(previousNode.getIdVertex())-1;
					size_t incl_end2 =  data.getAlignment(itP).getbegin2()-1;
					pw_alignment al(resultNeedleman.first,resultNeedleman.second, startOnPreviousNode,0, incl_end1,incl_end2, idx1, idx2);
					data.add_pw_alignment(al);
				//*	newGraph.addEdge(3,previousNode.getIdVertex(),it.getIdVertex());
					int start = 0;
					newGraph.setStart(previousNode.getIdVertex(),start);
					newGraph.setEnd(previousNode.getIdVertex(),incl_end2);
					std::cout << "5' end " << previousNode.getIdVertex() << " " << it.getIdVertex() << std::endl;

					}
				}
			}
			//Verify that the 3' end of the read fit to the graph
			if(data.getAlignment(itP).getend2() > (seqOfRead.length() - maxNumberOfBases) && !(data.getAlignment(itP).getend2() == seqOfRead.length()-1)){
				std::string partOfRead, partOfNextNode;
				std::cout << "Look at the end of the read from position "<< seqOfRead.length()<<" "<< data.getAlignment(itP).getend2()<< " node "<<data.getAlignment(itP).getreference1() <<std::endl;
				partOfRead = extractPartOfSeq(seqOfRead, data.getAlignment(itP).getend2(), seqOfRead.length());
				Graph::Vertex it = vertices.find(data.getAlignment(itP).getreference1())->second;
				size_t lengthToLook = partOfRead.size() + maxNumberOfBases; //TODO change that : statistiques ?
				for(unsigned int i = 0; i < it.getEndToStart().size(); ++i){
					Graph::Vertex nextNode = vertices.find(it.getEndToStart()[i])->second;
					std::cout << "next node "<< nextNode.getIdVertex()<< " length "<< partOfRead.size()<<std::endl;
					if(lengthToLook > data.get_seq_size(nextNode.getIdVertex())){
						rest = lengthToLook - data.get_seq_size(nextNode.getIdVertex());
						lengthToLook = data.get_seq_size(nextNode.getIdVertex());
					}
					int endOnNextNode =  lengthToLook ;
					partOfNextNode = extractPartOfSeq(data.getSequence(nextNode.getIdVertex()), 0,endOnNextNode );
					std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,partOfNextNode,3);
					size_t idx1 = nextNode.getIdVertex();
					size_t idx2 = data.getAlignment(itP).getreference2();
					size_t incl_end1 = data.get_seq_size(nextNode.getIdVertex())-1;
					size_t incl_end2 =  data.getAlignment(itP).getbegin2()-1;
					pw_alignment al(resultNeedleman.first,resultNeedleman.second, endOnNextNode,0, incl_end1,incl_end2, idx1, idx2);
					data.add_pw_alignment(al);
				//*	newGraph.addEdge(3,it.getIdVertex(),nextNode.getIdVertex());
					int start = data.getAlignment(itP).getend2();
					int end = seqOfRead.length()-1;
					std::cout << "start " << start << " end "<< end << std::endl;
					newGraph.setStart(nextNode.getIdVertex(),start);
					newGraph.setEnd(nextNode.getIdVertex(),end);
					std::cout << "3' end addEdge " << it.getIdVertex() <<" " <<nextNode.getIdVertex()<< std::endl;
				}
			}
			for(size_t itPOther = 0 ; itPOther < data.numAlignments(); ++itPOther){ //TODO double loop not so good...
				if(vectorAl[itP].getreference1() == vectorAl[itPOther].getreference1() && vectorAl[itP].getbegin2() == vectorAl[itPOther].getbegin2())
					continue;
				if(vectorAl[itPOther].getbegin2() - vectorAl[itP].getend2() > maxNumberOfBases)
					break; // if we reach the max distance between two node, we stop looking at itP and go to hte next one

				//check if the pw_alignment is good !
				//pw_alignment alignment =data.getAlignment(itP);
				//pw_alignment *tmp;
				//tmp = &alignment;
				//alignment.print();
				//std::cout << "rrrrrrr " << data.alignment_fits_ref(tmp)<< std::endl;

				// verify if there is overlap or gap
				//vector<int> tmpPath;
				std::map<int,Graph::Vertex>::iterator it = vertices.find(vectorAl[itP].getreference1());
				//std::map<int,Graph::Vertex>::iterator itNextNode = vertices.find(vectorAl[itPOther].getreference1());
				std::map<int,int>previous;


				// Case 1 : Perfect
				if(vectorAl[itP].getend2()+1 == vectorAl[itPOther].getbegin2()){
					std::cout << "life is perfect  ! pp[itP].getend2() "<< vectorAl[itP].getend2()<<" blop "<<vectorAl[itP].getbegin2() <<" vectorAl[itPOther].getbegin2() " <<vectorAl[itPOther].getbegin2()<<" " <<vectorAl[itP].getreference1() << " " << vectorAl[itPOther].getreference1() << std::endl;
					//TODO use dijkstra or just verify that i and i+1 are directly link and update their distance ?

					if( std::find(it->second.getEndToStart().begin(), it->second.getEndToStart().end(), vectorAl[itPOther].getreference1())!= it->second.getEndToStart().end()){
						std::cout << " there is a direct link =) "<<std::endl;
					//*	newGraph.addEdge(3,it->second.getIdVertex(),vectorAl[itPOther].getreference1());
						int start = vectorAl[itP].getbegin2();
						int end = vectorAl[itP].getend2();
						newGraph.setStart(vectorAl[itP].getreference1(),start);
						newGraph.setEnd(vectorAl[itP].getreference1(),end);
						int start2 = vectorAl[itPOther].getbegin2();
						int end2 = vectorAl[itPOther].getend2();
						newGraph.setStart(vectorAl[itPOther].getreference1(),start2);
						newGraph.setEnd(vectorAl[itPOther].getreference1(),end2);
						std::cout << "add edge Perfect " << it->second.getIdVertex() <<" " <<vectorAl[itPOther].getreference1() << std::endl;
					}
				}
				//Case 2 : Gap
				else if( ((vectorAl[itP].getend2()+1 < vectorAl[itPOther].getbegin2()) && (vectorAl[itPOther].getbegin2()- vectorAl[itP].getend2() < maxNumberOfBases))){
					std::cout << " gap between this two nodes ! vectorAl[itP).getend2() "<< vectorAl[itP].getend2()<<" blop "<<vectorAl[itP].getbegin2() <<" vectorAl[itPOther].getbegin2() " <<vectorAl[itPOther].getbegin2()<< data.get_seq_name(vectorAl[itP].getreference1()) << " " << data.get_seq_name(vectorAl[itPOther].getreference1()) << std::endl;
					// take the gap sequence of the read, take the node after vectorAl[itP] that is different from vectorAl[itPOther]
					// call needleman option 1
					// verify that the new node and pp[itPOther] are reachable
					std::string gapSeq = extractPartOfSeq(seqOfRead, vectorAl[itP].getend2(), vectorAl[itPOther].getbegin2());
					std::cout << " gapSeq " << gapSeq.size() << std::endl;
						//look at sons of pp[itP]
					/*	for(){
							needleman(gapSeq,nextNode,1);
							m.cost_function()
						}
							//take min of cost_function
							path = dijkstra(minNode,pp[itPOther]);
							if(path ==nul)
								//take the next min node
								//if all node not reachable error !
						}
			*/
				}
				//Case 3 : Overlap
				else if (vectorAl[itP].getend2()+1 > vectorAl[itPOther].getbegin2() )// do nothing
					std::cout << " overlap ! vectorAl[itP).getend2() " << vectorAl[itP].getend2()<< " vectorAl[itP).getbegin2() "<<vectorAl[itP].getbegin2()<<" vectorAl[itPOther).getbegin2() "<< vectorAl[itPOther].getbegin2() << " "<<data.get_seq_name(vectorAl[itP].getreference1()) << " " << data.get_seq_name(vectorAl[itPOther].getreference1()) << std::endl;
					//allPath.push_back(std::make_pair(tmpPath,itNextNode->second.getDistance()));
			}
		}
	//update cost score of the data
	std::cout << "number of alignment at the end "<< data.numAlignments()<< std::endl;
	initAllCostScore(data, newGraph);
}

std::string Graph::extractPartOfSeq(dnastring seq, int start, int end){
	std::string subSeq = "";
	for(int i= start; i < end ; i++){
		subSeq += seq.at(i);
	}
	return subSeq;
}

void Graph::findFinalPath(){
	std::cout << " Final Path "<<std::endl;
	for(std::map<int,Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it){
		std::cout << getScore(it->first) << std::endl;
	}
	//for(std::vector<pair<std::vector<int>,double>>::iterator itallPath = allPath.begin(); itallPath != allPath.end(); ++itallPath){
			//std::copy(itallPath->first.begin(), itallPath->first.end(), std::ostream_iterator<int>(std::cout, " "));
			//std::cout <<"distance " <<itallPath->second<<std::endl;

		//}
}
