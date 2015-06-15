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
void Graph::Vertex::setRtoF(const int& v){
	if((std::find(RtoF.begin(),RtoF.end(),v)) == RtoF.end())
		RtoF.push_back(v);
}

void Graph::Vertex::setRtoR(const int& v){
	if((std::find(RtoR.begin(),RtoR.end(),v) == RtoR.end()))
		RtoR.push_back(v);
}

void Graph::Vertex::setFtoF(const int& v){
	if((std::find(FtoF.begin(),FtoF.end(),v) == FtoF.end()))
		FtoF.push_back(v);
}

void Graph::Vertex::setFtoR(const int& v){
	if((std::find(FtoR.begin(),FtoR.end(),v) == FtoR.end()))
		FtoR.push_back(v);
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
const std::vector<int>& Graph::Vertex::getRtoF()const{return RtoF;}
const std::vector<int>& Graph::Vertex::getRtoR()const{return RtoR;}
const std::vector<int>& Graph::Vertex::getFtoF()const{return FtoF;}
const std::vector<int>& Graph::Vertex::getFtoR()const{return FtoR;}
const double& Graph::Vertex::getCostScore(){return costScore;}
const double& Graph::Vertex::getDistance(){return distance;}
const std::vector<int>& Graph::Vertex::getPrevious(){return previous;}
const int& Graph::Vertex::getStartOnRead(){	return startOnRead;}
const int& Graph::Vertex::getEndOnRead(){return endOnRead;}
const std::string& Graph::Vertex::getName(){return name;}

//-------------------------------
//		oriented_vertex
//------------------------------------------------------
//Graph::OrientedVertex::OrientedVertex();
Graph::OrientedVertex::OrientedVertex(const Graph::Vertex& v, bool b): vertex(v.getIdVertex()) { //not so good copy of object
//	Graph::Vertex vertex = getVertex(id);
	vertex = v;
	isForward = b;
}

bool Graph::OrientedVertex::vertexIsForward(pw_alignment p){
	bool b;
	size_t startOnRead, endOnRead;

	if(startOnRead < endOnRead)
		b = true;
	else
		b = false;
	return b;
}
void Graph::OrientedVertex::getSuccessors(std::vector<Graph::OrientedVertex>& oriented) const{
	if(isForward == true){
		Graph::Vertex tmp = vertex;
		vertex.getFtoF();

		for( uint it =0; it < vertex.getFtoF().size(); ++it){
			OrientedVertex orientedObject(vertex.getFtoF()[it],true);
			oriented.push_back(orientedObject);
		}
		for( uint it =0; it < vertex.getFtoR().size(); ++it){
			OrientedVertex orientedObject(vertex.getFtoR()[it],false);
			oriented.push_back(orientedObject);
		}
	}
	else{
		for( uint it =0; it < vertex.getRtoR().size(); ++it){
			OrientedVertex orientedObject(vertex.getRtoR()[it],false);
			oriented.push_back(orientedObject);
				}
		for( uint it =0; it < vertex.getRtoF().size(); ++it){
			OrientedVertex orientedObject(vertex.getRtoF()[it],false);
			oriented.push_back(orientedObject);
		}
	}
}


const Graph::Vertex& Graph::OrientedVertex::getVertex()const {return vertex;}
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
	if + + = FtoF
	if - + = RtoR
	if - - = RtoR
	if + - = FtoR
	*/
	std::string accession = "noAcc";
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
					std::map<std::string, size_t>::iterator findseq1 = longname2seqidx.find(tmpAccName1);
					if(findseq1==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence in dot File: " << tmpAccName1 << std::endl;
						exit(1);
							}
					origin = findseq1->second;
					std::string tmpAccName2 = accession + ":"+name2;
					std::map<std::string, size_t>::iterator findseq2 = longname2seqidx.find(tmpAccName2);
					if(findseq2==longname2seqidx.end()) {
						std::cerr << "Error: unknown sequence in dot File: " << tmpAccName2 << std::endl;
						exit(1);
					}
					destination = findseq2->second;
					//case 1 : RtoF
					if((tmp1 == "-") && (tmp2 == "+")){
						addEdge(1,origin,name1,destination,name2);
						continue;
					}
					//case 2 : RtoR
					if((tmp1 == "-") && (tmp2 == "-")){
						addEdge(2,origin,name1,destination,name2);
						continue;
					}
					//case 3 : FtoF
					if((tmp1 == "+") && (tmp2 == "+")){
						addEdge(3,origin,name1,destination,name2);
						continue;
					}
					//case 4 : FtoR
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
				v1.setRtoF(id2);
				break;
			case 2 :
				v1.setRtoR(id2);
				break;
			case 3 :
				v1.setFtoF(id2);
				break;
			case 4 :
				v1.setFtoR(id2);
				break;
		}
		v1.setName(name1);
		addVertex(v1);

	}
	else{ // if exist, update
		switch(typeOfEdge){
			case 1 :
				it->second.setRtoF(id2);
				break;
			case 2 :
				it->second.setRtoR(id2);
				break;
			case 3 :
				it->second.setFtoF(id2);
				break;
			case 4 :
				it->second.setFtoR(id2);
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
	os << " vertex : " << idVertex << " score "<< costScore << " pos["<<startOnRead<<"-"<<endOnRead<<"] name "<< name<<"\n";
	if(previous.size() != 0)
		os<<"previous node(s) " ;
	for(unsigned int i=0; i<previous.size(); ++i){
		os << previous[i] << " ";
	}
	os <<"\t";
	if(RtoF.size() != 0)
		os <<"RtoF: ";
	for(unsigned int i=0 ; i < RtoF.size(); ++i){
		 os << RtoF[i] <<" ";
	}
	os <<"\t";
	if(RtoR.size() != 0)
		os <<"RtoR: ";
	for(unsigned int i=0 ; i < RtoR.size(); ++i){
		os << RtoR[i] <<" ";
	}
	os <<"\t";
	if(FtoF.size() != 0)
		os <<"FtoF: ";
	for(unsigned int i=0 ; i < FtoF.size(); ++i){
		os << FtoF[i] <<" ";
	}
	os <<"\t";
	if(FtoR.size() != 0)
		os <<"FtoR: ";
	for(unsigned int i=0 ; i < FtoR.size(); ++i){
		os << FtoR[i] <<" ";
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
	// TODO Problem infinit run !!
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
			std::vector<int> tmpneighbor = notVisited.begin()->getFtoF(); // I need to sort neighbor to find the next one with the smaller cost
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
				if ( std::find(notVisited.begin(), notVisited.end(),*it2)== notVisited.end()) //TODO pb with circle  && it2->getIdVertex() != start
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
	std::ofstream outs("encode",std::ofstream::binary);
	m.train(outs);
	for(size_t i=0 ; i<data.numAlignments(); ++i){
		double c1;
		double c2;
		double m1;
		double m2;
		m.cost_function(data.getAlignment(i),c1,c2,m1,m2,outs);
		if(VERBOSE) std::cout << "cost function with "<<i << " "<<data.get_seq_name(data.getAlignment(i).getreference1())<<" c1 "<<c1<< std::endl;
		newGraph.setScore(i,c1);//c1 or c2 ?
	}
}

void Graph::lookPrevious(std::string partOfRead,Graph::Vertex it, all_data& data, Graph& newGraph){
	std::string partOfPreviousNode, newPartOfRead, restPartOfRead ;
	pw_alignment al;
	int start = 0;
	int end = partOfRead.size()-1;
	size_t lengthToLook = partOfRead.size(); //TODO change that : stats ?
	if(it.getPrevious().size() == 0 )
		std::cout << "no more Previous "<< lengthToLook-1<<std::endl;
	for(unsigned int i = 0; i < it.getPrevious().size(); ++i){
		Graph::Vertex previousNode = vertices.find(it.getPrevious()[i])->second;
		size_t lengthPreviousNode = data.get_seq_size(previousNode.getIdVertex());
		if(partOfRead.size() > lengthPreviousNode){
			std::cout << " Attention : node smaller than read " << std::endl;
			lengthToLook = lengthPreviousNode;
			end = it.getStartOnRead()-1;
			start = end - lengthToLook -1 ;
			restPartOfRead = extractPartOfSeq(partOfRead, 0,start-1); // take from the start of the read until the start of the part we look now
			partOfRead = extractPartOfSeq(partOfRead, start,end );
		}
		int startOnPreviousNode =  lengthPreviousNode - lengthToLook;
		size_t incl_end1 = startOnPreviousNode + lengthToLook -1;//-1 ?

		partOfPreviousNode = extractPartOfSeq(data.getSequence(previousNode.getIdVertex()), startOnPreviousNode,incl_end1 );
		std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,partOfPreviousNode,2);
		size_t idx1 = previousNode.getIdVertex();
		size_t idx2 = data.numSequences()-1;//the Read is the last sequence in fasta file

		size_t incl_end2 =  it.getStartOnRead() -1;// -1 ?
		pw_alignment al(resultNeedleman.first,resultNeedleman.second, startOnPreviousNode,start, incl_end1,incl_end2, idx1, idx2);
		// Does the alignment is ok ?
		pw_alignment *tmp;
		tmp = &al;
		data.alignment_fits_ref(tmp);
		data.add_pw_alignment(*tmp);
		int idOnNewGraph =  data.findIdAlignment(it.getName(),it.getStartOnRead(), it.getEndOnRead());
		previousNode.setStartOnRead(start);
		previousNode.setEndOnRead(end);
		newGraph.addEdge(3,data.numAlignments()-1,previousNode.getName(),idOnNewGraph,it.getName());
		newGraph.setStart(data.numAlignments()-1,start);
		newGraph.setEnd(data.numAlignments()-1,end);
		if (restPartOfRead.size() != 0 )
			lookPrevious(restPartOfRead,previousNode,data,newGraph);
		else{// Link node to startNode
			newGraph.addEdge(3,-1,"startNode",data.numAlignments()-1,previousNode.getName());
		}
	}
}
void Graph::lookNext(std::string partOfRead,Graph::Vertex it, all_data& data, Graph& newGraph){
	size_t lengthToLook = partOfRead.size();
	std::string newPartOfRead, restPartOfRead,partOfNextNode;
	int startToExtract, endToExtract;
	for(unsigned int i = 0; i < it.getFtoF().size(); ++i){
		Graph::Vertex nextNode = vertices.find(it.getFtoF()[i])->second;

		size_t lengthNextNode = data.get_seq_size(nextNode.getIdVertex());
		startToExtract = 0;
		if( lengthToLook > lengthNextNode ){
			lengthToLook = lengthNextNode;
			endToExtract = lengthToLook-1;
			restPartOfRead = extractPartOfSeq(partOfRead, endToExtract, partOfRead.size()-1);
			partOfRead = extractPartOfSeq(partOfRead, startToExtract,endToExtract);
		}
		else{
			endToExtract = lengthToLook -1;
		}
		partOfNextNode = extractPartOfSeq(data.getSequence(nextNode.getIdVertex()), startToExtract,endToExtract);
		std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,partOfNextNode,3);
		size_t idx1 = nextNode.getIdVertex();
		size_t idx2 = data.numSequences()-1; // the read is the last one in the fasta file
		size_t startOnRead = it.getEndOnRead()+1;
		size_t endOnRead =  startOnRead + lengthToLook-1;
		pw_alignment al(resultNeedleman.first,resultNeedleman.second, startToExtract,startOnRead,endToExtract,endOnRead, idx1, idx2);
		// Does the alignment is ok ?
		pw_alignment *tmp;
		tmp = &al;
	//	al.print();
		data.alignment_fits_ref(tmp);
		data.add_pw_alignment(*tmp);
		int idOnNewGraph =  data.findIdAlignment(it.getName(),it.getStartOnRead(), it.getEndOnRead());
		nextNode.setStartOnRead(startOnRead);
		nextNode.setEndOnRead(endOnRead);
		newGraph.addEdge(3,idOnNewGraph,it.getName(),data.numAlignments()-1,nextNode.getName());
		newGraph.setStart(data.numAlignments()-1,startOnRead);
		newGraph.setEnd(data.numAlignments()-1,endOnRead);
		if (restPartOfRead.size() != 0 )
			lookNext(restPartOfRead,nextNode,data,newGraph);
		else //link to endNode
			newGraph.addEdge(3,data.numAlignments()-1,nextNode.getName(),-2,"endNode");
	}
}
void Graph::DFS(Graph::Vertex v1, int target, int length, all_data data, std::vector<int> path){
	//find path with a corresponding length, for gap
	std::cout << "DFS "<< v1.getIdVertex() << " " << target << " " << length << std::endl;
		if(length <= 0 && std::find(v1.getFtoF().begin(),v1.getFtoF().end(),target) != v1.getFtoF().end()){
		std::cout << "path in DFS ";
//		path.insert(path.end(),target);
		std::copy(path.begin(), path.end(), std::ostream_iterator<int>(std::cout, " "));
		std::cout <<std::endl;
		return;
	}
	if(length <= 0 ){
		std::cout << "length <= 0 " << length << std::endl;
		return;
	}
	for(int i = 0 ; i < v1.getFtoF().size(); ++i){
//	for(std::vector<int>::iterator i = v1.getFtoF().begin(); i != v1.getFtoF().end(); ++i){

		size_t lengthNextNode = data.get_seq_size(v1.getFtoF()[i]);
		std::cout << " length "<<length << " " << lengthNextNode << std::endl;
	//	length = length - lengthNextNode;
		Graph::Vertex nextNode = vertices.find(v1.getFtoF()[i])->second;
		path.insert(path.begin(),v1.getFtoF()[i]);
		DFS(nextNode,target,length-lengthNextNode,data,path);
		path.pop_back();
	}


}
void Graph::lookGap(std::string partOfRead, Graph::Vertex v1, int idEnd, all_data data,Graph& newGraph){
	size_t endToExtract;
	size_t lengthToLook = partOfRead.size();
	std::string restPartOfRead, partOfNextNode;
	for(unsigned int i = 0; i < v1.getFtoF().size(); ++i){
		if(vertices.find(v1.getFtoF()[i])->first == idEnd ) // && gapSeq.size() != 0
			continue; //I want to look at other neighbor
		Graph::Vertex nextNode = vertices.find(v1.getFtoF()[i])->second;
		size_t lengthNextNode = data.get_seq_size(nextNode.getIdVertex());
		size_t startToExtract = 0;
		if( lengthToLook > lengthNextNode ){
			lengthToLook = lengthNextNode;
			endToExtract = lengthToLook-1;
			restPartOfRead = extractPartOfSeq(partOfRead, endToExtract, partOfRead.size()-1);
			partOfRead = extractPartOfSeq(partOfRead, startToExtract,endToExtract);
		}
		else{
			endToExtract = lengthToLook -1;
		}
		partOfNextNode = extractPartOfSeq(data.getSequence(nextNode.getIdVertex()), startToExtract,endToExtract);
		std::pair<std::string,std::string>resultNeedleman = runNeedleman(partOfRead,partOfNextNode,1); // Global or semi global ?
		size_t idx1 = nextNode.getIdVertex();
		size_t idx2 = data.numSequences()-1; // the read is the last one in the fasta file
		size_t startOnRead = v1.getEndOnRead()+1;
		size_t endOnRead =  startOnRead + lengthToLook-1;
		pw_alignment al(resultNeedleman.first,resultNeedleman.second, startToExtract,startOnRead,endToExtract,endOnRead, idx1, idx2);
		// Does the alignment is ok ?
		pw_alignment *tmp;
		tmp = &al;
		al.print();
		data.alignment_fits_ref(tmp);
		data.add_pw_alignment(*tmp);
		int idOnNewGraph =  data.findIdAlignment(v1.getName(),v1.getStartOnRead(), v1.getEndOnRead());
		nextNode.setStartOnRead(startOnRead);
		nextNode.setEndOnRead(endOnRead);
		newGraph.addEdge(3,idOnNewGraph,v1.getName(),data.numAlignments()-1,nextNode.getName());
		newGraph.setStart(data.numAlignments()-1,startOnRead);
		newGraph.setEnd(data.numAlignments()-1,endOnRead);
		if (restPartOfRead.size() != 0)
			lookGap(restPartOfRead, nextNode, idEnd, data,newGraph);

		else if (std::find(nextNode.getFtoF().begin(),nextNode.getFtoF().end(),idEnd) == nextNode.getFtoF().end())
			std::cout << "pb it dont reach the end of my gap, I should have not add this path to my newGraph "<< std::endl;
	}

}

void Graph::parseData(all_data& d, Graph& newGraph){

	//TODO build a graph to big, because if the node already exist with local alignment, it can rebuild it with global alignment, not a problem for score (score better with local alignment) but increase time of search in newGraph


	std::cout << "parse data "<< std::endl;
	std::string fastaFile = "/ebio/abt6_projects7/small_projects/mdubarry/Documents/graph/graph/output/all.fasta";
	std::string samFile = "/ebio/abt6_projects7/small_projects/mdubarry/Documents/graph/graph/output/referenceRead.sam";
	std::string dotFile = "/ebio/abt6_projects7/small_projects/mdubarry/Documents/data/graph.dot";
	all_data data = all_data();
	data.read_fasta_sam(fastaFile,samFile);

	Graph::Vertex v1= Graph::Vertex(-1);
	v1.setName("startNode");
	v1.setCostScore(0);
	newGraph.addVertex(v1);
	Graph::Vertex v2= Graph::Vertex(-2);
	v2.setName("endNode");
	v2.setCostScore(0);
	newGraph.addVertex(v2);
	std::map< std::string, size_t> longname2seqidx = data.getLongname2seqidx();
	readDotFile(dotFile,longname2seqidx);
	//init
	size_t maxNumberOfBases = 2000; // use a distance that change regarding read legnth : take 1/4 of read length ?

	vector<pw_alignment> vectorAl = data.getAlignments();
	std::sort(vectorAl.begin(),vectorAl.end(),sort_pw_alignment());
	dnastring seqOfRead= data.getSequence(data.getAlignment(0).getreference2()); // save read sequence TODO pointeur
	size_t numInitAl = data.numAlignments();
	for(size_t itP =0; itP < numInitAl; ++itP){
		std::map<int,Graph::Vertex>::iterator mapVertexItP = vertices.find(vectorAl[itP].getreference1());
		bool strand;
		if( vectorAl[itP].getbegin2() < vectorAl[itP].getend2())
			strand = true;
		else
			strand = false;
		OrientedVertex OV = OrientedVertex(mapVertexItP->second,strand);
		std::vector<OrientedVertex> oriented;
		std::cout << mapVertexItP->second;
		std::cout << OV.getSuccessors(oriented);
	//	OV.getSuccessors(oriented);
		exit(0);
		// Link node to startNode
		if(vectorAl[itP].getbegin2() ==0){
			int idOnNewGraph =  data.findIdAlignment(mapVertexItP->second.getName(),vectorAl[itP].getbegin2(), vectorAl[itP].getend2());
			newGraph.addEdge(3,-1,"startNode",idOnNewGraph,mapVertexItP->second.getName());
		}
		// Link node to endNode
		if(vectorAl[itP].getend2() == seqOfRead.length()-1){
			int idOnNewGraph =  data.findIdAlignment(mapVertexItP->second.getName(),vectorAl[itP].getbegin2(), vectorAl[itP].getend2());
			newGraph.addEdge(3,idOnNewGraph,mapVertexItP->second.getName(),-2,"endNode");
		}
		//Verify that the 5' end of the read fit to the graph
		if(vectorAl[itP].getbegin2() < maxNumberOfBases && vectorAl[itP].getbegin2() !=0){
			std::string partOfPreviousNode;
			std::cout << "Look at the start of the read from position "<< vectorAl[itP].getbegin2()<< " node "<<vectorAl[itP].getreference1() <<std::endl; //call needleman option 2
			std::string partOfRead = extractPartOfSeq(seqOfRead, 0, vectorAl[itP].getbegin2()-1);
			Graph::Vertex it = vertices.find(vectorAl[itP].getreference1())->second;
			it.setStartOnRead(vectorAl[itP].getbegin2());
			it.setEndOnRead(vectorAl[itP].getend2());
		//	lookPrevious(partOfRead,it,data,newGraph);
		}
		//Verify that the 3' end of the read fit to the graph
		size_t tmp ;
		tmp = seqOfRead.length() - maxNumberOfBases;
		if (vectorAl[itP].getend2() > tmp && vectorAl[itP].getend2() != seqOfRead.length() - 1 ){ //&& vectorAl[itP].getend2() != seqOfRead.length() - 1
			std::string partOfRead, partOfNextNode;
			std::cout << "Look at the end of the read from position "<< vectorAl[itP].getend2()<< " node "<<vectorAl[itP].getreference1() <<std::endl;
			partOfRead = extractPartOfSeq(seqOfRead, vectorAl[itP].getend2()+1, seqOfRead.length()-1);
			Graph::Vertex it = vertices.find(vectorAl[itP].getreference1())->second;
			it.setStartOnRead(vectorAl[itP].getbegin2());
			it.setEndOnRead(vectorAl[itP].getend2());
	//		lookNext(partOfRead,it,data,newGraph);
		}

		for(size_t itPOther = itP+1 ; itPOther < data.numAlignments(); ++itPOther){	//loop start at itP+1 because we dont need to look at node locate before and at the node[itP] with is the same
			if(vectorAl[itPOther].getbegin2() - vectorAl[itP].getend2() > maxNumberOfBases)
				break; // if we reach the max distance between two node, we stop looking at itP and go to the next one

			std::map<int,Graph::Vertex>::iterator mapVertexItPOther = vertices.find(vectorAl[itPOther].getreference1());
			//std::map<int,int>previous;

				// Case 1 : Perfect
			if(vectorAl[itP].getend2()+1 == vectorAl[itPOther].getbegin2()){
				std::cout << "life is perfect  ! pp[itP].getend2() "<< vectorAl[itP].getend2()<<" blop "<<vectorAl[itP].getbegin2() <<" vectorAl[itPOther].getbegin2() " <<vectorAl[itPOther].getbegin2()<<" " <<vectorAl[itP].getreference1() << " " << vectorAl[itPOther].getreference1() << std::endl;
				//TODO use dijkstra or just verify that i and i+1 are directly link ?

				if( std::find(mapVertexItP->second.getFtoF().begin(), mapVertexItP->second.getFtoF().end(), mapVertexItPOther->first)!= mapVertexItP->second.getFtoF().end()){
					int idNode1 =  data.findIdAlignment(mapVertexItP->second.getName(), vectorAl[itP].getbegin2(),vectorAl[itP].getend2() );
					int idNode2 =  data.findIdAlignment(mapVertexItPOther->second.getName(), vectorAl[itPOther].getbegin2(), vectorAl[itPOther].getend2());
					newGraph.addEdge(3,idNode1,mapVertexItP->second.getName(),idNode2,mapVertexItPOther->second.getName());
					int start = vectorAl[itP].getbegin2();
					int end = vectorAl[itP].getend2();
					newGraph.setStart(idNode1,start);
					newGraph.setEnd(idNode1,end);
					int start2 = vectorAl[itPOther].getbegin2();
					int end2 = vectorAl[itPOther].getend2();
					newGraph.setStart(idNode2,start2);
					newGraph.setEnd(idNode2,end2);
				}
			}
				//Case 2 : Gap
			else if( ((vectorAl[itP].getend2()+1 < vectorAl[itPOther].getbegin2()) && (vectorAl[itPOther].getbegin2()- vectorAl[itP].getend2() < maxNumberOfBases))){
				// take the gap sequence of the read, take the node after vectorAl[itP] that is different from vectorAl[itPOther]

				std::string gapSeq = extractPartOfSeq(seqOfRead, vectorAl[itP].getend2()+1, vectorAl[itPOther].getbegin2()-1);
				std::cout << " gapSeq " << gapSeq.size()<< " " <<  vectorAl[itP].getend2()+1 << " "<< vectorAl[itPOther].getbegin2()-1<< std::endl;
				std::vector<int> path;
				int lengthToLook = gapSeq.size();
				std::cout << " gap between this two alignment ! "<< mapVertexItP->first << " " << vectorAl[itPOther].getreference1()<< " length "<< lengthToLook << std::endl;
				vectorAl[itP].print();
				vectorAl[itPOther].print();
				DFS(mapVertexItP->second,  vectorAl[itPOther].getreference1(),lengthToLook, data,  path);
				//look at sons of pp[itP]
			//	Graph::Vertex v1 = vertices.find(vectorAl[itP].getreference1())->second;
			//	lookGap( gapSeq, v1, vectorAl[itPOther].getreference1(),data,newGraph);
			}
				//Case 3 : Overlap
			//	else if (vectorAl[itP].getend2()+1 > vectorAl[itPOther].getbegin2() )// do nothing
			//		std::cout << " overlap ! vectorAl[itP).getend2() " << vectorAl[itP].getend2()<< " vectorAl[itP).getbegin2() "<<vectorAl[itP].getbegin2()<<" vectorAl[itPOther).getbegin2() "<< vectorAl[itPOther].getbegin2() << " "<<data.get_seq_name(vectorAl[itP].getreference1()) << " " << data.get_seq_name(vectorAl[itPOther].getreference1()) << std::endl;
		}
	}
	//update cost score of the data
	std::cout << "number of alignment at the end "<< data.numAlignments()<< std::endl;
	initAllCostScore(data, newGraph);
//	newGraph.findFinalPath(seqOfRead.length());
//	newGraph.dijkstra(-1,-2);
}


std::string Graph::extractPartOfSeq(dnastring seq, int start, int end){
	if( start > end){
		std::cerr << "Problem to extract part of seq start " << start << " end "<< end<< std::endl;
		exit(1);
	}
	std::string subSeq = "";
	for(int i= start; i <= end ; ++i){
		subSeq += seq.at(i);
	}
	return subSeq;
}

/*
void Graph::findFinalPath(size_t lengthRead){ //TODO remove it !
	// TODO find path not between 0 and length but min and max
	std::cout << " Final Path "<<std::endl;
	for(std::map<int, Graph::Vertex>::iterator it = vertices.begin(); it!=vertices.end(); ++it) {
		if (it->second.getStartOnRead() == 0){
			for(std::map<int, Graph::Vertex>::iterator it2 = vertices.begin(); it2!=vertices.end(); ++it2){
				size_t end = it2->second.getEndOnRead();
				if( end == lengthRead-1){
					dijkstra(it->first,it2->first);
				}
			}
		}
	}
}


*/
