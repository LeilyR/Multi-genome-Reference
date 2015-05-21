#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <typeinfo>
#include <iterator>

#include "data.hpp"
#include "pw_alignment.hpp"

// ----------------------------
//		Graph
//-----------------------------
class Graph{

private:
	// ----------------------------
	//		Vertex
	//-----------------------------
	class Vertex{

	public:
		//Vertex();
		Vertex(std::string name);
		//Vertex(Vertex const& v);
		//~Vertex();
		const std::string& getIdVertex()const;
		friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
		void printVertex(std::ostream &os)const;
		void setStartToStart(const std::string& name);
		void setStartToEnd(const std::string& name);
		void setEndToStart(const std::string& name);
		void setEndToEnd(const std::string& name);
		void setCostScore(const int& score); //TODO double
		void setDistance(const int& distance);
		const std::vector<std::string>& getStartToStart();
		const std::vector<std::string>& getStartToEnd();
		const std::vector<std::string>& getEndToStart();
		const std::vector<std::string>& getEndToEnd();

		const int& getCostScore();
		const int& getDistance();
		bool operator < (const Graph::Vertex& v)const
				    {
				        return (distance < v.distance);
				    }

	private:
		std::string idVertex;
		bool visited;
		std::string previous;
		int costScore;
		int distance;
		std::vector<std::string> startToStart;
		std::vector<std::string> startToEnd;
		std::vector<std::string> endToStart;
		std::vector<std::string> endToEnd;

	};

public:
	//Graph();
	void readDotFile(std::string dotFile); // Graph from a dot file -> Transform to constructor ?
	//Graph(Graph const& g);
	//~Graph();
	void addVertex(Vertex v);
	void addEdge(int typeOfEdge,std::string name1, std::string name2);
	const std::map<std::string,Vertex>& getVertices()const;
	void setScore(const std::string& name, const double& score);
	friend std::ostream& operator<<(std::ostream &os, Graph & g);
	friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
	void printGraph(std::ostream &os);
	void updateDistance(Vertex v1, Vertex v2, std::map<std::string,int>& maxScore,std::map<std::string,std::string>& previous);
	void dijkstra(std::string start,std::string end);
	void readAlignment(std::string fastaFile, std::string samFile);
	void buildNewGraph();
	friend bool operator== ( const Vertex &v1, const Vertex &v2);


private:
	std::map<std::string,Vertex> vertices;
};

#endif // GRAPH_HPP


