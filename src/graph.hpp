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
		Vertex(int idVertex);
		//Vertex(Vertex const& v);
		//~Vertex();
		const int& getIdVertex()const;
		friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
		void printVertex(std::ostream &os)const;
		void setStartToStart(const int& id);
		void setStartToEnd(const int& id);
		void setEndToStart(const int& id);
		void setEndToEnd(const int& id);
		void setCostScore(const double& score);
		void setDistance(const double& distance);
		const std::vector<int>& getStartToStart();
		const std::vector<int>& getStartToEnd();
		const std::vector<int>& getEndToStart();
		const std::vector<int>& getEndToEnd();

		const double& getCostScore();
		const double& getDistance();
		bool operator < (const Graph::Vertex& v)const
				    {
				        return (distance < v.distance);
				    }

	private:
		int idVertex;
		std::string name;
		bool visited;
		std::string previous;
		int costScore;
		int distance;
		std::vector<int> startToStart;
		std::vector<int> startToEnd;
		std::vector<int> endToStart;
		std::vector<int> endToEnd;

	};

public:
	//Graph();
	//Graph(Graph const& g);
	//~Graph();
	void readDotFile(std::string dotFile); // Graph from a dot file -> Transform to constructor ?
	void addVertex(Vertex v);
	void addEdge(int typeOfEdge,int id1, int id2);
	const std::map<int,Vertex>& getVertices()const;
	void setScore(const int& id, const double& score);
	friend std::ostream& operator<<(std::ostream &os, Graph & g);
	friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
	void printGraph(std::ostream &os);
	//void updateDistance(Vertex v1, Vertex v2, std::map<std::string,int>& maxScore,std::map<std::string,std::string>& previous);
	void dijkstra(int start,int end);
	void readAlignment(std::string fastaFile, std::string samFile);
	friend bool operator== ( const Vertex &v1, const Vertex &v2);


private:
	std::map<int,Vertex> vertices;
};

#endif // GRAPH_HPP


