#include"stdafx.h"
#include"preprocess.h"

void generateMultiGraph() {
	sparserow graph, trigraph;
	string dir = "C:\\Users\\dul262\\Desktop\\DATA\\SNAPData\\data\\Orkut\\";
	string path = dir + "Orkut.txt";
	createGraph(path, graph);

	cout << "Create Graph complete!" << endl;

	//string tripath = dir + "trigraph.txt";
	//createTriGraph(graph, trigraph, tripath);

	//string fourpath = dir + "fourgraph.txt";
	//createFourGraph(graph, trigraph, fourpath);

	//string fivepath = dir + "fivegraph.txt";
	//createFiveGraph(graph, trigraph, fivepath);
}

