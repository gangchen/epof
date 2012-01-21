//graphLoader.cpp

#include "graphLoader.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
using namespace std;
using namespace std;

GraphLoader::GraphLoader()
{
	m_nArcs=0;
}
GraphLoader::GraphLoader(const char *szFileName)
{
	LoadGraphFromFile(szFileName);
}

void GraphLoader::LoadGraphFromFile(const char * szFileName)
{
	m_RawGraph.clear();
	ifstream InFile(szFileName);
	if(!InFile)
	{
		cout<<"Can't Open The Input File!"<<endl;
		return ;
	}
	GraphArc arc;
	for(;;)
	{
		if(!InFile)
			break;
		InFile>>arc.m_szFrom>>arc.m_szTo>>arc.m_fWeight;
		m_RawGraph.push_back(arc);
		InFile.ignore(100,'\n');
	}
	m_RawGraph.pop_back();
	InFile.close();
	m_nArcs=(int)m_RawGraph.size();
}
