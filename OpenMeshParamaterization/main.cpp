#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "GameTimer.h"
#include "OMmodel.h"

void printUsage()
{
	cout << "usage:" << endl;
	cout << "OpenMeshParamaterization.exe infile outfile solver save" << endl;
	cout << "solver:  " << "1.SparseLU, 2.BICGSTAB	" << endl;
	cout << "save:    " << "1.as 2D Vertex Coordinates, 2.as 2D Texture Coordinates" << endl;
	cout << endl << endl;
}

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		printUsage();
		return -1;
	}
	GameTimer timer;
	timer.Reset();
	timer.Start();
	OMPmodel model;
	model.Param(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
	timer.Stop();
	cout << "The time used is: " << timer.TotalTime() <<endl;
	return 0;
}