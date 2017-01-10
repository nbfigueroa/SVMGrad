/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 *
 *  main.cpp
 *
 *  Created on : Aug 14, 2012
 *  Author     : Ashwini Shukla
 *  Email      : ashwini.shukla@epfl.ch
 *  Website    : lasa.epfl.ch
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>

#include "ASVMLearning.h"

#include "signal.h"
using namespace std;
ASVM_SMO_Solver smo_solver;
#ifdef USE_NLOPT
ASVM_NLopt_Solver nlopt_solver;
#endif
//ASVM_ALGLIB_Solver solver;

//void exit_cleanly(int sig)
//{
//	cout<<"Force Exiting at current solution..."<<endl;
//	solver.force_stop();
//}

int main(int argc, char **argv)
{

	char datafile[1025]; strcpy(datafile, "data/data.txt");
	char outfile[1025]; strcpy(outfile, "models/learned_svm.txt");
	char paramfile[1025]; strcpy(paramfile, "./optparams.ini");
	char solvername[1025]; strcpy(solvername, "smo");
	double kernel_width = 1.0;
	double initial = 1e-8;
	int target_class = 0;

	if(argc == 1)
	{
		cout<<"==================== Usage ======================"<<endl;
		cout<<"./bin/train [OPTION1 VALUE1] [OPTION2 VALUE2] ..."<<endl<<endl;
		cout<<"OPTIONS  \t\t"<<"VALUES"<<endl;
		cout<<"--solver \t\t"<<"Which solver to use (smo/nlopt) [default: "<<solvername<<"]"<<endl;
		cout<<"--data   \t\t"<<"Data file name with extension [default: "<<datafile<<"]"<<endl;
		cout<<"--output \t\t"<<"Output file name with extension [default: "<<outfile<<"]"<<endl;
		cout<<"--param  \t\t"<<"Optimization parameter file [default: "<<paramfile<<"]"<<endl;
		cout<<"--sigma  \t\t"<<"Kernel Width [default: "<<kernel_width<<"]"<<endl;
		cout<<"--tclass \t\t"<<"Zero based target class [default: "<<target_class<<"]"<<endl;
		cout<<"--xinit \t\t"<<"Initial guess for optimization [default: 0]"<<endl;
		cout<<"        \t\t\to "<<"Positive number (p): All dimenstions are set equal to p"<<endl;
		cout<<"        \t\t\to "<<"-1: SVM Classifier solution"<<endl;
		cout<<"================================================="<<endl;
		return -6;
	}
	for(int i=0;i<argc;i++)
	{
		if(!strcmp(argv[i], "--data"))
			strcpy(datafile, argv[i+1]);
		else if(!strcmp(argv[i], "--output"))
			strcpy(outfile, argv[i+1]);
		else if(!strcmp(argv[i], "--param"))
			strcpy(paramfile, argv[i+1]);
		else if(!strcmp(argv[i], "--sigma"))
			kernel_width = atof(argv[i+1]);
		else if(!strcmp(argv[i], "--xinit"))
			initial = atof(argv[i+1]);
		else if(!strcmp(argv[i], "--tclass"))
			target_class = atoi(argv[i+1]);
		else if(!strcmp(argv[i], "--solver"))
			strcpy(solvername, argv[i+1]);
	}

	if(initial < 1e-8)
	{
		initial = 1e-8;
		cout<<"Overriding the initial value to 1e-8"<<endl;
	}
	asvmdata input;

	if(!input.loadFromFile(datafile))
		return -10;

	/* Example for filling data manually */
	asvmdata cpp_data;
	cpp_data.dim = input.dim;
	int num_targets = input.tar.size();
	int dim = input.dim;
	unsigned int i,j,k,l;
		for(i=0;i<num_targets;i++)
		{
			target t;
			t.targ = new double[dim];
			t.dim = dim;
			for(j=0;j<dim;j++)
				t.targ[j] = 0;

			int num_traj = input.tar[i].traj.size();

			for(j=0;j<num_traj;j++)
			{
				trajectory tr;

				int traj_len = input.tar[i].traj[j].nPoints;

				tr.dim = dim;
				tr.nPoints = traj_len;
				tr.y = new int[traj_len];
				tr.coords = new double*[traj_len];
				tr.vel = new double*[traj_len];
				for(k=0;k<traj_len;k++)
				{
					tr.coords[k] = new double[dim];
					tr.vel[k] = new double[dim];
				}

				for(k=0;k<dim;k++)
				{
					for(l=0;l<traj_len;l++)
						tr.coords[l][k] = input.tar[i].traj[j].coords[l][k];

				}
				for(l=0;l<traj_len;l++)
					tr.y[l] = i;

				for(l=0;l<dim;l++)
					t.targ[l] += tr.coords[traj_len-1][l];

				t.traj.push_back(tr);
			}

			for(l=0;l<dim;l++)
				t.targ[l] /= (double)num_traj;

			cpp_data.addTarget(t);
		}
		cpp_data.isOkay = true;
		/* End Example for filling data manually */


	input.setParams("rbf", kernel_width, initial);
	cpp_data.setParams("rbf", kernel_width, initial);

	int code=0;

	asvm svmobj;
	if(!strcmp(solvername, "nlopt"))
	{
#ifdef USE_NLOPT
		nlopt_solver.configure(paramfile);
		code = nlopt_solver.learn(input, target_class, &svmobj);
#else
		cout<<"NLOPT is not compiled!!"<<endl;
#endif
	}
	else if(!strcmp(solvername, "smo"))
	{
		smo_solver.configure(paramfile);
//		code = smo_solver.learn(input, target_class, &svmobj);
		code = smo_solver.learn(cpp_data, target_class, &svmobj);
	}


		svmobj.printinfo();
		svmobj.saveToFile(outfile);

	return code;

}



