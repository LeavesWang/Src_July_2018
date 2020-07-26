#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TChain.h"

using namespace std;

void MergRoot()
{
	string sSet[2]={"PS_270_382", "RS_270_382"};
	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	string sCombRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+to_string(runMin)+"--"+to_string(runMax)+".root";	
	ostringstream ssRun;
	string sRoot="";
	
	// TChain *chaiRoot = new TChain("tData", "combination of trees");
	// for(runNum=runMin; runNum<=runMax; runNum++)
		// for(int i=0; i<2; i++)
		// {
			// string sfSet="/home/kailong/ExpData/Jul2018/Src/runNum"+sSet[i]+".dat";
			// ifstream fSet(sfSet.c_str());
			// string sRead;
			// bool isFound=false;
			// while(getline(fSet, sRead))
				// if( sRead.find(to_string(runNum)) != string::npos )
				// {
					// isFound=true;
					// break;
				// }
			// if(isFound)
			// {
				// ssRun.str("");
				// ssRun<<setw(4)<<setfill('0')<<runNum;		
				// sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
				// printf("\n**********Now adding %s to %s!**********\n\n", sRoot.c_str(), sCombRoot.c_str());
				// chaiRoot->Add(sRoot.c_str());
				// break;
			// }
		// }
	// TFile *fCombRoot = new TFile(sCombRoot.c_str(), "RECREATE");	
	// chaiRoot->Merge(sCombRoot.c_str());
	// fCombRoot->Close();
	
	for(runNum=runMin; runNum<=runMax; runNum++)
		for(int i=0; i<2; i++)
		{
			string sfSet="/home/kailong/ExpData/Jul2018/Src/runNum"+sSet[i]+".dat";
			ifstream fSet(sfSet.c_str());
			string sRead;
			bool isFound=false;
			while(getline(fSet, sRead))
				if( sRead.find(to_string(runNum)) != string::npos )
				{
					isFound=true;
					break;
				}
			if(isFound)
			{
				ssRun.str("");
				ssRun<<setw(4)<<setfill('0')<<runNum;		
				sRoot+=" /home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root ";
				break;
			}
		}
	
	string sCmd1="rm -f "+sCombRoot;
	string sCmd2="hadd "+sCombRoot+sRoot;
	cout<<sCmd1<<endl;
	system(sCmd1.c_str());
	cout<<sCmd2<<endl;
	system(sCmd2.c_str());
}

#ifndef __CINT__
void StandaloneApplication(int argc, char** argv)
{
	MergRoot();
}

int main(int argc, char** argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	// app.Run();
	return 0;
}
#endif