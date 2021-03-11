#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <unistd.h>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TH1F.h"
#include "THStack.h"

using namespace std;

struct StrtMesytec
{
	int modID;
	int data[32];
	int modRes;
	int modEC_TS;
};

void Root2CalMCP()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetCanvasDefH(1080);
	gStyle->SetCanvasDefW(1920);

	const int QdcUp = 3840;
	// const int QdcMcpLow[8]={1000, 760, 810, 850, 790, 790, 780, 790}; //old pedestal values
	const int QdcMcpLow[8] = {759, 765, 802, 791, 752, 764, 765, 760}; //new pedestal values
	// const double CalQDC[8][2]={{0,1}, {0,1}, {0,1}, {0,1}, {982.999,20.2568}, {824.977,17.6707}, {-21.8547,12.8791}, {766.523,19.5068}};
	// const double mcpGainMat[8][2]={{0,1}, {0,1}, {0,1}, {0,1}, {155.116,21.0876}, {315.652,18.1773}, {27.9926,12.8659}, {217.036,19.7808}}; //from fitting #270--#385
	// const double mcpGainMat[8][2]={{0,1}, {0,1}, {0,1}, {0,1}, {223.877,21.5615}, {381.772,19.2887}, {56.4084,12.5894}, {237.447,23.2303}}; //from fitting #383--#385
	string sRoot, sCalMcp;
	StrtMesytec mqdcMCP;
	int run;
	double xMCP[3], yMCP[3];
	int nMCP[3];
	double calQdcMcp[12];

	long long iEnt, nEnt;
	int i, j;
	int iLow[4] = {6, 5, 4, 7};

	int runMin, runMax, runNum;
	cout << "Input minimum and maximum numbers of run: ";
	cin >> runMin >> runMax;

	sCalMcp = "/home/kailong/ExpData/Jul2018/CalMCP/calMCP-run-" + to_string(runMin) + "--" + to_string(runMax) + ".root";
	TFile *fCalMcp = new TFile(sCalMcp.c_str(), "RECREATE");
	TTree *tCalMcp = new TTree("tCalMcp", "tree for calibrating MCP");
	tCalMcp->Branch("run", &run, "run/I");
	tCalMcp->Branch("xMCP", xMCP, "xMCP[3]/D");
	tCalMcp->Branch("yMCP", yMCP, "yMCP[3]/D");

	ostringstream ssRun;
	for (runNum = runMin; runNum <= runMax; runNum++)
	{
		ssRun.str("");
		ssRun << setw(4) << setfill('0') << runNum;

		sRoot = "/home/kailong/ExpData/Jul2018/RootData/run-" + ssRun.str() + "-00.root";
		printf("\n**********Now converting %s to %s!**********\n\n", sRoot.c_str(), sCalMcp.c_str());

		TFile *fRoot = new TFile(sRoot.c_str());
		if (fRoot->IsZombie())
		{
			cout << "Error in opening " << sRoot << "!\n";
			continue;
		}

		TTree *tData;
		fRoot->GetObject("tData", tData);
		if (!tData)
		{
			cout << "Error read the tree of tData!\n";
			continue;
		}

		double mcpGainMat[8][2];
		for (i = 0; i < 8; i++)
		{
			mcpGainMat[i][0] = 0;
			mcpGainMat[i][1] = 1;
		}

		tData->SetEstimate(-1);
		nEnt = tData->GetEntries();
		double *errCh = new double[nEnt];
		for (iEnt = 0; iEnt < nEnt; iEnt++)
			errCh[iEnt] = 0.5;
		string strCanv = "gainMat_" + to_string(runNum);
		TCanvas *cGain = new TCanvas(strCanv.c_str(), strCanv.c_str());
		cGain->Divide(2, 2);
		bool goodRun = true;
		// TH1F *res[4];
		// TH1F *histMcp[8];
		TGraphErrors *gr[4];
		for (i = 0; i < 4; i++)
		{
			j = iLow[i];
			cGain->cd(i + 1);
			string sCut = "mqdcMCP.data[" + to_string(i) + "]>" + to_string(QdcMcpLow[i]) + "&&mqdcMCP.data[" + to_string(i) + "]<3840&&mqdcMCP.data[" + to_string(j) + "]>" + to_string(QdcMcpLow[j]) + "&&mqdcMCP.data[" + to_string(j) + "]<3840";

			string sDraw = "mqdcMCP.data[" + to_string(i) + "]-" + to_string(QdcMcpLow[i]) + ":mqdcMCP.data[" + to_string(j) + "]-" + to_string(QdcMcpLow[j]) + ">>h" + to_string(i) + "_" + to_string(j);

			printf("Now drawing {%s} of {%s} when {%s}\n\n", sDraw.c_str(), sRoot.c_str(), sCut.c_str());

			long long nData = tData->Draw(sDraw.c_str(), sCut.c_str(), "goff");
			if (nData < 2)
			{
				goodRun = false;
				break;
			}
			double *highGain = tData->GetV1();
			double *lowGain = tData->GetV2();
			gr[i] = new TGraphErrors(nData, lowGain, highGain, errCh, errCh);
			gr[i]->Draw("AP");
			gr[i]->SetTitle(("high" + to_string(i) + "_vs_low" + to_string(iLow[i])).c_str());
			TFitResultPtr fitRes = gr[i]->Fit("pol1", "MS");
			mcpGainMat[j][0] = fitRes->Parameter(0);
			mcpGainMat[j][1] = fitRes->Parameter(1);
			// mcpGainMat[j][0]=0;
			// mcpGainMat[j][1]=1;
			// TF1 *fitFcn=gr->GetFunction("pol1");
			// histMcp[i]=new TH1F(("hist_mqdc_mcp_"+to_string(i)).c_str(), ("hist_mqdc_mcp_"+to_string(i)).c_str(), 4000,0,4000);
			// histMcp[j]=new TH1F(("hist_mqdc_mcp_"+to_string(j)).c_str(), ("hist_mqdc_mcp_"+to_string(j)).c_str(), 2000,0,40000);

			// cGain->cd(i+5);
			// for(iEnt=0; iEnt<nEnt; iEnt++)
			// {
			// TRandom3 r1(0);
			// if(mqdcMCP.data[i]>QdcMcpLow[i]&&mqdcMCP.data[i]<QdcUp)
			// histMcp[i]->Fill(mqdcMCP.data[i]-QdcMcpLow[i]);
			// if(mqdcMCP.data[j]>QdcMcpLow[j]&&mqdcMCP.data[j]<QdcUp)
			// histMcp[j]->Fill(mcpGainMat[j][0]+mcpGainMat[j][1]*(mqdcMCP.data[j]-QdcMcpLow[j]+r1.Uniform(-0.5,0.5)));
			// }
			// histMcp[j]->SetLineColor(kRed);
			// THStack *hs = new THStack("hs","hs");
			// hs->Add(histMcp[i]);
			// hs->Add(histMcp[j]);
			// hs->Draw();
			// res[i]->Fill(highGain[iEnt]-fitFcn->Eval(lowGain[iEnt]));
			// res[i]->Draw(NOSTACK");
		}
		delete[] errCh;
		if (!goodRun)
			continue;
		cGain->SaveAs(("/home/kailong/ExpData/Jul2018/Graphs/Charts/" + strCanv + ".png").c_str());

		memset(&mqdcMCP, 0, sizeof(mqdcMCP));

		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);

		// nEnt=tData->GetEntries();
		for (iEnt = 0; iEnt < nEnt; iEnt++)
		{
			tData->GetEntry(iEnt);
			TRandom3 r(0);
			run = runNum;
			memset(xMCP, 0, sizeof(xMCP));
			memset(yMCP, 0, sizeof(yMCP));
			memset(nMCP, 0, sizeof(nMCP));
			memset(calQdcMcp, 0, sizeof(calQdcMcp));
			for (i = 0; i < 8; i++)
				if (mqdcMCP.data[i] > QdcMcpLow[i] && mqdcMCP.data[i] < QdcUp)
				{
					j = i / 4;
					nMCP[j]++;
					calQdcMcp[i] = mcpGainMat[i][0] + mcpGainMat[i][1] * (mqdcMCP.data[i] - QdcMcpLow[i] + r.Uniform(-0.5, 0.5));
				}
			for (i = 0; i < 4; i++)
			{
				if (mqdcMCP.data[i] > QdcMcpLow[i] && mqdcMCP.data[i] < QdcUp)
				{
					nMCP[2]++;
					calQdcMcp[i + 8] = mqdcMCP.data[i] - QdcMcpLow[i] + r.Uniform(-0.5, 0.5);
				}
				j = iLow[i];
				if (mqdcMCP.data[i] >= QdcUp && mqdcMCP.data[j] > QdcMcpLow[j] && mqdcMCP.data[j] < QdcUp && mcpGainMat[j][0] + mcpGainMat[j][1] * (mqdcMCP.data[j] - QdcMcpLow[j]) > QdcUp)
				// if(mqdcMCP.data[i]>=QdcUp&&mqdcMCP.data[j]>QdcMcpLow[j]&&mqdcMCP.data[j]<QdcUp)
				{
					nMCP[2]++;
					calQdcMcp[i + 8] = mcpGainMat[j][0] + mcpGainMat[j][1] * (mqdcMCP.data[j] - QdcMcpLow[j] + r.Uniform(-0.5, 0.5));
				}
			}

			if (nMCP[0] == 4)
			{
				xMCP[0] = (calQdcMcp[0] + calQdcMcp[3] - calQdcMcp[1] - calQdcMcp[2]) / (calQdcMcp[0] + calQdcMcp[1] + calQdcMcp[2] + calQdcMcp[3]);

				yMCP[0] = (calQdcMcp[1] + calQdcMcp[3] - calQdcMcp[0] - calQdcMcp[2]) / (calQdcMcp[0] + calQdcMcp[1] + calQdcMcp[2] + calQdcMcp[3]);
			}

			if (nMCP[1] == 4)
			{
				xMCP[1] = (calQdcMcp[6] + calQdcMcp[7] - calQdcMcp[5] - calQdcMcp[4]) / (calQdcMcp[6] + calQdcMcp[5] + calQdcMcp[4] + calQdcMcp[7]);

				yMCP[1] = (calQdcMcp[5] + calQdcMcp[7] - calQdcMcp[6] - calQdcMcp[4]) / (calQdcMcp[6] + calQdcMcp[5] + calQdcMcp[4] + calQdcMcp[7]);
			}

			if (nMCP[2] == 4)
			{
				xMCP[2] = (calQdcMcp[8] + calQdcMcp[11] - calQdcMcp[9] - calQdcMcp[10]) / (calQdcMcp[8] + calQdcMcp[9] + calQdcMcp[10] + calQdcMcp[11]);

				yMCP[2] = (calQdcMcp[9] + calQdcMcp[11] - calQdcMcp[8] - calQdcMcp[10]) / (calQdcMcp[8] + calQdcMcp[9] + calQdcMcp[10] + calQdcMcp[11]);
			}

			if (nMCP[0] == 4 || nMCP[1] == 4 || nMCP[2] == 4)
				tCalMcp->Fill();
		} //end of whole tree
		fRoot->Close();
	} //end of runs
	fCalMcp->cd();
	tCalMcp->Write();
	fCalMcp->Close();
	cout << "\n**********The program has been finished, press Ctr-C to quit!**********\n";
} //end of whole function

#ifndef __CINT__
void StandaloneApplication(int argc, char **argv)
{
	Root2CalMCP();
}

int main(int argc, char **argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	app.Run();
	return 0;
}
// int main()
// {
// 	Root2CalMCP();
// 	return 0;
// }
#endif