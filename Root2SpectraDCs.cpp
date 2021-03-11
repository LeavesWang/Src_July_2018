#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <cmath>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TDirectory.h"

using namespace std;

struct StrtMesytec
{
	int modID;
	int data[32];
	int modRes;
	int modEC_TS;
};

struct StrtS800
{
	int tS;
	int eC;
	int trig;				// =1: Coincidence (usual case, S800 scintillator provides); =16: Secondary (our ToF or MCP provides)
	int tof[16];			//S800 ToF packet
	int scin[2][2];			//S800 scintillator: first [2]: two PMTs, second [2]: [0] energy; [1] time
	int ionCham[16];		//S800 ion chamber; [16] means 16 segaments
	int crdcCath[2][5][64]; //[2]: two CRDC; [5]: [0]---sample; [1]-->[4]---energy
	int crdcAnode[2][2];	// [0]: energy; [1]: time
	int mesyTDC[16];		//[0]: time from S800; [2]: time from A1900
};

void Root2SpectraDCs()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	if (access("../Graphs/SpectraDcs", F_OK) != 0)
		system("mkdir ../Graphs/SpectraDcs");
	string sDirADC = "../Graphs/SpectraDcs/MADC32/";
	if (access(sDirADC.c_str(), F_OK) != 0)
		system(("mkdir " + sDirADC).c_str());
	string sDirTDC = "../Graphs/SpectraDcs/MTDC32/";
	if (access(sDirTDC.c_str(), F_OK) != 0)
		system(("mkdir " + sDirTDC).c_str());
	string sDirTofQDC = "../Graphs/SpectraDcs/V792QDCTOF/";
	if (access(sDirTofQDC.c_str(), F_OK) != 0)
		system(("mkdir " + sDirTofQDC).c_str());
	string sDirMcpQDC = "../Graphs/SpectraDcs/MQDC32MCP/";
	if (access(sDirMcpQDC.c_str(), F_OK) != 0)
		system(("mkdir " + sDirMcpQDC).c_str());
	string sDirIC = "../Graphs/SpectraDcs/S800IC/";
	if (access(sDirIC.c_str(), F_OK) != 0)
		system(("mkdir " + sDirIC).c_str());
	string sDirCRDC = "../Graphs/SpectraDcs/S800CRDC/";
	if (access(sDirCRDC.c_str(), F_OK) != 0)
		system(("mkdir " + sDirCRDC).c_str());

	const int IdAdc[12] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24};
	const int IdTdc[8] = {1, 3, 5, 7, 9, 11, 13, 15};
	const int IdQdcTof[8] = {1, 3, 5, 6, 8, 10, 12, 15};
	const int IdQdcMcp[8] = {1, 3, 5, 7, 9, 11, 13, 15};

	string sCanvAdc, sCanvTdc, sCanvQdcTof, sCanvQdcMcp, sCanvIC, sCanvCrdcAnod, sCanvCrdcCath;
	string sAnodeCRDC[2] = {"Energy", "Time"};

	int i = 0, j = 0, k = 0;

	StrtMesytec madc, mtdc, cqdcTOF, mqdcMCP;
	StrtS800 s800;

	int runMin, runMax, runNum;
	cout << "Input minimum and maximum numbers of run: ";
	cin >> runMin >> runMax;

	long long nEnt = 0, iEnt = 0;

	string strRun, sRoot, sDirRun;
	string sDraw, sCut, sCanv;

	ostringstream ssRun;
	for (runNum = runMin; runNum <= runMax; runNum++)
	{
		ssRun.str("");
		ssRun << setw(4) << setfill('0') << runNum;
		strRun = ssRun.str();
		sRoot = "../RootData/run-" + strRun + "-00.root";

		printf("\n**********Now plotting spectra of DCs of %s!**********\n", sRoot.c_str());

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
		sDirRun = "../Graphs/SpectraDcs/run" + strRun + "/";
		if (access(sDirRun.c_str(), F_OK) != 0)
			system(("mkdir " + sDirRun).c_str());
		memset(&madc, 0, sizeof(madc));
		memset(&mtdc, 0, sizeof(mtdc));
		memset(&cqdcTOF, 0, sizeof(cqdcTOF));
		memset(&mqdcMCP, 0, sizeof(mqdcMCP));
		memset(&s800, 0, sizeof(s800));

		tData->SetBranchAddress("madc", &madc);
		tData->SetBranchAddress("mtdc", &mtdc);
		tData->SetBranchAddress("cqdcTOF", &cqdcTOF);
		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);
		tData->SetBranchAddress("s800", &s800);

		TH1F *hAdc[12];
		for (i = 0; i < 12; i++)
			hAdc[i] = new TH1F(("hAdc" + to_string(i)).c_str(), ("Spectrum of ADC-" + to_string(IdAdc[i]) + ";Channel;Count/ch.").c_str(), 7670, 10, 7680);

		TH1F *hTdc[8];
		for (i = 0; i < 8; i++)
			hTdc[i] = new TH1F(("hTdc" + to_string(i)).c_str(), ("Spectrum of TDC-" + to_string(IdTdc[i]) + ";Channel;Count/ch.").c_str(), 65516, 10, 65526);

		TH1F *hQdcTof[8];
		for (i = 0; i < 8; i++)
			hQdcTof[i] = new TH1F(("hQdcTof" + to_string(i)).c_str(), ("Spectrum of ToF's QDC-" + to_string(IdQdcTof[i]) + ";Channel;Count/ch.").c_str(), 4076, 10, 4086);

		TH1F *hQdcMcp[8];
		for (i = 0; i < 8; i++)
			hQdcMcp[i] = new TH1F(("hQdcMcp" + to_string(i)).c_str(), ("Spectrum of MCP's QDC-" + to_string(IdQdcMcp[i]) + ";Channel;Count/ch.").c_str(), 3830, 10, 3840);

		TH1F *hAdcIC[16];
		for (i = 0; i < 16; i++)
			hAdcIC[i] = new TH1F(("hAdcIC" + to_string(i)).c_str(), ("Spectrum of S800 Ion Chamber's ADC-" + to_string(i) + ";Channel;Count/ch.").c_str(), 4076, 10, 4086);

		TH1F *hAdcCrdcAnod[2][2];
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				hAdcCrdcAnod[i][j] = new TH1F(("hAdcCrdcAnode_" + to_string(i) + "_" + to_string(j)).c_str(), ("Spectrum of S800 CRDC" + to_string(i + 1) + "-Anode's " + sAnodeCRDC[j] + ";Channel;Count/ch.").c_str(), 4076, 10, 4086);

		TH2F *hCrdcCath[2];
		for (i = 0; i < 2; i++)
			hCrdcCath[i] = new TH2F(("hCrdcCath_" + to_string(i) + "_" + to_string(j)).c_str(), ("Energy vs pad index of S800 CRDC" + to_string(i + 1) + "-Cathode;iPad;Energy [ch.]").c_str(), 256, 0, 256, 1004,10,1014);

		nEnt = tData->GetEntries();
		for (iEnt = 0; iEnt < nEnt; iEnt++)
		{
			tData->GetEntry(iEnt);

			for (i = 0; i < 12; i++)
				hAdc[i]->Fill(madc.data[IdAdc[i]]);

			for (i = 0; i < 8; i++)
				hTdc[i]->Fill(mtdc.data[IdTdc[i]]);

			for (i = 0; i < 8; i++)
				hQdcTof[i]->Fill(cqdcTOF.data[IdQdcTof[i]]);

			for (i = 0; i < 8; i++)
				hQdcMcp[i]->Fill(mqdcMCP.data[IdQdcMcp[i]]);

			for (i = 0; i < 16; i++)
				hAdcIC[i]->Fill(s800.ionCham[i]);

			for (i = 0; i < 2; i++)
				for (j = 0; j < 2; j++)
					hAdcCrdcAnod[i][j]->Fill(s800.crdcAnode[i][j]);

			for (i = 0; i < 2; i++)
				for (j = 1; j <= 4; j++)
					for (k = 0; k < 64; k++)
						hCrdcCath[i]->Fill((j - 1) * 64 + k, s800.crdcCath[i][j][k]);
		}
		sCanvAdc = "Spectra_MADC32_run" + strRun;
		TCanvas *canvAdc = new TCanvas(sCanvAdc.c_str(), sCanvAdc.c_str(), 4800, 2700);
		canvAdc->Divide(4, 3);
		for (i = 0; i < 12; i++)
		{
			canvAdc->cd(i + 1);
			hAdc[i]->Draw();
		}
		canvAdc->SaveAs((sDirADC + sCanvAdc + ".png").c_str());
		canvAdc->SaveAs((sDirRun + sCanvAdc + ".png").c_str());

		sCanvTdc = "Spectra_MTDC32_run" + strRun;
		TCanvas *canvTdc = new TCanvas(sCanvTdc.c_str(), sCanvTdc.c_str(), 4800, 1800);
		canvTdc->Divide(4, 2);
		for (i = 0; i < 8; i++)
		{
			canvTdc->cd(i + 1);
			hTdc[i]->Draw();
		}
		canvTdc->SaveAs((sDirTDC + sCanvTdc + ".png").c_str());
		canvTdc->SaveAs((sDirRun + sCanvTdc + ".png").c_str());

		sCanvQdcTof = "Spectra_V792QDCTOF_run" + strRun;
		TCanvas *canvQdcTof = new TCanvas(sCanvQdcTof.c_str(), sCanvQdcTof.c_str(), 4800, 1800);
		canvQdcTof->Divide(4, 2);
		for (i = 0; i < 8; i++)
		{
			canvQdcTof->cd(i + 1);
			hQdcTof[i]->Draw();
		}
		canvQdcTof->SaveAs((sDirTofQDC + sCanvQdcTof + ".png").c_str());
		canvQdcTof->SaveAs((sDirRun + sCanvQdcTof + ".png").c_str());

		sCanvQdcMcp = "Spectra_MQDC32MCP_run" + strRun;
		TCanvas *canvQdcMcp = new TCanvas(sCanvQdcMcp.c_str(), sCanvQdcMcp.c_str(), 4800, 1800);
		canvQdcMcp->Divide(4, 2);
		for (i = 0; i < 8; i++)
		{
			canvQdcMcp->cd(i + 1);
			hQdcMcp[i]->Draw();
		}
		canvQdcMcp->SaveAs((sDirMcpQDC + sCanvQdcMcp + ".png").c_str());
		canvQdcMcp->SaveAs((sDirRun + sCanvQdcMcp + ".png").c_str());

		sCanvIC = "Spectra_ADC7164hS800IC_run" + strRun;
		TCanvas *canvAdcIC = new TCanvas(sCanvIC.c_str(), sCanvIC.c_str(), 4800, 3600);
		canvAdcIC->Divide(4, 4);
		for (i = 0; i < 16; i++)
		{
			canvAdcIC->cd(i + 1);
			hAdcIC[i]->Draw();
		}
		canvAdcIC->SaveAs((sDirIC + sCanvIC + ".png").c_str());
		canvAdcIC->SaveAs((sDirRun + sCanvIC + ".png").c_str());

		sCanvCrdcAnod = "Spectra_ADC7164hS800CRDC_Anode_run" + strRun;
		TCanvas *canvAdcCrdcAnod = new TCanvas(sCanvCrdcAnod.c_str(), sCanvCrdcAnod.c_str(), 2400, 1800);
		canvAdcCrdcAnod->Divide(2, 2);
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
			{
				canvAdcCrdcAnod->cd(1 + j + 2 * i);
				hAdcCrdcAnod[i][j]->Draw();
			}
		canvAdcCrdcAnod->SaveAs((sDirCRDC + sCanvCrdcAnod + ".png").c_str());
		canvAdcCrdcAnod->SaveAs((sDirRun + sCanvCrdcAnod + ".png").c_str());

		sCanvCrdcCath = "Spectra_XLM72S800CRDC_Cathode_run" + strRun;
		TCanvas *canvCrdcCath = new TCanvas(sCanvCrdcCath.c_str(), sCanvCrdcCath.c_str(), 3600, 1800);
		canvCrdcCath->Divide(1, 2);
		for (i = 0; i < 2; i++)
		{
			canvCrdcCath->cd(i + 1);
			hCrdcCath[i]->Draw("colz");
		}
		canvCrdcCath->SaveAs((sDirCRDC + sCanvCrdcCath + ".png").c_str());
		canvCrdcCath->SaveAs((sDirRun + sCanvCrdcCath + ".png").c_str());

		delete canvCrdcCath;
		delete canvAdcCrdcAnod;
		delete canvAdcIC;
		delete canvQdcMcp;
		delete canvQdcTof;
		delete canvTdc;
		delete canvAdc;

		for (i = 0; i < 2; i++)
			delete hCrdcCath[i];

		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				delete hAdcCrdcAnod[i][j];

		for (i = 0; i < 16; i++)
			delete hAdcIC[i];

		for (i = 0; i < 8; i++)
			delete hQdcMcp[i];

		for (i = 0; i < 8; i++)
			delete hQdcTof[i];

		for (i = 0; i < 8; i++)
			delete hTdc[i];

		for (i = 0; i < 12; i++)
			delete hAdc[i];

		fRoot->Close();
		delete fRoot;
	}
}

#ifndef __CINT__
/*
void StandaloneApplication(int argc, char** argv)
{
	Root2SpectraDCs();
}

int main(int argc, char** argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	app.Run();
	return 0;
}
*/
int main()
{
	Root2SpectraDCs();
	return 0;
}
#endif