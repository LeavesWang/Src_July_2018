#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <omp.h>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TH1S.h"
#include "TH1C.h"
#include "TROOT.h"

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
	int trig;
	int tof[8];				//[0]-->[1]: ORTEC TAC+Phillips ADC; [2]-->[7]: Phillips TDC
	int crdcCath[2][5][64]; //[2]: two CRDC; [5]: [0]---sample; [1]-->[4]---energy
	int crdcAnode[2][2];	// [0]: energy; [1]: time
	int hodoEgy[32];		// 32 crystals
	int hodoTime;
	int pin[5];
	int mesyTDC[16];
};

struct StrtAna
{
	double tof[4]; //[4]: [0] TAC+ADC+clock; [1] TAC+ADC; [2]: regular CFD+TDC; [3] MCFD16+TDC
	double tD[4][8][8];
	double egyPMT[2][4]; //energies in 4 PMTs and each plastic
	double egyPla[2];	 //[2]: [0] is energy in plastic at S800, [1] is energy in plastic at A1900
	double xPlaT[4][2];	 //[2]: [0] for S800 plastic from time info, [1] for A1900 plastic from time info
	double yPlaT[4][2];
	double xPlaQ[2]; //[2]: [0] for S800 plastic from amp info, [1] for A1900 plastic from amp info
	double yPlaQ[2];
	double xMCP[2]; //two gain settings: [0] high gain setting; [1] comb-gain setting
	double yMCP[2]; //two gain settings: [0] high gain setting; [1] comb-gain setting
	double delE[5];
	double tke;
	double beta[4];
	double gamma[4];
	double Z[4];
	double dZ[4];
	double brho[2]; //two gain settings: [0] high gain setting; [1] comb-gain setting
	double AoQ[4][2];
	double Q[4][2];
	double ZmQ[4][2];
	double ZImQ[4][2];
	double A[4][2];
	double Araw[4][2];
	double Am2Q[4][2];
	double Am3Q[4][2];
	double Am2Z[4][2];
	double Am3Z[4][2];
	double dAm2Z[4][2];
	double dAm3Z[4][2];
	int numTime[4][2]; // [2]: [0] number of fired PMTs of S800 Plastic; [1] number of fired PMTs of A1900 Plastic
	int sigTime[4][2]; // [2]: [0] indexes (combination of 1,2,3,4) of fired PMTs of S800 plastic; [1] indexes (combination of 5,6,7,8) of fired PMTs of A1900 plastic;
	int Zi[4];
	int Am2Zi[4][2];
	int Am3Zi[4][2];
	double xCrdc[2][2];		 //first [2]: two CRDCs; second [2]: [0] is from gravity center; [1] is from Gaussian fit
	double xEgyCrdc[2][256]; //256 energies for 256 cathode pads
	double yCrdc[2];		 //y is only considered from electron's drift time
	int mulHodo;
	double egyHodo[32];
};

struct StrtPid
{
	double tof;
	double beta;
	double gamma;
	double Z;
	double dZ;
	double AoQ;
	double Q;
	double ZmQ;
	double ZImQ;
	double A;
	double Araw;
	double Am2Q;
	double Am3Q;
	double Am2Z;
	double Am3Z;
	double dAm2Z;
	double dAm3Z;
	int Zi;
	int Am2Zi;
	int Am3Zi;
};

const int AdcPinLow = 10, AdcPinUp = 4096;
const int AdcTofLow = 500, AdcTofUp = 7680;
const int TdcTofLow[16] = {12000, 12000, 12000, 12000, 35000, 29000, 36000, 36000, 10000, 10000, 10000, 10000, 32000, 12000, 1, 14500};
const int TdcTofUp[16] = {22000, 22000, 22000, 22000, 55000, 49000, 56000, 56000, 20000, 20000, 20000, 20000, 52000, 32000, 1, 34500};
const int QdcTofLow[8] = {729, 750, 754, 778, 770, 780, 760, 780};
const int QdcMcpLow[8] = {759, 765, 802, 791, 752, 764, 765, 760}; //new pedestal values
const int QdcUp = 3840;

const double CALADC[12] = {6.46209, 6.59645, 6.56230, 6.57185, 6.44156, 6.58265, 6.64827, 6.52219, 6.45537, 6.42844, 6.65406, 6.43436}; //unit: ps/ch
const double CALTDC = 3.90625;																											//ps/ch
const double CALPIN[5][2] = {{0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}};
// const double CALTKE[6]={81.071121, 1.071346, 0.662579, 3.013299, 2.826749, 0};  //unit: MeV/ch
const double CALTKE[6] = {-174.696239, 1.103438, 0.733263, 2.960807, 2.920445, 7.258505};

const double CALXMCP[2][10] = {{0.733129, 26.744387, -0.091781, 1.043661, 0.047598, 9.192684, 2.637526, -0.929438, 2.056948, 0.576781}, {0.802060, 26.063777, -0.897100, 1.296354, 1.163047, 11.688516, 3.208674, -1.230582, -2.736673, 3.004569}};	   //[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*x*x*x+[6]*y*y*y+[7]*x*y+[8]*x*x*y+[9]*x*y*y  //two gain settings: [0] high gain setting; [1] comb-gain setting
const double CALYMCP[2][10] = {{3.652901, 19.180574, 1.578795, -1.716251, 0.330541, 11.410052, -0.641449, -0.958885, 0.507911, 5.328422}, {3.727687, 18.762661, -0.510623, -1.588110, -0.511162, 10.227921, -1.138502, 0.227536, 0.858179, 4.114189}}; //[0]+[1]*y+[2]*y*y+[3]*x+[4]*x*x+[5]*y*y*y+[6]*x*x*x+[7]*x*y+[8]*y*y*x+[9]*x*x*y  //two gain settings: [0] high gain setting; [1] comb-gain setting
// const double CALXMCP[2][10]={{0,1,0,0,0,0,0,0,0,0}, {0,1,0,0,0,0,0,0,0,0}}; //[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*x*x*x+[6]*y*y*y+[7]*x*y+[8]*x*x*y+[9]*x*y*y //for raw pos
// const double CALYMCP[2][10]={{0,1,0,0,0,0,0,0,0,0}, {0,1,0,0,0,0,0,0,0,0}}; //[0]+[1]*y+[2]*y*y+[3]*x+[4]*x*x+[5]*y*y*y+[6]*x*x*x+[7]*x*y+[8]*y*y*x+[9]*x*x*y //for raw pos

const double CALTOF[4][2] = {{487.3099, -0.001}, {513.3008, -0.001}, {578.8793, -0.001}, {500, -0.001}};
const double BRHO0 = 3.7211; //Tm
const double DISP = 106.84;	 // mm/%
const double LOF = 60.763;	 //m
const double CALZ[4][2] = {{1.2768, 5.9166}, {1.1976, 5.9335}, {1.3354, 5.9030}, {0, 1}};
// const double CALZ[4][2]={{0,1},{0,1},{0,1},{0,1}};

const double CALTOF_PID[2] = {578.8793, -0.001};
const double CALZ_PID[2] = {1.3354, 5.9030};

const int iLow[4] = {6, 5, 4, 7};
const string sSet[2] = {"PS_270_382", "RS_270_382"};

void Root2Ana()
{
	int i, j, k, m, n, p;

	double tDifRang[8][8][3][3] = {0}; //First [3]: [0] target peak; [1] 1st other peak; [2] 2nd other peak.  Second [3]: [0] shift value; [1] lower limit of the peak; [2] upper limit of the peak

	tDifRang[0][4][0][1] = 91946.1615;
	tDifRang[0][4][0][2] = 92418.2385;
	tDifRang[0][4][0][0] = 0;
	tDifRang[0][4][2][1] = 137547.8215;
	tDifRang[0][4][2][2] = 137858.1785;
	tDifRang[0][4][2][0] = -45520.8;
	tDifRang[0][5][0][1] = 92098.117;
	tDifRang[0][5][0][2] = 92588.683;
	tDifRang[0][5][0][0] = 0;
	tDifRang[0][5][1][1] = 46656.1985;
	tDifRang[0][5][1][2] = 47095.0015;
	tDifRang[0][5][1][0] = 45467.8;
	tDifRang[0][5][2][1] = 137702.6235;
	tDifRang[0][5][2][2] = 137985.3765;
	tDifRang[0][5][2][0] = -45500.6;
	tDifRang[0][6][0][1] = 92299.8165;
	tDifRang[0][6][0][2] = 92683.1835;
	tDifRang[0][6][0][0] = 0;
	tDifRang[0][6][2][1] = 137816.3455;
	tDifRang[0][6][2][2] = 138141.6545;
	tDifRang[0][6][2][0] = -45487.5;
	tDifRang[0][7][0][1] = 92552.782;
	tDifRang[0][7][0][2] = 92928.818;
	tDifRang[0][7][0][0] = 0;
	tDifRang[0][7][2][1] = 138063.1035;
	tDifRang[0][7][2][2] = 138362.8965;
	tDifRang[0][7][2][0] = -45472.2;
	tDifRang[1][4][0][1] = 91030.6975;
	tDifRang[1][4][0][2] = 91482.3025;
	tDifRang[1][4][0][0] = 0;
	tDifRang[1][4][2][1] = 136646.0365;
	tDifRang[1][4][2][2] = 136935.9635;
	tDifRang[1][4][2][0] = -45534.5;
	tDifRang[1][5][0][1] = 91208.8365;
	tDifRang[1][5][0][2] = 91645.5635;
	tDifRang[1][5][0][0] = 0;
	tDifRang[1][5][1][1] = 45725.1385;
	tDifRang[1][5][1][2] = 46164.8615;
	tDifRang[1][5][1][0] = 45482.2;
	tDifRang[1][5][2][1] = 136805.889;
	tDifRang[1][5][2][2] = 137068.111;
	tDifRang[1][5][2][0] = -45509.8;
	tDifRang[1][6][0][1] = 91374.381;
	tDifRang[1][6][0][2] = 91755.619;
	tDifRang[1][6][0][0] = 0;
	tDifRang[1][6][2][1] = 136923.1805;
	tDifRang[1][6][2][2] = 137210.8195;
	tDifRang[1][6][2][0] = -45502;
	tDifRang[1][7][0][1] = 91622.6245;
	tDifRang[1][7][0][2] = 92006.3755;
	tDifRang[1][7][0][0] = 0;
	tDifRang[1][7][2][1] = 137163.732;
	tDifRang[1][7][2][2] = 137438.268;
	tDifRang[1][7][2][0] = -45486.5;
	tDifRang[2][4][0][1] = 91184.7185;
	tDifRang[2][4][0][2] = 91627.8815;
	tDifRang[2][4][0][0] = 0;
	tDifRang[2][4][2][1] = 136788.724;
	tDifRang[2][4][2][2] = 137079.276;
	tDifRang[2][4][2][0] = -45527.7;
	tDifRang[2][5][0][1] = 91344.4205;
	tDifRang[2][5][0][2] = 91799.9795;
	tDifRang[2][5][0][0] = 0;
	tDifRang[2][5][1][1] = 45878.2845;
	tDifRang[2][5][1][2] = 46317.7155;
	tDifRang[2][5][1][0] = 45474.2;
	tDifRang[2][5][2][1] = 136945.1495;
	tDifRang[2][5][2][2] = 137212.8505;
	tDifRang[2][5][2][0] = -45506.8;
	tDifRang[2][6][0][1] = 91528.282;
	tDifRang[2][6][0][2] = 91901.718;
	tDifRang[2][6][0][0] = 0;
	tDifRang[2][6][2][1] = 137069.6845;
	tDifRang[2][6][2][2] = 137350.3155;
	tDifRang[2][6][2][0] = -45495;
	tDifRang[2][7][0][1] = 91769.18;
	tDifRang[2][7][0][2] = 92158.42;
	tDifRang[2][7][0][0] = 0;
	tDifRang[2][7][2][1] = 137303.185;
	tDifRang[2][7][2][2] = 137586.815;
	tDifRang[2][7][2][0] = -45481.2;
	tDifRang[3][4][0][1] = 90039.978;
	tDifRang[3][4][0][2] = 90500.822;
	tDifRang[3][4][0][0] = 0;
	tDifRang[3][4][2][1] = 135654.7095;
	tDifRang[3][4][2][2] = 135969.2905;
	tDifRang[3][4][2][0] = -45541.6;
	tDifRang[3][5][0][1] = 90218.2185;
	tDifRang[3][5][0][2] = 90672.5815;
	tDifRang[3][5][0][0] = 0;
	tDifRang[3][5][1][1] = 44722.1455;
	tDifRang[3][5][1][2] = 45187.0545;
	tDifRang[3][5][1][0] = 45490.8;
	tDifRang[3][5][2][1] = 135828.6465;
	tDifRang[3][5][2][2] = 136103.3535;
	tDifRang[3][5][2][0] = -45520.6;
	tDifRang[3][6][0][1] = 90370.358;
	tDifRang[3][6][0][2] = 90784.842;
	tDifRang[3][6][0][0] = 0;
	tDifRang[3][6][2][1] = 135939.9985;
	tDifRang[3][6][2][2] = 136238.0015;
	tDifRang[3][6][2][0] = -45511.4;
	tDifRang[3][7][0][1] = 90609.9805;
	tDifRang[3][7][0][2] = 91044.4195;
	tDifRang[3][7][0][0] = 0;
	tDifRang[3][7][2][1] = 136171.7905;
	tDifRang[3][7][2][2] = 136474.2095;
	tDifRang[3][7][2][0] = -45495.8;
	tDifRang[4][5][0][1] = -66.9775;
	tDifRang[4][5][0][2] = 433.8495;
	tDifRang[4][5][0][0] = 0;
	tDifRang[4][5][1][1] = -45668.889;
	tDifRang[4][5][1][2] = -45029.711;
	tDifRang[4][5][1][0] = 45532.736;
	tDifRang[4][6][0][1] = 111.184;
	tDifRang[4][6][0][2] = 475.256;
	tDifRang[4][6][0][0] = 0;
	tDifRang[4][6][1][1] = 45479.5005;
	tDifRang[4][6][1][2] = 46148.0995;
	tDifRang[4][6][1][0] = -45520.58;
	tDifRang[4][7][0][1] = 351.784;
	tDifRang[4][7][0][2] = 717.014;
	tDifRang[4][7][0][0] = 0;
	tDifRang[5][6][0][1] = -45.1135;
	tDifRang[5][6][0][2] = 281.8895;
	tDifRang[5][6][0][0] = 0;
	tDifRang[5][6][1][1] = 45370.395;
	tDifRang[5][6][1][2] = 45912.205;
	tDifRang[5][6][1][0] = -45522.912;
	tDifRang[5][7][0][1] = 198.884;
	tDifRang[5][7][0][2] = 541.504;
	tDifRang[5][7][0][0] = 0;
	tDifRang[5][7][1][1] = 45623.7555;
	tDifRang[5][7][1][2] = 46138.8445;
	tDifRang[5][7][1][0] = -45511.106;
	tDifRang[6][7][0][1] = 125.1625;
	tDifRang[6][7][0][2] = 360.3755;
	tDifRang[6][7][0][0] = 0;
	tDifRang[0][1][0][1] = 820.682;
	tDifRang[0][1][0][2] = 1020.174;
	tDifRang[0][1][0][0] = 0;
	tDifRang[0][2][0][1] = 677.7995;
	tDifRang[0][2][0][2] = 868.7385;
	tDifRang[0][2][0][0] = 0;
	tDifRang[0][3][0][1] = 1753.806;
	tDifRang[0][3][0][2] = 2054.154;
	tDifRang[0][3][0][0] = 0;
	tDifRang[1][2][0][1] = -222.349;
	tDifRang[1][2][0][2] = -71.923;
	tDifRang[1][2][0][0] = 0;
	tDifRang[1][3][0][1] = 877.234;
	tDifRang[1][3][0][2] = 1089.562;
	tDifRang[1][3][0][0] = 0;
	tDifRang[2][3][0][1] = 1014.0495;
	tDifRang[2][3][0][2] = 1247.9505;
	tDifRang[2][3][0][0] = 0;

	double calHodo[32][2] = {0};
	int nH = 0;
	string strRead;
	ifstream fCalHodo("fCalHodo.dat");
	while (!fCalHodo.eof() && fCalHodo.peek() != EOF)
	{
		if (fCalHodo.peek() != '#')
		{
			fCalHodo >> calHodo[nH][1] >> calHodo[nH][0];
			getline(fCalHodo, strRead);
			nH++;
		}
		else
			getline(fCalHodo, strRead);
	}
	fCalHodo.close();

	// TFile *fNonLin=new TFile("fNonLin.root");
	// TCanvas *cINLTacAdc;
	// TGraph *grINL[12];
	// vector<vector<double> > inl(12, vector<double>(8192, 0));
	// fNonLin->GetObject("cINLTacAdc", cINLTacAdc);
	// cINLTacAdc->Draw();
	// for(i=0; i<12; i++)
	// {
	// 	cINLTacAdc->cd(i+1);
	// 	grINL[i]=(TGraph*)(gPad->GetPrimitive(("INL_TacAdc_"+to_string(i)).c_str()));

	// 	for(j=0; j<grINL[i]->GetN(); j++)
	// 	{
	// 		double vINL=grINL[i]->GetY()[j];
	// 		if(abs(vINL)>0)
	// 		{
	// 			k=(int)(grINL[i]->GetX()[j]);
	// 			inl[i][k]=vINL;
	// 		}
	// 	}
	// }
	// fNonLin->Close();
	// delete fNonLin;

	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	string sRoot, sAna;
	StrtMesytec madc, mtdc, mqdcTOF, mqdcMCP;
	StrtS800 s800;
	StrtAna ana;
	StrtPid pid;

	double tPMT[8], timeDet[2];
	double calQdcMcp[2][4];

	long long iEnt, nEnt, jEnt, lEnt, rEnt;
	const int LjEnt = 0, UjEnt = 1;
	int ts[UjEnt + LjEnt] = {0};
	int tsSi = 0, tsHodo = 0;
	int iAna;
	double eHodo[32];
	int nHodo;
	int nQdcTof[2];
	int nGoodEvt, nGoodPin, nGoodMcp[2];
	bool goodPin[5];
	double b;
	string setting;
	int run;
	double dTsHodoSi;
	double xMcpRaw = 0, yMcpRaw = 0;

	int runMin, runMax, runNum;
	cout << "Input minimum and maximum numbers of run: ";
	cin >> runMin >> runMax;

	ostringstream ssRun;
	for (runNum = runMin; runNum <= runMax; runNum++)
	{
		sAna = "/home/kailong/ExpData/Jul2018/AnaData/ana-run-" + to_string(runNum) + ".root";
		TFile *fAna = new TFile(sAna.c_str(), "RECREATE");
		TTree *tAna = new TTree("tAna", "tree for data analysis");

		tAna->Branch("setting", &setting);
		tAna->Branch("run", &run, "run/I");
		// tAna->Branch("dTsHodoSi", &dTsHodoSi, "dTsHodoSi/D");
		tAna->Branch("ana", &ana, "tof[4]/D:tD[4][8][8]/D:egyPMT[2][4]/D:egyPla[2]/D:xPlaT[4][2]/D:yPlaT[4][2]/D:xPlaQ[2]/D:yPlaQ[2]/D:xMCP[2]/D:yMCP[2]/D:delE[5]/D:tke/D:beta[4]/D:gamma[4]/D:Z[4]/D:dZ[4]/D:brho[2]/D:AoQ[4][2]/D:Q[4][2]/D:ZmQ[4][2]/D:ZImQ[4][2]/D:A[4][2]/D:Araw[4][2]/D:Am2Q[4][2]/D:Am3Q[4][2]/D:Am2Z[4][2]/D:Am3Z[4][2]/D:dAm2Z[4][2]/D:dAm3Z[4][2]/D:numTime[4][2]/I:sigTime[4][2]/I:Zi[4]/I:Am2Zi[4][2]/I:Am3Zi[4][2]/I:xCrdc[2][2]/D:xEgyCrdc[2][256]/D:yCrdc[2]/D:mulHodo/I:egyHodo[32]/D");
		tAna->Branch("pid", &pid, "tof/D:beta/D:gamma/D:Z/D:dZ/D:AoQ/D:Q/D:ZmQ/D:ZImQ/D:A/D:Araw/D:Am2Q/D:Am3Q/D:Am2Z/D:Am3Z/D:dAm2Z/D:dAm3Z/D:Zi/I:Am2Zi/I:Am3Zi/I");

		run = runNum;

		ssRun.str("");
		ssRun << setw(4) << setfill('0') << runNum;

		sRoot = "/home/kailong/ExpData/Jul2018/RootData/run-" + ssRun.str() + "-00.root";
		if (access("/home/kailong/ExpData/Jul2018/AnaData", F_OK) != 0)
			system("mkdir /home/kailong/ExpData/Jul2018/AnaData");
		printf("\n**********Now converting %s to %s!**********\n\n", sRoot.c_str(), sAna.c_str());

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

		if (runNum >= 270)
		{
			tData->SetEstimate(-1);
			nEnt = tData->GetEntries();
			double *errCh = new double[nEnt];
			for (iEnt = 0; iEnt < nEnt; iEnt++)
				errCh[iEnt] = 0.5;
			// string strCanv="gainMat_"+to_string(runNum);
			// TCanvas *cGain=new TCanvas(strCanv.c_str(), strCanv.c_str());
			// cGain->Divide(2,2);
			TGraphErrors *gr[4];
			for (i = 0; i < 4; i++)
			{
				j = iLow[i];
				// cGain->cd(i+1);
				string sCut = "mqdcMCP.data[" + to_string(i) + "]>" + to_string(QdcMcpLow[i]) + "&&mqdcMCP.data[" + to_string(i) + "]<3840&&mqdcMCP.data[" + to_string(j) + "]>" + to_string(QdcMcpLow[j]) + "&&mqdcMCP.data[" + to_string(j) + "]<3840";

				string sDraw = "mqdcMCP.data[" + to_string(i) + "]-" + to_string(QdcMcpLow[i]) + ":mqdcMCP.data[" + to_string(j) + "]-" + to_string(QdcMcpLow[j]);

				// printf("Now drawing {%s} of {%s} when {%s}\n\n", sDraw.c_str(), sRoot.c_str(), sCut.c_str());

				long long nData = tData->Draw(sDraw.c_str(), sCut.c_str(), "goff");
				if (nData < 2)
					continue;
				double *highGain = tData->GetV1();
				double *lowGain = tData->GetV2();
				gr[i] = new TGraphErrors(nData, lowGain, highGain, errCh, errCh);
				// gr[i]->Draw("AP");
				// gr[i]->SetTitle(("high"+to_string(i)+"_vs_low"+to_string(iLow[i])).c_str());
				TFitResultPtr fitRes = gr[i]->Fit("pol1", "SQ");
				int fitSt = fitRes;
				if (fitSt != 0 && fitSt != 4000)
					continue;
				mcpGainMat[j][0] = fitRes->Parameter(0);
				mcpGainMat[j][1] = fitRes->Parameter(1);
				// printf("%f %f\n",mcpGainMat[j][0],mcpGainMat[j][1]);
			}
			delete[] errCh;
			// cGain->SaveAs(("/home/kailong/ExpData/Jul2018/Graphs/Charts/"+strCanv+".png").c_str());
			// cGain->Close();
			for (i = 0; i < 4; i++)
				if (!gr[i])
					delete gr[i];
			// delete cGain;
		}

		memset(&madc, 0, sizeof(madc));
		memset(&mtdc, 0, sizeof(mtdc));
		memset(&mqdcTOF, 0, sizeof(mqdcTOF));
		memset(&mqdcMCP, 0, sizeof(mqdcMCP));
		memset(&s800, 0, sizeof(s800));

		tData->SetBranchAddress("madc", &madc);
		tData->SetBranchAddress("mtdc", &mtdc);
		tData->SetBranchAddress("mqdcTOF", &mqdcTOF);
		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);
		tData->SetBranchAddress("s800", &s800);

		for (i = 0; i < 2; i++)
		{
			string sfSet = "/home/kailong/ExpData/Jul2018/Src/runNum" + sSet[i] + ".dat";
			ifstream fSet(sfSet.c_str());
			string sRead;
			bool isFound = false;
			while (getline(fSet, sRead))
				if (sRead.find(to_string(runNum)) != string::npos)
				{
					isFound = true;
					break;
				}
			if (isFound)
			{
				setting = sSet[i];
				break;
			}
			if (i == 1 && !isFound)
				setting = "other_270_382";
		}

		nEnt = tData->GetEntries();
		for (iEnt = 0; iEnt < nEnt; iEnt++)
		// for(iEnt=0; iEnt<1000; iEnt++)
		{
			tData->GetEntry(iEnt);
			tsHodo = 0;
			if (s800.trig == 1 || s800.trig == 16)
			{
				TRandom3 r1(0);

				tsHodo = s800.tS;
				nHodo = 0;
				memset(eHodo, 0, sizeof(eHodo));
				for (m = 0; m < 32; m++)
					if (s800.hodoEgy[m] > 10 && s800.hodoEgy[m] < 3000 && s800.hodoTime > 50 && s800.hodoTime < 2200)
					{
						nHodo++;
						eHodo[m] = calHodo[m][0] + calHodo[m][1] * (s800.hodoEgy[m] + r1.Uniform(-0.5, 0.5));
					}
				// if(nHodo<2)
				// continue;
				memset(ts, 0, sizeof(ts));
				n = 0;
				lEnt = (iEnt - LjEnt) > 0 ? (iEnt - LjEnt) : 0;
				rEnt = (iEnt + UjEnt) < nEnt ? (iEnt + UjEnt) : nEnt;
				for (jEnt = lEnt; jEnt < rEnt; jEnt++)
					if (s800.trig == 1 || s800.trig == 16)
					{
						tData->GetEntry(jEnt);
						if (s800.trig == 1 || s800.trig == 16)
							ts[n++] = s800.tS;
					}
				bool goodTS = false;
				for (j = 0; j < n - 1; j++)
				{
					goodTS = ts[j] < ts[j + 1];
					if (!goodTS)
						break;
				}
				if (!goodTS)
				{
					lEnt = iEnt;
					rEnt = iEnt + 1;
				}

				tsSi = 0;
				for (jEnt = lEnt; jEnt < rEnt; jEnt++)
				{
					tData->GetEntry(jEnt);
					if (s800.trig != 1 && s800.trig != 16)
						continue;
					TRandom3 r(0);

					dTsHodoSi = 0;
					memset(&ana, 0, sizeof(ana));
					memset(&pid, 0, sizeof(pid));

					tsSi = s800.tS;
					dTsHodoSi = 0.1 * (tsHodo - tsSi);
					ana.mulHodo = nHodo;
					for (i = 0; i < 32; i++)
						ana.egyHodo[i] = eHodo[i];

					ana.tke = CALTKE[0];
					nGoodPin = 0;
					memset(goodPin, 0, sizeof(goodPin));
					for (i = 0; i < 5; i++)
						if (s800.pin[i] > AdcPinLow && s800.pin[i] < AdcPinUp)
						{
							goodPin[i] = true;
							nGoodPin++;

							ana.delE[i] = CALPIN[i][0] + CALPIN[i][1] * (s800.pin[i] + r.Uniform(-0.5, 0.5));
							ana.tke += CALTKE[i + 1] * (s800.pin[i] + r.Uniform(-0.5, 0.5));
						}

					for (m = 0; m < 2; m++)
						if (s800.crdcAnode[m][0] > 0 && s800.crdcAnode[m][1] > 0)
						{
							//calculate y from the drift time to anode
							ana.yCrdc[m] = s800.crdcAnode[m][1];

							TH1S hx("hx", "hx", 256, 0, 256);
							TH1C hc("hc", "hc", 256, 0, 256);
							for (i = 1; i <= 4; i++)
								for (j = 0; j < 64; j++)
								{
									k = (i - 1) * 64 + j;
									if (k > 90 && k < 130 && s800.crdcCath[m][i][j] > 100 && s800.crdcCath[m][i][j] < 990)
									{
										ana.xEgyCrdc[m][k] = s800.crdcCath[m][i][j];
										hx.Fill(k, s800.crdcCath[m][i][j]);
										hc.Fill(k);
									}
								}
							if (hx.GetEntries() > 2)
							{
								double maxPad = hx.GetBinCenter(hx.GetMaximumBin());
								hx.SetAxisRange(maxPad - 10, maxPad + 10, "X");
								hc.SetAxisRange(maxPad - 10, maxPad + 10, "X");
								if (hc.Integral() > 2)
								{
									//calculate x from gravity center
									ana.xCrdc[m][0] = hx.GetMean();
									//calculate x from Gaussian fit
									TFitResultPtr fitX = hx.Fit("gaus", "S0Q");
									int fitSt = fitX;
									if (fitSt != 0 && fitX->Parameter(1) > 0 && fitX->Parameter(1) < 256)
										ana.xCrdc[m][1] = fitX->Parameter(1);
								}
							}
						}
					memset(nQdcTof, 0, sizeof(nQdcTof));
					for (i = 0; i < 8; i++)
						if (mqdcTOF.data[2 * i + 1] > QdcTofLow[i] && mqdcTOF.data[2 * i + 1] < QdcUp)
						{
							j = i - i / 4 * 4;
							nQdcTof[i / 4]++;
							ana.egyPMT[i / 4][j] = mqdcTOF.data[2 * i + 1] + r.Uniform(-0.5, 0.5) - QdcTofLow[i];
						}
					if (nQdcTof[0] == 4)
					{
						ana.xPlaQ[0] = log(ana.egyPMT[0][1] / ana.egyPMT[0][3]);
						ana.yPlaQ[0] = log(ana.egyPMT[0][2] / ana.egyPMT[0][0]);
						ana.egyPla[0] = (ana.egyPMT[0][0] + ana.egyPMT[0][1] + ana.egyPMT[0][2] + ana.egyPMT[0][3]) / 4;
					}
					if (nQdcTof[1] == 4)
					{
						ana.xPlaQ[1] = log(ana.egyPMT[1][1] / ana.egyPMT[1][3]);
						ana.yPlaQ[1] = log(ana.egyPMT[1][0] / ana.egyPMT[1][2]);
						ana.egyPla[1] = (ana.egyPMT[1][0] + ana.egyPMT[1][1] + ana.egyPMT[1][2] + ana.egyPMT[1][3]) / 4;
					}
					if (nGoodPin > 1 && goodPin[0])
					{
						memset(calQdcMcp, 0, sizeof(calQdcMcp));
						memset(nGoodMcp, 0, sizeof(nGoodMcp));
						if (runNum >= 270)
						{
							for (i = 0; i < 4; i++)
							{
								if (mqdcMCP.data[i] > QdcMcpLow[i] && mqdcMCP.data[i] < QdcUp)
								{
									nGoodMcp[0]++;
									calQdcMcp[0][i] = mqdcMCP.data[i] - QdcMcpLow[i] + r.Uniform(-0.5, 0.5);

									nGoodMcp[1]++;
									calQdcMcp[1][i] = calQdcMcp[0][i];
								}
								m = iLow[i];
								if (mqdcMCP.data[i] >= QdcUp && mqdcMCP.data[m] > QdcMcpLow[m] && mqdcMCP.data[m] < QdcUp)
								{

									calQdcMcp[1][i] = mcpGainMat[m][0] + mcpGainMat[m][1] * (mqdcMCP.data[m] - QdcMcpLow[m] + r.Uniform(-0.5, 0.5));
									if (calQdcMcp[1][i] > QdcUp - QdcMcpLow[i])
										nGoodMcp[1]++;
								}
							}

							for (i = 0; i < 2; i++)
								if (nGoodMcp[i] == 4)
								{
									xMcpRaw = (calQdcMcp[i][0] + calQdcMcp[i][3] - calQdcMcp[i][1] - calQdcMcp[i][2]) / (calQdcMcp[i][0] + calQdcMcp[i][1] + calQdcMcp[i][2] + calQdcMcp[i][3]);

									yMcpRaw = (calQdcMcp[i][1] + calQdcMcp[i][3] - calQdcMcp[i][0] - calQdcMcp[i][2]) / (calQdcMcp[i][0] + calQdcMcp[i][1] + calQdcMcp[i][2] + calQdcMcp[i][3]);

									ana.xMCP[i] = CALXMCP[i][0] + CALXMCP[i][1] * xMcpRaw + CALXMCP[i][2] * pow(xMcpRaw, 2) + CALXMCP[i][3] * yMcpRaw + CALXMCP[i][4] * pow(yMcpRaw, 2) + CALXMCP[i][5] * pow(xMcpRaw, 3) + CALXMCP[i][6] * pow(yMcpRaw, 3) + CALXMCP[i][7] * xMcpRaw * yMcpRaw + CALXMCP[i][8] * pow(xMcpRaw, 2) * yMcpRaw + CALXMCP[i][9] * xMcpRaw * pow(yMcpRaw, 2);

									ana.yMCP[i] = CALYMCP[i][0] + CALYMCP[i][1] * yMcpRaw + CALYMCP[i][2] * pow(yMcpRaw, 2) + CALYMCP[i][3] * xMcpRaw + CALYMCP[i][4] * pow(xMcpRaw, 2) + CALYMCP[i][5] * pow(yMcpRaw, 3) + CALYMCP[i][6] * pow(xMcpRaw, 3) + CALYMCP[i][7] * yMcpRaw * xMcpRaw + CALYMCP[i][8] * pow(yMcpRaw, 2) * xMcpRaw + CALYMCP[i][9] * yMcpRaw * pow(xMcpRaw, 2);
								}
						}
						for (i = 0; i < 2; i++)
							ana.brho[i] = BRHO0 * (1 + ana.xMCP[i] / DISP / 100);

						nGoodEvt = 0;
						for (iAna = 0; iAna < 4; iAna++)
						{
							memset(tPMT, 0, sizeof(tPMT));
							memset(timeDet, 0, sizeof(timeDet));

							if (iAna == 0) //For TAC+ADC+clock
							{
								for (i = 0; i < 8; i++)
								{
									k = i / 4;
									if (madc.data[i] > AdcTofLow && madc.data[i] < AdcTofUp)
									{
										ana.numTime[0][k]++;
										ana.sigTime[0][k] = 10 * ana.sigTime[0][k] + (i + 1);
										tPMT[i] = (CALADC[i] + 0) * (-madc.data[i] + r.Uniform(-0.5, 0.5));
										timeDet[k] += tPMT[i];
									}
								}
								for (i = 0; i < 7; i++)
									for (j = i + 1; j < 8; j++)
										if (abs(tPMT[i]) > 0 && abs(tPMT[j]) > 0)
										{
											ana.tD[0][i][j] = tPMT[i] - tPMT[j];
											ana.tD[0][j][i] = -ana.tD[0][i][j];
										}
								if (ana.numTime[0][0] > 0 && ana.numTime[0][1] > 0)
									ana.tof[0] = timeDet[1] / ana.numTime[0][1] - timeDet[0] / ana.numTime[0][0];
							}

							if (iAna == 1) //For TAC+ADC
							{
								for (i = 8; i < 12; i++)
								{
									m = i - 8;
									p = i - 4;
									if (madc.data[i] > AdcTofLow && madc.data[i] < AdcTofUp)
									{
										ana.numTime[1][0]++;
										ana.numTime[1][1]++;
										ana.sigTime[1][0] = 10 * ana.sigTime[1][0] + (m + 1);
										ana.sigTime[1][1] = 10 * ana.sigTime[1][1] + (p + 1);
										ana.tD[1][p][m] = (CALADC[i] + 0) * (madc.data[i] + r.Uniform(-0.5, 0.5));
										ana.tD[1][m][p] = -ana.tD[1][p][m];
										ana.tof[1] += ana.tD[1][p][m];
									}
								}
								if (ana.numTime[1][0] == 4)
									ana.tof[1] /= ana.numTime[1][0];
							}

							if (iAna >= 2)
							{
								for (j = 0; j < 8; j++)
								{
									if (iAna == 2)
									{
										k = j + 1;
										p = j;
									}
									if (iAna == 3)
									{
										k = 2 * j + 17;
										p = j + 8;
									}
									n = j / 4;
									if (mtdc.data[k] > TdcTofLow[p] && mtdc.data[k] < TdcTofUp[p])
									{
										ana.numTime[iAna][n]++;
										ana.sigTime[iAna][n] = 10 * ana.sigTime[iAna][n] + (j + 1);
										tPMT[j] = CALTDC * (mtdc.data[k] + r.Uniform(-0.5, 0.5));
										timeDet[n] += tPMT[j];
									}
								}
								for (i = 0; i < 7; i++)
									for (j = i + 1; j < 8; j++)
										if (abs(tPMT[i]) > 0 && abs(tPMT[j]) > 0)
										{
											ana.tD[iAna][i][j] = tPMT[i] - tPMT[j];
											ana.tD[iAna][j][i] = -ana.tD[iAna][i][j];
										}
								if (ana.numTime[iAna][0] == 4 && ana.numTime[iAna][1] == 4)
									ana.tof[iAna] = timeDet[1] / ana.numTime[iAna][1] - timeDet[0] / ana.numTime[iAna][0];
							}
						}

						// // For combining 2 or 3 peaks of TAC (NoClk) + ADC
						if (ana.numTime[0][0] == 4 && ana.numTime[0][1] == 4 && ana.numTime[2][0] == 4 && ana.numTime[2][1] == 4)
						{
							bool isShift[8][8] = {0};
							for (i = 0; i < 7; i++)
								for (j = i + 1; j < 8; j++)
								{
									isShift[i][j] = false;
									if (abs(ana.tD[0][i][j]) > 0 && abs(ana.tD[2][i][j]) > 0)
									{
										// for (k = 0; k < 1; k++) // For the most statistic peak (k=0)
										for (k = 0; k < 3; k++)  // For all peaks
											if (abs(tDifRang[i][j][k][1]) > 0 && abs(tDifRang[i][j][k][2]) > 0 && (ana.tD[0][i][j] - ana.tD[2][i][j]) > tDifRang[i][j][k][1] && (ana.tD[0][i][j] - ana.tD[2][i][j]) < tDifRang[i][j][k][2])
											{
												isShift[i][j] = true;
												ana.tD[0][i][j] += tDifRang[i][j][k][0];
												break; //because only one "k" is possible for one event
											}
									}
									if (!isShift[i][j])
										ana.tD[0][i][j] = 0;
									ana.tD[0][j][i] = -ana.tD[0][i][j];
								}
							ana.tof[0] = 0;
							k = 0;
							for (i = 0; i < 4; i++)
								for (j = 4; j < 8; j++)
									if (abs(ana.tD[0][j][i]) > 0 && isShift[i][j])
									{
										ana.tof[0] += ana.tD[0][j][i];
										k++;
									}
							if(k==16)
								ana.tof[0] /= k;
							else
								ana.tof[0] = 0;
						}

						for (iAna = 0; iAna < 4; iAna++)
							if (ana.numTime[iAna][0] == 4 && ana.numTime[iAna][1] == 4 && abs(ana.tof[iAna]) > 0)
							{
								ana.tof[iAna] = CALTOF[iAna][0] + CALTOF[iAna][1] * ana.tof[iAna];

								if (iAna != 1 && ana.numTime[iAna][0] == 4)
								{
									// ana.xPlaT[iAna][0]=(ana.tD[iAna][2][0]+ana.tD[iAna][3][1]+ana.tD[iAna][2][1]+ana.tD[iAna][3][0])/4;
									// ana.yPlaT[iAna][0]=(ana.tD[iAna][0][2]+ana.tD[iAna][3][1]+ana.tD[iAna][0][1]+ana.tD[iAna][3][2])/4;
									ana.xPlaT[iAna][0] = ana.tD[iAna][3][1];
									ana.yPlaT[iAna][0] = ana.tD[iAna][0][2];
								}
								if (iAna != 1 && ana.numTime[iAna][1] == 4)
								{
									// ana.xPlaT[iAna][1]=(ana.tD[iAna][4][6]+ana.tD[iAna][7][5]+ana.tD[iAna][4][5]+ana.tD[iAna][7][6])/4;
									// ana.yPlaT[iAna][1]=(ana.tD[iAna][6][4]+ana.tD[iAna][7][5]+ana.tD[iAna][6][5]+ana.tD[iAna][7][4])/4;
									ana.xPlaT[iAna][1] = ana.tD[iAna][7][5];
									ana.yPlaT[iAna][1] = ana.tD[iAna][6][4];
								}

								b = LOF / ana.tof[iAna] / 0.299792458;
								if (b > 0 && b < 1)
								{
									ana.beta[iAna] = b;
									ana.gamma[iAna] = 1 / sqrt(1 - b * b);

									ana.Z[iAna] = sqrt(ana.delE[0] / (log(5930 / (1 / b / b - 1)) / b / b - 1));
									ana.Z[iAna] = CALZ[iAna][0] + CALZ[iAna][1] * ana.Z[iAna];
									ana.Zi[iAna] = TMath::Nint(ana.Z[iAna]);
									ana.dZ[iAna] = ana.Z[iAna] - ana.Zi[iAna];

									if (runNum == 150 || runNum == 152 || runNum == 153 || (nGoodMcp[1] == 4 && nQdcTof[0] > 0 && nQdcTof[1] > 0))
									{
										for (k = 0; k < 2; k++)
										{
											ana.AoQ[iAna][k] = ana.brho[k] / ana.beta[iAna] / ana.gamma[iAna] * 0.32184043;
											ana.Q[iAna][k] = ana.tke / (931.4940954 * (ana.gamma[iAna] - 1) * ana.AoQ[iAna][k]);
											ana.ZmQ[iAna][k] = ana.Z[iAna] - ana.Q[iAna][k];
											ana.ZImQ[iAna][k] = ana.Zi[iAna] - ana.Q[iAna][k];
											ana.A[iAna][k] = ana.AoQ[iAna][k] * ana.Q[iAna][k];
											ana.Araw[iAna][k] = ana.AoQ[iAna][k] * ana.Zi[iAna];
											ana.Am2Q[iAna][k] = ana.A[iAna][k] - 2 * ana.Q[iAna][k];
											ana.Am3Q[iAna][k] = ana.A[iAna][k] - 3 * ana.Q[iAna][k];
											ana.Am2Z[iAna][k] = ana.Araw[iAna][k] - 2 * ana.Zi[iAna];
											ana.Am3Z[iAna][k] = ana.Araw[iAna][k] - 3 * ana.Zi[iAna];
											ana.Am2Zi[iAna][k] = TMath::Nint(ana.Am2Z[iAna][k]);
											ana.Am3Zi[iAna][k] = TMath::Nint(ana.Am3Z[iAna][k]);
											ana.dAm2Z[iAna][k] = ana.Am2Z[iAna][k] - ana.Am2Zi[iAna][k];
											ana.dAm3Z[iAna][k] = ana.Am3Z[iAna][k] - ana.Am3Zi[iAna][k];
											// if(ana.Am3Z[iAna][k]<0)
											// ana.dAm3Z[iAna][k]+=1;
										}
										nGoodEvt++;
									}
								}
							}

						if (nGoodEvt > 1 && ana.numTime[2][0] == 4 && ana.numTime[2][1] == 4 && (runNum == 150 || runNum == 152 || runNum == 153 || (nGoodMcp[1] == 4 && nQdcTof[0] == 4 && nQdcTof[1] == 4))) //standard condition  //begin to fill branch of pid
						{
							pid.tof = CALTOF_PID[0] + CALTOF_PID[1] * (ana.tD[2][4][0] + ana.tD[2][5][1] + ana.tD[2][6][2] + ana.tD[2][7][3]) / 4;

							b = LOF / pid.tof / 0.299792458;
							if (b > 0 && b < 1)
							{
								pid.beta = b;
								pid.gamma = 1 / sqrt(1 - b * b);

								pid.Z = sqrt(ana.delE[0] / (log(5930 / (1 / b / b - 1)) / b / b - 1));
								pid.Z = CALZ_PID[0] + CALZ_PID[1] * pid.Z;
								pid.Zi = TMath::Nint(pid.Z);
								pid.dZ = pid.Z - pid.Zi;
								pid.AoQ = ana.brho[1] / pid.beta / pid.gamma * 0.32184;
								pid.Q = ana.tke / (931.4940954 * (pid.gamma - 1) * pid.AoQ);
								pid.ZmQ = pid.Z - pid.Q;
								pid.ZImQ = pid.Zi - pid.Q;
								pid.A = pid.AoQ * pid.Q;
								pid.Araw = pid.AoQ * pid.Zi;
								pid.Am2Q = pid.A - 2 * pid.Q;
								pid.Am3Q = pid.A - 3 * pid.Q;
								pid.Am2Z = pid.Araw - 2 * pid.Zi;
								pid.Am3Z = pid.Araw - 3 * pid.Zi;
								pid.Am2Zi = TMath::Nint(pid.Am2Z);
								pid.Am3Zi = TMath::Nint(pid.Am3Z);
								pid.dAm2Z = pid.Am2Z - pid.Am2Zi;
								pid.dAm3Z = pid.Am3Z - pid.Am3Zi;
								// if(pid.Am3Z<0)
								// pid.dAm3Z+=1;
								tAna->Fill();
							}
						}
					}
				}
			}
		} //end of whole tree
		fRoot->Close();
		delete fRoot;

		fAna->cd();
		tAna->Write();
		fAna->Close();
		delete fAna;
	}
} //end of whole function

void StandaloneApplication(int argc, char **argv)
{
	Root2Ana();
}

int main(int argc, char **argv)
{
	ROOT::EnableThreadSafety();
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	// app.Run();
	return 0;
}
