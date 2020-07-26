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
#include "TCanvas.h"
#include "TH1F.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h" 

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
	int tof[8]; //[0]-->[1]: ORTEC TAC+Phillips ADC; [2]-->[7]: Phillips TDC
	int crdcCath[2][5][64]; //[2]: two CRDC; [5]: [0]---sample; [1]-->[4]---energy
	int crdcAnode[2][2]; // [0]: energy; [1]: time	
	int hodoEgy[32]; // 32 crystals
	int hodoTime;
	int pin[5];
	int mesyTDC[16];
};

struct StrtBasAlign
{
	double mean[32];
	double meanErr[32];
};

struct StrtAdvAlign
{
	double ave[32];
	double aveErr[32];
};

void Root2Align()
{
	const int MODNUM=5;
	const int NCHAN[MODNUM]={12, 16, 8, 8, 5};
	const int RNGMOD[MODNUM][2]={{100,7680}, {1,65536}, {800,3840}, {700,3840}, {10,4010}};
	const int DIVMOD[MODNUM][2]={{4,3}, {4,4}, {4,2}, {4,2}, {3,2}};
	const int SETNUM=2;
	
	string sModule[MODNUM]={"Madc", "Mtdc", "MqdcTof", "MqdcMcp", "Pin"};
	string sRawBrh[5]={"madc", "mtdc", "mqdcTOF", "mqdcMCP", "s800"};
	string sRawLeaf[5]={"madc.data", "mtdc.data", "mqdcTOF.data", "mqdcMCP.data", "s800.pin"};
	string sAlignRawLeaf[5]={"alignMadc.data", "alignMtdc.data", "alignMqdcTof.data", "alignMqdcMcp.data", "alignS800.pin"};
	
	
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1111);
	
	int i, j, k;
	int iEntry;
	string sName, sTitle;
	ostringstream ssRun;
	TCut cutRunNum;
	string sCutRunNum;
	string sRoot;
	
	if(access("/home/kailong/ExpData/Jul2018/AlignData", F_OK)!=0)
		system("mkdir /home/kailong/ExpData/Jul2018/AlignData");
	
	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	
	string sBasAlign="/home/kailong/ExpData/Jul2018/AlignData/BasAlignRun_"+to_string(runMin)+"-"+to_string(runMax)+".root";
	StrtMesytec mstc[4];
	StrtS800 s800;
		
	TFile *fRoot, *fBasAlign, *fAdvAlign;
	TTree *tData, *tBasAlign, *tAdvAlign;
	char yn1;
	printf("Do you want to generate the total basic root files of alignibrating information from run#%d to run#%d? (y/n)--->", runMin, runMax);
	cin>>yn1;
	if(yn1=='y'||yn1=='Y')
	{			
		printf("\n**********Now generating %s!**********\n\n", sBasAlign.c_str());
		
		fBasAlign=new TFile(sBasAlign.c_str(), "RECREATE");		
		tBasAlign=new TTree("tBasAlign", "Mean, Sigma and N of histogram of data of each channel for every run");
		tBasAlign->Branch("runNum", &runNum, "runNum/I");
		StrtBasAlign basAlign[MODNUM];
		for(i=0; i<MODNUM; i++)
		{
			sName="basAlign"+sModule[i];
			tBasAlign->Branch(sName.c_str(), &basAlign[i], "mean[32]/D:meanErr[32]/D");
		}
		
		TCanvas *cBas[MODNUM];
		for(i=0; i<MODNUM; i++)
		{
			sName="cBas"+sModule[i];
			sTitle="canvas of histograms for "+sModule[i];
			cBas[i]=new TCanvas(sName.c_str(), sTitle.c_str());
			cBas[i]->Divide(DIVMOD[i][0], DIVMOD[i][1]);
			tBasAlign->Branch(sName.c_str(), "TCanvas", &cBas[i]);			
		}
		
		TH1F *hMod[MODNUM][32];
		Option_t *OPhMod;
		TCut cutMod;
		string shMod, sOPhMod, sCutMod;
		
		for(i=0; i<MODNUM; i++)
			for(j=0; j<NCHAN[i]; j++)
			{
				k=-1;
				switch(i)
				{
					case 1:
						if(j<8)
							k=j+1;
						else
							k=2*j+1;
						break;
					case 2:
						k=2*j+1;
						break;
					default:
						k=j;
				}
				shMod="h"+sModule[i]+"_"+to_string(k);
				hMod[i][j]=new TH1F(shMod.c_str(), shMod.c_str(), RNGMOD[i][1]-RNGMOD[i][0], RNGMOD[i][0], RNGMOD[i][1]);
				tBasAlign->Branch(shMod.c_str(), "TH1F", &hMod[i][j]);
			}
			
		for(runNum=runMin; runNum<=runMax; runNum++)
		{
			ssRun.str("");
			ssRun<<setw(4)<<setfill('0')<<runNum;
			
			sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
			printf("\n**********Now processing %s!**********\n\n", sRoot.c_str());
			
			memset(basAlign, 0, sizeof(basAlign));
			for(i=0; i<MODNUM; i++)
			{
				cBas[i]->Clear("D");
				for(j=0; j<NCHAN[i]; j++)
					hMod[i][j]->Reset();
			}
			
			fRoot=new TFile(sRoot.c_str());
			if(fRoot->IsZombie())
			{
				cout<<"Error in opening "<<sRoot<<"!\n";
				continue;
			}
			fRoot->GetObject("tData",tData);
			if(!tData)
			{
				cout<<"Error read the tree of tData!\n";
				continue;
			}
			
			fBasAlign->cd();
			
			for(i=0; i<MODNUM; i++)
				for(j=0; j<NCHAN[i]; j++)
				{
					cBas[i]->cd(j+1);
					k=-1;
					switch(i)
					{
						case 1:
							if(j<8)
								k=j+1;
							else
								k=2*j+1;
							break;
						case 2:
							k=2*j+1;
							break;
						default:
							k=j;
					}
					shMod="h"+sModule[i]+"_"+to_string(k);
					sOPhMod=sRawLeaf[i]+"["+to_string(k)+"]>>"+shMod;
					OPhMod=sOPhMod.c_str();
					sCutMod=sRawLeaf[i]+"["+to_string(k)+"]>"+to_string(RNGMOD[i][0])+"&&"+sRawLeaf[i]+"["+to_string(k)+"]<"+to_string(RNGMOD[i][1]);
					cutMod=sCutMod.c_str();
					tData->Draw(OPhMod, cutMod);
					basAlign[i].mean[j]=hMod[i][j]->GetMean();
					basAlign[i].meanErr[j]=hMod[i][j]->GetMeanError();
				}
			
			tBasAlign->Fill();
			tData->Delete();			
			fRoot->Close();
		}
		fBasAlign->cd();
		tBasAlign->Write();
		tBasAlign->Delete();
		fBasAlign->Close();
	}
	
	string sAdvAlign="/home/kailong/ExpData/Jul2018/AlignData/AdvAlignRun_"+to_string(runMin)+"-"+to_string(runMax)+".root";
	string sSet[SETNUM]={"PS_270_382", "RS_270_382"};
	StrtAdvAlign advAlign[MODNUM];
	string runSet, sfSet;
	int iSet;
	ifstream fSet;
	int iRun, nSelect;
	char yn2;
	printf("Do you want to generate the advanced root files for aligning information from run#%d to run#%d? (y/n)--->", runMin, runMax);
	cin>>yn2;
	if(yn2=='y'||yn2=='Y')
	{
		fBasAlign = new TFile(sBasAlign.c_str());
		if(!(fBasAlign->IsZombie()))
		{	
			fBasAlign->GetObject("tBasAlign",tBasAlign);
			if(!(tBasAlign->IsZombie()))
			{
				fAdvAlign=new TFile(sAdvAlign.c_str(), "RECREATE");
				tAdvAlign=new TTree("tAdvAlign", "Average of each channel for particular runs");
				tAdvAlign->Branch("runSet", &runSet);
				for(i=0; i<MODNUM; i++)
				{
					sName="advAlign"+sModule[i];
					tAdvAlign->Branch(sName.c_str(), &advAlign[i], "ave[32]/D:aveErr[32]/D");
				}
				
				TCanvas *cAdv[MODNUM];
				for(i=0; i<MODNUM; i++)
				{
					sName="cAdv"+sModule[i];
					sTitle="canvas of graphs for "+sModule[i];
					cAdv[i]=new TCanvas(sName.c_str(), sTitle.c_str());
					cAdv[i]->Divide(DIVMOD[i][0], DIVMOD[i][1]);
					tAdvAlign->Branch(sName.c_str(), "TCanvas", &cAdv[i]);
				}
				
				TGraphErrors *grMod[MODNUM][32];
				Option_t *opGrMod;
				string sGrMod, sOpGrMod;	
				for(i=0; i<MODNUM; i++)
					for(j=0; j<NCHAN[i]; j++)
					{
						sGrMod="grAve"+sModule[i]+"_"+to_string(j);
						grMod[i][j]=new TGraphErrors();
						grMod[i][j]->SetNameTitle(sGrMod.c_str(), sGrMod.c_str());
						tAdvAlign->Branch(sGrMod.c_str(), "TGraphErrors", &grMod[i][j]);
					}
				
				for(iSet=0; iSet<SETNUM; iSet++)
				{
					sfSet="/home/kailong/ExpData/Jul2018/Src/runNum"+sSet[iSet]+".dat";
					printf("\n**********Now selecting the data from %s to %s! [according to %s]\n\n", sBasAlign.c_str(), sAdvAlign.c_str(), sfSet.c_str());
					fSet.open(sfSet.c_str());
					fSet>>iRun;
					sCutRunNum="runNum=="+to_string(iRun);
					while(!fSet.eof())
					{
						fSet>>iRun;
						sCutRunNum=sCutRunNum+"||runNum=="+to_string(iRun);
					}
					fSet.close();
					cutRunNum=sCutRunNum.c_str();
					
					memset(advAlign, 0, sizeof(advAlign));
					for(i=0; i<MODNUM; i++)
						cAdv[i]->Clear("D");
					
					runSet=sSet[iSet];
					
					for(i=0; i<MODNUM; i++)
						for(j=0; j<NCHAN[i]; j++)
						{
							cAdv[i]->cd(j+1);
							sOpGrMod="runNum:basAlign"+sModule[i]+".mean["+to_string(j)+"]:basAlign"+sModule[i]+".meanErr["+to_string(j)+"]";
							opGrMod=sOpGrMod.c_str();
							nSelect=tBasAlign->Draw(opGrMod, cutRunNum, "goff");
							grMod[i][j]->Set(nSelect);
							for(k=0; k<nSelect; k++)
							{
								grMod[i][j]->SetPoint(k, *(tBasAlign->GetV1()+k), *(tBasAlign->GetV2()+k));
								grMod[i][j]->SetPointError(k, 0, *(tBasAlign->GetV3()+k));
							}
							grMod[i][j]->Draw("ALP*");
							TFitResultPtr fitInfo=grMod[i][j]->Fit("pol0","SQM");
							advAlign[i].ave[j]=fitInfo->Parameter(0);
							advAlign[i].aveErr[j]=fitInfo->ParError(0);
						}
					tAdvAlign->Fill();
				}				
				tAdvAlign->Write();
				tAdvAlign->Delete();
				fAdvAlign->Close();
			}
			tBasAlign->Delete();
		}
		fBasAlign->Close();
	}
	
	string sAlign="/home/kailong/ExpData/Jul2018/AlignData/AlignRun_"+to_string(runMin)+"-"+to_string(runMax)+".root";
	char yn3;
	printf("Do you want to generate the final aligning root files from run#%d to run#%d? (y/n)--->", runMin, runMax);
	cin>>yn3;
	int runNum2;
	if(yn3=='y'||yn3=='Y')
	{
		fAdvAlign=new TFile(sAdvAlign.c_str());
		if(!(fAdvAlign->IsZombie()))
		{	
			fAdvAlign->GetObject("tAdvAlign",tAdvAlign);
			if(!(tAdvAlign->IsZombie()))
			{
				double aveMod[SETNUM][MODNUM][32];
				TCut cutSet;
				string sCutSet;
				for(iSet=0; iSet<SETNUM; iSet++)
				{
					sCutSet="runSet==\""+sSet[iSet]+"\"";
					cutSet=sCutSet.c_str();
					tAdvAlign->Draw("advAlignMadc.ave:advAlignMtdc.ave:advAlignMqdcTof.ave:advAlignMqdcMcp.ave:advAlignPin.ave", cutSet, "goff");
					for(i=0; i<MODNUM; i++)
						for(j=0; j<32; j++)
							aveMod[iSet][i][j]=*(tAdvAlign->GetVal(i)+j);
				}
				
				fBasAlign=new TFile(sBasAlign.c_str());
				fBasAlign->GetObject("tBasAlign",tBasAlign);
				
				StrtMesytec alignMstc[4];
				StrtS800 alignS800;
				string runSet2;
				TFile *fAlign=new TFile(sAlign.c_str(), "RECREATE");
				TTree *tAlign=new TTree("tAlign", "tree for aligning runs");
				tAlign->Branch("runNum", &runNum, "runNum/I");
				tAlign->Branch("runSet", &runSet2);
				for(i=0; i<4; i++)
				{
					sName="align"+sModule[i];
					tAlign->Branch(sName.c_str(), &alignMstc[i], "modID/I:data[32]/I:modRes/I:modEC_TS/I");
				}
				tAlign->Branch("alignS800", &alignS800, "tS/I:eC/I:trig/I:tof[8]/I:crdcCath[2][5][64]/I:crdcAnode[2][2]/I:hodoEgy[32]/I:hodoTime/I:pin[5]/I:mesyTDC[16]/I");
				TTree *tDrawHist=new TTree("tDrawHist", "tree for histograms of each run after alignment");
				
				tDrawHist->Branch("runNum", &runNum2, "runNum/I");
				StrtBasAlign histInfo[MODNUM];
				for(i=0; i<MODNUM; i++)
				{
					sName="hist"+sModule[i]+"Info";
					tDrawHist->Branch(sName.c_str(), &histInfo[i], "mean[32]/D:meanErr[32]/D");
				}
				
				TCanvas *cHist[MODNUM];
				for(i=0; i<MODNUM; i++)
				{
					sName="cHist"+sModule[i];
					sTitle="canvas of histograms for "+sModule[i]+" after alignment";
					cHist[i]=new TCanvas(sName.c_str(), sTitle.c_str());
					cHist[i]->Divide(DIVMOD[i][0], DIVMOD[i][1]);
					tDrawHist->Branch(sName.c_str(), "TCanvas", &cHist[i]);			
				}
				
				TH1F *hist[MODNUM][32];
				Option_t *opHist;
				TCut cutHist;
				string sHist, sOpHist, sCutHist;
				
				for(i=0; i<MODNUM; i++)
					for(j=0; j<NCHAN[i]; j++)
					{
						sHist="hist"+sModule[i]+"_"+to_string(j);
						hist[i][j]=new TH1F(sHist.c_str(), sHist.c_str(), RNGMOD[i][1]-RNGMOD[i][0], RNGMOD[i][0], RNGMOD[i][1]);
						tDrawHist->Branch(sHist.c_str(), "TH1F", &hist[i][j]);
					}
					
				TTree *tDrawGraph=new TTree("tDrawGraph", "tree for average information for particular runs after alignment");
				tDrawGraph->Branch("runSet", &runSet2);
				StrtAdvAlign graphInfo[MODNUM];
				for(i=0; i<MODNUM; i++)
				{
					sName="graph"+sModule[i]+"Info";
					tDrawGraph->Branch(sName.c_str(), &graphInfo[i], "ave[32]/D:aveErr[32]/D");
				}
				
				TCanvas *cGraph[MODNUM];
				for(i=0; i<MODNUM; i++)
				{
					sName="cGraph"+sModule[i];
					sTitle="canvas of graphs for "+sModule[i]+" after alignment";
					cGraph[i]=new TCanvas(sName.c_str(), sTitle.c_str());
					cGraph[i]->Divide(DIVMOD[i][0], DIVMOD[i][1]);
					tDrawGraph->Branch(sName.c_str(), "TCanvas", &cGraph[i]);
				}
				
				TGraphErrors *graph[MODNUM][32];
				Option_t *opGraph;
				string sGraph, sOpGraph;	
				for(i=0; i<MODNUM; i++)
					for(j=0; j<NCHAN[i]; j++)
					{
						sGraph="graph"+sModule[i]+"_"+to_string(j);
						graph[i][j]=new TGraphErrors();
						graph[i][j]->SetNameTitle(sGraph.c_str(), sGraph.c_str());
						tDrawGraph->Branch(sGraph.c_str(), "TGraphErrors", &graph[i][j]);
					}
								
				double meanModRun[MODNUM][32];
				double aveModTot[MODNUM][32];
				int nRun=runMax-runMin+1;
				long long *nEntries=new long long[nRun];
				long long *firstEntry=new long long[nRun];
				for(i=0; i<nRun; i++)
				{
					nEntries[i]=0;
					firstEntry[i]=0;
				}
				int runIndex=-1;
				for(runNum=runMin; runNum<=runMax; runNum++)
				{
					runIndex++;
					firstEntry[runIndex]=tAlign->GetEntries();
					ssRun.str("");
					ssRun<<setw(4)<<setfill('0')<<runNum;
					sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
					printf("\n**********Now converting %s to %s with tree of tAlign!**********\n\n", sRoot.c_str(), sAlign.c_str());
					
					sCutRunNum="runNum=="+to_string(runNum);
					cutRunNum=sCutRunNum.c_str();
					tBasAlign->Draw("basAlignMadc.mean:basAlignMtdc.mean:basAlignMqdcTof.mean:basAlignMqdcMcp.mean:basAlignPin.mean", cutRunNum, "goff");
					for(i=0; i<MODNUM; i++)
						for(j=0; j<32; j++)
							meanModRun[i][j]=*(tBasAlign->GetVal(i)+j);
					
					for(iSet=0; iSet<SETNUM; iSet++)
					{
						sfSet="/home/kailong/ExpData/Jul2018/Src/runNum"+sSet[iSet]+".dat";
						fSet.open(sfSet.c_str());
						while(!fSet.eof())
						{
							fSet>>iRun;
							if(iRun==runNum)
								break;
						}
						fSet.close();
						if(iRun==runNum)
						{
							runSet2=sSet[iSet];
							memmove(aveModTot, aveMod[iSet], sizeof(aveMod[iSet]));
							break;
						}
						else
						{
							runSet2="other_270_382";
							memset(aveModTot, 0, sizeof(aveModTot));	
						}
						fSet.close();
					}
					
					fRoot=new TFile(sRoot.c_str());
					if(fRoot->IsZombie())
					{
						cout<<"Error in opening "<<sRoot<<"!\n";
						continue;
					}		
					fRoot->GetObject("tData",tData);
					if(!tData)
					{
						cout<<"Error read the tree of tData!\n";
						continue;
					}
					memset(mstc, 0, sizeof(mstc));
					memset(&s800, 0, sizeof(s800));
					
					for(i=0; i<4; i++)
						tData->SetBranchAddress(sRawBrh[i].c_str(), &mstc[i]);
					tData->SetBranchAddress("s800", &s800);
					
					for(iEntry=0; iEntry<tData->GetEntries(); iEntry++)
					{
						tData->GetEntry(iEntry);
						for(i=0; i<4; i++)
							alignMstc[i]=mstc[i];
						alignS800=s800;
						
						for(i=0; i<4; i++)
							for(j=0; j<NCHAN[i]; j++)
							{
								switch(i)
								{
									case 1:
										if(j<8)
											k=j+1;
										else
											k=2*j+1;
										break;
									case 2:
										k=2*j+1;
										break;
									default:
										k=j;
								}
								
								if(mstc[i].data[k]>RNGMOD[i][0]&&mstc[i].data[k]<RNGMOD[i][1])
									alignMstc[i].data[k]=alignMstc[i].data[k]+(aveModTot[i][j]-meanModRun[i][j]);
							}
						for(i=0; i<5; i++)
							if(s800.pin[i]>RNGMOD[4][0]&&s800.pin[i]<RNGMOD[4][1])
								alignS800.pin[i]=alignS800.pin[i]+(aveModTot[4][i]-meanModRun[4][i]);
						tAlign->Fill();
					}
					fRoot->Close();
					
					nEntries[runIndex]=tAlign->GetEntries()-firstEntry[runIndex];
				}
				fAlign->cd();
				tAlign->Write();
				
				runIndex=0;
				for(runNum2=runMin; runNum2<=runMax; runNum2++)
				{
					cout<<"nEntries[runIndex]="<<nEntries[runIndex]<<", firstEntry[runIndex]="<<firstEntry[runIndex]<<endl;
					ssRun.str("");
					ssRun<<setw(4)<<setfill('0')<<runNum2;
					sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
					printf("\n**********Now converting %s to %s with tree of tDrawHist!**********\n\n", sRoot.c_str(), sAlign.c_str());
					memset(histInfo, 0, sizeof(histInfo));
					for(i=0; i<MODNUM; i++)
					{
						cHist[i]->Clear("D");
						for(j=0; j<NCHAN[i]; j++)
							hist[i][j]->Reset();
					}
					
					for(i=0; i<MODNUM; i++)
						for(j=0; j<NCHAN[i]; j++)
						{
							cHist[i]->cd(j+1);
							k=-1;
							switch(i)
							{
								case 1:
									if(j<8)
										k=j+1;
									else
										k=2*j+1;
									break;
								case 2:
									k=2*j+1;
									break;
								default:
									k=j;
							}
							sHist="hist"+sModule[i]+"_"+to_string(j);
							sOpHist=sAlignRawLeaf[i]+"["+to_string(k)+"]>>"+sHist;
							opHist=sOpHist.c_str();
							sCutHist="runNum=="+to_string(runNum2);
							cutHist=sCutHist.c_str();
							tAlign->Draw(opHist, cutHist, "", nEntries[runIndex], firstEntry[runIndex]);
							histInfo[i].mean[j]=hist[i][j]->GetMean();
							histInfo[i].meanErr[j]=hist[i][j]->GetMeanError();
						}
					tDrawHist->Fill();
					runIndex++;
				}
				tDrawHist->Write();
				
				for(iSet=0; iSet<SETNUM; iSet++)
				{
					sfSet="/home/kailong/ExpData/Jul2018/Src/runNum"+sSet[iSet]+".dat";
					printf("\n**********Now generating the tree of tDrawGraph in %s! [according to %s]\n\n", sAlign.c_str(), sfSet.c_str());
					fSet.open(sfSet.c_str());
					fSet>>iRun;
					sCutRunNum="runNum=="+to_string(iRun);
					while(!fSet.eof())
					{
						fSet>>iRun;
						sCutRunNum=sCutRunNum+"||runNum=="+to_string(iRun);
					}
					fSet.close();
					cutRunNum=sCutRunNum.c_str();
					
					memset(graphInfo, 0, sizeof(graphInfo));
					for(i=0; i<MODNUM; i++)
						cGraph[i]->Clear("D");
					
					runSet2=sSet[iSet];
					
					for(i=0; i<MODNUM; i++)
						for(j=0; j<NCHAN[i]; j++)
						{
							cGraph[i]->cd(j+1);
							sOpGraph="runNum:hist"+sModule[i]+"Info.mean["+to_string(j)+"]:hist"+sModule[i]+"Info.meanErr["+to_string(j)+"]";
							opGraph=sOpGraph.c_str();
							nSelect=tDrawHist->Draw(opGraph, cutRunNum, "goff");
							graph[i][j]->Set(nSelect);
							for(k=0; k<nSelect; k++)
							{
								graph[i][j]->SetPoint(k, *(tDrawHist->GetV1()+k), *(tDrawHist->GetV2()+k));
								graph[i][j]->SetPointError(k, 0, *(tDrawHist->GetV3()+k));
							}
							graph[i][j]->Draw("ALP*");
							TFitResultPtr fitGraph=graph[i][j]->Fit("pol0","SQM");
							graphInfo[i].ave[j]=fitGraph->Parameter(0);
							graphInfo[i].aveErr[j]=fitGraph->ParError(0);
						}
					tDrawGraph->Fill();
				}
				
				tDrawGraph->Write();
				tDrawGraph->Delete();
				tDrawHist->Delete();
				tAlign->Delete();
				fAlign->Close();
				tBasAlign->Delete();
				fBasAlign->Close();
			}
			tAdvAlign->Delete();
		}
		fAdvAlign->Close();
	}
	
}

#ifndef __CINT__

// void StandaloneApplication(int argc, char** argv)
// {
	// Root2Align();
// }

// int main(int argc, char** argv)
// {
	// TApplication app("ROOT Application", &argc, argv);
	// StandaloneApplication(app.Argc(), app.Argv());
	// app.Run();
	// return 0;
// }

int main()
{
	Root2Align();
	return 0;
}
#endif
