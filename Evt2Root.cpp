#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <unistd.h>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"

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
	int trig;  // =1: Coincidence (usual case); =16: Secondary (MCP provides)
	int tof[8]; //[0]-->[1]: ORTEC TAC+Phillips ADC; [2]-->[7]: Phillips TDC
	int crdcCath[2][5][64]; //[2]: two CRDC; [5]: [0]---sample; [1]-->[4]---energy
	int crdcAnode[2][2]; // [0]: energy; [1]: time	
	int hodoEgy[32]; // 32 crystals
	int hodoTime;
	int pin[5];
	int mesyTDC[16];  //[0]: time from S800; [2]: time from A1900
};

void Evt2Root()
{	
	string sEvt, sRoot;
	unsigned char ringHead[4];
	unsigned int evtWord;
	int run[2];
	StrtMesytec madc, mtdc, mqdcTOF, mqdcMCP;
	StrtS800 s800;
	int ringSize, ringType, bHSig, srcID;
	int nWords, iCh, iConnect, iCRDC, iHodo;
	int i, iCtrlCRDC, iEvt;
	unsigned char buf[8192];
	unsigned char *pBuf;
	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	ostringstream ssRun1, ssRun2, ssIEvt;	
	for(runNum=runMin; runNum<=runMax; runNum++)
	{
		ssRun1.str("");
		ssRun1<<runNum;
		ssRun2.str("");
		ssRun2<<setw(4)<<setfill('0')<<runNum;
		for(iEvt=0; iEvt<=15; iEvt++)
		{
			ssIEvt.str("");
			ssIEvt<<setw(2)<<setfill('0')<<iEvt;
			
			// sEvt="/disk2/e12022/data/experiment/run"+ssRun1.str()+"/run-"+ssRun2.str()+"-"+ssIEvt.str()+".evt";
			sEvt="/home/kailong/ExpData/Jul2018/EvtData/run"+ssRun1.str()+"/run-"+ssRun2.str()+"-"+ssIEvt.str()+".evt";
			if(access("/home/kailong/ExpData/Jul2018/RootData", F_OK)!=0)
				system("mkdir /home/kailong/ExpData/Jul2018/RootData");
			sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun2.str()+"-"+ssIEvt.str()+".root";
			printf("\n**********Now converting %s to %s!**********\n\n", sEvt.c_str(), sRoot.c_str());
			
			ifstream fEvt(sEvt.c_str(), ios::binary);
			if(!fEvt.is_open())
			{
				cout<<"Error in Opening "<<sEvt<<endl;
				continue;
			}
			TFile *fRoot = new TFile(sRoot.c_str(), "RECREATE");
			TTree *tData = new TTree("tData", "tree of data");
			tData->Branch("run", run, "run[2]/I");  //[0]: run number;  [1]: sub-number of evt
			tData->Branch("madc", &madc, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("mtdc", &mtdc, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("mqdcTOF", &mqdcTOF, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("mqdcMCP", &mqdcMCP, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("s800", &s800, "tS/I:eC/I:trig/I:tof[8]/I:crdcCath[2][5][64]/I:crdcAnode[2][2]/I:hodoEgy[32]/I:hodoTime/I:pin[5]/I:mesyTDC[16]/I");
			int j=0;
			run[0]=runNum;
			run[1]=iEvt;
			while(!fEvt.eof())
			{
				j++;
				//cout<<"j= "<<j<<", ";
				memset(&madc, 0, sizeof(madc));
				memset(&mtdc, 0, sizeof(mtdc));
				memset(&mqdcTOF, 0, sizeof(mqdcTOF));
				memset(&mqdcMCP, 0, sizeof(mqdcMCP));
				memset(&s800, 0, sizeof(s800));
				
				fEvt.read((char*)ringHead, sizeof(ringHead));
				ringSize=ringHead[0] | ringHead[1]<<8 | ringHead[2]<<16 | ringHead[3]<<24;
				if(ringSize>8192)
					continue;
				
				fEvt.read((char*)buf, ringSize-4);
				pBuf=buf;			
				
				ringType=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
				pBuf+=4;
				
				if(ringType==30)
				{
					// cout<<"ringSize= "<<ringSize<<endl;
					// cout<<"ringType= "<<ringType<<endl;
					bHSig=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
					pBuf+=4;
					if(bHSig==20)
						pBuf+=16;
					pBuf+=12;
					srcID=*pBuf | *(pBuf+1)<<8;
					if(srcID==2)
						continue;				
					
					pBuf+=42;
				
					//Decode ADC data
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
								
					if(evtWord!=0xFFFFFFFF)
					{				
						pBuf+=4;
						madc.modID=(evtWord>>16)&0x1F;
					
						madc.modRes=(evtWord>>12)&0x7;
						nWords=evtWord&0xFFF;
						
						//printf("ADC Header: %x\n", evtWord);
						
						for(i=1; i<=nWords-1; i++)
						{
							evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
							pBuf+=4;
							iCh=(evtWord>>16)&0x1F;
							madc.data[iCh]=evtWord&0x1FFF;
							// printf("adc[%d]=%d\n", iCh, madc.data[iCh]);
						}
						
						evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
						pBuf+=4;
						madc.modEC_TS=evtWord&0x3FFFFFFF;
					}
						
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;				
					if(evtWord!=0xFFFFFFFF)
						continue;				
					pBuf+=4;
					
					//Decode TDC data
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
					if(evtWord!=0xFFFFFFFF)
					{	
						pBuf+=4;				

						mtdc.modID=(evtWord>>16)&0xFF;
						
						mtdc.modRes=(evtWord>>12)&0xF;
						nWords=evtWord&0xFFF;
						
						//printf("TDC Header: %x\n", evtWord);
						
						for(i=1; i<=nWords-1; i++)
						{
							evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
							pBuf+=4;
							iCh=(evtWord>>16)&0x3F;
							if(iCh>=0&&iCh<32)
							{
								mtdc.data[iCh]=evtWord&0xFFFF;
								// printf("tdc[%d]=%d\n", iCh, mtdc.data[iCh]);
							}
						}
						
						evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
						pBuf+=4;
						mtdc.modEC_TS=evtWord&0x3FFFFFFF;
						//printf("TDC End: %x\n", evtWord);
					}
					
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;				
					if(evtWord!=0xFFFFFFFF)
						continue;
					pBuf+=4;
					
					//Decode QDC data of TOF
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;				
					if(evtWord!=0xFFFFFFFF)
					{
						pBuf+=4;

						mqdcTOF.modID=(evtWord>>16)&0xFF;
						
						nWords=evtWord&0xFFF;
						
						//printf("QDC Header: %x\n", evtWord);
						
						for(i=1; i<=nWords-1; i++)
						{
							evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
							pBuf+=4;
							iCh=(evtWord>>16)&0x1F;
							mqdcTOF.data[iCh]=evtWord&0xFFF;
							//printf("qdcTOF[%d]=%d\n", iCh, mqdcTOF.data[iCh]);
						}
						
						evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
						pBuf+=4;
						mqdcTOF.modEC_TS=evtWord&0x3FFFFFFF;
						//printf("QDC End: %x\n", evtWord);
					}
					
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;				
					if(evtWord!=0xFFFFFFFF)
						continue;
					pBuf+=4;
					
					//Decode QDC data of MCP
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
					if(evtWord!=0xFFFFFFFF)
					{
						pBuf+=4;

						mqdcMCP.modID=(evtWord>>16)&0xFF;
						
						nWords=evtWord&0xFFF;
						
						//printf("QDC Header: %x\n", evtWord);
						
						for(i=1; i<=nWords-1; i++)
						{
							evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
							pBuf+=4;
							iCh=(evtWord>>16)&0x1F;
							mqdcMCP.data[iCh]=evtWord&0xFFF;
							//printf("qdcMCP[%d]=%d\n", iCh, mqdcMCP.data[iCh]);
						}
						
						evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
						pBuf+=4;
						mqdcMCP.modEC_TS=evtWord&0x3FFFFFFF;
						//printf("QDC End: %x\n", evtWord);
					}
					
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;				
					if(evtWord!=0xFFFFFFFF)
						continue;

					pBuf+=4;
					
					if(pBuf-buf>=ringSize-4)
						continue;
					
					pBuf+=60;
					if(pBuf-buf>=ringSize-4)
						continue;
					evtWord=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
					s800.tS=evtWord&0x3FFFFFFF;
					
					pBuf+=12;
					if(pBuf-buf>=ringSize-4)
						continue;
					s800.eC=*pBuf | *(pBuf+1)<<8 | *(pBuf+2)<<16 | *(pBuf+3)<<24;
					// cout<<"s800.eC="<<s800.eC<<endl;
					
					pBuf+=6;
					if(pBuf-buf>=ringSize-4)
						continue;
					nWords=*pBuf | *(pBuf+1)<<8;
					// cout<<"triggernWords="<<nWords<<endl;
					pBuf+=4;
					if(pBuf-buf>=ringSize-4)
						continue;
					s800.trig=*pBuf | *(pBuf+1)<<8;
					pBuf+=(2*nWords-4);  //skip trigger packet 0X5801
					if(pBuf-buf>=ringSize-4)
						continue;
					
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=4;  //skip 0X000* and 0X5802
					// cout<<"nWords="<<nWords<<endl;
					if(pBuf-buf>=ringSize-4)
						continue;
					
					for(i=0; i<nWords-2; i++)
					{
						evtWord=*pBuf | *(pBuf+1)<<8;
						pBuf+=2; //go to next word
						iCh=(evtWord>>12)&0xF;
						s800.tof[iCh-4]=evtWord&0xFFF;
						// cout<<"s800.tof[iCh-4]="<<s800.tof[iCh-4]<<endl;
					}
					
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=(2*nWords); //skip scintillator packet 0X5810
					if(pBuf-buf>=ringSize-4)
						continue;
					
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=(2*nWords); //skip ion-chamber packet 0x5820
					// cout<<"ionchambernWords="<<nWords<<endl;
					if(pBuf-buf>=ringSize-4)
						continue;
					
					for(iCRDC=0; iCRDC<2; iCRDC++)
					{
						pBuf+=6; //skip head of CRDC packet
						if(pBuf-buf>=ringSize-4)
							continue;
						nWords=*pBuf | *(pBuf+1)<<8;
						pBuf+=6;
						if(pBuf-buf>=ringSize-4)
							continue;
						// cout<<"CRDCnWords="<<nWords<<endl;
						iCh=-1;
						for(i=0; i<nWords-3; i++)
						{
							evtWord=*pBuf | *(pBuf+1)<<8;
							pBuf+=2;
							iCtrlCRDC=(evtWord>>15)&0x1;
							if(iCtrlCRDC==1)
							{
								iCh=evtWord&0x3F;
								s800.crdcCath[iCRDC][0][iCh]=(evtWord>>6)&0x1FF;
							}
							if(iCtrlCRDC==0&&iCh!=-1)
							{
								iConnect=1+((evtWord>>10)&0x3);
								s800.crdcCath[iCRDC][iConnect][iCh]=evtWord&0x3FF;
							}
						}
						
						pBuf+=4; //skip the head of CRDC Anode sub-packet
						if(pBuf-buf>=ringSize-4)
							continue;
						s800.crdcAnode[iCRDC][0]=*pBuf | *(pBuf+1)<<8;
						pBuf+=2;
						if(pBuf-buf>=ringSize-4)
							continue;
						s800.crdcAnode[iCRDC][1]=*pBuf | *(pBuf+1)<<8;
						pBuf+=2;
						if(pBuf-buf>=ringSize-4)
							continue;
					}
					
					for(iHodo=0; iHodo<2; iHodo++)
					{
						nWords=*pBuf | *(pBuf+1)<<8;
						pBuf+=6;
						if(pBuf-buf>=ringSize-4)
							continue;
						// cout<<"HODOWords="<<nWords<<endl;
						for(i=0; i<nWords-3; i++)
						{
							evtWord=*pBuf | *(pBuf+1)<<8;
							pBuf+=2;
							iCh=(evtWord>>12)&0xF;
							s800.hodoEgy[iCh+16*iHodo]=evtWord&0xFFF;
						}
					}
					pBuf+=10; //skip head of hit-pattern sub-packet
					if(pBuf-buf>=ringSize-4)
						continue;
					evtWord=*pBuf | *(pBuf+1)<<8;
					s800.hodoTime=evtWord&0xFFF;
					pBuf+=2;
					if(pBuf-buf>=ringSize-4)
						continue;
					
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=(2*nWords); //skip TPPACs packet 0X5870
					if(pBuf-buf>=ringSize-4)
						continue;
					
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=(2*nWords); //skip OBJECT PIN packet 0X58a0
					if(pBuf-buf>=ringSize-4)
						continue;
					
					//process pin data: 0X5805
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=4;
					if(pBuf-buf>=ringSize-4)
						continue;
					for(i=0; i<nWords-2; i++)
					{
						evtWord=*pBuf | *(pBuf+1)<<8;
						pBuf+=2;
						iCh=(evtWord>>12)&0xF;
						s800.pin[iCh]=evtWord&0xFFF;
						// cout<<"s800.pin[iCh]="<<s800.pin[iCh]<<endl;
					}
					
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=(2*nWords); //skip Galotte packet 0X58d0
					if(pBuf-buf>=ringSize-4)
						continue;
					
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=(2*nWords); //skip LaBr packet 0X58E0
					if(pBuf-buf>=ringSize-4)
						continue;
					
					//process MesytecTDC data: 0X58F0
					nWords=*pBuf | *(pBuf+1)<<8;
					pBuf+=4;
					if(pBuf-buf>=ringSize-4)
						continue;
					// cout<<"MTDCnWords="<<nWords<<endl;
					for(i=0; i<nWords-2; i++)
					{
						evtWord=*pBuf | *(pBuf+1)<<8;
						pBuf+=2;
						iCh=evtWord&0xF;
						s800.mesyTDC[iCh]=*pBuf | *(pBuf+1)<<8;
						pBuf+=2;
					}
					
					cout<<runNum<<"--"<<j<<":  madc.modEC_TS="<<madc.modEC_TS<<",  s800.tS="<<s800.tS<<",  s800.trig="<<s800.trig<<endl;
					tData->Fill();						
				}
			}
			fRoot->Write();
			fRoot->Close();
			fEvt.close();	
		}		
	}			
}

#ifndef __CINT__
/*
void StandaloneApplication(int argc, char** argv)
{
	Evt2Root();
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
	Evt2Root();
	return 0;
}
#endif
