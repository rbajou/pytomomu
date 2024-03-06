#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TBranch.h"
#include "TMath.h"
#include <string>

// string strrun = "run9";
// string analysis_filename = "/local/home/rb277364/Projects/tomomu/inb72/mnt/INB72_2024_"+strrun+"_analyse.root";
// string light_analysis_filename = "/local/home/rb277364/Data/tomomu/inb72/"+strrun+"/INB72_2024_"+strrun+"_analyse1.root";

int shorten_analysis1_inb72(int runid)
{
	string strrun = "run"+to_string(runid) ; 
	string analysis_filename = "/local/home/rb277364/Projects/tomomu/inb72/mnt/INB72_2024_"+strrun+"_analyse.root";
	string light_analysis_filename = "/local/home/rb277364/Projects/tomomu/Data/inb72/" + strrun + "/INB72_2024_"+strrun+"_analyse1.root";

	//TFile *fout = new TFile(light_analysis_filename.c_str(),"RECREATE");
	//TFile fout(light_analysis_filename.c_str(), "");
    //TTree *Tout = new TTree("T","event");

	TFile f(analysis_filename.c_str(), "");
    TTree *T = (TTree*)f.Get("T");

	int evn;
	double evttime;
	int nlay=8;
	int INB72_NClus[nlay];
	double INB72_StripMaxAmpl[nlay][61];

	double INB72_ClusAmpl[nlay][300];
	double INB72_ClusSize[nlay][300];
	double INB72_ClusPos[nlay][300];
	double INB72_ClusMaxStripAmpl[nlay][300];
	double INB72_ClusMaxSample[nlay][300];
	double INB72_ClusTOT[nlay][300];
	double INB72_ClusT[nlay][300];
	int INB72_ClusMaxStrip[nlay][300];

	double INB72_ClusAmpl_s[nlay][8];
	double INB72_ClusSize_s[nlay][8];
	double INB72_ClusPos_s[nlay][8];
	double INB72_ClusMaxStripAmpl_s[nlay][8];
	double INB72_ClusMaxSample_s[nlay][8];
	double INB72_ClusTOT_s[nlay][8];
	double INB72_ClusT_s[nlay][8];
	int INB72_ClusMaxStrip_s[nlay][8];



	// Link branches element to variables
	T->SetBranchAddress("evn",&evn);
	T->SetBranchAddress("evttime",&evttime);
	T->SetBranchAddress("INB72_NClus",INB72_NClus);
	T->SetBranchAddress("INB72_ClusAmpl",INB72_ClusAmpl);
	T->SetBranchAddress("INB72_ClusSize",INB72_ClusSize);
	T->SetBranchAddress("INB72_ClusPos",INB72_ClusPos);
	T->SetBranchAddress("INB72_ClusMaxStripAmpl",INB72_ClusMaxStripAmpl);
	T->SetBranchAddress("INB72_ClusMaxSample",INB72_ClusMaxSample);
	T->SetBranchAddress("INB72_ClusTOT",INB72_ClusTOT);
	T->SetBranchAddress("INB72_ClusT",INB72_ClusT);
	T->SetBranchAddress("INB72_StripMaxAmpl",INB72_StripMaxAmpl);
	T->SetBranchAddress("INB72_ClusMaxStrip",INB72_ClusMaxStrip);


	// Link branches elements to variables

	TFile *fout = new TFile(light_analysis_filename.c_str(),"RECREATE");
    TTree *Tout = new TTree("T","event");
    
	Tout->Branch("evn", &evn, "evn/I");
	
	Tout->Branch("evttime", &evttime, "evttime/D");
	Tout->Branch("INB72_NClus", INB72_NClus, "INB72_NClus[8]/I");
	Tout->Branch("INB72_StripMaxAmpl", INB72_StripMaxAmpl, "INB72_StripMaxAmpl[8][61]/D");

	Tout->Branch("INB72_ClusAmpl", INB72_ClusAmpl_s, "INB72_ClusAmpl[8][8]/D");
	Tout->Branch("INB72_ClusSize", INB72_ClusSize_s, "INB72_ClusSize[8][8]/D");
	Tout->Branch("INB72_ClusPos", INB72_ClusPos_s, "INB72_ClusPos[8][8]/D");
	Tout->Branch("INB72_ClusMaxStripAmpl", INB72_ClusMaxStripAmpl_s, "INB72_ClusMaxStripAmpl[8][8]/D");
	Tout->Branch("INB72_ClusMaxSample", INB72_ClusMaxSample_s, "INB72_ClusMaxSample[8][8]/D");
	Tout->Branch("INB72_ClusTOT", INB72_ClusTOT_s, "INB72_ClusTOT[8][8]/D");
	Tout->Branch("INB72_ClusT", INB72_ClusT_s, "INB72_ClusT[8][8]/D");
	Tout->Branch("INB72_ClusMaxStrip", INB72_ClusMaxStrip, "INB72_ClusMaxStrip[8][8]/I");


	Long64_t nentries = T->GetEntries();


	for(Long64_t i=0; i<nentries;i++)
	{
		if(i%1000==0) std::cout << "event " << i << "/"<<nentries << std::endl;
		if(i%50000==0) Tout->Write();
		

		T->GetEntry(i);

		for(int j=0; j<nlay; j++)
		{
			for(int k=0; k<8; k++)
			{
				INB72_ClusAmpl_s[j][k] = INB72_ClusAmpl[j][k];
				INB72_ClusSize_s[j][k] = INB72_ClusSize[j][k];
				INB72_ClusPos_s[j][k] = INB72_ClusPos[j][k];
				INB72_ClusMaxStripAmpl_s[j][k] = INB72_ClusMaxStripAmpl[j][k];
				INB72_ClusMaxSample_s[j][k] = INB72_ClusMaxSample[j][k];
				INB72_ClusTOT_s[j][k] = INB72_ClusTOT[j][k];
				INB72_ClusT_s[j][k] = INB72_ClusT[j][k]; 
			}
		}

		Tout->Fill();

	}

	Tout->Write();
	
	fout->Close();

	return 0;
}
