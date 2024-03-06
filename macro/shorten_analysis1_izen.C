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

int shorten_analysis1_izen(int runid)
{
	string strrun = "run"+to_string(runid) ; 
	string prefix = "IZEN_COLIS_MID";
	string analysis_filename = "/local/home/rb277364/Projects/tomomu/izen/mnt/"+prefix+"_"+strrun+"_analyse.root";
	string light_analysis_filename = "/local/home/rb277364/Projects/tomomu/Data/izen/" + strrun + "/"+prefix+"_"+strrun+"_analyse1.root";

	//TFile *fout = new TFile(light_analysis_filename.c_str(),"RECREATE");
	//TFile fout(light_analysis_filename.c_str(), "");
    //TTree *Tout = new TTree("T","event");

	TFile f(analysis_filename.c_str(), "");
    TTree *T = (TTree*)f.Get("T");

	int evn;
	double evttime;
	int nlay=24;
	int MGv3_NClus[nlay];
	double MGv3_StripMaxAmpl[nlay][61];

	double MGv3_ClusAmpl[nlay][300];
	double MGv3_ClusSize[nlay][300];
	double MGv3_ClusPos[nlay][300];
	double MGv3_ClusMaxStripAmpl[nlay][300];
	double MGv3_ClusMaxSample[nlay][300];
	double MGv3_ClusTOT[nlay][300];
	double MGv3_ClusT[nlay][300];
	int MGv3_ClusMaxStrip[nlay][300];

	double MGv3_ClusAmpl_s[nlay][24];
	double MGv3_ClusSize_s[nlay][24];
	double MGv3_ClusPos_s[nlay][24];
	double MGv3_ClusMaxStripAmpl_s[nlay][24];
	double MGv3_ClusMaxSample_s[nlay][24];
	double MGv3_ClusTOT_s[nlay][24];
	double MGv3_ClusT_s[nlay][24];
	int MGv3_ClusMaxStrip_s[nlay][24];



	// Link branches element to variables
	T->SetBranchAddress("evn",&evn);
	T->SetBranchAddress("evttime",&evttime);
	T->SetBranchAddress("MGv3_NClus",MGv3_NClus);
	T->SetBranchAddress("MGv3_ClusAmpl",MGv3_ClusAmpl);
	T->SetBranchAddress("MGv3_ClusSize",MGv3_ClusSize);
	T->SetBranchAddress("MGv3_ClusPos",MGv3_ClusPos);
	T->SetBranchAddress("MGv3_ClusMaxStripAmpl",MGv3_ClusMaxStripAmpl);
	T->SetBranchAddress("MGv3_ClusMaxSample",MGv3_ClusMaxSample);
	T->SetBranchAddress("MGv3_ClusTOT",MGv3_ClusTOT);
	T->SetBranchAddress("MGv3_ClusT",MGv3_ClusT);
	T->SetBranchAddress("MGv3_StripMaxAmpl",MGv3_StripMaxAmpl);
	T->SetBranchAddress("MGv3_ClusMaxStrip",MGv3_ClusMaxStrip);


	// Link branches elements to variables

	TFile *fout = new TFile(light_analysis_filename.c_str(),"RECREATE");
    TTree *Tout = new TTree("T","event");
    
	Tout->Branch("evn", &evn, "evn/I");
	
	Tout->Branch("evttime", &evttime, "evttime/D");
	Tout->Branch("MGv3_NClus", MGv3_NClus, "MGv3_NClus[24]/I");
	Tout->Branch("MGv3_StripMaxAmpl", MGv3_StripMaxAmpl, "MGv3_StripMaxAmpl[24][61]/D");

	Tout->Branch("MGv3_ClusAmpl", MGv3_ClusAmpl_s, "MGv3_ClusAmpl[24][8]/D");
	Tout->Branch("MGv3_ClusSize", MGv3_ClusSize_s, "MGv3_ClusSize[24][8]/D");
	Tout->Branch("MGv3_ClusPos", MGv3_ClusPos_s, "MGv3_ClusPos[24][8]/D");
	Tout->Branch("MGv3_ClusMaxStripAmpl", MGv3_ClusMaxStripAmpl_s, "MGv3_ClusMaxStripAmpl[24][8]/D");
	Tout->Branch("MGv3_ClusMaxSample", MGv3_ClusMaxSample_s, "MGv3_ClusMaxSample[24][8]/D");
	Tout->Branch("MGv3_ClusTOT", MGv3_ClusTOT_s, "MGv3_ClusTOT[24][8]/D");
	Tout->Branch("MGv3_ClusT", MGv3_ClusT_s, "MGv3_ClusT[24][8]/D");
	Tout->Branch("MGv3_ClusMaxStrip", MGv3_ClusMaxStrip, "MGv3_ClusMaxStrip[24][8]/I");


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
				MGv3_ClusAmpl_s[j][k] = MGv3_ClusAmpl[j][k];
				MGv3_ClusSize_s[j][k] = MGv3_ClusSize[j][k];
				MGv3_ClusPos_s[j][k] = MGv3_ClusPos[j][k];
				MGv3_ClusMaxStripAmpl_s[j][k] = MGv3_ClusMaxStripAmpl[j][k];
				MGv3_ClusMaxSample_s[j][k] = MGv3_ClusMaxSample[j][k];
				MGv3_ClusTOT_s[j][k] = MGv3_ClusTOT[j][k];
				MGv3_ClusT_s[j][k] = MGv3_ClusT[j][k]; 
			}
		}

		Tout->Fill();

	}

	Tout->Write();
	
	fout->Close();

	return 0;
}
