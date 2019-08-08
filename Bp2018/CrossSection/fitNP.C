#include "uti.h"
//#include "parameters.h"
#include "TF1.h"

const int nBins=1;
double ptBins[nBins+1] = {5.,100.};

Double_t minhisto=5.0;
Double_t maxhisto=6.0;
Double_t nbinsmasshisto=50;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

TString weight;
TString weightgen;
TString seldata;
TString selmc;
TString selmceff;
TString selmcgen;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;

Double_t npyield;
Double_t npyieldErr;

void fitNP(int usePbPb=1, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=1, Double_t luminosity=1., int doweight=0, TString collsyst="", TString outputfile="ROOTfiles/NPFitPbPb.root", Float_t centmin=0., Float_t centmax=0.)
{
	collisionsystem=collsyst;
	hiBinMin = centmin*2;
	hiBinMax = centmax*2;
	centMin = centmin;
	centMax = centmax;

	double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing);

	if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
	bool isPbPb=(bool)(usePbPb);

	if(!isPbPb)
	{
		seldata = Form("%s&&%s",trgselection.Data(),cut.Data());
		selmceff = Form("%s&&%s",trgselection.Data(),cut.Data());
		selmcgen = Form("%s",cutmcgen.Data());
	}
	else
	{
		seldata = Form("%s&&%s&&hiBin>%f&&hiBin<%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);
		selmceff = Form("%s&&%s&&hiBin>%f&&hiBin<%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);
		selmcgen = Form("%s&&hiBin>%f&&hiBin<%f",cutmcgen.Data(),hiBinMin,hiBinMax);
	}

	selmc = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);

	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadRightMargin(0.043);
	gStyle->SetPadLeftMargin(0.18);
	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadBottomMargin(0.145);
	gStyle->SetTitleX(.0f);

	void clean0 (TH1D* h);
	TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centmin, Float_t centmax);

	doweight=0;

	if(doweight==0) {
		weightgen="1";
		weight="1";
	}
	if(doweight<0 || doweight>=1) std::cout<<"ERROR, this weighting option is not defined"<<std::endl;

	std::cout<<"we are using weight="<<weight<<std::endl;

	TFile* inf = new TFile(inputdata.Data());

	TTree* nt = (TTree*)inf->Get("Bfinder/ntKp");
	nt->AddFriend("hltanalysis/HltTree");
	nt->AddFriend("hiEvtAnalyzer/HiTree");
	nt->AddFriend("skimanalysis/HltTree");
	nt->AddFriend("Bfinder/ntGen");

	TTree* ntMC = (TTree*)inf->Get("Bfinder/ntKp");
	ntMC->AddFriend("hltanalysis/HltTree");
	ntMC->AddFriend("hiEvtAnalyzer/HiTree");
	ntMC->AddFriend("skimanalysis/HltTree");
	ntMC->AddFriend("Bfinder/ntGen");

	TF1 *totalmass;

	TFile* outf = new TFile(outputfile.Data(),"recreate");
	for(int i=0;i<nBins;i++)
	{
		TF1* f = fit(nt,nt,ptBins[i],ptBins[i+1],isMC,isPbPb, totalmass,centmin, centmax);
		//double yield = f->Integral(minhisto,maxhisto)/binwidthmass;
		//double yieldErr = f->Integral(minhisto,maxhisto)/binwidthmass*f->GetParError(0)/f->GetParameter(0);
	}  
}

void clean0(TH1D* h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		if(h->GetBinContent(i)==0) h->SetBinError(i,1);
	}
}

TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC,bool isPbPb,TF1* &total,Float_t centmin, Float_t centmax)
{
	//cout<<cut.Data()<<endl;
	static Int_t count=0;
	count++;
	TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
	TH1D* h = new TH1D(Form("h%d",count),"",nbinsmasshisto,minhisto,maxhisto);

	//TF1* f = new TF1(Form("f%d",count),"[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [3]*Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5]) + [6]*Gaus(x,[7],[8])/(sqrt(2*3.14159)*[8]) + [9]+[10]*x ");
	//TF1* f = new TF1(Form("f%d",count),"[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [3]*Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5]) + [9]+[10]*x ");
	//TF1* f = new TF1(Form("f%d",count),"[0]*TMath::Erf((x-[1])/[2])+[0]+[9]+[10]*x");
	TF1* f = new TF1(Form("f%d",count),"[0]*TMath::Erf((x-[1])/[2])+[0]+[6]*TMath::Gaus(x,[7],[8])/(sqrt(2*3.14159)*[8])+[9]+[10]*x");
	
	if(isMC==1) nt->Project(Form("h%d",count),"Bmass",Form("%s*(%s&&Bpt>%f&&Bpt<%f)","1",seldata.Data(),ptmin,ptmax));   
	else nt->Project(Form("h%d",count),"Bmass",Form("(%s&&Bpt>%f&&Bpt<%f)",seldata.Data(),ptmin,ptmax));   

	clean0(h);

	h->Draw();

	//error fn
	f->SetParLimits(0, 1, 1e3);
	f->SetParLimits(1, 4.9, 5.3);
	f->SetParLimits(2, -0.5, 0);
	f->SetParameter(0,20);
	f->SetParameter(1,5.14);
	f->SetParameter(2,-0.1);
	
	f->SetParLimits(3, 1e-2, 1e4);
	f->SetParLimits(4, 5.06, 5.10);
	f->SetParLimits(5, 0.001, 0.1);
	f->SetParameter(3,1e2);
	f->SetParameter(4,5.07);
	f->SetParameter(5,0.05);

	f->SetParLimits(6, 1e-1, 1e2);
	f->SetParLimits(7, 5.32, 5.36);
	f->SetParLimits(8, 0.005, 0.04);
	f->SetParameter(6,10);
	f->SetParameter(7,5.35);
	f->SetParameter(8,0.04);

	f->SetParLimits(9, 0, 1e5);
	f->SetParLimits(10, -500,  100);
	f->SetParameter(9,1e3);
	f->SetParameter(10,-1);

	h->GetEntries();

	h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
	h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
	h->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);
	h->SetMarkerSize(0.8);
	h->SetMarkerStyle(20);

	TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x");
	background->SetParameter(0,f->GetParameter(9));
	background->SetParameter(1,f->GetParameter(10));
	background->SetLineColor(4);
	background->SetRange(minhisto,maxhisto);
	background->SetLineStyle(2);

	//TF1 *mass = new TF1(Form("fmass%d",count),"[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [3]*Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5]) + [6]*Gaus(x,[7],[8])/(sqrt(2*3.14159)*[8])");
	//mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4),f->GetParameter(5),f->GetParameter(6),f->GetParameter(7),f->GetParameter(8));
	//TF1 *mass = new TF1(Form("fmass%d",count),"[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [3]*Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5])");
	//mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4),f->GetParameter(5));
	//TF1 *mass = new TF1(Form("fmass%d",count),"[0]*TMath::Erf((x-[1])/[2]) + [0]");
	//mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2));
	TF1 *mass = new TF1(Form("fmass%d",count),"[0]*TMath::Erf((x-[1])/[2])+[0]+[3]*TMath::Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5])");
	mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(6),f->GetParameter(7),f->GetParameter(8));
	mass->SetParError(0,f->GetParError(0));
	mass->SetParError(1,f->GetParError(1));
	mass->SetParError(2,f->GetParError(2));
	mass->SetParError(3,f->GetParError(6));
	mass->SetParError(4,f->GetParError(7));
	mass->SetParError(5,f->GetParError(8));
	mass->SetLineColor(2);

	h->SetXTitle("m_{#mu#muK} (GeV/c^{2})");
	double yunit = 1000./nbinsmasshisto;
	h->SetYTitle(Form("Events / (%.0f MeV/c^{2})",yunit));//basically, the unit is equivalent to binwidthmass in MeV
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
	h->GetXaxis()->SetTitleOffset(1.3);
	h->GetYaxis()->SetTitleOffset(1.8);
	h->GetXaxis()->SetLabelOffset(0.007);
	h->GetYaxis()->SetLabelOffset(0.007);
	h->GetXaxis()->SetTitleSize(0.045);
	h->GetYaxis()->SetTitleSize(0.045);
	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleFont(42);
	h->GetXaxis()->SetLabelFont(42);
	h->GetYaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetYaxis()->SetLabelSize(0.04);
	h->SetMarkerSize(0.8);
	h->SetMarkerStyle(20);
	h->SetStats(0);
	h->Draw("e");
	background->Draw("same");   
	mass->SetRange(minhisto,maxhisto);
	mass->Draw("same");
	mass->SetLineStyle(2);
	mass->SetFillStyle(3004);
	mass->SetFillColor(2);
	f->Draw("same");
	
	npyield = mass->Integral(minhisto,maxhisto)/binwidthmass;
	//npyieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);

	TLegend* leg = new TLegend(0.65,0.58,0.82,0.88,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(h,"Data","pl");
	//leg->AddEntry(h,"Bs #rightarrow J/#psi + X","pl");
	leg->AddEntry(f,"Fit","l");
        leg->AddEntry(mass,"Peaking BG","f");
	leg->AddEntry(background,"Combinatorial","l");
	leg->Draw("same");

	TLatex* texYield = new TLatex(0.55,0.51,Form("Yield:%.2f", npyield));
	texYield->SetNDC();
	texYield->SetTextAlign(12);
	texYield->SetTextSize(0.035);
	texYield->SetTextFont(42);
	//texYield->Draw();

	TLatex* texChi = new TLatex(0.55,0.475,Form("#chi^{2}/nDOF:%.2f/%d=%.2f", f->GetChisquare(), f->GetNDF(), f->GetChisquare()/f->GetNDF()));
	texChi->SetNDC();
	texChi->SetTextAlign(12);
	texChi->SetTextSize(0.035);
	texChi->SetTextFont(42);
        texChi->Draw();

	TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
	texCms->SetNDC();
	texCms->SetTextAlign(12);
	texCms->SetTextSize(0.04);
	texCms->SetTextFont(42);
	texCms->Draw();

	TLatex* texCol;
	if(collisionsystem=="pp"||collisionsystem=="PP") texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","pp"));
	else texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","PbPb"));
	texCol->SetNDC();
	texCol->SetTextAlign(32);
	texCol->SetTextSize(0.04);
	texCol->SetTextFont(42);
	texCol->Draw();

	TLatex* tex;

	tex = new TLatex(0.22,0.78,Form("%.0f < p_{T} < %.0f GeV/c",ptmin,ptmax));
	tex->SetNDC();
	tex->SetTextFont(42);
	tex->SetTextSize(0.04);
	tex->SetLineWidth(2);
	tex->Draw();

	if(centMax>0){
		TString texper="%";
		tex = new TLatex(0.22,0.71,Form("Cent. %.0f-%.0f%s",centMin,centMax,texper.Data()));
		tex->SetNDC();
		tex->SetTextColor(1);
		tex->SetTextFont(42);
		tex->SetTextSize(0.045);
		tex->SetLineWidth(2);
		tex->Draw();
	}

	tex = new TLatex(0.22,0.83,"|y| < 2.4");
	tex->SetNDC();
	tex->SetTextFont(42);
	tex->SetTextSize(0.04);
	tex->SetLineWidth(2);
	tex->Draw();

	total=f;

	h->Write();
	f->Write();
	if(!isPbPb) c->SaveAs(Form("plotNP/BMass%s.pdf",collisionsystem.Data()));
	else c->SaveAs(Form("plotNP/BMass%s_%.0f_%.0f.pdf",collisionsystem.Data(),centMin,centMax));
	return mass;
}


int main(int argc, char *argv[])
{
	if(argc==14)
	{
		fitNP(atoi(argv[1]),argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11],atof(argv[12]),atof(argv[13]));
		return 0;
	}
	else if(argc==12)
	{
		fitNP(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11]);
		return 0;
	}
	else
	{
		std::cout << "Wrong number of inputs" << std::endl;
		return 1;
	}
}



double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing){
	double squarederroronsigma=(1+smear)*(1+smear)*errwidth*errwidth+width*width*errsmearing*errsmearing;
	double erroronsigma=TMath::Sqrt(squarederroronsigma);
	return erroronsigma;

}

