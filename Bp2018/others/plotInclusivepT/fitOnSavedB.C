#include "uti.h"
#include "parameters.h"
#include "TF1.h"

Double_t setparam0=100.;
Double_t setparam1=5.28;
Double_t setparam2=0.05;
Double_t setparam3=0.03;
Double_t fixparam1=5.279;

Double_t minhisto=5.0;
Double_t maxhisto=6.0;
Double_t nbinsmasshisto=50;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

TString weight;
TString weightgen;
TString weightdata;
TString seldata;
TString selmc;
TString selmceff;
TString selmcgen;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;
double _ErrCor=1;

int _nBins = nBins;
double *_ptBins = ptBins;

void fitOnSavedB(int usePbPb=0, TString inputhist="ROOTfiles/hPtSpectrumSaveHistBplusPP.root", TString inputmc="/data/HeavyFlavourRun2/MC2015/Bntuple/pp/Bntuple20160606_pp_Pythia8_BuToJpsiK_Bpt5p0_Pthat5.root", TString trgselection="1",  TString cut="TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&Bmass>5&&Bmass<6&&Btrk1Pt>0.9&&Bchi2cl>1.32e-02&&(Bd0/Bd0Err)>3.41&&cos(Bdtheta)>-3.46e-01&&Bmu1pt>1.5&&Bmu2pt>1.5&&Blxy>0.025", TString cutmcgen="TMath::Abs(Gy)<2.4&&abs(GpdgId)==521&&GisSignal==1", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=100.)
{
	collisionsystem=collsyst;
	if(collisionsystem=="ppInc"||collisionsystem=="PbPbInc"){
		_nBins = nBinsInc;
		_ptBins = ptBinsInc;
	}

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
		selmceff = Form("%s",cut.Data());
		selmcgen = Form("%s",cutmcgen.Data());
	}
	else
	{
		seldata = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);
		selmceff = Form("%s&&hiBin>=%f&&hiBin<=%f",cut.Data(),hiBinMin,hiBinMax);
		selmcgen = Form("%s&&hiBin>=%f&&hiBin<=%f",cutmcgen.Data(),hiBinMin,hiBinMax);
	}

	selmc = Form("%s",cut.Data());

	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadRightMargin(0.043);
	gStyle->SetPadLeftMargin(0.18);
	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadBottomMargin(0.145);
	gStyle->SetTitleX(.0f);

	void clean0 (TH1D* h);
	TF1* fit (TFile* inf, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centmin, Float_t centmax, TString npfit);

    weightgen="1";
    weight="1";
    weightdata="1";
    if(doDataCor == 1){
    	weightdata="0";
    }

	TFile* inf = new TFile(inputhist.Data());

	TH1D* hPt = new TH1D("hPt","",_nBins,_ptBins);

	TH1D* hMean = new TH1D("hMean","",_nBins,_ptBins);                       
	TF1 *totalmass;

	TString outputf;
	outputf = Form("%s",outputfile.Data());
	TFile* outf = new TFile(outputf.Data(),"recreate");
	outf->cd();

	for(int i=0;i<_nBins;i++)
	{
		TF1* f = fit(inf,_ptBins[i],_ptBins[i+1],isMC,isPbPb, totalmass,centmin, centmax, npfit);
		hMean->SetBinContent(i+1,f->GetParameter(1));
		hMean->SetBinError(i+1,f->GetParError(1));  
		double yield = f->Integral(minhisto,maxhisto)/binwidthmass;
		double yieldErr = f->Integral(minhisto,maxhisto)/binwidthmass*f->GetParError(0)/f->GetParameter(0);
		printf("yield: %f, yieldErr: %f\n", yield, yieldErr);
		yieldErr = yieldErr*_ErrCor;
		hPt->SetBinContent(i+1,yield/(_ptBins[i+1]-_ptBins[i]));
		hPt->SetBinError(i+1,yieldErr/(_ptBins[i+1]-_ptBins[i]));
	}  


	TCanvas* cPt =  new TCanvas("cPt","",600,600);
	cPt->SetLogy();
	hPt->SetXTitle("D^{0} p_{T} (GeV/c)");
	hPt->SetYTitle("Uncorrected dN(D^{0})/dp_{T}");
	hPt->Sumw2();
	hPt->Draw();
	hPt->Write();
	hMean->Write();
	outf->Close();

}

void clean0(TH1D* h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		if(h->GetBinContent(i)==0) h->SetBinError(i,1);
	}
}

TF1 *fit(TFile *inf, Double_t ptmin, Double_t ptmax, int isMC,bool isPbPb,TF1* &total,Float_t centmin, Float_t centmax, TString npfit)
{
	static Int_t count=0;
	count++;
	TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
	TH1D* h = (TH1D*)inf->Get(Form("h-%d",count));
	TH1D* hMCSignal = (TH1D*)inf->Get(Form("hMCSignal-%d",count));

	TString iNP = npfit;
	TF1* f = new TF1(Form("f%d",count),"[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*("+iNP+")");
	f->SetNpx(5000);
	f->SetLineWidth(5);

	clean0(h);

	//h->Draw();
	f->SetParLimits(4,-1000,1000);
	f->SetParLimits(2,0.01,0.05);
	f->SetParLimits(8,0.01,0.05);
	f->SetParLimits(7,0,1);
	f->SetParLimits(5,0,1000);

	f->SetParameter(0,setparam0);
	f->SetParameter(1,setparam1);
	f->SetParameter(2,setparam2);
	f->SetParameter(8,setparam3);
	f->FixParameter(1,fixparam1);
	f->FixParameter(5,0);
	if(weightdata != "1"){
		int maxb = h->GetMaximumBin();
		double _max = h->GetBinContent(maxb);
		double _maxE = h->GetBinError(maxb);
		_ErrCor = (_maxE/_max)/(1/sqrt(_max));
		f->SetParLimits(0,0,1e5);
		f->SetParLimits(4,-1e5,1e5);
		f->SetParLimits(5,0,1e4);
	}
	h->GetEntries();

	hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
	hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
	f->ReleaseParameter(1);
	hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
	hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
	hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
	hMCSignal->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);

	f->FixParameter(1,f->GetParameter(1));
	f->FixParameter(2,f->GetParameter(2));
	f->FixParameter(7,f->GetParameter(7));
	f->FixParameter(8,f->GetParameter(8));
	f->ReleaseParameter(5);
	printf("Fixed para.:\n");
	printf("%f, %f, %f\n", f->GetParameter(2), f->GetParameter(7), f->GetParameter(8));
	h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
	h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
	f->ReleaseParameter(1);
	h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
	h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
	h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
	h->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);
	//accessing fir results
	//https://root.cern.ch/root/html/ROOT__Fit__FitResult.html
	//TFitResultPtr fr = h->Fit(Form("f%d",count),"L m s e","",minhisto,maxhisto);
	//Int_t fitStatus = fr;
	//printf("fit result status: %d\n", fitStatus);
	//printf("Central val: %f\n", fr->GetParams()[0]);
	//printf("HESSE err: %f (%.2f%)\n" , fr->GetErrors()[0], fr->GetErrors()[0]/fr->GetParams()[0]*100);
	//printf("Minos err: %f (%.2f%), %f (%.2f%)\n", fr->LowerError(0), -fr->LowerError(0)/fr->GetParams()[0]*100, fr->UpperError(0), fr->UpperError(0)/fr->GetParams()[0]*100);
	//printf("diff in %: (%.2f%), (%.2f%)\n", -fr->LowerError(0)/fr->GetParams()[0]*100-fr->GetErrors()[0]/fr->GetParams()[0]*100, fr->UpperError(0)/fr->GetParams()[0]*100-fr->GetErrors()[0]/fr->GetParams()[0]*100);
	
	if(weightdata != "1"){
		//      h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
		//      h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
		//      h->Fit(Form("f%d",count),"m","",minhisto,maxhisto);
	}

	TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x");
	background->SetParameter(0,f->GetParameter(3));
	background->SetParameter(1,f->GetParameter(4));
	background->SetLineColor(4);
	background->SetRange(minhisto,maxhisto);
	//background->SetLineStyle(2);//PAS
	background->SetLineStyle(7);//paper
	//background->SetLineWidth(5);
	background->SetLineWidth(9);

	TF1 *Bkpi = new TF1(Form("fBkpi%d",count),"[0]*("+iNP+")");
	Bkpi->SetParameter(0,f->GetParameter(5));
	Bkpi->SetRange(minhisto,maxhisto);
	Bkpi->SetLineStyle(1);
	//Bkpi->SetFillStyle(3004);//PAS
	Bkpi->SetFillStyle(3005);//paper
	//Bkpi->SetLineColor(kGreen+1);//PAS
	//Bkpi->SetFillColor(kGreen+1);//PAS
	Bkpi->SetLineColor(kGreen+4);//paper
	Bkpi->SetFillColor(kGreen+4);//paper
	//Bkpi->SetLineWidth(5);
	Bkpi->SetLineWidth(9);

	TF1 *mass = new TF1(Form("fmass%d",count),"[0]*([3]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[3])*Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
	mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(7),f->GetParameter(8));
	mass->SetParError(0,f->GetParError(0));
	mass->SetParError(1,f->GetParError(1));
	mass->SetParError(2,f->GetParError(2));
	mass->SetParError(7,f->GetParError(7));
	mass->SetParError(8,f->GetParError(8));
	//mass->SetLineColor(2);//PAS
	//mass->SetFillColor(2);//PAS
	mass->SetFillColor(kOrange-3);//paper
	mass->SetLineColor(kOrange-3);//paper
	//mass->SetFillStyle(3004);//PAS
	mass->SetFillStyle(3002);//paper
	//mass->SetLineWidth(5);
	mass->SetLineWidth(9);
	//mass->SetLineStyle(2);//PAS
	mass->SetLineStyle(7);//paper

	//h->SetXTitle("m_{#mu#muK} (GeV/c^{2})");
	h->SetXTitle("m_{B} (GeV/c^{2})");
	h->SetYTitle("Events / (20 MeV/c^{2})");
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	//h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
	//h->GetXaxis()->SetTitleOffset(1.3);
	//h->GetYaxis()->SetTitleOffset(1.8);
	h->GetXaxis()->SetTitleOffset(0.9);
	h->GetYaxis()->SetTitleOffset(1.3);
	//h->GetXaxis()->SetLabelOffset(0.007);
	//h->GetYaxis()->SetLabelOffset(0.007);
	//h->GetXaxis()->SetTitleSize(0.045);
	//h->GetYaxis()->SetTitleSize(0.045);
	h->GetXaxis()->SetTitleSize(0.07);
	h->GetYaxis()->SetTitleSize(0.07);
	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleFont(42);
	h->GetXaxis()->SetLabelFont(42);
	h->GetYaxis()->SetLabelFont(42);
	//h->GetXaxis()->SetLabelSize(0.04);
	//h->GetYaxis()->SetLabelSize(0.04);
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetYaxis()->SetLabelSize(0.06);
	//h->SetMarkerSize(1.0);
	h->SetMarkerSize(1.55);
	h->SetMarkerStyle(20);
	h->SetLineColor(1);
	h->SetLineWidth(5);
	h->SetStats(0);
	h->GetXaxis()->SetNdivisions(-50205);
	h->Draw("e");
	Bkpi->SetRange(5.00,5.60);
	Bkpi->Draw("same");
	background->Draw("same");   
	mass->SetRange(5.16,5.40);
	mass->Draw("same");
	f->Draw("same");
	c->RedrawAxis();

	Double_t yield = mass->Integral(minhisto,maxhisto)/binwidthmass;
	Double_t yieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);

	//TLegend* leg = new TLegend(0.62,0.58,0.82,0.88,NULL,"brNDC");
    TLegend *leg = new TLegend(0.525,0.38,0.85,0.70,NULL,"brNDC");//paper
	leg->SetBorderSize(0);
	//leg->SetTextSize(0.04);
	leg->SetTextSize(0.06);
	leg->SetTextFont(42);
	//leg->SetTextFont(62);
	leg->SetFillStyle(0);
	leg->AddEntry(h,"Data","pl");
	leg->AddEntry(f,"Fit","l");
	//leg->AddEntry(mass,"B^{+} signal","f");
	leg->AddEntry(mass,"Signal","f");
	leg->AddEntry(background,"Combinatorial","l");
	leg->AddEntry(Bkpi,"B #rightarrow J/#psi X","f");
	leg->Draw("same");

	TLatex* texChi = new TLatex(0.58,0.55, Form("#chi^{2}/nDOF: %.2f/%d = %.2f", f->GetChisquare(), f->GetNDF(), f->GetChisquare()/f->GetNDF()));
	texChi->SetNDC();
	texChi->SetTextAlign(12);
	texChi->SetTextSize(0.03);
	texChi->SetTextFont(42);
	//texChi->Draw();
	printf("NDF: %d, chi: %f, prob: %f, recal: %f\n", f->GetNDF(), f->GetChisquare(), f->GetProb(), TMath::Prob(f->GetChisquare(), f->GetNDF()));
	float _chi2 = 0;
	for(int i = 0; i < 50; i++){
		float _data = h->GetBinContent(i+1);
		float _fit = f->Eval(h->GetBinCenter(i+1));
		float _err = h->GetBinError(i+1);
		float _chi = (_data-_fit)*(_data-_fit)/_err/_err;
		//printf("%f, %f, %f, %f\n", _data, _fit, _err, _chi);
		_chi2 += _chi;
	}
	printf("self cal chi2: %f\n", _chi2);

	//TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
	//TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS}");
	//texCms->SetNDC();
	//texCms->SetTextAlign(12);
	//texCms->SetTextSize(0.04);
	//texCms->SetTextFont(42);
	//texCms->Draw();
	TLatex* texcms = new TLatex(0.22,0.87,"CMS");
	texcms->SetNDC();
	texcms->SetTextAlign(13);
	texcms->SetTextFont(62);//61
	texcms->SetTextSize(0.08);
	texcms->SetLineWidth(2);
	texcms->Draw();

    //TLatex* texB = new TLatex(0.81,0.30,"B^{+}");
    TLatex* texB = new TLatex(0.22,0.73,"B^{#plus}+B^{#minus}");
    texB->SetNDC();
    texB->SetTextFont(42);
    texB->SetTextSize(0.07);
    texB->SetLineWidth(2);
    texB->Draw();

	TLatex* texCol;
	//if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc"||collisionsystem=="PbPbInc") texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s} = 5.02 TeV","pp"));
	//else texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","PbPb"));
	//if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc"||collisionsystem=="PbPbInc") texCol= new TLatex(0.935,0.93, Form("28.0 pb^{-1} (#sqrt{s} = 5.02 TeV %s)","pp"));
	//else texCol= new TLatex(0.93,0.93, Form("350.1 #mub^{-1} (#sqrt{s_{NN}} = 5.02 TeV %s)","PbPb"));
	if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc"||collisionsystem=="PbPbInc") texCol= new TLatex(0.945,0.94, Form("28.0 pb^{-1} (%s 5.02 TeV)","pp"));
	else texCol= new TLatex(0.94,0.94, Form("351 #mub^{-1} (%s 5.02 TeV)","PbPb"));
	texCol->SetNDC();
	texCol->SetTextAlign(32);
	//texCol->SetTextSize(0.04);
	texCol->SetTextSize(0.06);
	texCol->SetTextFont(42);
	texCol->Draw();

	TLatex* tex;

	//tex = new TLatex(0.22,0.73,Form("%.0f < p_{T} < %.0f GeV/c",ptmin,ptmax));
	//tex = new TLatex(0.488,0.84,Form("%.0f < p_{T} < %.0f GeV/c",ptmin,ptmax));
	tex = new TLatex(0.488,0.82,Form("%.0f < p_{T} < %.0f GeV/c",ptmin,ptmax));
	tex->SetNDC();
	tex->SetTextFont(42);
	//tex->SetTextSize(0.04);
	tex->SetTextSize(0.06);
	tex->SetLineWidth(2);
	tex->Draw();

	if(centMax>0){
		TString texper="%";
		tex = new TLatex(0.22,0.71,Form("Centrality %.0f-%.0f%s",centMin,centMax,texper.Data()));//0.2612903,0.8425793
		tex->SetNDC();
		tex->SetTextColor(1);
		tex->SetTextFont(42);
		tex->SetTextSize(0.045);
		tex->SetLineWidth(2);
		//if(collisionsystem!="pp"&&collisionsystem!="PP"&&collisionsystem!="ppInc"&&collisionsystem!="PbPbInc") tex->Draw();
	}

	//tex = new TLatex(0.22,0.78,"|y| < 2.4");
	//tex = new TLatex(0.700,0.77,"|y| < 2.4");
	tex = new TLatex(0.725,0.75,"|y| < 2.4");
	tex->SetNDC();
	tex->SetTextFont(42);
	//tex->SetTextSize(0.04);
	tex->SetTextSize(0.06);
	tex->SetLineWidth(2);
	tex->Draw();

	total=f;

    TF1* t = (TF1*)h->GetFunction(Form("f%d",count))->Clone();
	h->GetFunction(Form("f%d",count))->Delete();
    t->Draw("same");
	h->Draw("e same");
	h->Write();

	//if(!isPbPb) c->SaveAs(Form("plotFitsOnSaved/BMass%s_%d.pdf",collisionsystem.Data(),count));
	//else c->SaveAs(Form("plotFitsOnSaved/BMass%s_%.0f_%.0f_%d.pdf",collisionsystem.Data(),centMin,centMax,count));
	TString _postfix = "";
	if(weightdata!="1") _postfix = "_EFFCOR";
	if(isPbPb && isMC==0) 
		c->SaveAs(Form("plotFitsOnSaved/data_PbPb_%.0f_%.0f%s.pdf",ptmin,ptmax,_postfix.Data()));
	else if(isPbPb && isMC==1) 
		c->SaveAs(Form("plotFitsOnSaved/mc_PbPb_%.0f_%.0f%s.pdf",ptmin,ptmax,_postfix.Data()));
	else if(!isPbPb && isMC==0) 
		c->SaveAs(Form("plotFitsOnSaved/data_pp_%.0f_%.0f%s.pdf",ptmin,ptmax,_postfix.Data()));
	else 
		c->SaveAs(Form("plotFitsOnSaved/mc_pp_%.0f_%.0f%s.pdf",ptmin,ptmax,_postfix.Data()));

	return mass;
}

int main(int argc, char *argv[])
{
	if(argc==16)
	{
		fitOnSavedB(atoi(argv[1]),argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11],argv[12],atoi(argv[13]),atof(argv[14]),atof(argv[15]));
		return 0;
	}
	else if(argc==14)
	{
		fitOnSavedB(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11], argv[12], atoi(argv[13]));
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
