const int nBins=1;
double ptBins[nBins+1] = {5., 50.};
//double ptBins[nBins+1] = {5., 10., 15., 20., 50.};
//double ptBins[nBins+1] = {5., 10., 15., 20., 30., 40., 50.};
//double ptBins[nBins+1] = {5., 7., 10., 15., 20., 30., 40., 50., 60.};

//For Centrality RAA Analysis, designate p_t bins.
const int nBinsInc=1;
double ptBinsInc[nBinsInc+1] = {5., 50.};
//double ptBinsInc[nBinsInc+1] = {5., 60.};

const int nBinsFine=45;
double ptBinsFine[nBinsFine+1] = {5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50.};

const int nBinsY=4;
double ptBinsY[nBinsY+1] = {0.0, 0.5, 1.0, 1.5, 2.4};

//Use appropriate values for the designated centrality bins!!
/*
const int nBinsCent=2;
double ptBinsCent[nBinsCent+1] = {0.*2, 30.*2, 90.*2};
double TAA[nBinsCent] = {15.41, 1.705};
double npart[nBinsCent] = {269.1, 54.44};
*/


const int nBinsCent=1;
double ptBinsCent[nBinsCent+1] = {0.*2, 90.*2};
double TAA[nBinsCent] = {6.274};
double npart[nBinsCent] = {126};


/*
const int nBinsCent=1;
double ptBinsCent[nBinsCent+1] = {0.*2, 30.*2};
double TAA[nBinsCent] = {15.41};
double npart[nBinsCent] = {269.1};
*/

/*
const int nBinsCent=1;
double ptBinsCent[nBinsCent+1] = {30.*2, 90.*2};
double TAA[nBinsCent] = {1.705};
double npart[nBinsCent] = {54.44};
*/

/*
const int nBinsCent=4;
double ptBinsCent[nBinsCent+1] = {0.*2, 10*2., 30.*2, 50.*2, 90.*2};
double TAA[nBinsCent] = {23.09, 11.57, 3.939, 0.5888};
double npart[nBinsCent] = {357.4, 224.9, 108.8, 27.28};
*/

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/Glauber5TeVPbPbNewParameters                                                                      


const int nBinsReweight=7;
double ptBinsReweight[nBinsReweight+1]={5., 7., 10., 15., 20., 30., 50., 100.};

const int BIN_NUM=241;
const int HMIN=0;
const int HMAX=120;

const double binsize=((double)(HMAX)-(double)(HMIN))/(double)(BIN_NUM);
Double_t BRchain=6.02061e-5; //previous: 6.09604e-5


