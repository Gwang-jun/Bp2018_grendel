//HltInfo
int           Bf_HLT_Run;
ULong64_t     Bf_HLT_Event;
int           Bf_HLT_LumiBlock;
void setHltBranch(TTree* hltroot)
{
  hltroot->SetBranchAddress("Run",&Bf_HLT_Run);
  hltroot->SetBranchAddress("Event",&Bf_HLT_Event);
  hltroot->SetBranchAddress("LumiBlock",&Bf_HLT_LumiBlock);
}

//hiEvtInfo
unsigned int       Bf_HiTree_Run;
unsigned long long Bf_HiTree_Evt;
unsigned int       Bf_HiTree_Lumi;
void setHiTreeBranch(TTree* hiroot)
{
  hiroot->SetBranchAddress("run",&Bf_HiTree_Run);
  hiroot->SetBranchAddress("evt",&Bf_HiTree_Evt);
  hiroot->SetBranchAddress("lumi",&Bf_HiTree_Lumi);
}
