void BuildTrees(const char* path){

  std::vector<std::string> energies;
  energies.push_back("0p5GeV");
  energies.push_back("1GeV");
  energies.push_back("2GeV");
  energies.push_back("3GeV");

  std::vector<std::string> particles;
  particles.push_back("muon");
  particles.push_back("electron");
  particles.push_back("piplus");
  particles.push_back("proton");

  for(std::vector<std::string>::iterator eIter=energies.begin();eIter!=energies.end();++eIter){
    for(std::vector<std::string>::iterator pIter=particles.begin();pIter!=particles.end();++pIter){
      std::cout<<"Building tree for "<<*pIter<<"@"<<*eIter<<"..."<<std::endl;
      BuildTree(std::string(path)+"ntuple_"+*pIter+"_"+*eIter+".root","var_"+*pIter+"_"+*eIter+".root");
    }
  }
}

void BuildTree(std::string inFile,std::string outFile){
  TFile* fIn=new TFile(inFile.c_str());

  TTree* tr=(TTree*)fIn->Get("valrec/valrec");
  LArPID* pidObj=0;
  tr->SetBranchAddress("PID",&pidObj);
 
 TFile* fOut=new TFile(outFile.c_str(),"RECREATE");
 TTree* mvaTree = new TTree("mvaTree","mvaTree");
 

  double   AvgedEdxAmpStart;   
  double   AvgedEdxAmpEnd;     
  double   AvgedEdxAreaStart;  
  double   AvgedEdxAreaEnd;    
  double   EvalRatio;	       
  double   ChargeRatioCoreHalo;
  double   Conicalness;	       
  double   Concentration;      
  double   ChargeLongRatio;
  double   ChargeLongRatioHalf;
  double   ChargeEndRatio;
  double   ChargeEndRatio10;


  mvaTree->Branch("AvgedEdxAmpStart",&AvgedEdxAmpStart,"AvgedEdxAmpStart/D");   
  mvaTree->Branch("AvgedEdxAmpEnd",&AvgedEdxAmpEnd,"AvgedEdxAmpEnd/D");     
  mvaTree->Branch("AvgedEdxAreaStart",&AvgedEdxAreaStart,"AvgedEdxAreaStart/D");  
  mvaTree->Branch("AvgedEdxAreaEnd",&AvgedEdxAreaEnd,"AvgedEdxAreaEnd/D");    
  mvaTree->Branch("EvalRatio",&EvalRatio,"EvalRatio/D");	       
  mvaTree->Branch("ChargeRatioCoreHalo",&ChargeRatioCoreHalo,"ChargeRatioCoreHalo/D");
  mvaTree->Branch("Conicalness",&Conicalness,"Conicalness/D");	       
  mvaTree->Branch("Concentration",&Concentration,"Concentration/D");      
  mvaTree->Branch("ChargeLongRatio",&ChargeLongRatio,"ChargeLongRatio/D");      
  mvaTree->Branch("ChargeLongRatioHalf",&ChargeLongRatioHalf,"ChargeLongRatioHalf/D");      
  mvaTree->Branch("ChageEndRatio",&ChargeEndRatio,"ChargeEndRatio/D");      
  mvaTree->Branch("ChageEndRatio10",&ChargeEndRatio10,"ChargeEndRatio10/D");      


 for(int iEntry=0;iEntry<tr->GetEntries();++iEntry){
   tr->GetEntry(iEntry);
   if(pidObj->NTracks&&pidObj->ChargeRatioCoreHalo[0]<100){

     AvgedEdxAmpStart=max(pidObj->AvgedEdxAmpStart[0],0.);   
     AvgedEdxAmpEnd=max(pidObj->AvgedEdxAmpEnd[0],0.);     
     AvgedEdxAreaStart=max(pidObj->AvgedEdxAreaStart[0],0.);  
     AvgedEdxAreaEnd=max(pidObj->AvgedEdxAreaEnd[0],0.);    
     EvalRatio=max(pidObj->EvalRatio[0],0.);	           
     ChargeRatioCoreHalo=max(pidObj->ChargeRatioCoreHalo[0],0.);
     Conicalness=max(pidObj->Conicalness[0],0.);	           
     Concentration=max(pidObj->Concentration[0],0.);      
     ChargeLongRatio=max(pidObj->ChargeLongRatio[0],0.);      
     ChargeLongRatioHalf=max(pidObj->ChargeLongRatioHalf[0],0.);      
     ChargeEndRatio=max(pidObj->ChargeEndRatio[0],0.);      
     ChargeEndRatio10=max(pidObj->ChargeEndRatio10[0],0.);      

     mvaTree->Fill();
   }
 }

 fIn.Close();
 fOut.Write();
 fOut.Close();
}

void TrainMVA(std::vector<std::string> signalFiles,std::vector<std::string> backgroundFiles,std::string outputFile,std::string jobName){

  TFile* fOut = new TFile(outputFile.c_str(),"RECREATE");
  TMVA::Factory* factory = new TMVA::Factory( jobName.c_str(), fOut, "" );

  std::vector<TTree*> sigTrees;

  for(std::vector<std::string>::iterator fIter=signalFiles.begin();fIter!=signalFiles.end();++fIter){
    TFile* fIn=new TFile(fIter->c_str());
    factory->AddSignalTree((TTree*)fIn->Get("mvaTree"));
  }
  for(std::vector<std::string>::iterator fIter=backgroundFiles.begin();fIter!=backgroundFiles.end();++fIter){
    TFile* fIn=new TFile(fIter->c_str());
    factory->AddBackgroundTree((TTree*)fIn->Get("mvaTree"));
  }

  factory->AddVariable("AvgedEdxAmpStart",'F');   
  factory->AddVariable("AvgedEdxAmpEnd",'F');     
  //factory->AddVariable("AvgedEdxAreaStart",'F');  
  //factory->AddVariable("AvgedEdxAreaEnd",'F');    
  factory->AddVariable("EvalRatio",'F');	       
  factory->AddVariable("ChargeRatioCoreHalo",'F');
  factory->AddVariable("ChargeLongRatio",'F');
  factory->AddVariable("ChargeEndRatio",'F');
  //factory->AddVariable("Conicalness",'F');	       
  //factory->AddVariable("Concentration",'F');      

  factory->BookMethod( TMVA::Types::kTMlpANN, "TMlp_ANN", "" );
  factory->BookMethod( TMVA::Types::kBDT, "BDT", "" );
  //factory->BookMethod( TMVA::Types::kKNN, "kNN", "" );
  //factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", "" );
  //factory->BookMethod( TMVA::Types::kCuts, "Cuts", "FitMethod=SA" );
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  fOut.Write();
  fOut.Close();
}

void TrainAll(){

  std::vector<std::string> energies;
  energies.push_back("0p5GeV");
  energies.push_back("1GeV");
  energies.push_back("2GeV");
  energies.push_back("3GeV");

  std::vector<std::string> particles;
  particles.push_back("muon");
  particles.push_back("electron");
  particles.push_back("piplus");
  particles.push_back("proton");

  for(std::vector<std::string>::iterator energy=energies.begin();energy!=energies.end();++energy){
    std::vector<std::string> particleFilenames[4];
    
    int iFilename=0;
    for(std::vector<std::string>::iterator pIter=particles.begin();pIter!=particles.end();++pIter){
      
      particleFilenames[iFilename++].push_back(std::string("/data/t2k/phsmaj/lbne/data/scratch/haigh/mva/var_")+*pIter+"_"+*energy+".root");
    }
    
    for(int iPart=0;iPart<particles.size();++iPart){
      for(int jPart=iPart+1;jPart<particles.size();++jPart){
	std::string mvaName("mva_");
	mvaName=mvaName+particles[iPart]+"_"+particles[jPart]+"_"+*energy;
	TrainMVA(particleFilenames[iPart],particleFilenames[jPart],mvaName+".root",mvaName);
      }
      
      std::vector<std::string> bkgFiles;
      for(int jPart=0;jPart<particles.size();++jPart){
	if(jPart==iPart) continue;
	for(std::vector<std::string>::iterator fIter=particleFilenames[jPart].begin();
	    fIter!=particleFilenames[jPart].end();++fIter){
	  bkgFiles.push_back(*fIter);
	}
      }
      std::string mvaName("mva_");
      mvaName=mvaName+particles[iPart]+"_All_"+*energy;
      TrainMVA(particleFilenames[iPart],bkgFiles,mvaName+".root",mvaName);
    }
  }
}

void MakePlots(std::string outFileName){

  std::vector<std::string> methods;
  methods.push_back("TMlp_ANN");
  methods.push_back("BDT");
  //methods.push_back("Cuts");
  //methods.push_back("kNN");
  //methods.push_back("Likelihood");
  gStyle->SetOptStat("0000");

  TFile fOut(outFileName.c_str(),"RECREATE");

  TSystemDirectory dir("./","./");
  TIter nextFile(dir.GetListOfFiles());
  TSystemFile* file;
  while(file = (TSystemFile*)nextFile.Next()){
    std::string fName=file->GetName();
    std::cout<<fName<<std::endl;
    if(fName.find("mva_")!=std::string::npos){
      fOut.cd();
      std::string caName=fName.substr(0,fName.size()-5);
      TCanvas* ca=new TCanvas(caName.c_str(),caName.c_str());
      TLegend* le=new TLegend(0.2,0.1,0.4,0.3);
      TFile f(fName.c_str());
      int iMethod=0;
      for(std::vector<std::string>::iterator mIter=methods.begin();mIter!=methods.end();++mIter){
	std::string methodFolder;
	if(mIter->find("ANN")!=std::string::npos){
	  methodFolder="TMlpANN";
	}
	else{
	  methodFolder=*mIter;
	}
	std::cout<<(std::string("Method_")+methodFolder+"/"+*mIter+"/MVA_"+*mIter+"_rejBvsS").c_str()<<std::endl;
	TH1D* h=(TH1D*)f.Get((std::string("Method_")+methodFolder+"/"+*mIter+"/MVA_"+*mIter+"_rejBvsS").c_str());
	le->AddEntry(h);
	if(iMethod==0){
	  h->SetTitle(caName.c_str());
	  h->SetLineColor(kRed);
	  h->Draw();
	  h->GetYaxis().SetRangeUser(0,1.2);
	  firstMethod=false;
	}
	else{
	  h->SetLineColor(iMethod==1?kBlue:kGreen);
	  h->Draw("same");
	}
	++iMethod;
      }
      le->Draw();
      fOut.cd();
      std::cout<<caName.c_str()<<std::endl;
      ca->Write(caName.c_str());
      ca->SaveAs((std::string("./plots/")+caName+".png").c_str());
    }
  }
  fOut.Close();
}
