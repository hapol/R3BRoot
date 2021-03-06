void eventDisplay()
{
  FairRunAna *fRun= new FairRunAna();
  
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo*  parIo1 = new FairParRootFileIo();
  parIo1->open("r3bpar.root");
  rtdb->setFirstInput(parIo1);
  rtdb->print();
  
  fRun->SetInputFile("r3bsim.root");
  fRun->SetOutputFile("test.root");
  
  FairEventManager *fMan= new FairEventManager();
  FairMCTracks *Track =  new FairMCTracks ("Monte-Carlo Tracks");
  FairMCPointDraw *LandPoints =   new FairMCPointDraw ("LandPoint",kOrange,  kFullSquare);
  
  fMan->AddTask(Track);
  
  fMan->AddTask(LandPoints);
  
  fMan->Init();
}
