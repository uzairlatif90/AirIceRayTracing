#include "MultiRayAirIceRefraction.cc"

void RunMultiRayCode_loop(){
  
  ////All variables are in m here
  double AntennaDepth=-200;////Depth of antenna in the ice
  double IceLayerHeight=3000;////Height where the ice layer starts off
  double AntennaNumber=0;
  
  MultiRayAirIceRefraction::MakeRayTracingTable(AntennaDepth*100, IceLayerHeight*100, 0);

  double GridStepSizeHb=123;
  double GridStepSizeThb=0.23;

  double GridStartThb=90.2;
  double GridStopThb=179.8;
  double GridWidthThb=GridStopThb-GridStartThb;

  double GridStartHb=IceLayerHeight+1;////just set a non-zeronumber for now
  double GridStopHb=100000;////just set a non-zeronumber for now
  if(AntennaDepth>=0){
    GridStartHb=AntennaDepth+IceLayerHeight+1;
  }
  double GridWidthHb=GridStopHb-GridStartHb;

  int GridPointsb=100;////just set a non-zeronumber for now
  int TotalStepsHb=100;////just set a non-zeronumber for now
  int TotalStepsThb=100;////just set a non-zeronumber for now
  
  TotalStepsHb=(GridWidthHb/GridStepSizeHb)+1;
  TotalStepsThb=(GridWidthThb/GridStepSizeThb)+1;
  GridPointsb=TotalStepsHb*TotalStepsThb;
  
  TH2D *h2b=new TH2D("","",100,(GridStartHb-1),(GridStopHb+1),100,GridStartThb-1,GridStopThb+1);
  TH2D *h2c=new TH2D("","",100,(GridStartHb-1),(GridStopHb+1),100,GridStartThb-1,GridStopThb+1);
  TH2D *h2corr=new TH2D("","",100,(GridStartHb-1),(GridStopHb+1),100,GridStartThb-1,GridStopThb+1);
  //TH2D *launchtest=new TH2D("","",100,(GridStartHb-1),(GridStopHb+1),100,GridStartThb-1,GridStopThb+1);

  TH1D * h1=new TH1D("","",200,0,20000);
  TH1D * h1error_dRR=new TH1D("","",200,0,100);
  TH1D * h1error=new TH1D("","",200,0,100);

  TH1D * h1RTduration=new TH1D("","",200,0,0.5);
  
  TGraph2D* gr2b=new TGraph2D();
  TGraph2D* gr2c=new TGraph2D();
  TGraph2D* gr2corr=new TGraph2D();
  
  int count1=0;
  int count2=0;
  int count3=0;
  int count4=0;
  
  double minb=1e100;
  double minc=1e100;
  
  for(int ih=0;ih<TotalStepsHb;ih++){
    for(int ith=0;ith<TotalStepsThb;ith++){

      double hR=GridStartHb+GridStepSizeHb*ih;
      double thR=GridStartThb+GridStepSizeThb*ith;
      double TotalHorizontalDistance=0;
      if(AntennaDepth<0){
	TotalHorizontalDistance=(hR-IceLayerHeight-AntennaDepth)*tan((180-thR)*(MultiRayAirIceRefraction::pi/180.0));
      }
      if(AntennaDepth>=0){
	TotalHorizontalDistance=(hR-(IceLayerHeight+AntennaDepth))*tan((180-thR)*(MultiRayAirIceRefraction::pi/180.0));
      }
      
      
      double dummy[20];

      double opticalPathLengthInIce;
      double opticalPathLengthInAir;
      double geometricalPathLengthInIce;
      double geometricalPathLengthInAir;
      double launchAngle;
      double horizontalDistanceToIntersectionPoint;
      double reflectionCoefficientS;
      double reflectionCoefficientP;
      bool checksol=true;

      auto t1 = std::chrono::high_resolution_clock::now();
      checksol=MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint( hR*100,  TotalHorizontalDistance*100, AntennaDepth*100,  IceLayerHeight*100, opticalPathLengthInIce, opticalPathLengthInAir, geometricalPathLengthInIce, geometricalPathLengthInAir,launchAngle,horizontalDistanceToIntersectionPoint,reflectionCoefficientS,reflectionCoefficientP);
      auto t2 = std::chrono::high_resolution_clock::now();
      double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
      
      h1RTduration->Fill(duration/1000);

      double rtresult=opticalPathLengthInIce/100;//horizontalDistanceToIntersectionPoint/100;//launchAngleB*(180.0/MultiRayAirIceRefraction::pi);
      if(checksol==false || std::isnan(rtresult)==true || rtresult>1e10 ){
	rtresult=-1000;
      }
      
      ////For recording how much time the process took
      double opticalPathLengthInIceB;
      double opticalPathLengthInAirB;
      double geometricalPathLengthInIceB;
      double geometricalPathLengthInAirB;
      double launchAngleB;
      double horizontalDistanceToIntersectionPointB;
      double reflectionCoefficientSB;
      double reflectionCoefficientPB;
      
      auto t1b = std::chrono::high_resolution_clock::now();
      checksol=MultiRayAirIceRefraction::GetHorizontalDistanceToIntersectionPoint_Table( hR*100,  TotalHorizontalDistance*100, AntennaDepth*100,  IceLayerHeight*100, 0, opticalPathLengthInIceB, opticalPathLengthInAirB, geometricalPathLengthInIceB, geometricalPathLengthInAirB, launchAngleB,horizontalDistanceToIntersectionPointB,reflectionCoefficientSB,reflectionCoefficientPB);
     
      auto t2b = std::chrono::high_resolution_clock::now();                                      
      double Duration = std::chrono::duration_cast<std::chrono::nanoseconds>( t2b - t1b ).count(); 
      double interresult=opticalPathLengthInIceB/100;//horizontalDistanceToIntersectionPointB/100;//launchAngleB*(180.0/MultiRayAirIceRefraction::pi);

      
      if(checksol==false || std::isnan(interresult)==true || interresult>1e10 ){
	interresult=-1000;
      }

      if(rtresult!=-1000){
	h2b->Fill((hR),thR,rtresult);
	gr2b->SetPoint(count1,(hR),thR,rtresult);
	count1++;
      }
      
      if(rtresult<minb && rtresult!=-1000){
        minb=rtresult;
      }
      
      if(interresult<minc && interresult!=-1000){
        minc=interresult;
      }
      
      if(interresult!=-1000){
	h2c->Fill((hR),thR,interresult);
	gr2c->SetPoint(count2,(hR),thR,interresult);
	count2++;
      }
      
      if(rtresult!=-1000 && interresult!=-1000){
	h1->Fill(Duration);
	h2corr->Fill((hR),thR,((rtresult-interresult)/rtresult)*100);
	gr2corr->SetPoint(count3,(hR),thR,((rtresult-interresult)/rtresult)*100);
	count3++;
	h1error_dRR->Fill((fabs(rtresult-interresult)/rtresult) *100);
	h1error->Fill(fabs(rtresult-interresult));
      }
	
	
      if(rtresult==-1000 && interresult!=-1000){
	//cout<<count4<<" we have a problem! "<<hR<<" "<<thR<<" "<<rtresult<<" "<<interresult<<endl;
	count4++;
      }
      
      //if(fabs(rtresult-interresult)>100){
      //cout<<ih<<" "<<ith<<" "<<hR<<" "<<thR<<" "<<rtresult<<" "<<interresult<<" THD "<<TotalHorizontalDistance<<endl;
      //}
      
    }
  }

  cout<<minb<<" "<<h2b->GetMaximum()<<endl;
  cout<<minc<<" "<<h2c->GetMaximum()<<endl;
  
  h2b->GetZaxis()->SetRangeUser(minb,h2b->GetMaximum());
  h2c->GetZaxis()->SetRangeUser(minc,h2c->GetMaximum());
  
  gr2b->GetZaxis()->SetRangeUser(minb,gr2b->GetMaximum());
  gr2c->GetZaxis()->SetRangeUser(minc,gr2c->GetMaximum());

  gr2b->SetTitle("RayTrace results for OP_Ice; Tx Height in Air (m); Straight Line Angle (deg); OP_Ice (m)");
  gr2c->SetTitle("Interpolated results for OP_Ice; Tx Height in Air (m); Straight Line Angle (deg); OP_Ice (m)");
  gr2corr->SetTitle("Percentage Error for OP_Ice; Tx Height in Air (m); Straight Line Angle (deg); |#DeltaOP_Ice|/OP_Ice x 100");
  h1error->SetTitle("Absolute Error for OP_Ice; |#DeltaOP_Ice| (m) ; Cumulative Fraction of Tx Positions;");
  h1error_dRR->SetTitle("Percentage Error for OP_Ice; |#DeltaOP_Ice|/OP_Ice x 100 ; Cumulative Fraction of Tx Positions;");
  h1->SetTitle("Time taken to do interpolation; Duration (ns) ; Interpolation calls;");
  
  TCanvas *c2a=new TCanvas("c2a","c2a");
  c2a->cd();
  //c2a->SetLogz();
  c2a->SetGridx();
  c2a->SetGridy();
  gr2b->Draw("cont4z");
  c2a->SaveAs("OP_Ice_RT.png");
  
  TCanvas *c2b=new TCanvas("c2b","c2b");
  c2b->cd();
  //c2b->SetLogz();
  c2b->SetGridx();
  c2b->SetGridy();
  gr2c->Draw("cont4z");
  c2b->SaveAs("OP_Ice_Inter.png");

  TCanvas *c2c=new TCanvas("c2c","c2c");
  c2c->cd();
  ////c2c->SetLogz();
  c2c->SetGridx();
  c2c->SetGridy();
  gr2corr->Draw("cont4z");
  c2c->SaveAs("OP_Ice_RT_IT.png");

  TH1* h1error_cum = h1error->GetCumulative();
  TH1* h1error_dRR_cum = h1error_dRR->GetCumulative();

  h1error_cum->Scale(1./count3);
  h1error_dRR_cum->Scale(1./count3);

  TCanvas *c2d=new TCanvas("c2d","c2d");
  c2d->Divide(2,1);
  c2d->cd(1);
  ////c2d->cd(1)->SetLogy();
  c2d->cd(1)->SetGridx();
  c2d->cd(1)->SetGridy();
  h1error_cum->SetStats(kFALSE);
  h1error_cum->Draw("");

  c2d->cd(2);
  ////c2d->cd(2)->SetLogy();
  c2d->cd(2)->SetGridx();
  c2d->cd(2)->SetGridy();
  h1error_dRR_cum->SetStats(kFALSE);
  h1error_dRR_cum->Draw("");
  c2d->SaveAs("OP_Ice_RT_IT_hist.png");

  TCanvas *c2e=new TCanvas("c2e","c2e");
  c2e->cd();
  c2e->SetLogy();
  c2e->SetGridx();
  c2e->SetGridy();
  h1->Draw("");
  c2e->SaveAs("OP_Ice_Duration_hist.png");

  // TCanvas *c2f=new TCanvas("c2f","c2f");
  // c2f->cd();
  // //c2e->SetLogy();
  // c2f->SetGridx();
  // c2f->SetGridy();
  // //c2e->SaveAs("OP_Ice_Duration_hist.png");
  
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(2,2);
  c1->cd(1);
  //c1->cd(1)->SetLogz();
  c1->cd(1)->SetGridx();
  c1->cd(1)->SetGridy();
  gr2b->Draw("cont4z");
  c1->cd(2);
  //c1->cd(2)->SetLogz();
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  gr2c->Draw("cont4z");
  c1->cd(3);
  ////c1->cd(3)->SetLogz();
  c1->cd(3)->SetGridx();
  c1->cd(3)->SetGridy();
  gr2corr->Draw("cont4z");
  c1->cd(4);
  ////c1->cd(4)->SetLogy();
  c1->cd(4)->SetGridx();
  c1->cd(4)->SetGridy();
  h1error_dRR_cum->Draw();
  c1->SaveAs("OP_Ice_All.png");

  // TCanvas *c3=new TCanvas("c3","c3");
  // c3->cd();
  // c3->SetLogy();
  // c3->SetGridx();
  // c3->SetGridy();
  // h1RTduration->Draw("");
  
}
