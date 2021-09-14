#include "variable_function16.h"
#include <TCanvas.h>

void analyze(const char* path, const char* output, int ch1, int ch2, int ch3, int ch4, int nBoards =4, int npoints=1024, double tscale = 1.0)
{
    TString file_name;
    TString file_ext = ".root";
    TString file;
   // std::cout<< "file name"<< std::endl;
   // std::cin >> file_name ;
    //file = ouptut;
    double a1,a2,a3,a4;
    TFile* loadFile = TFile::Open(path);
    TTree* loadTree = (TTree*) loadFile->Get("p");
    TFile* new_file = new TFile(output,"recreate","8");
    TTree* new_tree = new TTree("wfm","waveform_analysis");
    TH1F* amplitude1 = new TH1F("amp1",";Pulse Amplitude;Entries;Amplitude (mV)", 100, 0, 500);
    TH1F* amplitude2 = new TH1F("amp2",";Pulse Amplitude;Entries;Amplitude (mV)", 100, 0, 500);
    TH1F* amplitude3 = new TH1F("amp3","Pulse Amplitude;Entries;Amplitude (mV)", 100, 0, 500);
    TH1F* amplitude4 = new TH1F("amp4",";Pulse Amplitude;Entries;Amplitude (mV)", 100, 0, 500);
    int tracker = 0;
    int Entries = loadTree -> GetEntries();
    std::cout << "total events = " << Entries << endl;
 //   int npoints; //loadTree->GetLeaf(w2)->GetLen();
  //  std::cout <<"npoints = ?" << std::endl;
   // std::cin >> npoints;
   // std::cout <<"npoints = " << npoints << std::endl;
    m_inoise = ceil(npoints*0.25);
    double w1f[npoints], w2f[npoints], w3f[npoints], w4f[npoints];
    double t1f[npoints], t2f[npoints], t3f[npoints], t4f[npoints];
    int boardf{0};
    //branches to load from input file
    loadTree->SetBranchAddress("board",&boardf, &b_boardf);
    loadTree->SetBranchAddress("t1",&t1f, &b_t1f);
    loadTree->SetBranchAddress("t2",&t2f, &b_t2f);
    loadTree->SetBranchAddress("t3",&t3f, &b_t3f);
    loadTree->SetBranchAddress("t4",&t4f, &b_t4f);
    
    //channel 1 branches_____________________________________
    TBranch *bboard0 = new_tree->Branch("board0", &board0);
    if(ch1 == 1)
    {
        cout << "ch 1 is created." << endl;
  //      TBranch *bwp1 = new_tree->Branch("w1", &w1);
  //      TBranch *bwbase1 = new_tree->Branch("wbase1", &wbase1);
    //    TBranch *btime1 = new_tree->Branch("t1",&t1);
        TBranch *bcfd1 = new_tree->Branch("cfd1",&cfd1);
        TBranch *brms1 = new_tree->Branch("rms1",&rms1);
        TBranch *bpmax1 =  new_tree->Branch("pmax1",&pmax1);
        TBranch *btmax1 =  new_tree->Branch("tmax1",&tmax1);
        TBranch *bpulse_area1 = new_tree->Branch("pulse_area1",&pulse_area1);
        TBranch *brise1_1040 = new_tree->Branch("rise1_1040",&rise1_1040);
        TBranch *brise1_1090 = new_tree->Branch("rise1_1090",&rise1_1090);
        TBranch *bdvdt1 = new_tree->Branch("dvdt1",&dvdt1);
        TBranch *bt01 = new_tree->Branch("t01",&t01);
        loadTree->SetBranchAddress("c1",&w1f, &b_w1f);
    }
    //channel 2 branches_____________________________________
    if(ch2 == 2)
    {
        cout << "ch 2 is created." << endl;
 //       TBranch *bwp2 = new_tree->Branch("w2", &w2);
//        TBranch *bwbase2 = new_tree->Branch("wbase2", &wbase2);
   //     TBranch *btime2 = new_tree->Branch("t2",&t2);
        TBranch *bcfd2 = new_tree->Branch("cfd2",&cfd2);
        TBranch *brms2 = new_tree->Branch("rms2",&rms2);
        TBranch *bpmax2 =  new_tree->Branch("pmax2",&pmax2);
        TBranch *btmax2 =  new_tree->Branch("tmax2",&tmax2);
        TBranch *bpulse_area2 = new_tree->Branch("pulse_area2",&pulse_area2);
        TBranch *brise2_1040 = new_tree->Branch("rise2_1040",&rise2_1040);
        TBranch *brise2_1090 = new_tree->Branch("rise2_1090",&rise2_1090);
        TBranch *bdvdt2 = new_tree->Branch("dvdt2",&dvdt2);
        TBranch *bt02 = new_tree->Branch("t02",&t02);
        loadTree->SetBranchAddress("c2",&w2f, &b_w2f);
    }
    //channel 3 branches_____________________________________
    if(ch3 == 3)
    {
        cout << "ch 3 is created." << endl;
     //   TBranch *bwp3 = new_tree->Branch("w3", &w3);
     //   TBranch *bwbase3 = new_tree->Branch("wbase3", &wbase3);
      //  TBranch *btime3 = new_tree->Branch("t3", &t3);
        TBranch *bcfd3  = new_tree->Branch("cfd3", &cfd3);
        TBranch *brms3 = new_tree->Branch("rms3", &rms3);
        TBranch *bpmax3 =  new_tree->Branch("pmax3", &pmax3);
        TBranch *btmax3 =  new_tree->Branch("tmax3",&tmax3);
        TBranch *bpulse_area3 = new_tree->Branch("pulse_area3", &pulse_area3);
        TBranch *brise3_1040 = new_tree->Branch("rise3_1040", &rise3_1040);
        TBranch *brise3_1090 = new_tree->Branch("rise3_1090", &rise3_1090);
        TBranch *bdvdt3 = new_tree->Branch("dvdt3", &dvdt3);
        TBranch *bt03 = new_tree->Branch("t03",&t03);
        loadTree->SetBranchAddress("c3",&w3f, &b_w3f);
    }
    //channel 4 branches_____________________________________
    if(ch4 == 4)
    {
        cout << "ch 4 is created." << endl;
  //      TBranch *bwp4 = new_tree->Branch("w4", &w4);
//        TBranch *bwbase4 = new_tree->Branch("wbase4", &wbase4);
    //    TBranch *btime4 = new_tree->Branch("t4", &t4);
        TBranch *bcfd4  = new_tree->Branch("cfd4", &cfd4);
        TBranch *brms4 = new_tree->Branch("rms4", &rms4);
        TBranch *bpmax4 =  new_tree->Branch("pmax4", &pmax4);
        TBranch *btmax4 =  new_tree->Branch("tmax4",&tmax4);
        TBranch *bpulse_area4 = new_tree->Branch("pulse_area4", &pulse_area4);
        TBranch *brise4_1040 = new_tree->Branch("rise4_1040", &rise4_1040);
        TBranch *brise4_1090 = new_tree->Branch("rise4_1090", &rise4_1090);
        TBranch *bdvdt4 = new_tree->Branch("dvdt4", &dvdt4);
        TBranch *bt04 = new_tree->Branch("t04",&t04);
        loadTree->SetBranchAddress("c4",&w4f, &b_w4f);
    }
    
    std::cout << "analysis section---------------------" << endl;
    
    for(int fil = 0; fil < Entries/nBoards; fil++){
        for (int j{0}; j<nBoards; j++){
            loadTree->GetEntry(fil*nBoards +j);
            for(int i = 0; i < npoints; i++)
                   {
                  //loop for averaging windows, could do something better like sgolay?
                       //if (i%1==0) {cout << "Processing point  " << i << endl;}
                    /*   if (i<10){
                           if(ch1 == 1){a1 = 0; }
                           if(ch2 == 2){a2 = 0 ;}
                           if(ch3 == 3){a3 = 0;}
                           if(ch4 == 4){a4 = 0 ;}
                       }else{
                           if(ch1 == 1){ a1 = w1f[i]+ w1f[i-1]+ w1f[i-2]+ w1f[i-3]+w1f[i-4]+w1f[i-5]+ w1f[i-6]+ w1f[i-7]+ w1f[i-8]+w1f[i-9]; a1*=0.1;}
                           if(ch2 == 2){ a2 = w2f[i]+ w2f[i-1]+ w2f[i-2]+ w2f[i-3]+w2f[i-4]+w2f[i-5]+ w2f[i-6]+ w2f[i-7]+ w2f[i-8]+w2f[i-9]; a2*=0.1;}
                           if(ch3 == 3){ a3 = w3f[i]+ w3f[i-1]+ w3f[i-2]+ w3f[i-3]+w3f[i-4]+w3f[i-5]+ w3f[i-6]+ w3f[i-7]+ w3f[i-8]+w3f[i-9]; a3*=0.1;}
                           if(ch4 == 4){ a4 = w4f[i]+ w4f[i-1]+ w4f[i-2]+ w4f[i-3]+w4f[i-4]+w4f[i-5]+ w4f[i-6]+ w4f[i-7]+ w4f[i-8]+w4f[i-9]; a4*=0.1;}
                       }*/
                       if(ch1 == 1){ a1 = -1.0*w1f[i];}//-1 to invert pulses
                       if(ch2 == 2){ a2 = -1.0*w2f[i];}
                       if(ch3 == 3){ a3 = -1.0*w3f[i];}
                       if(ch4 == 4){ a4 = -1.0*w4f[i];}
                       
                       if(ch1 == 1){w1V.push_back( a1 ); t1V.push_back( t1f[i] ); }
                       if(ch2 == 2){w2V.push_back( a2); t2V.push_back(  t2f[i] ); }
                       if(ch3 == 3){w3V.push_back( a3 ); t3V.push_back(  t3f[i] ); }
                       if(ch4 == 4){w4V.push_back( a4); t4V.push_back(  t4f[i] ); }
            }
            if (fil%100==0){
                        if (j==0 )cout << "Processing event  " << fil << endl;
                        cout<< " The board number is "<<boardf<<endl;
                   }
            board0.push_back(boardf);
            if(ch1 == 1)
            {
               int points = w1V.size();
               base_line(points, w1V, m_inoise);
               pmax1.push_back( pulse_max(points, w1V, m_inoise) );
               amplitude1->Fill(pulse_max(points, w1V, m_inoise));
               for(int k =0; k<101; k++){ cfd1.push_back( cfd_index(points, k, t1V, w1V) ); }
             //  cfd1.push_back(cfd1V);
               tmax1.push_back( time_max(points, t1V, w1V, m_inoise) );
               rms1.push_back( noise_rms(points, w1V, m_inoise) );
               pulse_area1.push_back( pulse_area(points, t1V, w1V, m_inoise) );
               rise1_1040.push_back( rise_time(points, t1V, w1V, 0.1, 0.4) );
               rise1_1090.push_back( rise_time(points, t1V, w1V, 0.1, 0.9) );
               dvdt1.push_back( pulse_dvdt_cfd(points, 20, 1, t1V, w1V) );
               t01.push_back(th_t0(points, t1V, w1V,5) );
            }
            if(ch2 == 2)
            {
               int points = w2V.size();
               base_line(points, w2V, m_inoise);
               pmax2.push_back( pulse_max(points, w2V, m_inoise));
               amplitude2->Fill(pulse_max(points, w2V, m_inoise));
               for(int k =0; k<101; k++){ cfd2.push_back( cfd_index(points, k, t2V, w2V) ); }
            //   cfd2.push_back(cfd2V);
               tmax2.push_back( time_max(points, t2V, w2V, m_inoise) );
               rms2.push_back( noise_rms(points, w2V, m_inoise) );
               pulse_area2.push_back( pulse_area(points, t2V, w2V, m_inoise) );
               rise2_1040.push_back( rise_time(points, t2V, w2V, 0.1, 0.4) );
               rise2_1090.push_back(rise_time(points, t2V, w2V,0.10,0.9 ) );
               dvdt2.push_back( pulse_dvdt_cfd(points, 20, 1, t2V, w2V) );
               t02.push_back(th_t0(points, t2V, w2V,5) );
            }
            if(ch3 == 3)
            {
               int points = w3V.size();
               base_line(points, w3V, m_inoise);
               pmax3.push_back( pulse_max(points, w3V, m_inoise) );
               for(int k =0; k<101; k++){ cfd3.push_back( cfd_index(points, k, t3V, w3V) ); }
          //     cfd3.push_back(cfd3V);
               tmax3.push_back( time_max(points, t3V, w3V, m_inoise) );
               rms3.push_back( noise_rms(points, w3V, m_inoise) );
               pulse_area3.push_back( pulse_area(points, t3V, w3V, m_inoise) );
               rise3_1040.push_back( rise_time(points, t3V, w3V, 0.1, 0.4) );
               rise3_1090.push_back( rise_time(points, t3V, w3V, 0.1, 0.9) );
               dvdt3.push_back( pulse_dvdt_cfd(points, 20, 1, t3V, w3V) );
               t03.push_back(th_t0(points, t3V, w3V,5) );
            }
            if(ch4 == 4)
            {
               int points = w4V.size();
               base_line(points, w4V, m_inoise);
            //   for(int i = 0; i < points; i++){w4V[i] = -1.0*w4V[i]; wbase4[i]=-1.0*wbase4[i];}
               pmax4.push_back( pulse_max(points, w4V, m_inoise) );
               for(int k =0; k<101; k++){ cfd4.push_back( cfd_index(points, k, t4V, w4V) ); }
          //     cfd4.push_back(cfd4V);
               tmax4.push_back( time_max(points, t4V, w4V, m_inoise) );
               rms4.push_back( noise_rms(points, w4V, m_inoise) );
               pulse_area4.push_back( pulse_area(points, t4V, w4V, m_inoise) );
               rise4_1040.push_back( rise_time(points, t4V, w4V, 0.1, 0.4) );
               rise4_1090.push_back( rise_time(points, t4V, w4V,0.10,0.9 ) );
               dvdt4.push_back( pulse_dvdt_cfd(points, 20, 1, t4V, w4V) );
               t04.push_back(th_t0(points,t4V, w4V,5) );
           //    for(int i = 0; i < points; i++){w4[i] = -1.0*w4[i]; wbase4[i]=-1.0*wbase4[i];}
           //     if(ch1 == 1){w1.push_back( w1V ); t1.push_back( t1V ); wbase1.push_back( wbase1V );}
            //    if(ch2 == 2){w2.push_back( w2V); t2.push_back(  t2V); wbase2.push_back( wbase2V);}
             //   if(ch3 == 3){w3.push_back( w3V ); t3.push_back(  t3V ); wbase3.push_back( wbase3V);}
              //  if(ch4 == 4){w4.push_back( w4V); t4.push_back(  t4V ); wbase4.push_back( wbase4V );}
            }
                w1V.clear();
                w2V.clear();
                w3V.clear();
                w4V.clear();
                t1V.clear();
                t2V.clear();
                t3V.clear();
                t4V.clear();
         //       cfd1V.clear();
          //      cfd2V.clear();
          //      cfd3V.clear();
           //     cfd4V.clear();
        }///boards loop
        //fill new tree and clear vectors
        new_tree->Fill();
        board0.clear();
        w1.clear();
        t1.clear();
        pmax1.clear();
        tmax1.clear();
        cfd1.clear();
        rms1.clear();
        pulse_area1.clear();
        rise1_1040.clear();
        rise1_1090.clear();
        dvdt1.clear();
        t01.clear();
    //    wbase1.clear();
        
        w2.clear();
        t2.clear();
        pmax2.clear();
        tmax2.clear();
        cfd2.clear();
        rms2.clear();
        pulse_area2.clear();
        rise2_1040.clear();
        rise2_1090.clear();
        dvdt2.clear();
        t02.clear();
//        wbase2.clear();
        
        w3.clear();
        t3.clear();
        pmax3.clear();
        tmax3.clear();
        cfd3.clear();
        rms3.clear();
        pulse_area3.clear();
        rise3_1040.clear();
        rise3_1090.clear();
        dvdt3.clear();
        t03.clear();
   //     wbase3.clear();
        
        w4.clear();
        t4.clear();
        pmax4.clear();
        tmax4.clear();
        cfd4.clear();
        rms4.clear();
        pulse_area4.clear();
        rise4_1040.clear();
        rise4_1090.clear();
        dvdt4.clear();
        t04.clear();
    //    wbase4.clear();

    }//entries/nBoards loop
    new_file->Write();
    /*
    TCanvas* c1 = new TCanvas("");
    bool first = true;
    if(ch1 == 1){
        amplitude1->Draw("");
    }
    if(ch1 == 2){
        amplitude2->SetLineWidth(3);
        amplitude2->SetLineColor(2);

        if (first) {
            amplitude2->Draw("");
            first = false;
        }else amplitude2->Draw("same");
    }
     */
  //  if(ch1 == 3)  amplitude3->Draw("same");
  // if(ch1 == 4)  amplitude4->Draw("same");
    // c1->SaveAs("");
 //   new_file->Close();
    delete new_file, new_tree, amplitude1, amplitude2, amplitude3, amplitude4;
}




