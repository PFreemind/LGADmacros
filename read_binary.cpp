/*
   Name:           read_binary.cpp
   Created by:     Stefan Ritt <stefan.ritt@psi.ch>
   Date:           July 30th, 2014

   Purpose:        Example file to read binary data saved by DRSOsc.
 
   Compile and run it with:
 
      gcc -o read_binary read_binary.cpp
 
      ./read_binary <filename>

   This program assumes that a pulse from a signal generator is split
   and fed into channels #1 and #2. It then calculates the time difference
   between these two pulses to show the performance of the DRS board
   for time measurements.
p
   $Id: read_binary.cpp 22290 2016-04-27 14:51:37Z ritt $
*/

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

//for ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

typedef struct {
   char           tag[3];
   char           version;
} FHEADER;

typedef struct {
   char           time_header[4];
} THEADER;

typedef struct {
   char           bn[2];
   unsigned short board_serial_number;
} BHEADER;

typedef struct {
   char           event_header[4];
   unsigned int   event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short range;	// range center in mV
} EHEADER;

typedef struct {
   char           tc[2];
   unsigned short trigger_cell;
} TCHEADER;

typedef struct {
   char           c[1];
   char           cn[3];
} CHEADER;

/*-----------------------------------------------------------------------------*/

int main(int argc, const char * argv[])
{
   FHEADER  fh;
   THEADER  th;
   BHEADER  bh;
   EHEADER  eh;
   TCHEADER tch;
   CHEADER  ch;
   
   unsigned int scaler;
   unsigned short voltage[1024];
   double waveform[16][4][1024], time[16][4][1024];
   float bin_width[16][4][1024];
   int i, j, b, chn, n, chn_index, n_boards;
   double tt1, tt2, tt3, tt4, dt, dt34;
   char filename[256];	//for input binary file
   char outroot[256];	//for output root file

   int ndt;
   double threshold, sumdt, sumdt2;

   if (argc == 3) {
      strcpy(filename, argv[1]);
      strcpy(outroot, argv[2]);
   }
   else if (argc == 2) {
      printf("Error: both input binary file and output root file should be specified!\n");
      return -1;
   }
   else {
      printf("Error: input binary file and output root file should be specified!\n");
      return -1;
   }

// ---------------for ROOT

   TFile* rfile = new TFile(outroot, "RECREATE");
   TTree* rtree = new TTree("p", "tree for drs4 analysis");
   //rtree->Branch("t1", &t1, "t1/D"); //br for time of threshold crossing signal in 1 ch
  // rtree->Branch("t2", &t2, "t2/D"); //br for time of threshold crossing signal in 2 ch
   int ncell;
   const int ncellMax = 1030;
   double c1[ncellMax], t1[ncellMax];	//variable size array
   double c2[ncellMax], t2[ncellMax];
   double c3[ncellMax], t3[ncellMax], clock;
   int board;
//------for other channels 
   double c4[ncellMax], t4[ncellMax];

   rtree->Branch("ncell", &ncell, "ncell/I"); 
   rtree->Branch("c1", c1, "c1[ncell]/D");
   rtree->Branch("t1", t1, "t1[ncell]/D");

//------for other channels 
   rtree->Branch("c2", c2, "c2[ncell]/D");
   rtree->Branch("t2", t2, "t2[ncell]/D");
    
   rtree->Branch("c3", c3, "c3[ncell]/D");
   rtree->Branch("t3", t3, "t3[ncell]/D");
    
   rtree->Branch("c4", c4, "c4[ncell]/D");
   rtree->Branch("t4", t4, "t4[ncell]/D");
   rtree->Branch("clock", &clock, "clock/D");

   rtree->Branch("board", &board, "board/I");

//----------------

   // open the binary waveform file
   FILE *f = fopen(filename, "r");
   if (f == NULL) {
      printf("Cannot find file \'%s\'\n", filename);
      return 0;
   }

   // read file header
   fread(&fh, sizeof(fh), 1, f);
   if (fh.tag[0] != 'D' || fh.tag[1] != 'R' || fh.tag[2] != 'S') {
      printf("Found invalid file header in file \'%s\', aborting.\n", filename);
      return 0;
   }
   
   if (fh.version != '2') {
      printf("Found invalid file version \'%c\' in file \'%s\', should be \'2\', aborting.\n", fh.version, filename);
      return 0;
   }

   // read time header
   fread(&th, sizeof(th), 1, f);
   if (memcmp(th.time_header, "TIME", 4) != 0) {
      printf("Invalid time header in file \'%s\', aborting.\n", filename);
      return 0;
   }

   for (b = 0 ; ; b++) {
      // read board header
      fread(&bh, sizeof(bh), 1, f);
      if (memcmp(bh.bn, "B#", 2) != 0) {
         // probably event header found
         fseek(f, -4, SEEK_CUR);
         break;
      }
      
//    printf("Found data for board #%d\n", bh.board_serial_number);
      
      // read time bin widths
      memset(bin_width[b], sizeof(bin_width[0]), 0);
      for (chn=0 ; chn<5 ; chn++) {
         fread(&ch, sizeof(ch), 1, f);
         if (ch.c[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;
         }
         i = ch.cn[2] - '0' - 1;
         printf("Found timing calibration for channel #%d\n", i+1);
         fread(&bin_width[b][i][0], sizeof(float), 1024, f);
	/*my printf
		printf("bin width %d \n", bin_width[b][i][10]); */
         // fix for 2048 bin mode: double channel
         if (bin_width[b][i][1023] > 10 || bin_width[b][i][1023] < 0.01) {
            for (j=0 ; j<512 ; j++)
               bin_width[b][i][j+512] = bin_width[b][i][j];
	/*my printf
		printf("bin width %d \n", bin_width[b][i][j+512]); */
         }
      }
   }
   n_boards = b;
   
   // initialize statistics
   ndt = 0;
   sumdt = sumdt2 = 0;
   
   // loop over all events in the data file
   for (n=0 ; ; n++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;
      
      if ( !(eh.event_serial_number%100) ) {
	 printf("Found event #%d\n", eh.event_serial_number);
      }
      clock = eh.year * 365.*3600.*24.+ eh.month*30.*3600.*24.+ eh.day*3600.*24.+eh.hour*3600+eh.minute*60+eh.second+eh.millisecond*0.001;
      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {
         
         // read board header
         fread(&bh, sizeof(bh), 1, f);
         if (memcmp(bh.bn, "B#", 2) != 0) {
            printf("Invalid board header in file \'%s\', aborting.\n", filename);
            return 0;
         }
         
         // read trigger cell
         fread(&tch, sizeof(tch), 1, f);
         if (memcmp(tch.tc, "T#", 2) != 0) {
            printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
            return 0;
         }

         if (n_boards > 1&& !(eh.event_serial_number%100) )
            printf("Found data for board #%d\n", bh.board_serial_number);
         board = bh.board_serial_number;
          
         // reach channel data
         for (chn=0 ; chn<4 ; chn++) {
            
            // read channel header
            fread(&ch, sizeof(ch), 1, f);
            if (ch.c[0] != 'C') {
               // event header found
               fseek(f, -4, SEEK_CUR);
               break;
            }
            chn_index = ch.cn[2] - '0' - 1;
	//	printf("print channel %d \n",chn_index);
            fread(&scaler, sizeof(int), 1, f);
            fread(voltage, sizeof(short), 1024, f);
            
            for (i=0 ; i<1024 ; i++) {
               // convert data to volts
               waveform[b][chn_index][i] = (voltage[i] / 65536. + eh.range/1000.0 - 0.5); //calculation of amplitudes values for each cell
	
		//for ROOT
                ncell = i;
                if(chn_index == 0) {c1[i] = waveform[b][chn_index][i];}
                if(chn_index == 1) {c2[i] = waveform[b][chn_index][i];}
                if(chn_index == 2) {c3[i] = waveform[b][chn_index][i];}
                if(chn_index == 3) {c4[i] = waveform[b][chn_index][i];}
                board = bh.board_serial_number;

               // calculate time for this cell
               for (j=0,time[b][chn_index][i]=0 ; j<i ; j++){
                  time[b][chn_index][i] += bin_width[b][chn_index][(j+tch.trigger_cell) % 1024];
               }
            }
         } // end of the channel loop (chn)
         // align cell #0 of all channels
         tt1 = time[b][0][(1024-tch.trigger_cell) % 1024];
	//my print;
	// printf("t1 %1.6lf \n",time[b][0][(1024-tch.trigger_cell) % 1024]);
         for (chn=1 ; chn<4 ; chn++) {
            tt2 = time[b][chn][(1024-tch.trigger_cell) % 1024];
//adding channels 3 and 4
            tt3 = time[b][chn][(1024-tch.trigger_cell) % 1024];
            tt4 = time[b][chn][(1024-tch.trigger_cell) % 1024];
	//my prinf
	//printf("t2 %f for %d %d %d \n",time[b][chn][(1024-tch.trigger_cell) % 1024], b, chn, i);
            dt = tt1 - tt2;
            dt34 = tt3 - tt4;
            for (i=0 ; i<1024 ; i++) {
               time[b][chn][i] += dt;  //each element of time gets dt correction
		//my print;
	      // printf("time %1.6lf for %d %d %d \n",time[b][chn][i], b, chn, i);
	    }
		
         }
         tt1 = tt2 = tt3 = tt4 = 0;
         threshold = -0.045; //my threshold, used to be 0.3


	    //for ROOT
	 for(i=0 ; i<1024 ; i++) {
	    t1[i] = time[b][0][i];
	    t2[i] = time[b][1][i];
        t3[i] = time[b][0][i];
        t4[i] = time[b][1][i];
	 }

             // find peak in channel 1 above threshold
         for (i=0 ; i<1022 ; i++) {

            if (waveform[b][0][i] < threshold && waveform[b][0][i+1] >= threshold) {
               tt1 = (threshold-waveform[b][0][i])/(waveform[b][0][i+1]-waveform[b][0][i])*(time[b][0][i+1]-time[b][0][i])+time[b][0][i];
		//my prinf
		//printf("t1 recalc %1.6lf %d \n",t1, i);
               break;
            }

	 }
         
         // find peak in channel 2 above threshold
         for (i=0 ; i<1022 ; i++) {
            if (waveform[b][1][i] < threshold && waveform[b][1][i+1] >= threshold) {
               tt2 = (threshold-waveform[b][1][i])/(waveform[b][1][i+1]-waveform[b][1][i])*(time[b][1][i+1]-time[b][1][i])+time[b][1][i];
		//my prinf
		//printf("t2 recalc %1.6lf \n",t2);
               break;
            }
	 }

         // calculate distance of peaks with statistics
         if (tt1 > 0 && tt2 > 0) {
            ndt++;
            dt = tt2 - tt1;
            sumdt += dt;
            sumdt2 += dt*dt;
         }
          rtree->Fill();
      } //end of the boards loop
   } // end of the events loop
   
   // print statistics
   printf("dT = %1.3lfns +- %1.1lfps\n", sumdt/ndt, 1000*sqrt(1.0/(ndt-1)*(sumdt2-1.0/ndt*sumdt*sumdt)));

   rfile->Write();
   rfile->Close();
   return 1;
}

