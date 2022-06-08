import DESYLGADs as dl
import argparse

#source root before running:  source /opt/root/bin/thisroot.sh 
#written for ROOT with python 3.8, earlier python versions may not work (though likely will)

# main part of code, to be executed
#could do argparsing here for run, path, and time bin
parser = argparse.ArgumentParser(description='convert .txt files from CAEN D5742')
parser.add_argument('-r', '--run', type=str, help='the run number, written as 4 digits with leading zeros (ex. 0117)', default = '0117')
args = parser.parse_args()

run = args.run#'0117'
dpath = '../data/TB/CAEN//Run'+run+'/'  # '/data/LGADwaveforms/TB/CAEN/Run'+run+'/'
timeBin = 0.4

# convert files to root files 

#read_single_caen(dpath+'TR_0_0.txt', dpath+'/run'+run+'_16.root', 16, polarity[16], timeBin) #trigger file

dl.read_caen_binary(dpath+'TR_0_0.dat', dpath+'/run'+run+'_16.root', 16, dl.polarity[16], timeBin)
for i in range(16):
  dl.read_caen_binary(dpath+'wave_'+str(i)+'.dat', dpath+'/run'+run+'_'+str(i)+'.root', i, dl.polarity[i], timeBin)

dl.mergeCAEN(run,'../data/TB/CAEN/' )

dl.addTrigID(dpath+'CAENmergerd_'+run+'.root', 16, 7 )


  






