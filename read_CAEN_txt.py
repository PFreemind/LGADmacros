import ROOT
from array import array
import statistics
import argparse
#source root before running:  source /opt/root/bin/thisroot.sh 
#written for ROOT with python 3.8, earlier python versions may not work (though likely will)

def geti0(c, imax, pmax, cfd = 50., timeBin = 0.4):
  i0 = imax
  deltai=0
  epsilon = 0.0000000001
  for i in range(imax):
     if c[imax -i] < float(cfd/100.)*pmax:
       i0 = imax -i 
       delta = c[i0+1] - c[i0] +epsilon #voltage difference between successive points, epsilon added to avoid division by 0 (though this should be impossile given the conditional?)
       delta50 = (pmax*cfd/100. - c[i0]) #votlage difference between 50% mark
       deltai =   delta50/delta
       break 
  i0 = float(i0) +deltai # interpolation is a little buggy, fix this later
  return i0

def read_single_caen(inF, output,  iteration, polarity, timeBin = 0.4, ):
    #use arparser here...
    #hard-coded for now...
    outfile = ROOT.TFile(output, 'RECREATE')
    tree = ROOT.TTree('tree'+str( iteration),'tree'+str( iteration))
    c = array('f', [0.]*1024)
    temp = array('f', [0.]*1024)
    t =  array('f', [0.]*1024)
    pmax = array('f', [0.])
    tmax = array('f', [0.])
    t0 = array('f', [0.])
    i0 = array('f', [0.])
    
    tree.Branch("t", t, "t[1024]/F")
    tree.Branch("c", c, "c[1024]/F")
    tree.Branch("pmax", pmax, "pmax/F")
    tree.Branch("tmax", tmax, "tmax/F")
    tree.Branch("t0", t0, "t0/F")

    f = open(inF, 'r')
    lines = f.readlines()

    for i in range(1024):
        t[i]=( (timeBin *i))
    i = 0
    evt = 0
    for line in lines:
 #     try:   #try statement if headers are included, as of now they typically are not
  #      float(line)
  #    except:
  #       # print(line)
  #      continue
      temp[i] = float(line) *0.25 *polarity#*-1 #negative one to change pulse polarity, 0.25 for ADV to mV
      i = i+1
      if i == 1024:
 #       print(c[-1])
        pedestal = statistics.mean(temp[:50])
        #print ("the pedestal is: "+str(pedestal))
        for j in range(1024): c[j]= float(temp[j] - pedestal)
        pmax[0] = max(c)
        imax = c.index(pmax[0])
        tmax[0] = timeBin* imax
        i0 = geti0(c, imax, pmax[0])
        t0[0] = i0*timeBin
        #fill the output tree
        tree.Fill()
        i = 0
        #reset voltage array
        #c = array('f', [0.]*1024)
        #temp = array('f', [0.]*1024)
        evt= evt+1
        if evt%1000==0:
          print(str(evt)+" events processed for channel "+str(iteration))
          print("pmax: "+str(pmax[0])+", tmax: "+str(tmax[0])+" baseline correction:"+str(pedestal)+" t0:"+str(t0[0]))

       #   print(temp)
      #write tree of time and voltage values to file
    print ("processed "+str(evt)+" events for channel "+str(iteration)) 
    f.close()  
    outfile.Write()
    outfile.Close()

def mergeCAEN( run = '0000', hdir='/data/LGADwaveforms/TB/CAEN/' ):
  dpath = hdir+'Run'+run+'/'
  out = ROOT.TFile(dpath+'CAENmergerd'+'_'+run+'.root', 'RECREATE')
  print('merging root files')
  for i in range(17):
    #read the root file
    f = ROOT.TFile(dpath+'run'+run+'_'+str(i)+'.root', 'READ')
    oldtree = f.Get('tree'+str(i))
    out.cd()
#    folder = ROOT.TDirectory('c'+str(i),'c'+str(i))
    newtree = oldtree.CloneTree()
    out.Write()
    print('tree written for channel '+str(i))
    
def getTrigID(c, t0, level = 350., timeBin = 0.4): #funciton to read bitstream of TLU, 7 bits at 40 MHz (25 ns period) assumed
    trigID = 0
    for i in range(7): # read 7 bits for now
       i0 = int( (t0 + 25./2. + 25. *i)/timeBin) # int( (t0 + 25./2.)/timeBin) +63*i#
       if i0 >1021: i0 = 1021
   #    print("i0 is "+str(i0)+", the time value is "+str(c[i0]) )
       bit = round ((c[i0-1]+c[i0]+c[i0+1])/3/level) #statistics.mean(c[i0-2:i0+2])  / level  ) # mean of the center of the bitstream Vpp
       trigID = trigID + bit * pow(2, i)
   # print ("the triggerID is:"+str(trigID))
    return int(trigID)
       
def addTrigID(merged, trigch, trigIDch, timeBin = 0.4):#function to add trigger IDs to merged trees
  f = ROOT.TFile(merged, 'READ')
  out = ROOT.TFile(merged+'.trigID.root', 'RECREATE')
  for i in range(17):
    oldtree = f.Get('tree'+str(i))
    out.cd()
    tree= oldtree.CloneTree()
    out.Write()
  trigID = array('i', [0])

  newtree = ROOT.TTree('treeTrigID', 'treeTrigID')
  newtree.Branch("trigID",trigID, "trigID/I")
  trigTree = f.Get('tree'+str(trigch))
  IDTree = f.Get('tree'+str(trigIDch))
  #get t0 branch from trigger, bitstream begins 25 ns after this
  nEvt = IDTree.GetEntries()
  for i in range(nEvt):
    trigTree.GetEntry(i)
    t0 = trigTree.t0 + 25
    IDTree.GetEntry(i)
    c = IDTree.c
  #  print(c[640])
    trigID[0] = int(getTrigID(c,t0, 350., timeBin))
    #print(trigID)
    newtree.Fill()
  out.cd()
  out.Write()

  out.Close()
  f.Close()

#def cleanTrees() #function to read merged trees and combine in a single tree with new names for branches and c and t arrays removed
#open file
#loop over trees
#delete old trees

# main part of code, to be executed
#could do argparsing here for run, path, and time bin
parser = argparse.ArgumentParser(description='convert .txt files from CAEN D5742')
parser.add_argument('-r', '--run', type=str, help='the run number', default = '0117')
args = parser.parse_args()

run = args.run#'0117'
dpath = '/data/LGADwaveforms/TB/CAEN/Run'+run+'/'
timeBin = 0.4
polarity = [-1,  #polarity of pulses 
-1,
-1,
-1,
-1,
-1, 
-1, 
1, #triggerID
-1,
-1,
-1,
-1,
-1,
-1,
-1,
-1,
1  #TLU triggerpulse
]

# convert files to root files 

read_single_caen(dpath+'TR_0_0.txt', dpath+'/run'+run+'_16.root', 16, polarity[16], timeBin) #trigger file
for i in range(16):
  read_single_caen(dpath+'wave_'+str(i)+'.txt', dpath+'/run'+run+'_'+str(i)+'.root', i, polarity[i], timeBin)

mergeCAEN(run)
addTrigID(dpath+'CAENmergerd_'+run+'.root', 16, 7 )

#merge root files, renaming the branch c to c# or putting them in folders

#open ROOT file for writing



#run analysis on each folder to get t0, tmax, pmax
#get trig ID from a specific branch/channel
'''
def read_CAEN_binary(inF, output,  iteration, polarity, timeBin = 0.4, ):
    #use arparser here...
    #hard-coded for now...
    outfile = ROOT.TFile(output, 'RECREATE')
    tree = ROOT.TTree('tree'+str( iteration),'tree'+str( iteration))
    c = array('f', [0.]*1024)
    temp = array('f', [0.]*1024)
    t =  array('f', [0.]*1024)
    pmax = array('f', [0.])
    tmax = array('f', [0.])
    t0 = array('f', [0.])
    i0 = array('f', [0.])
    
    tree.Branch("t", t, "t[1024]/F")
    tree.Branch("c", c, "c[1024]/F")
    tree.Branch("pmax", pmax, "pmax/F")
    tree.Branch("tmax", tmax, "tmax/F")
    tree.Branch("t0", t0, "t0/F")

    f = open(inF, 'rb')
    lines = f.readlines()
    
    try:
    byte = f.read(1)
    while byte != "":
        # Do stuff with byte.
        byte = f.read(1)
    finally:
      f.close()
'''
'''

    for i in range(1024):
        t[i]=( (timeBin *i))
    i = 0
    evt = 0
    for line in lines:
 #     try:   #try statement if headers are included, as of now they typically are not
  #      float(line)
  #    except:
  #       # print(line)
  #      continue
      temp[i] = float(line) *0.25 *polarity#*-1 #negative one to change pulse polarity, 0.25 for ADV to mV
      i = i+1
      if i == 1024:
 #       print(c[-1])
        pedestal = statistics.mean(temp[:50])
        #print ("the pedestal is: "+str(pedestal))
        for j in range(1024): c[j]= float(temp[j] - pedestal)
        pmax[0] = max(c)
        imax = c.index(pmax[0])
        tmax[0] = timeBin* imax
        i0 = geti0(c, imax, pmax[0])
        t0[0] = i0*timeBin
        #fill the output tree
        tree.Fill()
        i = 0
        #reset voltage array
        #c = array('f', [0.]*1024)
        #temp = array('f', [0.]*1024)
        evt= evt+1
        if evt%1000==0:
          print(str(evt)+" events processed for channel "+str(iteration))
          print("pmax: "+str(pmax[0])+", tmax: "+str(tmax[0])+" baseline correction:"+str(pedestal)+" t0:"+str(t0[0]))

       #   print(temp)
      #write tree of time and voltage values to file
    print ("processed "+str(evt)+" events for channel "+str(iteration)) 
    f.close()  
    outfile.Write()
    outfile.Close()
f = open("myfile", "rb")
try:
    byte = f.read(1)
    while byte != "":
        # Do stuff with byte.
        byte = f.read(1)
finally:
    f.close()

'''
 
  
read_CAEN
  






