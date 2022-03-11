import ROOT
import numpy as np
from array import array
import statistics



def geti0(c, imax, pmax, cfd = 50., timeBin = 0.4):
  i0 = int(0)
  epsilon = 0.0000000001
  for i in range(imax):
     if c[imax -i] < cfd/100.*pmax:
       i0 = imax -i
  delta = c[i0+1] - c[i0] +epsilon #voltage difference between successive points, epsilon added to avoid division by 0
  delta50 = (pmax*cfd/100. - c[i0]) #votlage difference between 50% mark
  deltai =   delta50/delta
  i0 = float(i0) +deltai
  return i0


  
def read_single_caen(input, output,  iteration, timeBin = 0.4, ):
    #use arparser here...
    #hard-coded for now...
    outfile = ROOT.TFile(output, 'RECREATE')
    tree = ROOT.TTree('tree'+str( iteration),'tree'+str( iteration))

    w, h = 1, 1024
    c = array('f', [0.]*1024)#[][ [ 0 for y in range( 1024 ) ] for x in range( nch ) ]
    t =  array('f', [0.]*1024)#[]
    pmax = array('f', [0.])
    tmax = array('f', [0.])
    t0 = array('f', [0.])
    i0 = array('f', [0.])
    
    tree.Branch("t", t, "t[1024]/F")
    tree.Branch("c", c, "c[1024]/F")
    tree.Branch("pmax", pmax, "pmax/F")
    tree.Branch("tmax", tmax, "tmax/F")
    tree.Branch("t0", t0, "t0/F")
  #  tree.Branch("trigID", trigID, "trigID/I")

    f = open(input, 'r')
    lines = f.readlines()

    for i in range(1024):
        t[i]=( (timeBin *i))
    i = 0
    evt = 0
    for line in lines:
      try:
        float(line)
      except:
         # print(line)
          continue
      c[i] = float(line) *0.25 #*-1 #negative one to change pulse polarity, 0.25 for ADV to mV
      i = i+1
      if i == 1024:
       # print(c)
        pedestal = statistics.mean(c[:50])
        for j in range(1024): c[j]= c[j] - pedestal
        pmax[0] = max(c)
        imax = c.index(pmax[0])
        tmax[0] = timeBin* imax
        i0 = geti0(c, imax, pmax[0])
        t0[0] = i0*timeBin
        #fill the output tree
        tree.Fill()
        i = 0
        #reset voltage array
        c = array('f', [0.]*1024)
        evt= evt+1
        if evt%100==0:
          print(str(evt)+" events processed")
          print("pmax: "+str(pmax[0])+", tmax: "+str(tmax[0])+" baseline correction:"+str(pedestal)+" t0:"+str(t0[0]))

       #   print(temp)
      #write tree of time and voltage values to file
    outfile.Write()
    f.close()


def mergeCAEN( run = '0000'):
  hdir =  '/Users/patfreeman/Desktop/UCSB/DESY_LGADs/wavedump/CAEN/run'+run+'/'

  out = ROOT.TFile(hdir+'CAENmergerd'+'_'+run+'.root', 'RECREATE')
  
 
  for i in range(16):
    #read the root file
    f = ROOT.TFile(hdir+'run'+run+'_'+str(i)+'.root', 'READ')
    oldtree = f.Get('tree'+str(i))
    out.cd()
#    folder = ROOT.TDirectory('c'+str(i),'c'+str(i))
    newtree = oldtree.CloneTree()
    out.Write()
    print('tree written')
    
    
def getTrigID(c, t0, level = 250., timeBin = 0.4):
    trigID = 0
    for i in range(6): # read 7 bits for now
       i0 = int( (t0 + 25./2. + 25. *i)/timeBin)
       bit = round ( statistics.mean(c[i0-2:i0+2])  / level  ) # mean of the center of the bitstream Vpp
       trigID = trigID + bit * pow(2, i)
       
     
def addTrigID(merged, trigch, trigIDch):

  trigID = array('i', [0])

  out = ROOT.TFile('trigID'+merged, RECREATE)
  newtree = ROOT.TTree('treeTrigID')
  newtree.Branch("trigID",trigID, "trigID/I")
  
  f = ROOT.TFile(merged, 'READ')
  trigTree = f.Get('tree'+str(trigch))
  IDTree = f.Get('tree'+str(trigIDch))
  #get t0 branch from trigger, bitstream begins 25 ns after this
  for evt in trigTree:
    t0 = trigTree.Get("t0") + 25
    c = IDtree.Get("c")
    trigID = getTrigID(c,t0)
    newtree.Fill()
  out.cd()
  newtree.Write()
  trigTree.Write()
  IDtree.Write()

  #copy old trees to new file

# main part of code, to be executed
# convert files to root files
for i in range(16):
  read_single_caen('/Users/patfreeman/Desktop/UCSB/DESY_LGADs/wavedump/CAEN/run0000/wave_'+str(i)+'.txt', '/Users/patfreeman/Desktop/UCSB/DESY_LGADs/wavedump/CAEN/run0000/run0000_'+str(i)+'.root', i, 0.4)
  
mergeCAEN('0000')

addTrigID(CAENmergerd_0000.root, 6, 7 )

#merge root files, renaming the branch c to c# or putting them in folders

#open ROOT file for writing



#run analysis on each folder to get t0, tmax, pmax
#get trig ID from a specific branch/channel
  
  
  






