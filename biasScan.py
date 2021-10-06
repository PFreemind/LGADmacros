#! /usr/bin/env  python
import ROOT
 
#read csv files
files = [
"run_0012.csv",
"run_0013.csv",
"run_0014.csv",
"run_0015.csv",
"run_0016.csv",
"run_0017.csv",
"run_0018.csv",
"run_0019.csv",
"run_0020.csv",
"run_0021.csv"
]

bias = [150, 140, 130, 120, 110, 100, 90, 80, 70, 60]

dir = "/Users/patrick/DESY_LGADs/data/TB/plots/"

gT = ROOT.TGraphErrors()
gRT = ROOT.TGraphErrors()
gP = ROOT.TGraphErrors()

gT.SetTitle(";Bias [V]; Time resolution [ps]")
gRT.SetTitle(";Bias [V]; Rise time [ns]")
gP.SetTitle(";Bias [V]; Landau MPV [mV]")
i=0

for fi in files:
  print "the file is"+fi
  f = open(dir+fi,"r")
  lines = f.readlines()
  f.close()
  gT.SetPoint(i, bias[i], float(lines[0].split(",")[1])*1000/pow(2,0.5))
  gT.SetPointError(i, 0, float(lines[0].split(",")[2])*1000/pow(2,0.5))
  gRT.SetPoint(i, bias[i], float(lines[1].split(",")[1]))
  gRT.SetPointError(i, 0, float(lines[1].split(",")[2]))
  gP.SetPoint(i, bias[i], float(lines[2].split(",")[1])*1000)
  gP.SetPointError(i, 0, float(lines[2].split(",")[2])*1000)
  print "the time resolution is "+lines[0].split(",")[1]
  
  i=i+1


c1 = ROOT.TCanvas("c1")
gT.SetMarkerStyle(26)
gT.Draw()
c1.SaveAs(dir+"timeResvBias.C")
c1.SaveAs(dir+"timeResvBias.pdf")

c2 = ROOT.TCanvas("c2")
gRT.SetMarkerStyle(26)
gRT.Draw()
c2.SaveAs(dir+"riseTimevBias.C")
c2.SaveAs(dir+"riseTimevBias.pdf")

c3 = ROOT.TCanvas("c3")
gP.SetMarkerStyle(26)
gP.Draw()
c3.SaveAs(dir+"riseTimevBias.C")
c3.SaveAs(dir+"ampvBias.pdf")
