#! /usr/bin/env  python
import ROOT
 
#read csv files
files = [
"run_0048.csv",
"run_0047.csv",
"run_0039.csv",
"run_0040.csv",
"run_0041.csv",
"run_0042.csv",
"run_0043.csv",
"run_0044.csv",
"run_0045.csv",
"run_0046.csv",
]

bias = [160, 155, 150, 145, 140, 135, 130, 120, 110, 100]

dir = "/data/LGADwaveforms/TB/plots/"

gT = ROOT.TGraphErrors()
gRT = ROOT.TGraphErrors()
gP = ROOT.TGraphErrors()
gC = ROOT.TGraphErrors()

gT.SetTitle(";Bias [V]; Time resolution [ps]")
gRT.SetTitle(";Bias [V]; Rise time [ns]")
gP.SetTitle(";Bias [V]; Landau MPV [mV]")
gC.SetTitle(";Bias [V]; Landau MPV [C]")
i=0

for fi in files:
  print "file: "+fi
  f = open(dir+fi,"r")
  lines = f.readlines()
  f.close()
  gT.SetPoint(i, bias[i], float(lines[0].split(",")[1])*1000/pow(2,0.5))
  gT.SetPointError(i, 0, float(lines[0].split(",")[2])*1000/pow(2,0.5))
  gRT.SetPoint(i, bias[i], float(lines[1].split(",")[1]))
  gRT.SetPointError(i, 0, float(lines[1].split(",")[2]))
  gP.SetPoint(i, bias[i], float(lines[2].split(",")[1])*1000)
  gP.SetPointError(i, 0, float(lines[2].split(",")[2])*1000)
  gC.SetPoint(i, bias[i], float(lines[3].split(",")[1]))
  gC.SetPointError(i, 0, float(lines[3].split(",")[2]))
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
c3.SaveAs(dir+"ampvBias.C")
c3.SaveAs(dir+"ampvBias.pdf")

c4 = ROOT.TCanvas("c4")
gC.SetMarkerStyle(26)
gC.Draw()
c3.SaveAs(dir+"chargevBias.C")
c3.SaveAs(dir+"chargevBias.pdf")
