#!/bin/env python
###########################################
#  Author: Jean Jacquemier                #
#  Contact: jean.jacquemier@lapp.in2p3.fr #
#  last modified: March 18 2011           #
###########################################


import pprint
import os , sys
import getopt
from datetime import date


 
from ROOT import TCanvas,TH1F

#############################################################
def getUserTimeInFile(path,file):
  """return float value corresonding to User time (second)
  print in Geant4 application output"""

  file = open(path+"/"+file,"r")
  for line in file.readlines():
    if ( line.find("User") != -1 ):
      start = line.find("User") + 5
      end = line.find("Real")-2
      user = line[start:end]
  
  file.close()
  return user

#############################################################
def getReleaseInFile(path,file):
  """return string value corresonding to geant4 release 
  print in Geant4 application output"""

  file = open(path+"/"+file,"r")
  for line in file.readlines():
    if ( line.find("Geant4 version") != -1 ):  #Geant4 version geant4-09-04-ref-0
      start = line.find("geant4")
      end = line.find("(")-1
      release = line[start:end]
      release = release.strip()
  
  
  file.close()
  return release
#############################################################
def getReleaseResult(release):
  """return new dictionary containing file name without 
  extension as key and User time as value for a specific
  directory"""

  list = []
  try :
    # get directory listing
    list = os.listdir(release)

    print list
    # fill map with release and user time
    map = {}
    for file in list :
      user = getUserTimeInFile(release,file)
      file = file[0:len(file)-4]
      map[file] = float(user)

    return map

  except:
    print ( "Directory %s directory does not exit") % (release)
    sys.exit(0)

#############################################################
def addNewResult(release,dictionary,new_map):
  """Add new benches results or remplace if already exists,
   for a specific release in global dictionary"""

  for bench,result in new_map.iteritems():
    try:
      dictionary[bench][release]=result
    except:
      # first time dictonary is empty
      dictionary[bench]={}
      dictionary[bench][release]=result
  return dictionary

#############################################################
def printAllResult(base,dbfilename):
  """"print all benches results"""
  print getTextResult(base,dbfilename)

#############################################################
def getTextResult(base,dbfilename):
  """print all benches results"""


  d = date.today()
  map = loadMap(dbfilename)
  text =       "      ==============================\n"
  text=text +  "      Results of EM CPU benchmarking\n"
  text=text +  "      Created  "+ d.strftime("%d %B %Y") + "\n"
  text=text +  "      ==============================\n"
  text=text +  "\nApplications\n"
  text=text +  "  em1: 10 GeV e- in matrix 5x5 of PbWO4 crystals (CMS-type) cut = 0.7 mm, 100 events\n"
  text=text +  "  em2: 10 GeV e- in ATLAS barrel type sampling calorimeter cut = 0.7 mm 100 events\n"
  text=text +  "  em3: 10 GeV e- in ATLAS barrel type sampling calorimeter cut = 0.02 mm 10 events\n"
  
  text=text +  "\nPhysics lists\n"
  text=text +  "  standard: emstandard PhysicList\n"
  text=text +  "  emv: emstandard_opt1 PhysicList\n"
  text=text +  "  emx: emstandard_opt2 PhysicList\n"


  try:
    resultModule = __import__(dbfilename)
    operating = getattr(resultModule,"operating")
  except(ImportError,AttributeError):
    error = "file "+ file + ".py containing database does not exist"
    sys.exit(error)

  text=text +  "\n\n----------------------Operating system informations----------------------\n"
  text=text +  "Benches runnuing on " + operating
  text=text +  "\n--------------------------------------------------------------------------"



  for bench,results in map.iteritems():
    ref = results[base]
    text=text +  "\n\n--- Bench %s" % ( bench)
    for g4version in sorted(results.keys()):
      normalized = results[g4version] / ref
      user = results[g4version]
      percent = normalized * 100 - 100
      text=text +  "\n    -- %s => Normalized=%f(%f" % ( g4version,normalized,percent ) + "%) "
      text=text +  ", User time=%f sec "  % ( user  )

  return text

#############################################################
def constructReleaseName(major,minor,type,devNum):
  """return a string containing major,minor,type and devNum of
  a release"""
  name = "release_" + major + "." + minor + "." + type + devNum
  return name

#############################################################
def constructBenchName(benchType,physicList):
  """return a stringa benchType and physics list for a bench"""
  name = benchType + "_" + physicList
  return name

#############################################################
def rootPlot(dictionary,bench):
  """create a TH1F root histogram for a specific bench.
  Save histogram to gif file"""

  # create TCanvas
  c1 = TCanvas( 'c1','', 200, 10, 700, 500 )
  # add grid on y axe
  c1.SetGridy();
  # get number of dictionary entry 
  nx = len(dictionary)
  # create TH1F histogram
  histo = TH1F(bench,  bench + ': time normalization versus Geant4.9.4-ref00 ', nx,0,nx)
  histo.SetLineStyle(1)
  histo.SetMarkerStyle(3)


  # fill histogram with sorted release keys
  i = 1
  for release in sorted(dictionary.keys()):
    normalized = dictionary[release]
    histo.SetBinContent(i,float(normalized))
    label = release.rsplit("geant4-09-")[1]
    histo.GetXaxis().SetBinLabel(i,label)
    histo.LabelsOption("v")
    histo.GetXaxis().SetLabelSize(0.03)
    i = i + 1

  # remove stat panel
  histo.SetStats(0)

  # draw histogram with line through bin contents and poly markers
  histo.Draw("PL")
  c1.Print(bench+".gif");

#############################################################
def plotHisto(base,dictionary,bench):
  """"compute normalized User time then launch rootPlot() for a specific bench"""
  res = {}
  benchmap = dictionary[bench]

  # get reference User timer for the base bench
  ref = benchmap[base]
  
  # compute normalized User time
  for release,time in benchmap.iteritems():
    normalized = time/ref
    res[release]=normalized

  # launch root plot
  rootPlot(res,bench)

#############################################################
def loadMap(file):
  """load dictionary containing benches result from a file"""

  mapName = 'benchmark'
# try to load dictionary from previous test results
  try:
    resultModule = __import__(file)
  except(ImportError,AttributeError):
    print  "file "+ file + ".py containing database does not exist,create it now"
    map={}
    return map

  try:
    map = getattr(resultModule,mapName)
  except(ImportError,AttributeError):
    sys.exit("Attribute "+ mapName +" not found")

  try:
    op = getattr(resultModule, "operating")
  except(ImportError,AttributeError):
    sys.exit("Attribute operating not found")


  if ( operating != op):
    sys.exit("Warning other benches was not ran on this system")

  return map

#############################################################
def saveTestDb(map,database):
  """save dictionary contaning benches result to a file"""

  mapString = "benchmark = " + pprint.pformat(map)
  fullPath = database + ".py"
  file = open(fullPath, 'w')

  system = "operating =\"" + operating +"\"\n"
  file.write(system)
  file.write(mapString)
  file.close()

#############################################################
if __name__ == '__main__':
 
#constants definition

  # define 3 physics list which composed result name
  standard  = 'standard' # emstandard
  emv       = 'emv'      # emstandard_opt1
  emx       = 'emx'      # emstandard_opt2

  #define 3 Genat4 applications  which composed result name
  em1       = 'em1'  # cms10gev
  em2       = 'em2'  # atlasbar
  em3       = 'em3'  # atlasbar_20um

  #database file name
  dbfilename = 'EMdatabase'

#release used to normalize results
  base = "geant4-09-04-ref-00"

#operating system 
  system = os.uname()  #(system, node, release, version, machine, processor) 
  operating = system[1] + ": " + system[0] + " " + system[2] + ", processor " + system[4]
  
  #create list contsining the nine benches
  benches = ( constructBenchName(em1,standard) , constructBenchName(em1,emx), constructBenchName(em1,emv ),
constructBenchName(em2,standard) , constructBenchName(em2,emx), constructBenchName(em2,emv ),
constructBenchName(em3,standard) , constructBenchName(em3,emx), constructBenchName(em3,emv ))

  showGif = True
 # get command arguments
  if len(sys.argv) < 2 :
    if ( len(sys.argv) == 2):
      arg1 = sys.argv[1]
      if ( arg1 == "-h" ):
        print ""
        print "------------------ benchmark Help Menu --------------------"
        print ""
        print " compute normalised User time, save results to database file and create histograms:"
        print "   benchmark.py dir    :Example for geant4.9.ref00: benchmark.py geant4.9.ref00"
        print ""
        print " Show all existing results"
        print "   benchmark.py -s"
        print "---------------------------------------------------------"
        sys.exit()
      elif ( arg1 == "-s" ):
        printAllResult(base,dbfilename)
        sys.exit()
      else:
        sys.exit("Help : benchmark.py -h")
    else:
      sys.exit("Help : benchmark.py -h")
  else:
    dir  = sys.argv[1]
    if len(sys.argv) > 2:
      if sys.argv[2] == "-w":
        showGif = False


# load global dictionary from file
    dictionary=loadMap(dbfilename)

# get benchmark result from new results files
    release = getReleaseInFile(dir,"em1_standard.out")#constructReleaseName(maj,min,type,dev)
    new_map =  getReleaseResult(release)
    dictionary = addNewResult(release,dictionary,new_map)
# save new result to database file    
    saveTestDb(dictionary,dbfilename)


# plot with pyroot and save gif files
  for bench in benches:
    if showGif: 
      plotHisto(base,dictionary,bench)

    text = getTextResult(base,dbfilename)
    file = open("cpu_result.txt" , 'w')
    file.write(text)
    file.close()



