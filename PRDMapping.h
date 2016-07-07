// Mapping class is implemented as a singleton

#ifndef __PRDMAPPING__
#define __PRDMAPPING__
/*******************************************************************************
*  AMORE FOR PRD - PRD                                                         *
*  PRDMapping                                                                  *
*  PRD Module Class                                                            *
*  Author: Kondo GNANVO 18/08/2010                                             *
*******************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "TList.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#endif

using namespace std;

class PRDMapping /*: public TObject*/ {

public:

  ~PRDMapping(){
    Clear() ;
    ClearMapOfList(fDetectorListFromDetectorTypeMap) ;
    ClearMapOfList(fDetectorListFromReadoutBoardMap) ;
    ClearMapOfList(fPlaneListFromDetectorMap) ;
    ClearMapOfList(fFECIDListFromPlaneMap) ;
    ClearMapOfList(fAPVIDListFromFECIDMap) ;
    ClearMapOfList(fAPVIDListFromPlaneMap) ;
    ClearMapOfList(fAPVIDListFromDetectorMap) ;
    ClearMapOfList(fAPVToPadChannelMap) ;
    ClearMapOfList(fCartesianPlaneMap) ;
    ClearMapOfList(fActiveADCchannelsMap);
  }

  template <typename M> void ClearMapOfList( M & amap ) ;

  static PRDMapping * GetInstance() {
    if (instance == 0) instance = new PRDMapping() ;
    return instance;
  }

  void PrintMapping();
  void SaveMapping(const char * mappingCfgFilename);

  void LoadMapping(const char * mappingCfgFilename);
  //void LoadAPVtoPadMapping(const char * mappingCfgFilename);

  void SetAPVMap(TString detPlane, Int_t fecId, Int_t adcCh, Int_t apvNo, Int_t apvOrient, Int_t apvIndex, Int_t apvHdr, TString defaultAPV) ;
  //void SetAPVtoPadMapping(Int_t fecId, Int_t adcCh, Int_t padId, Int_t apvCh) ;
  void SetCartesianStripsReadoutMap(TString readoutBoard, TString detectorType, TString detector, Int_t detID, TString planeX,  Float_t sizeX, Int_t connectorsX, Int_t orientX, TString planeY,  Float_t sizeY, Int_t connectorsY, Int_t orientY) ;

  Int_t GetPRadStripMapping(Int_t apvID, Int_t chNo) ;
  //This function disabled by Xb

  void ComputeActiveADCchannelsMap() ;
  vector <Int_t> GetActiveADCchannels(Int_t fecID) {return fActiveADCchannelsMap[fecID];}

  void Clear();

  map <Int_t, Int_t>  GetAPVNoFromIDMap()            {return fAPVNoFromIDMap ;}
  map <Int_t, Int_t>  GetAPVIDFromAPVNoMap()         {return fAPVIDFromAPVNoMap ;}
  map <Int_t, Int_t>  GetAPVGainFromIDMap()          {return fAPVGainFromIDMap ;}
  map <Int_t, Int_t>  GetAPVOrientationFromIDMap()   {return fAPVOrientationFromIDMap;}
  map <Int_t, UInt_t> GetAPVHeaderLevelFromIDMap()   {return fAPVHeaderLevelFromIDMap;}
  map <Int_t, Int_t>  GetAPVIndexOnPlaneFromIDMap()  {return fAPVIndexOnPlaneFromIDMap;}

  map <Int_t, TString>   GetAPVFromIDMap()           {return fAPVFromIDMap ;} 
  map <Int_t, TString>   GetDetectorFromIDMap()      {return fDetectorFromIDMap ;}

  map <TString, Int_t>   GetPlaneIDFromPlaneMap()      {return fPlaneIDFromPlaneMap ;}
  map<TString, vector<Float_t> > GetCartesianPlaneMap() {return fCartesianPlaneMap ; }
  map <TString, list <Int_t> > GetAPVIDListFromDetectorMap()      {return fAPVIDListFromDetectorMap;}
  map <TString, list <Int_t> > GetAPVIDListFromPlaneMap() {return fAPVIDListFromPlaneMap;}

  list <Int_t> GetAPVIDListFromDetector(TString detName) {
    fAPVIDListFromDetectorMap[detName].unique() ;
    return fAPVIDListFromDetectorMap[detName];
  }

  list <Int_t> GetAPVIDListFromPlane(TString planeName) {
    fAPVIDListFromPlaneMap[planeName].unique() ;
    return fAPVIDListFromPlaneMap[planeName];
  }

  list <Int_t> GetFECIDListFromPlane(TString planeName) {
    fFECIDListFromPlaneMap[planeName].unique();
    return fFECIDListFromPlaneMap[planeName];
  }

  list <Int_t> GetAPVIDListFromFECID(Int_t fecId) {
    fAPVIDListFromFECIDMap[fecId].unique();
    return fAPVIDListFromFECIDMap[fecId];
  }

  list <TString> GetPlaneListFromDetector(TString detName) {
    fPlaneListFromDetectorMap[detName].unique();
    return  fPlaneListFromDetectorMap[detName] ;
  }

  Int_t GetFECIDFromAPVID(Int_t apvID)       {return (apvID >> 4 ) & 0xF;}
  Int_t GetADCChannelFromAPVID(Int_t apvID)  {return apvID & 0xF;}

  Int_t GetAPVNoFromID(Int_t apvID)          {return fAPVNoFromIDMap[apvID];}
  Int_t GetAPVIDFromAPVNo(Int_t apvID)       {return fAPVIDFromAPVNoMap[apvID];}
  Int_t GetAPVIndexOnPlane(Int_t apvID)      {return fAPVIndexOnPlaneFromIDMap[apvID];}
  Int_t GetAPVOrientation(Int_t apvID)       {return fAPVOrientationFromIDMap[apvID];}
  Int_t GetAPVIDFromName(TString apvName)    {return fAPVIDFromNameMap[apvName];}
  UInt_t GetAPVHeaderLevelFromID(Int_t apvID) {return fAPVHeaderLevelFromIDMap[apvID];}

  Int_t GetDetectorID(TString detName) {return fDetectorIDFromDetectorMap[detName];}
  TString GetDetectorFromID(Int_t detID)   {return fDetectorFromIDMap[detID];}

  TString GetDetectorTypeFromDetector(TString detectorName)   {return fDetectorTypeFromDetectorMap[detectorName];}
  TString GetReadoutBoardFromDetector(TString detectorName)   {return fReadoutBoardFromDetectorMap[detectorName];}
  TString GetDetectorFromPlane(TString planeName)             {return fDetectorFromPlaneMap[planeName];}

  vector< Float_t > GetCartesianReadoutMap(TString plane)   {return fCartesianPlaneMap[plane] ;}

  Int_t GetPlaneID(TString planeName) {
    return (Int_t) (fCartesianPlaneMap[planeName])[0] ;
  }

  Int_t GetEtaSector(TString planeName) {
    return (Int_t) (fCartesianPlaneMap[planeName])[0] ;
  }

  Float_t GetPlaneSize(TString planeName) {
    return (fCartesianPlaneMap[planeName])[1] ;
  }

  Int_t GetNbOfAPVsOnPlane(TString planeName)  {
    return  (Int_t) (fCartesianPlaneMap[planeName])[2] ;
  } 

  Int_t GetPlaneOrientation(TString planeName)  {
    return  (Int_t) (fCartesianPlaneMap[planeName])[3] ;
  } 
 
  TString GetPlaneFromAPVID(Int_t apvID) {return fPlaneFromAPVIDMap[apvID];}

  TString GetAPVstatus(Int_t apvID) {return fAPVstatusMap[apvID];}
  TString GetAPVFromID(Int_t apvID) {return fAPVFromIDMap[apvID];}
  TString GetAPV(TString detPlane, Int_t fecId, Int_t adcChannel, Int_t apvNo, Int_t apvIndex, Int_t apvID) ;

  vector<Int_t> GetPadChannelsMapping(Int_t apvID) {return fAPVToPadChannelMap[apvID]; }
  map < Int_t, vector<Int_t> > GetPadChannelsMapping() {return fAPVToPadChannelMap;}

  Int_t GetNbOfAPVs() {return fAPVNoFromIDMap.size();}
  Int_t GetNbOfFECs() {return fAPVIDListFromFECIDMap.size();}
  Int_t GetNbOfPlane() {return fDetectorFromPlaneMap.size();}
  Int_t GetNbOfDetectors() {return fReadoutBoardFromDetectorMap.size();}

  Bool_t IsAPVIDMapped(Int_t apvID) {
    map<Int_t, TString>::iterator itr ;
    itr = fAPVFromIDMap.find(apvID) ;
    if(itr != fAPVFromIDMap.end())
      return kTRUE ;
    else
      return kFALSE ;
  }

private:

  PRDMapping() {fNbOfAPVs = 0;}

  static PRDMapping * instance;
  Int_t fNbOfAPVs ;

  map<Int_t, UInt_t>   fAPVHeaderLevelFromIDMap;
  map<Int_t, Int_t>   fAPVNoFromIDMap, fAPVIDFromAPVNoMap, fAPVIndexOnPlaneFromIDMap,fAPVOrientationFromIDMap;
  map<TString, Int_t> fNbOfAPVsFromDetectorMap ;          

  map<Int_t, Int_t>      fAPVGainFromIDMap ; 
  map<TString, Int_t>    fAPVIDFromNameMap ; 
  map<Int_t, TString>    fAPVFromIDMap ;  
  map<Int_t, TString>    fAPVstatusMap ;  

  map<TString, Int_t>    fPlaneIDFromPlaneMap ;
  map<Int_t, TString>    fPlaneFromAPVIDMap;

  map<Int_t, TString>    fReadoutBoardFromIDMap ;  

  map<TString, TString>  fDetectorTypeFromDetectorMap ;
  map<TString, TString>  fReadoutBoardFromDetectorMap ;
  map<TString, Int_t>    fDetectorIDFromDetectorMap ;  
  map<Int_t, TString>    fDetectorFromIDMap ;   
  map<Int_t, TString>    fDetectorFromAPVIDMap; 
  map<TString, TString > fDetectorFromPlaneMap ;

  map<Int_t, list<Int_t> >   fAPVIDListFromFECIDMap;       

  map<TString, list<Int_t> > fFECIDListFromPlaneMap; 
  map<TString, list<Int_t> > fAPVIDListFromDetectorMap;     
  map<TString, list<Int_t> > fAPVIDListFromPlaneMap;

  map<TString, list<TString> > fDetectorListFromDetectorTypeMap ;
  map<TString, list<TString> > fDetectorListFromReadoutBoardMap;
  map<TString, list<TString> > fPlaneListFromDetectorMap ;

  map<Int_t, vector<Int_t> >     fAPVToPadChannelMap; 
  map<Int_t, vector <Int_t> >    fActiveADCchannelsMap ;
  map<TString, vector<Float_t> > fCartesianPlaneMap ;

  //ClassDef(PRDMapping, 1)

  //Xb: Get Bank ID from Mapping
public:
  set<int> GetBankIDSet();

  TString GetFECIPFromFECID(Int_t id);

  Int_t GetFECIDFromFECIP(TString);
private:
  map<Int_t, TString> fFECIPFromFECID;
  map<TString, Int_t> fFECIDFromFECIP;
};

#endif



