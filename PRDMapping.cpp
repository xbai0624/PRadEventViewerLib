#include "PRDMapping.h"

PRDMapping * PRDMapping::instance = 0 ;

void PRDMapping::SetCartesianStripsReadoutMap(TString readoutBoard, TString detectorType, TString detector, Int_t detID, TString planeX,  Float_t sizeX, Int_t connectorsX, Int_t orientX, TString planeY,  Float_t sizeY, Int_t connectorsY, Int_t orientY) {

  printf("   PRDMapping::SetDetectorMap() => readout=%s, detType=%s, det=%s, detID=%d, planeX=%s, SizeX=%f, connectorsX=%d, orientationX=%d, planeY=%s, SizeY=%f, connectorsY=%d, orientationY=%d \n", readoutBoard.Data(), detectorType.Data(), detector.Data(), detID, planeX.Data(), sizeX, connectorsX, orientX, planeY.Data(), sizeY, connectorsY, orientY) ;

  fDetectorFromIDMap[detID]     = detector ;
  fReadoutBoardFromIDMap[detID] = readoutBoard ; 

  fDetectorIDFromDetectorMap[detector] = detID ;
  fReadoutBoardFromDetectorMap[detector] = readoutBoard ;
  fDetectorTypeFromDetectorMap[detector] = detectorType ;

  fDetectorListFromDetectorTypeMap[detectorType].push_back(detector) ;
  fDetectorListFromReadoutBoardMap[readoutBoard].push_back(detector) ;

  fPlaneIDFromPlaneMap[planeX] = 0 ;
  fPlaneIDFromPlaneMap[planeY] = 1 ;

  fDetectorFromPlaneMap[planeX] = detector ;
  fDetectorFromPlaneMap[planeY] = detector ;

  fPlaneListFromDetectorMap[detector].push_back(planeX) ;
  fPlaneListFromDetectorMap[detector].push_back(planeY) ;

  fCartesianPlaneMap[planeX].push_back(0) ;
  fCartesianPlaneMap[planeX].push_back(sizeX) ;
  fCartesianPlaneMap[planeX].push_back(connectorsX) ;
  fCartesianPlaneMap[planeX].push_back(orientX) ;
  fCartesianPlaneMap[planeY].push_back(1) ;
  fCartesianPlaneMap[planeY].push_back(sizeY) ;
  fCartesianPlaneMap[planeY].push_back(connectorsY) ;
  fCartesianPlaneMap[planeY].push_back(orientY) ;
}

TString PRDMapping::GetAPV(TString detPlane, Int_t fecId, Int_t adcCh, Int_t apvNo, Int_t apvIndex, Int_t apvID) {
  stringstream out ;
  out << apvID ;
  TString apvIDStr = out.str();
  out.str("") ;
  out << fecId ;
  TString fecIDStr = out.str();
  out.str("") ;
  out <<  adcCh;
  TString adcChStr = out.str();
  out.str("") ;
  out <<  apvNo ;
  TString apvNoStr = out.str();
  out.str("") ;
  out <<  apvIndex ;
  TString apvIndexStr = out.str();
  out.str("") ;
  TString apvName = "apv" + apvNoStr + "_Id" + apvIDStr + "_adcCh" + adcChStr + "_FecId" + fecIDStr + "_" + detPlane ;
  return apvName ;
}

void  PRDMapping::SetAPVMap(TString detPlane, Int_t fecId, Int_t adcCh, Int_t apvNo, Int_t apvOrient, Int_t apvIndex, Int_t apvHdr, TString defaultAPV) {

  Int_t apvID = (fecId << 4) | adcCh ;

  TString apvName = GetAPV(detPlane, fecId, adcCh, apvNo, apvIndex, apvID) ;

  fAPVNoFromIDMap[apvID]           = apvNo ;
  fAPVIDFromAPVNoMap[apvNo]        = apvID ;
  fAPVFromIDMap[apvID]             = apvName ;
  fAPVstatusMap[apvID]             = defaultAPV ;
  fAPVHeaderLevelFromIDMap[apvID]  = apvHdr ;
  fAPVOrientationFromIDMap[apvID]  = apvOrient ;
  fAPVIndexOnPlaneFromIDMap[apvID] = apvIndex ;

  fAPVIDFromNameMap[apvName] = apvID ;
  fPlaneFromAPVIDMap[apvID]  = detPlane ;

  fAPVIDListFromFECIDMap[fecId].push_back(apvID);
  fFECIDListFromPlaneMap[detPlane].push_back(fecId);
  fAPVIDListFromPlaneMap[detPlane].push_back(apvID);

  TString detector = GetDetectorFromPlane(detPlane) ;
  fAPVIDListFromDetectorMap[detector].push_back(apvID);
}

void PRDMapping::PrintMapping() {

  map<TString, list<TString> >::const_iterator det_itr ;
  for(det_itr = fPlaneListFromDetectorMap.begin(); det_itr != fPlaneListFromDetectorMap.end(); ++det_itr) {
    TString detector = det_itr->first ;
    printf("  ======================================================================================================================\n") ;
    printf("  PRDMapping::PrintMapping() ==> Detector = %s \n",detector.Data()) ;

    list<TString> detPlaneList = det_itr->second ;
    list<TString>::const_iterator plane_itr ;
    for(plane_itr = detPlaneList.begin(); plane_itr != detPlaneList.end(); ++plane_itr) {
      TString detPlane = * plane_itr ;
      list <Int_t> fecList = GetFECIDListFromPlane(detPlane) ;
      list<Int_t>::const_iterator fec_itr ;
      for(fec_itr = fecList.begin(); fec_itr != fecList.end(); ++fec_itr) { 
	Int_t fecId = * fec_itr ;
	printf("  PRDMapping::PrintMapping()     Plane = %s,        FEC = %d \n", detPlane.Data(), fecId) ;

	list <Int_t> apvList = GetAPVIDListFromPlane(detPlane) ;
	list<Int_t>::const_iterator apv_itr ;
	for(apv_itr = apvList.begin(); apv_itr != apvList.end(); ++apv_itr) { 
	  Int_t apvID       = * apv_itr ;
	  Int_t apvNo       = GetAPVNoFromID(apvID);
	  Int_t apvIndex    = GetAPVIndexOnPlane(apvID);
	  Int_t apvOrient   = GetAPVOrientation(apvID);
	  Int_t fecID       = GetFECIDFromAPVID(apvID);
	  Int_t adcCh       = GetADCChannelFromAPVID(apvID);
	  Int_t apvHdrLevel = GetAPVHeaderLevelFromID(apvID);
	  TString  apvName  = GetAPVFromID(apvID) ;
	  TString defautAPV = GetAPVstatus(apvID);

	  if(fecID == fecId) printf("  PRDMapping::PrintMapping() ==> adcCh=%d,  apvName=%s,  apvID=%d, apvNo=%d,  apvIndex=%d,  apvOrientation=%d,  apvHdr=%d, apvStatus=%s\n", adcCh, apvName.Data(), apvID, apvNo, apvIndex, apvOrient, apvHdrLevel, defautAPV.Data()) ;
	}
      }
      printf("\n") ;
    } 
  }
  printf("======================================================================================================================\n") ;
  printf("  PRDMapping::PrintMapping() ==> Mapping of %d detectors, %d planes, %d FECs, %d APVs\n", GetNbOfDetectors(), GetNbOfPlane(), GetNbOfFECs(), GetNbOfAPVs());
  printf("======================================================================================================================\n") ;
}

void PRDMapping::SaveMapping(const char * file) {
  printf("  PRDMapping::SaveMapping() ==> Saving PRD Mapping to file [%s],\n", file) ;
  FILE * f = fopen(file, "w");  

  fprintf(f,"#################################################################################################\n") ;
  fprintf(f,"         readoutType  Detector    Plane  DetNo   Plane   size (mm)  connectors  orientation\n");
  fprintf(f,"#################################################################################################\n") ;

  map<TString, list<TString> >::const_iterator det_itr ;
  for(det_itr = fPlaneListFromDetectorMap.begin(); det_itr != fPlaneListFromDetectorMap.end(); ++det_itr) {
    TString detector      = det_itr->first ;
    TString readoutBoard = GetReadoutBoardFromDetector(detector) ;
    TString detectorType = GetDetectorTypeFromDetector(detector) ;

    if ( (readoutBoard == "CARTESIAN") ||  (readoutBoard == "UV_ANGLE_OLD") ){
      list<TString> detPlaneList = det_itr->second ;
      TString planeX    = detPlaneList.front() ;
      vector <Float_t> cartesianPlaneX = GetCartesianReadoutMap(planeX) ;
       Float_t sizeX     = cartesianPlaneX[1] ;
      Int_t connectorsX = (Int_t) (cartesianPlaneX[2]) ;
      Int_t orientX     = (Int_t) (cartesianPlaneX[3]) ;

      TString planeY    = detPlaneList.back() ;
      vector <Float_t> cartesianPlaneY = GetCartesianReadoutMap(planeY) ;
      Float_t sizeY     = cartesianPlaneY[1] ;
      Int_t connectorsY = (Int_t) (cartesianPlaneY[2]) ;
      Int_t orientY     = (Int_t) (cartesianPlaneY[3]) ;
      fprintf(f,"DET,  %s,   %s,   %s,   %s,  %f,   %d,   %d,   %s,   %f,   %d,   %d \n", readoutBoard.Data(), detectorType.Data(), detector.Data(), planeX.Data(), sizeX, connectorsX, orientX, planeY.Data(), sizeY, connectorsY, orientY ) ;
    }

    else {
      printf("  PRDMapping::SaveMapping() ==> detector readout board type %s is not yet implemented ==> PLEASE MOVE ON \n", readoutBoard.Data()) ;
      continue ;
    }
  }

  fprintf(f,"###############################################################\n") ;
  fprintf(f,"#     fecId   adcCh   detPlane  apvOrient  apvIndex    apvHdr #\n");
  fprintf(f,"###############################################################\n") ;
  map<Int_t, TString>::const_iterator apv_itr;
  for(apv_itr = fAPVFromIDMap.begin(); apv_itr != fAPVFromIDMap.end(); ++apv_itr){
    Int_t apvID       = apv_itr->first;
    Int_t fecId       = GetFECIDFromAPVID(apvID);
    Int_t adcCh       = GetADCChannelFromAPVID(apvID);
    TString detPlane  = GetPlaneFromAPVID(apvID) ;
    Int_t apvOrient   = GetAPVOrientation(apvID);
    Int_t apvIndex    = GetAPVIndexOnPlane(apvID);
    Int_t apvHdr      = GetAPVHeaderLevelFromID(apvID);
    TString defautAPV = GetAPVstatus(apvID);
    fprintf(f,"APV,   %d,     %d,     %s,     %d,    %d,   %d,  %s, \n", fecId, adcCh, detPlane.Data(), apvOrient, apvIndex, apvHdr, defautAPV.Data());
  }
  fclose(f);
}

void PRDMapping::LoadMapping(const char * mappingCfgFilename) {

  Clear() ;
  printf("  PRDMapping::LoadDefaultMapping() ==> Loading Mapping from %s \n", mappingCfgFilename) ;
  Int_t apvNo = 0 ;
  Int_t detID = 0 ;

  ifstream filestream (mappingCfgFilename, ifstream::in); 
  TString line;
  while (line.ReadLine(filestream)) {

    line.Remove(TString::kLeading, ' ');   // strip leading spaces
    if (line.BeginsWith("#")) continue ;   // and skip comments
    //printf("   PRDMapping::LoadDefaultMapping() ==> Scanning the mapping cfg file %s\n",line.Data()) ;

    //=== Create an array of the tokens separated by "," in the line;
    TObjArray * tokens = line.Tokenize(","); 

    //=== iterator on the tokens array
    TIter myiter(tokens); 
    while (TObjString * st = (TObjString*) myiter.Next()) {

      //== Remove leading and trailer spaces
      TString s = st->GetString().Remove(TString::kLeading, ' ' );
      s.Remove(TString::kTrailing, ' ' );                         

      //      printf("    PRDMapping::LoadDefaultMapping() ==> Data ==> %s\n",s.Data()) ;
      if(s == "DET") {
	TString readoutBoard = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );
	TString detectorType = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );
	TString detector     = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );

	if (readoutBoard == "CARTESIAN")  {
	  TString planeX           = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );
	  Float_t sizeX            = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atof();
	  Int_t   nbOfConnectorsX  = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
	  Int_t   orientationX     = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();

	  TString planeY           = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );
	  Float_t sizeY            = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atof();
	  Int_t   nbOfConnectorsY  = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
	  Int_t   orientationY     = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
	  SetCartesianStripsReadoutMap(readoutBoard, detectorType, detector, detID, planeX, sizeX, nbOfConnectorsX, orientationX, planeY, sizeY, nbOfConnectorsY, orientationY) ;
	}
	else {
	  printf("XXXXXXX PRDMapping::LoadDefaultMapping()==> detector with this readout board type %s is not yet implemented ==> PLEASE MOVE ON XXXXXXXXXXX \n", readoutBoard.Data()) ;
	  continue ;
	}
	detID++ ;
      }

      if(s == "APV") {
	Int_t   fecId     = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
	Int_t   adcCh     = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
        TString detPlane  = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );
        Int_t   apvOrient = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
	Int_t   apvIndex  = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
	Int_t   apvheader = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
        TString defautAPV = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );

	if (detPlane == "NULL") continue ;
	SetAPVMap(detPlane, fecId, adcCh, apvNo, apvOrient, apvIndex, apvheader, defautAPV) ;
	apvNo++ ;
      }
      if(s == "FEC") {
	Int_t   fecId     = (((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' )).Atoi();
        TString ip = ((TObjString*) myiter.Next())->GetString().Remove(TString::kLeading, ' ' );
	fFECIPFromFECID[fecId] = ip;
	fFECIDFromFECIP[ip] = fecId;
      }
    }
   tokens->Delete();
  }
  printf("  PRDMapping::LoadDefaultMapping() ==> ComputeActiveADCchannelsMap(  \n") ;
  ComputeActiveADCchannelsMap() ;
  printf("\n  ======================================================================================================================\n") ;
}

void PRDMapping::Clear() {
  printf("  PRDMapping::Clear() ==> Clearing Previous Mapping \n") ;

  fAPVIDFromAPVNoMap.clear() ; 
  fAPVIDFromNameMap.clear() ; 
  fAPVIDListFromDetectorMap.clear() ; 
  fAPVIDListFromPlaneMap.clear() ; 
  fAPVNoFromIDMap.clear() ; 
  fAPVFromIDMap.clear() ; 
  fAPVstatusMap.clear() ; 
  fAPVHeaderLevelFromIDMap.clear() ;

  fPlaneIDFromPlaneMap.clear() ;
  fDetectorIDFromDetectorMap.clear() ; 
  fDetectorFromIDMap.clear() ; 
  fDetectorFromAPVIDMap.clear() ;  
  fDetectorFromPlaneMap.clear() ;

  fPlaneFromAPVIDMap.clear() ;
  fReadoutBoardFromIDMap.clear() ;
  fReadoutBoardFromDetectorMap.clear() ;
  fNbOfAPVsFromDetectorMap.clear() ;
  fAPVOrientationFromIDMap.clear() ;
  fAPVIndexOnPlaneFromIDMap.clear() ;
  printf("  PRDMapping::Clear() ==> Previous Mapping cleared \n") ;
}

template <typename M> void ClearMapOfList( M & amap ) {
  for ( typename M::iterator it = amap.begin(); it != amap.end(); ++it ) {
    ((*it).second).clear();
  }
  amap.clear() ;
}

void PRDMapping::ComputeActiveADCchannelsMap() {
  fActiveADCchannelsMap.clear() ;
  map <Int_t, Int_t >::const_iterator adcChannel_itr ;
  for(adcChannel_itr = fAPVNoFromIDMap.begin(); adcChannel_itr != fAPVNoFromIDMap.end(); ++adcChannel_itr) {
    Int_t apvid = (* adcChannel_itr).first ;
    Int_t activeChannel = apvid & 0xF;
    Int_t fecId = (apvid >> 4 ) & 0xF;
    fActiveADCchannelsMap[fecId].push_back(activeChannel) ;
  }
}

Int_t PRDMapping::GetPRadStripMapping(Int_t apvID, Int_t chNo) {
  //------------ APV25 Internal Channel Mapping
  chNo = (32 * (chNo%4)) + (8 * (Int_t)(chNo/4)) - (31 * (Int_t)(chNo/16)) ;

  //------------ APV25 Channel to readout strip Mapping
  if ( (GetPlaneFromAPVID(apvID).Contains("X")) && (GetAPVIndexOnPlane(apvID)  == 11) ) {
    if (chNo % 2 == 0) 
      chNo = ( chNo / 2) + 48 ;
    else
      if (chNo < 96) 
	chNo = (95 - chNo) / 2 ;
      else
	chNo = 127 + (97 - chNo) / 2 ;
  }
  else { // NON (fDetectorType == "PRADGEM") && (fPlane.Contains("Y")) && (fAPVIndex == 11)
    if (chNo % 2 == 0) 
      chNo = ( chNo / 2) + 32 ;
    else 
      if (chNo < 64) 
	chNo = (63 - chNo) / 2 ;
      else
	chNo = 127 + (65 - chNo) / 2 ;
  }
  //   printf("PRDPedestal::PRadStripsFMapping ==>  APVID=%d, chNo=%d, stripNo=%d, \n",fAPVID, chno, chNo) ;
  return chNo ;
}

// Get Bank IDs 
set<int> PRDMapping::GetBankIDSet(){
  set<int> fec;
  fec.clear();
  map<Int_t, list<Int_t> >::const_iterator it;
  for(it=fAPVIDListFromFECIDMap.begin(); it!=fAPVIDListFromFECIDMap.end(); ++it)
    {
      fec.insert(it->first);
    }
  return fec;
}

TString PRDMapping::GetFECIPFromFECID(Int_t id)
{
  return fFECIPFromFECID[id];
}

Int_t PRDMapping::GetFECIDFromFECIP(TString ip)
{
  string str1 = (const char *)ip;

  str1.erase( remove(str1.begin(), str1.end(), ' '), str1.end() );
  map<TString, Int_t>::iterator it = fFECIDFromFECIP.begin();
  for(;it!=fFECIDFromFECIP.end();++it)
    {
      string str2 = (const char *) it->first;
      str2.erase( remove(str2.begin(), str2.end(), ' '), str2.end() );
      if( str1.compare(str2) == 0 ) 
	{
	  return it->second;
	}
    }
  return -1;
}
