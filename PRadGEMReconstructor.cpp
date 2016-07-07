//================================================//
//                                                //
// Xinzhan Bai                                    //
// 07/05/2016                                     //
//================================================//

#include <cassert>
#include <cmath>

#include <TH1F.h>

#include "PRadGEMReconstructor.h"
#include "PRDMapping.h"

#include "PRadDataHandler.h"
#include "PRadEventStruct.h"
#include "PRadGEMSystem.h"

// GEMHit
PRadDataHandler * GEMHit::fHandler = nullptr;

GEMHit::GEMHit(int hitID, int apvID, int chNo, 
               int zeroSupCut, TString isHitMaxOrTotalADCs)

    : fAPVID(apvID), fHitID(hitID), fAPVChNo(chNo),
      fSignalPeakTimeBin(0), fZeroSupCut(zeroSupCut),  
      fHitADCs(-10000), fPeakADCs(0.0), fIntegratedADCs(0.0),
      fIsHitMaxOrTotalADCs(isHitMaxOrTotalADCs), NCH(128)
{
  gem_srs = fHandler -> GetSRS();
  fTimeBinADCs.clear();

  mapping            = PRDMapping::GetInstance();
  fAPVIndexOnPlane   = mapping->GetAPVIndexOnPlane(fAPVID);
  fAPVOrientation    = mapping->GetAPVOrientation(fAPVID);
  fPlane             = mapping->GetPlaneFromAPVID(fAPVID);
  fPlaneID           = mapping->GetPlaneID(fPlane);
  fDetector          = mapping->GetDetectorFromPlane(fPlane);
  fDetectorID        = mapping->GetDetectorID(fDetector);
  fDetectorType      = mapping->GetDetectorTypeFromDetector(fDetector);
  fReadoutBoard      = mapping->GetReadoutBoardFromDetector(fDetector);
  fPlaneSize         = mapping->GetPlaneSize(fPlane);
  fNbOfAPVsOnPlane   = mapping->GetNbOfAPVsOnPlane(fPlane);

  SetStripNo();
  ComputePosition();
}

GEMHit::~GEMHit()
{
  fTimeBinADCs.clear();
  gem_srs -> Clear();
}

void GEMHit::TimingFindPeakTimeBin() 
{
  float currentMax = 0.0 ;
  map <int, float>::const_iterator  max_itr ;
  for(max_itr = fTimeBinADCs.begin(); max_itr != fTimeBinADCs.end(); ++max_itr) 
  {
    if (max_itr->second > currentMax) 
    {
      currentMax = max_itr->second ;
      fSignalPeakTimeBin = max_itr->first;
    }
  }
}

void GEMHit::AddTimeBinADCs(int timebin, float charges) 
{
  assert(timebin>=0);

  fTimeBinADCs[timebin] = charges;

  if(fZeroSupCut > 0 ) 
  {
    if (charges >  fPeakADCs) // bug 
        fPeakADCs = charges ;
    fIntegratedADCs += charges ;

    if (fIsHitMaxOrTotalADCs == "integratedADCs") 
        fHitADCs = fIntegratedADCs ;
    else                                          
        fHitADCs =  fPeakADCs;
  }
  else 
 {
    if (timebin ==0)
      fHitADCs = charges ;
    else
    {
      fHitADCs += charges;
      fHitADCs/=(timebin+1);
    }
  }
}

void GEMHit::ComputePosition() 
{
  float pitch = 0.4; //0.4mm
  if ( (fDetectorType == "PRADGEM") && (fPlane.Contains("X")) ) 
  {
      //mapping fix by xb
      fStripPosition =  -0.5 * (fPlaneSize + 12.8 /* 32 strips*0.4 pitch */-pitch ) 
                        + (pitch * fStripNo) ;
  }
  else 
  {
      fStripPosition =  -0.5 * (fPlaneSize - pitch) + (pitch * fStripNo) ;
  }
}

void GEMHit::SetStripNo() 
{

  fAbsoluteStripNo = StripMapping(fAPVChNo) ;

  if ( (fDetectorType == "PRADGEM")  && (fPlane.Contains("X")) )
    {
      int nbAPVsOnPlane = fNbOfAPVsOnPlane - 1 ;
      if(fAPVIndexOnPlane != 11)
      {
        if (fAPVIndexOnPlane > fNbOfAPVsOnPlane) fStripNo = -100000000 ;
	if(fAPVOrientation == 0) fAbsoluteStripNo = 127 - fAbsoluteStripNo ;
	//fStripNo = (fAbsoluteStripNo-16) + (NCH * (fAPVIndexOnPlane % (nbAPVsOnPlane))) ; 
	//the above line will cause chamber center to shift, do not use
	fStripNo = (fAbsoluteStripNo) + (NCH * (fAPVIndexOnPlane % (nbAPVsOnPlane))) ; 
	// take chamber center as origin, should not shift
      }
      else if (fAPVIndexOnPlane == 11)
      {
        int apvIndexOnPlane = fAPVIndexOnPlane - 1 ;
	if (apvIndexOnPlane > fNbOfAPVsOnPlane) fStripNo = -100000000 ;
	if(fAPVOrientation == 0) fAbsoluteStripNo = 127 - fAbsoluteStripNo ;
	//fStripNo = (fAbsoluteStripNo - 32 ) + (NCH * (apvIndexOnPlane % nbAPVsOnPlane)) ;
	fStripNo = (fAbsoluteStripNo - 16) + (NCH * (apvIndexOnPlane % nbAPVsOnPlane)) ;
      }
    }
  else
  {
    if (fAPVIndexOnPlane > fNbOfAPVsOnPlane) fStripNo = -100000000 ;
    if(fAPVOrientation == 0) fAbsoluteStripNo = 127 - fAbsoluteStripNo ;
    fStripNo = fAbsoluteStripNo + (NCH * (fAPVIndexOnPlane % fNbOfAPVsOnPlane)) ;
  }
}

int GEMHit::StripMapping(int chNo) 
{
  chNo = APVchannelCorrection(chNo) ;
  if (fDetectorType == "PRADGEM") chNo = PRadStripMapping(chNo) ;
  return chNo ;
}

int GEMHit::APVchannelCorrection(int chNo) 
{
  chNo = (32 * (chNo%4)) + (8 * (int)(chNo/4)) - (31 * (int)(chNo/16)) ;
  return chNo ;
}

int GEMHit::PRadStripMapping(int chNo) 
{
  if ( (fDetectorType == "PRADGEM") && (fPlane.Contains("X")) && (fAPVIndexOnPlane == 11) ) 
  {
    if (chNo % 2 == 0) chNo = ( chNo / 2) + 48 ; // originally was 48 xb: important
    else
      if (chNo < 96) chNo = (95 - chNo) / 2 ;
      else           chNo = 127 + (97 - chNo) / 2 ;
  }
  else 
  {
    if (chNo % 2 == 0) chNo = ( chNo / 2) + 32 ;
    else
      if (chNo < 64) chNo = (63 - chNo) / 2 ;
      else           chNo = 127 + (65 - chNo) / 2 ;
  }
  return chNo ;
}

// GEMCluster
PRadDataHandler * GEMCluster::fHandler = nullptr;

GEMCluster::GEMCluster(int minClusterSize, int maxClusterSize, TString isMaxorTotalADCs)

    : fNbOfHits(0), fClusterPeakTimeBin(0), fClusterTimeBin(0), 
      fClusterPeakADCs(0), fClusterTimeBinADC(0), 
      fClusterSumADCs(0), fposition(0), 
      fclusterCentralStrip(0), fstrip(0), fPlaneSize(512.), fNbAPVsOnPlane(10), 
      fMinClusterSize(minClusterSize), fMaxClusterSize(maxClusterSize),
      fIsClusterMaxOrSumADCs(isMaxorTotalADCs), fPlane("GEM1X"),
      fIsGoodCluster(true)
{
  gem_srs = fHandler -> GetSRS();

  fArrayOfHits = new TObjArray(maxClusterSize);
}

GEMCluster::~GEMCluster() 
{
  fArrayOfHits->Clear();
  delete fArrayOfHits;
  fClusterTimeBinADCs.clear() ;
  gem_srs -> Clear();
}

void GEMCluster::Timing() 
{
  float q = 0.0;
  int nBins = fClusterTimeBinADCs.size() ;
  for (int k = 0; k < nBins; k++) 
  {
    if (fClusterTimeBinADCs[k] > q) 
    {
      q = fClusterTimeBinADCs[k];
      fClusterTimeBin = k;
    }
  }
}

int GEMCluster::GetClusterTimeBin() 
{
  TObjArray &temp = *fArrayOfHits;

  int nbofhits =  GetNbOfHits() ;

  for (int i = 0; i < nbofhits; i++) 
  {
    map <int, float>  timeBinADCs = ((GEMHit*)temp[i])->GetTimeBinADCs() ;
    int nbOfTimeBins = timeBinADCs.size() ;
    fClusterTimeBinADCs.resize(nbOfTimeBins) ;
    for (int k = 0; k < nbOfTimeBins; k++)  
        fClusterTimeBinADCs[k] += timeBinADCs[k] ;
    timeBinADCs.clear() ;
  }
  Timing() ;
  return fClusterTimeBin ;
}

void GEMCluster::ClusterPositionPulseHeghtWeight() 
{  
  // Calculate the fposition and the total fClusterSumADCs
  float hitposition, q;
  TObjArray &temp = *fArrayOfHits;
  int nbofhits =  GetNbOfHits() ;
  for (int i = 0; i < nbofhits; i++) 
  {
    q  = ((GEMHit*)temp[i])->GetHitADCs() ;
    hitposition = ((GEMHit*)temp[i])->GetStripPosition() ;
    fClusterSumADCs += q ;
    fposition += q * hitposition ;
    if (q > fClusterPeakADCs) 
    {
      fClusterPeakTimeBin = ((GEMHit*)temp[i])->GetSignalPeakTimeBin() ;
      fClusterPeakADCs = q ;
    }
  }
  fposition /= fClusterSumADCs;
}

void GEMCluster::ClusterCentralStrip() 
{
  float p, dp ;
  float dpmin = 99;
  TObjArray &temp = *fArrayOfHits;
  int nbofhits =  GetNbOfHits() ;
  for (int i = 0; i < nbofhits; i++) 
  {
    p  = ((GEMHit*)temp[i])->GetStripPosition();
    dp = fabs(fposition - p);
    if (dp <= dpmin) 
    {
      fclusterCentralStrip = p;
      dpmin = dp;
    }
  }
}

float GEMCluster::GetClusterADCs() 
{
  float adcs = 0 ;
  if (fIsClusterMaxOrSumADCs == "maximumADCs") 
  {
    adcs = fClusterPeakADCs ;
  }
  else 
  {
    adcs = fClusterSumADCs ;
  }
  return adcs ;
}

void GEMCluster::AddHit(GEMHit *h) 
{
  fArrayOfHits->AddLast(h);
}

void GEMCluster::ClearArrayOfHits() 
{
  fArrayOfHits->Clear();
}

bool GEMCluster::IsGoodCluster() 
{
  fIsGoodCluster = true ;
  fNbOfHits = fArrayOfHits->GetEntries() ;
  if ( (fNbOfHits > fMaxClusterSize) || (fNbOfHits < fMinClusterSize) ) 
  {
    ClearArrayOfHits() ;
    fIsGoodCluster = false ;
    fNbOfHits = fArrayOfHits->GetEntries() ;
  }
  return fIsGoodCluster ;
}

int GEMCluster::Compare(const TObject *obj) const 
{
  int compare = (fClusterSumADCs < ((GEMCluster*)obj)->GetClusterADCs()) ? 1 : -1;
  return compare ;
}

void GEMCluster::ComputeClusterPosition() 
{
  ClusterPositionPulseHeghtWeight() ;
  ClusterCentralStrip() ;
}

// GEMZeroHitDecoder
PRadDataHandler * GEMZeroHitDecoder::fHandler = nullptr;

GEMZeroHitDecoder::GEMZeroHitDecoder(vector<GEM_Data> * gemdata)

    : fIsGoodClusterEvent(false), fMinClusterSize(1), fMaxClusterSize(20),
      NCH(128), fZeroSupCut(5), fFECID(0), fADCChannel(0), fAPVID(0),
      fAPVKey(0), fIsHitMaxOrTotalADCs("signalPeak"), 
      fIsClusterMaxOrTotalADCs("totalADCs"), gem_data(gemdata),
      nTimeBin(3), Zgem1(5300.0), Zgem2(5260.)
{
  gem_srs = fHandler -> GetSRS();

  fMapping = PRDMapping::GetInstance();
  fListOfHitsZero.clear();
  fListOfHitsZeroFromPlane.clear();
  fListOfClustersZeroFromPlane.clear();

  ProcessEvent();
}

GEMZeroHitDecoder::~GEMZeroHitDecoder()
{
  Clear();
  gem_srs -> Clear();
}

void GEMZeroHitDecoder::Clear()
{

  if(fListOfHitsZero.size() > 0)
  {
    for(auto &it : fListOfHitsZero)
    {
      delete it.second;
    }
  }
  fListOfHitsZero.clear();

  if( fListOfClustersZeroFromPlane.size() > 0)
  {
    for( auto & itt : fListOfClustersZeroFromPlane)
    {
      if( itt.second.size() <=0 ) 
          continue;
      for(auto & itc : itt.second)
      {
        delete itc;
      }
      itt.second.clear();
    }
  }
  fListOfClustersZeroFromPlane.clear();

  if( fListOfHitsZeroFromPlane.size() > 0)
  {
    for(auto &i : fListOfHitsZeroFromPlane)
    {
      if(i.second.size() <= 0) 
          continue;
      i.second.clear();
    }
  }
  fListOfHitsZeroFromPlane.clear(); 
}

void GEMZeroHitDecoder::ProcessEvent()
{
  EventHandler();
}

void GEMZeroHitDecoder::EventHandler()
{
  if(gem_data->size() == 0)
     return;
  int fSize = gem_data->size();
  int chNo = 0;
  int adc = 0;
  int time_sample = 0;
  for(int i = 0;i<fSize;++i)
  {
    fFECID = gem_data->at(i).addr.fec;
    fADCChannel = gem_data->at(i).addr.adc;
    chNo = gem_data->at(i).addr.strip;
    for(int ts = 0;ts<nTimeBin;++ts)
    {
      time_sample = ts;
      adc = gem_data->at(i).values[ts];

      fAPVID = (fFECID<<4)|fADCChannel;
      fAPVKey = fMapping->GetAPVNoFromID(fAPVID);
      int hitID = (fAPVKey << 8) | chNo ;
      /*
      cout<<"fec: "<<fFECID
          <<" adc: "<<fADCChannel
          <<" strip: "<<chNo
	  <<" timebin: "<<time_sample
	  <<" adc: "<<adc<<endl;
      */
      if( !fListOfHitsZero[hitID])
      {
        GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut, fIsHitMaxOrTotalADCs);
	fListOfHitsZero[hitID] = hit;
      }
      fListOfHitsZero[hitID] -> AddTimeBinADCs(time_sample, adc);
    }
  }

  GetListOfHitsZeroFromPlanes();
  ComputeClusters();
}

void GEMZeroHitDecoder::GetListOfHitsZeroFromPlanes() 
{
  map < Int_t, GEMHit * >::const_iterator hit_itr ;
  for(hit_itr = fListOfHitsZero.begin(); hit_itr != fListOfHitsZero.end(); ++hit_itr) 
  {
    GEMHit * hit = (* hit_itr).second ;
    TString planename = hit->GetPlane() ;
    fListOfHitsZeroFromPlane[planename].push_back(hit) ;
  }
}

TH1F* GEMZeroHitDecoder::GetZeroHit(TString str)
{
  TH1F * h1;
  int nbDetector = fMapping->GetNbOfDetectors();
  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);

    list<TString>::iterator it;
    for(it=planeList.begin();it!=planeList.end();++it)
    {
      if(*it == str)
      {
        TString hh = detectorName+"_"+(*it)+"_hit_distribution_zero_suppression";
        h1 = new TH1F(hh, hh, 2000, -fMapping->GetPlaneSize(*it)/2-100, fMapping->GetPlaneSize(*it)/2+100 );
        list< GEMHit* > hitList = fListOfHitsZeroFromPlane[ *it  ];
        list< GEMHit* >::iterator hit_it;
        for(hit_it=hitList.begin(); hit_it!=hitList.end();++hit_it)
        {
          Float_t pos = (*hit_it) -> GetStripPosition();
          Float_t adc = (*hit_it) -> GetHitADCs();
          h1 -> Fill(pos, adc);
        }
      }
    }
  }

  if(h1)
      return h1;
  else 
      return nullptr;
}

//-------------------------
//GEMZeroHitDecoder Cluster
//-------------------------
static bool CompareStripNo( TObject *obj1, TObject *obj2) 
{
  bool compare ;
  if ( ( (GEMHit*) obj1 )->GetStripNo() < ( ( GEMHit*) obj2 )->GetStripNo() ) 
      compare = true ;
  else 
      compare = false ;
  return compare ;
}

static bool CompareHitADCs( TObject *obj1, TObject *obj2) 
{
  bool compare ;
  if ( ( (GEMHit*) obj1 )->GetHitADCs() > ( ( GEMHit*) obj2 )->GetHitADCs()) 
      compare = true ;
  else 
      compare = false ;
  return compare ;
}

static bool CompareClusterADCs( TObject *obj1, TObject *obj2) 
{
  bool compare ;
  if ( ( (GEMCluster*) obj1 )->GetClusterADCs() > ( ( GEMCluster*) obj2 )->GetClusterADCs()) 
      compare = true ;
  else 
      compare = false ;
  return compare ;
}

static bool CompareClusterSize( TObject *obj1, TObject *obj2) 
{
  bool compare ;
  if ( ( (GEMCluster*) obj1 )->GetNbOfHits() > ( ( GEMCluster*) obj2 )->GetNbOfHits()) 
      compare = true ;
  else 
      compare = false ;
  return compare ;
}

void GEMZeroHitDecoder::ComputeClusters()
{
  map < TString, list <GEMHit*> >::const_iterator  hitsFromPlane_itr ;

  for (hitsFromPlane_itr = fListOfHitsZeroFromPlane.begin(); hitsFromPlane_itr != fListOfHitsZeroFromPlane.end(); ++hitsFromPlane_itr)
    {
      TString plane =  (*hitsFromPlane_itr).first ;
      list <GEMHit*> hitsFromPlane = (*hitsFromPlane_itr).second ;
      hitsFromPlane.sort(CompareStripNo) ;
      Int_t listSize = hitsFromPlane.size() ;

      if (listSize < fMinClusterSize)
        {
          fIsGoodClusterEvent = false ;
          continue ;
        }

      Int_t previousStrip = -2 ;
      Int_t clusterNo = -1 ;
      map<Int_t, GEMCluster *> clustersMap ;
      list <GEMHit *>::const_iterator hit_itr ;

      for (hit_itr = hitsFromPlane.begin(); hit_itr != hitsFromPlane.end(); hit_itr++)
        {
          GEMHit * hit =  * hit_itr ;

          Int_t currentStrip = hit->GetStripNo() ;
          if( plane.Contains("X") && ( (currentStrip<16) || (currentStrip > 1391) )  ) 
	      continue;

          if(currentStrip != (previousStrip + 1))
            {
              clusterNo++ ;
            }
          if(!clustersMap[clusterNo])
            {
              clustersMap[clusterNo] = new GEMCluster(fMinClusterSize, fMaxClusterSize, fIsClusterMaxOrTotalADCs) ;
              clustersMap[clusterNo]->SetNbAPVsFromPlane(hit->GetNbAPVsFromPlane());
              clustersMap[clusterNo]->SetAPVIndexOnPlane(hit->GetAPVIndexOnPlane());
              clustersMap[clusterNo]->SetPlaneSize(hit->GetPlaneSize());
              clustersMap[clusterNo]->SetPlane(hit->GetPlane());
            }
          clustersMap[clusterNo]->AddHit(hit) ;
          previousStrip = currentStrip;
        }

      map<Int_t, GEMCluster *>::const_iterator  cluster_itr ;
      for (cluster_itr = clustersMap.begin(); cluster_itr != clustersMap.end(); cluster_itr++)
        {
          GEMCluster * cluster = ( * cluster_itr ).second ;
          if (!cluster->IsGoodCluster())
            {
              delete cluster ;
              continue ;
            }
          cluster->ComputeClusterPosition() ;
          fListOfClustersZeroFromPlane[plane].push_back(cluster) ;
        }

      fListOfClustersZeroFromPlane[plane].sort(CompareClusterADCs) ;
      hitsFromPlane.clear() ;
      clustersMap.clear() ;
    }
}

TH1F* GEMZeroHitDecoder::GetCluster(TString str)
{
  TH1F * hc1;

  int nbDetector = fMapping->GetNbOfDetectors();

  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);

    list<TString>::iterator it;
    for(it=planeList.begin();it!=planeList.end();++it)
    {
      if(*it == str)
      {
        TString hh = detectorName+"_"+(*it)+"_Cluster_Distribution_zero_suppression";
        hc1 = new TH1F(hh, hh, 2000, -fMapping->GetPlaneSize(*it)/2-100, fMapping->GetPlaneSize(*it)/2+100 );
        list< GEMCluster* > clusterList = fListOfClustersZeroFromPlane[ *it  ];
        list< GEMCluster* >::iterator cluster_it;
        for(cluster_it=clusterList.begin(); cluster_it!=clusterList.end();++cluster_it)
        {
          Float_t pos = (*cluster_it) -> GetClusterPosition();
          Float_t adc = (*cluster_it) -> GetClusterADCs();
          hc1 -> Fill(pos, adc);
        }
      }
    }
  }
  if(hc1)
      return hc1;
  else 
      return nullptr;
}

void GEMZeroHitDecoder::GetClusterHyCal(vector<PRadGEMCluster> &gem1, 
                                        vector<PRadGEMCluster> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  /*
   * Convert GEM coordinate to HyCal coordinate
   *
   *    Move GEM coordinate to HyCal Hole center
   *    overlapping area: hole diameter (44mm)
   *    Chamber X side length: 550.4mm
   *
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must be pointing downward
   *
   *    beam downstream is z axis direction
   */

  float O_Transfer = 253.2;
  //float OverlapLength = 44;
  float z_gem1 = Zgem1; //mm
  float z_gem2 = Zgem2; //mm

  // cutting edge
  float edge1 = 0;
  float edge2 = 0;

  // parameters got from production data
  float xoffset = -0.3722;
  float yoffset = 0.1681;
  
  //the above offsets are projected on GEM2 z-plane
  //real GEM1 offsets should be projected back
  xoffset = xoffset*z_gem1/z_gem2;
  yoffset = yoffset*z_gem1/z_gem2;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      /*
       * overlapping area: according to frame design, overlapping area is in fact 44mm
       */
      if(((*itx)->GetClusterPosition() -xoffset -O_Transfer) <= edge1) 
      // remove overlapping area on GEM1, use the corresponding area on GEM2
      // do not use chamber edge as the cut line. 
      // choose some arbitrary line (eg.=0) and find the corresponding line on GEM2
      {
        float c_x = (*itx)->GetClusterADCs();
	float c_y = (*ity)->GetClusterADCs();
        float x = (*itx++)->GetClusterPosition() -O_Transfer - xoffset;
	float y = (*ity++)->GetClusterPosition() -yoffset; 
        gem1.push_back( PRadGEMCluster(x, y, z_gem1, c_x, c_y) ) ;
      }
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      if( ( O_Transfer - ((*itx2)->GetClusterPosition())  ) < edge2) 
          continue;
      {
        float c_x = (*itx2)->GetClusterADCs();
	float c_y = (*ity2)->GetClusterADCs();
        float x =  O_Transfer - ((*itx2++)->GetClusterPosition());
	float y =  -(*ity2++)->GetClusterPosition() ; 
	gem2.push_back(PRadGEMCluster(x, y, z_gem2, c_x, c_y));
      }
    }
  }

}

void GEMZeroHitDecoder::GetClusterBeamLine(vector<PRadGEMCluster> &gem1, 
                                           vector<PRadGEMCluster> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  /*
   * Convert GEM coordinate to BeamLine coordinate
   *
   *    Move GEM coordinate to Beam center
   *    Chamber X plane size: 550.4 mm
   *    overlapping area: hole diameter (44mm)
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must be pointing downward
   *
   *    beam downstream is z axis direction
   */

  float O_Transfer = 253.2;
  //float OverlapLength = 44;
 
  // cutting edge
  float edge1 = 0;
  float edge2 = 0;
  float z_gem1 = Zgem1; //mm
  float z_gem2 = Zgem2; //mm

  float xoffset_gem = -0.3722;
  float yoffset_gem = 0.1681;

  //the above offsets are projected on GEM2 z-plane
  //real GEM1 offsets should be projected back
  xoffset_gem = xoffset_gem*z_gem1/z_gem2;
  yoffset_gem = yoffset_gem*z_gem1/z_gem2;

  // got from production data
  float xoffset_beam = 1.804;
  float yoffset_beam = -0.1568;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      /*
       * overlapping area: according to frame design, overlapping area is in fact 44mm
       */
      if(((*itx)->GetClusterPosition() -xoffset_gem -O_Transfer) <= edge1) 
      // remove overlapping area on GEM1, use the corresponding area on GEM2
      // do not use chamber edge as the cut line. 
      // choose some arbitrary line (eg. = 0) and find the corresponding line on GEM2
      {
        float c_x = (*itx)->GetClusterADCs();
	float c_y = (*ity)->GetClusterADCs();
        float x = (*itx++)->GetClusterPosition() -O_Transfer - xoffset_gem - xoffset_beam ;
	float y =  (*ity++)->GetClusterPosition() -yoffset_gem - yoffset_beam;
	gem1.push_back(PRadGEMCluster(x, y, z_gem1, c_x, c_y));
      }
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      if( ( O_Transfer - ((*itx2)->GetClusterPosition())  ) < edge2) 
          continue;
      {
        float c_x = (*itx2)->GetClusterADCs();
	float c_y = (*ity2)->GetClusterADCs();
        float x =  O_Transfer - ((*itx2++)->GetClusterPosition()) - xoffset_beam;
	float y = -(*ity2++)->GetClusterPosition()  - yoffset_beam;
	gem2.push_back(PRadGEMCluster(x, y, z_gem2, c_x, c_y) );
      }
    }
  }

}

void GEMZeroHitDecoder::GetClusterGEM(vector<PRadGEMCluster> &gem1, 
                                      vector<PRadGEMCluster> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersZeroFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersZeroFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersZeroFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersZeroFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  float z_gem1 = Zgem1;//mm
  float z_gem2 = Zgem2;//mm

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      float c_x = (*itx)->GetClusterADCs();
      float c_y = (*ity)->GetClusterADCs();
      float x = (*itx++)->GetClusterPosition();
      float y = (*ity++)->GetClusterPosition();
      gem1.push_back(PRadGEMCluster(x, y, z_gem1, c_x, c_y));
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      float c_x = (*itx2)->GetClusterADCs();
      float c_y = (*ity2)->GetClusterADCs();
      float x =  (*itx2++)->GetClusterPosition(); 
      float y =  (*ity2++)->GetClusterPosition(); 
      gem2.push_back( PRadGEMCluster(x, y, z_gem2, c_x, c_y) );
    }
  }

}

// PRadGEMReconstructor
PRadDataHandler * PRadGEMReconstructor::fHandler = nullptr;

void PRadGEMReconstructor::gSetHandler( PRadDataHandler * handler)
{
  GEMHit::fHandler = handler;
  GEMCluster::fHandler = handler;
  GEMZeroHitDecoder::fHandler = handler;
  PRadGEMReconstructor::fHandler = handler;
}

PRadGEMReconstructor::PRadGEMReconstructor( PRadDataHandler * handler)
    : fPRadGEMSystem(nullptr), event(nullptr),
      pDecode(nullptr), Zgem1(5300.0), Zgem2(5260), 
      Zhycal(5820.0), Match_Criteria(40.)
{
  gSetHandler(handler);
  fPRadGEMSystem = fHandler -> GetSRS();

  fMapping = PRDMapping::GetInstance();

  fPRadGEMCluster.clear();
  gem1_local.clear();
  gem2_local.clear();
  gem1_hycal.clear();
  gem2_hycal.clear();
  gem1_beaml.clear();
  gem2_beaml.clear();

}

PRadGEMReconstructor::~PRadGEMReconstructor()
{
  Clear();
  fPRadGEMSystem->Clear();
}

void PRadGEMReconstructor::Clear() 
{
  fPRadGEMCluster.clear();
  gem1_local.clear();
  gem2_local.clear();
  gem1_hycal.clear();
  gem2_hycal.clear();
  gem1_beaml.clear();
  gem2_beaml.clear();

  //fPRadGEMSystem->Clear();
}

vector<PRadGEMCluster> &PRadGEMReconstructor::CoarseGEMReconstruct(const int &event_index) 
{
  if(fHandler) 
  {
    EventData &ev = fHandler->GetEvent(event_index);
    return CoarseGEMReconstruct(ev);
  }
  else 
    return CoarseGEMReconstruct(*event);
}

vector<PRadGEMCluster> &PRadGEMReconstructor::CoarseGEMReconstruct(EventData &ev) 
{
  if(ev.gem_data.size() != 0) 
  {
    event = &ev;
    return Reconstruct();
  }
  else 
  {
    Clear();
    return fPRadGEMCluster;
  }
}

vector<PRadGEMCluster> &PRadGEMReconstructor::Reconstruct()
{
  Clear();

  vector<GEM_Data> gem_data = event->gem_data;
  pDecode = new GEMZeroHitDecoder(&gem_data);
  PackClusters();
  pDecode -> Clear();
  event->clear();
  return fPRadGEMCluster;
}

void PRadGEMReconstructor::PackClusters()
{
  // pack data
  pDecode->GetClusterHyCal(gem1_hycal, gem2_hycal);
  pDecode->GetClusterGEM(gem1_local, gem2_local);
  pDecode->GetClusterBeamLine(gem1_beaml, gem2_beaml);

  HyCalGEMPosMatch();
}

vector<PRadGEMCluster> &PRadGEMReconstructor::GEMClusteringLocal(const int det)
{
  if(det == 0)
      return gem1_local;
  else if( det == 1)
      return gem2_local;
  else
  {
      cout<<"PRadGEMReconstructor::GEMClusteringLocal ERROR..."
          <<endl;
      exit(-1);
  }
}

vector<PRadGEMCluster> &PRadGEMReconstructor::GEMClusteringHyCal(const int det)
{
  if(det == 0)
      return gem1_hycal;
  else if( det == 1)
      return gem2_hycal;
  else
  {
      cout<<"PRadGEMReconstructor::GEMClusteringHyCal ERROR..."
          <<endl;
      exit(-1);
  }
}

vector<PRadGEMCluster> &PRadGEMReconstructor::GEMClusteringBeamL(const int det)
{
  if(det == 0)
      return gem1_beaml;
  else if( det == 1)
      return gem2_beaml;
  else
  {
      cout<<"PRadGEMReconstructor::GEMClusteringBeamL ERROR..."
          <<endl;
      exit(-1);
  }
}

void PRadGEMReconstructor::HyCalGEMPosMatch()
{
  // do the match
  auto hycal_hit = fHandler->GetHyCalCluster(*event);

  Match(gem1_beaml, gem2_beaml, &hycal_hit);

  // after match
  if(gem1_beaml.size() > 0)
  {
    for(auto &i : gem1_beaml)
    {
      fPRadGEMCluster.push_back(i);
    }
  }
  if(gem2_beaml.size() > 0)
  {
    for(auto &j : gem2_beaml)
    {
      fPRadGEMCluster.push_back(j);
    }
  }
}

int PRadGEMReconstructor::Match(vector<PRadGEMCluster> &gem1, 
                                 vector<PRadGEMCluster> &gem2, 
                                 vector<HyCalHit> *pHHit)
{

  if( (gem1.size() == 0) && (gem2.size() == 0))
    {

      pHHit->clear();
      return 0;
    }
  if(pHHit->size() == 0)
    {
      gem1.clear();
      gem2.clear();

      return 0;
    }

  double z_gem1 = Zgem1; 
  double z_gem2 = Zgem2;
  double z_hycal = Zhycal; 
  double res = Match_Criteria; 
  // a larger range, 60mm

  vector<PRadGEMCluster> res_gem1;
  vector<PRadGEMCluster> res_gem2;

  int nh = pHHit->size();
  for(int i=0;i<nh;i++)
    {

      int match_gem1 = 0;
      int match_gem2 = 0;
      float m_index = 0;
      float m_e = 0;

      // search GEM1
      int n = gem1.size();
      if(n > 0)
        {
          double x_hycal = (pHHit->at(i).x) *z_gem1/z_hycal;
          double y_hycal = (pHHit->at(i).y) *z_gem1/z_hycal;

          for(int j=0;j<n;j++)
            {
              double delta = dd(gem1[j].x, gem1[j].y, x_hycal, y_hycal); 
              if( delta < res )
                {
                  res = delta;
                  m_index = j;
                  m_e = pHHit->at(i).E;
                  match_gem1 = 1;
                  match_gem2 = 0;
                }
            }
        }

      //continue search GEM2
      n = gem2.size();
      if(n>0)
        {
          double x_hycal = (pHHit->at(i).x) *z_gem2/z_hycal;
          double y_hycal = (pHHit->at(i).y) *z_gem2/z_hycal;

          for(int j=0;j<n;j++)
            {
              double delta = dd(gem2[j].x, gem2[j].y, x_hycal, y_hycal); 
	      if( delta < res )
                {
                  res = delta;
                  m_index = j;
                  m_e = pHHit->at(i).E;
                  match_gem2 = 1;
                  match_gem1 = 0;
                }
            }
        }

      //after searching
      if( (match_gem1 == 1) && (match_gem2 == 1) ) 
           cout<<"GEMPhysHandler::HyCalGEMPosMatch(): error...."
               <<endl;
      if( match_gem1 == 1 ) 
      { 
          gem1[m_index].energy = m_e; 
	  res_gem1.push_back(gem1[m_index]);
      }
      if( match_gem2 == 1 ) 
      { 
          gem2[m_index].energy = m_e; 
	  res_gem2.push_back(gem2[m_index]);
      }
    }

  gem1.clear();
  gem2.clear();

  if(res_gem1.size() > 0)
    {
      gem1 = res_gem1;
      res_gem1.clear();
    }
  if(res_gem2.size() > 0)
    {
      gem2 = res_gem2;
      res_gem2.clear();
    }

  return gem1.size() + gem2.size();
}

float PRadGEMReconstructor::dd(const float &x1, const float &y1, 
                               const float &x2, const float &y2)
{
  return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ); 
}
