#ifndef _PRADGEMRECONSTRUCTOR_H
#define _PRADGEMRECONSTRUCTOR_H

#include <iostream>
#include <map>
#include <list>

#include <TObject.h>
#include <TObjArray.h>
#include <TString.h>

#include "PRadEventStruct.h"

using namespace std;

class PRDMapping;
class PRadGEMSystem;
class PRadGEMAPV;
class PRadEventStruct;
class PRadDataHandler;
class PRadDSTParser;

class GEMHit : public TObject
{
public:
  GEMHit();
  GEMHit(int hitID, int apvID, 
         int chNo, int zeroSupCut, TString isHitMaxOrTotalADCs);
  ~GEMHit();

  void TimingFindPeakTimeBin();
  int Compare(const TObject *obj) const {
       return (fStripNo > ((GEMHit*)obj)->GetStripNo()) ? 1: -1; 
  }
  void ComputePosition();

  void AddTimeBinADCs(int timebin, float charges);
  void ClearTimeBinADCs() {
      fTimeBinADCs.clear();
  }

  int GetAPVID() { 
      return fAPVID;
  }
  int GetAPVOrientation() {
      return fAPVOrientation;
  }
  int GetAPVIndexOnPlane() {
      return fAPVIndexOnPlane;
  }
  int GetNbAPVsFromPlane() {
      return fNbOfAPVsOnPlane;
  }
  float GetHitADCs() {
      return fHitADCs;
  }
  map<int, float> & GetTimeBinADCs() {
      return fTimeBinADCs;
  }

  int StripMapping(int chNo);
  int APVchannelCorrection(int chNo);
  int PRadStripMapping(int chNo);

  int GetSignalPeakTimeBin() {
      TimingFindPeakTimeBin();
      return fSignalPeakTimeBin;
  }

  TString GetPlane() {
      return fPlane;
  }
  float GetPlaneSize() {
      return fPlaneSize;
  }
  TString GetDetector() {
      return fDetector;
  }
  TString GetDetectorType() {
      return fDetectorType;
  }
  TString GetReadoutBoard() {
      return fReadoutBoard;
  }
  TString GetHitMaxOrTotalADCs() {
      return fIsHitMaxOrTotalADCs;
  }
  
  void SetStripNo();
  int GetStripNo() {
      return fStripNo;
  }
  int GetAbsoluteStripNo() {
      return fAbsoluteStripNo;
  }

  float GetStripPosition() {
      ComputePosition();
      return fStripPosition;
  }

private:
  map<int, float> fTimeBinADCs;
  int fDetectorID, fPlaneID, fAPVID, fHitID, fAPVChNo, fStripNo, fAbsoluteStripNo;
  int fAPVIndexOnPlane, fNbOfAPVsOnPlane, fAPVOrientation;
  int fSignalPeakTimeBin;
  int fZeroSupCut;

  float fHitADCs, fPeakADCs, fIntegratedADCs, fStripPosition, fPlaneSize,  fHitPedestalNoise;
  TString fPlane, fReadoutBoard, fDetectorType, fDetector, fIsHitMaxOrTotalADCs;

  int NCH;

  PRDMapping * mapping;

};

class GEMCluster : public TObject
{
public:
  GEMCluster(int minClusterSize, int maxClusterSize, TString isMaximumOrTotalCharges);
  ~GEMCluster();

  TObjArray* GetArrayOfHits(){
      return fArrayOfHits;
  }

  GEMHit *GetHit(int i) {
      TObjArray &temp = *fArrayOfHits;
      return (GEMHit*)temp[i];
  }

  void SetMinClusterSize(int min){
      fMinClusterSize = min;
  }
  void SetMaxClusterSize(int max){
      fMaxClusterSize = max;
  }

  void AddHit( GEMHit* h);
  int Compare(const TObject *obj) const;

  int &GetNbOfHits(){
      return fNbOfHits;
  }
  TString GetPlane(){
      return fPlane;
  }
  void SetPlane(TString planename) {
      fPlane = planename;
  }
  int GetNbAPVsFromPlane(){
      return fNbAPVsOnPlane;
  }
  void SetNbAPVsFromPlane(int nb) {
      fNbAPVsOnPlane = nb;
  }
  int GetAPVIndexOnPlane() {
      return fapvIndexOnPlane;
  }
  void SetAPVIndexOnPlane(int nb) {
      fapvIndexOnPlane = nb;
  }
  float GetPlaneSize() {
      return fPlaneSize;
  }
  void SetPlaneSize(float planesize) {
      fPlaneSize = planesize;
  }

  float &GetClusterPosition() {
      return fposition;
  }
  float &GetClusterCentralStrip(){
      return fclusterCentralStrip;
  }

  void Timing();
  int GetClusterTimeBin();
  int GetClusterPeakTimeBin() {
      return fClusterPeakTimeBin;
  }
  float GetClusterADCs();
  void SetClusterADCs(float adc){
      fClusterSumADCs = adc;
  }

  void ClearArrayOfHits();
  bool IsGoodCluster();
  void ClusterCentralStrip();
  void ClusterPositionPulseHeghtWeight();
  vector<float> GetClusterTimeBinADCs(){
      return fClusterTimeBinADCs;
  }
  void ComputeClusterPosition();

private:
  int fNbOfHits;
  TObjArray *fArrayOfHits;
  int fClusterPeakTimeBin, fClusterTimeBin;
  float fClusterPeakADCs, fClusterTimeBinADC, 
        fClusterSumADCs, fposition, 
	fclusterCentralStrip, fstrip, 
	fPlaneSize;
  int fapvID, fStripNo, fAbsoluteStripNo, 
      fapvIndexOnPlane, fNbAPVsOnPlane, 
      fMinClusterSize, fMaxClusterSize;
  TString fIsClusterMaxOrSumADCs, fPlane;
  bool fIsGoodCluster;
  vector<float> fClusterTimeBinADCs;
};

struct PRadGEMCluster
{
    float x;
    float y;
    float z;
    float x_charge;
    float y_charge;
    float energy;

    PRadGEMCluster(float xi, float yi, float zi,
               float cix = 0., float ciy = 0.,
               float ei = 0)
    : x(xi), y(yi), z(zi), x_charge(cix), y_charge(ciy), energy(ei) {}

    void SetEnergy(float e) {energy = e;}
    void SetX(float xp) {x = xp;}
    void SetY(float xp) {y = xp;}
    void SetZ(float zp) {z = zp;}
    void SetXCharge(float xp) {x_charge = xp;}
    void SetYCharge(float xp) {y_charge = xp;}

    PRadGEMCluster  operator = (const PRadGEMCluster &c)
    {
      return PRadGEMCluster(c.x, c.y, c.z, c.x_charge, c.y_charge, c.energy);
    }

};

class GEMZeroHitDecoder
{
public:
  GEMZeroHitDecoder(vector<GEM_Data> * gemdata);
  ~GEMZeroHitDecoder();

  void ProcessEvent();
  void EventHandler();
  
  void GetListOfHitsZeroFromPlanes();
  TH1F* GetZeroHit(TString plane);

  void ComputeClusters();
  map<TString, list<GEMCluster*> > GetListOfClustersFromPlanes() {
      return fListOfClustersZeroFromPlane;
  }
  TH1F* GetCluster(TString str);
  void GetClusterGEM(vector<PRadGEMCluster> &gem1, 
                     vector<PRadGEMCluster> &gem2);
  void GetClusterHyCal(vector<PRadGEMCluster> &gem1, 
                       vector<PRadGEMCluster> &gem2);
  void GetClusterBeamLine(vector<PRadGEMCluster> &gem1, 
                          vector<PRadGEMCluster> &gem2);
  void Clear();
private:
  bool fIsGoodClusterEvent;
  int fMinClusterSize, fMaxClusterSize;

  map<TString, list<GEMCluster*> > fListOfClustersZeroFromPlane;

  int NCH;

  int fZeroSupCut;
  int fFECID;
  int fADCChannel;
  int fAPVID;
  int fAPVKey;
  TString fIsHitMaxOrTotalADCs;
  TString fIsClusterMaxOrTotalADCs;
  PRDMapping * fMapping;
  map<int, GEMHit*> fListOfHitsZero;
  map<TString, list<GEMHit*> > fListOfHitsZeroFromPlane;

  vector<GEM_Data> *gem_data;
  int nTimeBin; 
  double Zgem1;
  double Zgem2;
};

class PRadGEMReconstructor
{
public:
  PRadGEMReconstructor(PRadDataHandler *h = nullptr);
  virtual ~PRadGEMReconstructor();
  //void InitConfig(const string &path); // to be implemented...
  void Clear();
  void SetHandler(PRadDataHandler* h) {
      fHandler = h; // for GEM HyCal Match
  }

  vector<PRadGEMCluster> &CoarseGEMReconstruct(const int &event_index = -1);
  vector<PRadGEMCluster> &CoarseGEMReconstruct(EventData &event);
  vector<PRadGEMCluster> &Reconstruct();
  void PackClusters();

  void HyCalGEMPosMatch();
  int Match(vector<PRadGEMCluster> &gem1, 
             vector<PRadGEMCluster> &gem2, 
             vector<HyCalHit> *pHHit);
  float dd(const float &x1, const float &y1, const float &x2, const float &y2);

  // GEM charactorization
  vector<PRadGEMCluster> &GEMClusteringLocal(const int det);
  vector<PRadGEMCluster> &GEMClusteringHyCal(const int det);
  vector<PRadGEMCluster> &GEMClusteringBeamL(const int det);

private:
  PRDMapping *fMapping;
  PRadDataHandler *fHandler;
  PRadGEMSystem *fPRadGEMSystem;
  EventData *event;
  GEMZeroHitDecoder *pDecode;

  double Zgem1;
  double Zgem2;
  double Zhycal;
  double Match_Criteria;

  vector<PRadGEMCluster> fPRadGEMCluster;
  vector<PRadGEMCluster> gem1_local;
  vector<PRadGEMCluster> gem2_local;
  vector<PRadGEMCluster> gem1_hycal;
  vector<PRadGEMCluster> gem2_hycal;
  vector<PRadGEMCluster> gem1_beaml;
  vector<PRadGEMCluster> gem2_beaml;

};

#endif
