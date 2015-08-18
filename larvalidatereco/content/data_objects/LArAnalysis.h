
#ifndef __LAR_ANALYSIS__
#define __LAR_ANALYSIS__

///
/// The LArAnalysis object is a container class designed to be written
/// to a ROOT TTree.
///
/// Passing this object by address to TTree::Branch() allows the
/// addressing of all the member data in a single function call
/// (rather than laboriously adding each variables branch address by hand)
/// - which is the method preferred  by ROOT.
///
///  This object should be as close to a struct as possible, as such there
///  are no  associated methods, with the exception of the required
///  c'tor(s) and d'tors and  Reset() - to be called at the end of an
///  event,after TTree::Fill().  It resets all STL vectors (invoking
///  std::vector::clear()): all the elements of the vector are dropped
///  their destructors are called, and  then they are removed from the
///  vector container, leaving the container with a  size of 0.
///
///  The contents of the TTree::Branch associated with this object are
///  flattened at  the event level. i.e. for a single event, the branch
///  contains only: basic data  types; STL vectors of basic data types;
///  physics vectors and STL vectors of physics  vectors.
///
///  Additional modularity may be required in the future, one may choose
///  whether to define a new analysis level object in this header, derive
///  a new object from  LArAnalysis or simply extend the existing data in
///  LArAnalysis.
///
///  A point about filling STL vectors : std::vector::push_back() calls
///  std::vector::resize()  which is unnecessarily innefficient if the
///  size of the vector required is known a priori  (often the case if
///  you're reading an event record). In such a case, the vector should
///  be resized once and the values filled by index. If you inadvertently
///  exceed the bounds  of your vector, it will be resized automatically
///  (unlike a fixed length array).
///
/// One writes the object in the usual way, first compile in CINT
///   root [1] .x LArAnalysis.hxx+
///
/// or from the command line
///   $ root -b -q LArAnalysis.hxx+
///
/// Then, create an instance and put it in a TTree/TFile
///   root [2] anal = new LArAnalysis()
///   root [3] f = TFile::Open("test.root","RECREATE")
/// (class TFile*)0x1a9dcc0
///
///   root [4] t = new TTree("LArTree","LArTree")
/// (class TTree*)0x1902a70
///
///   root [5] t->Branch("LArAnalysis",&anal)
/// (class TBranch*)0x1be6490
///
///   root [6] anal->NVtxMC = 123
/// (const int)123
///
///   root [7] anal->Vtx4PosMC.push_back(TLorentzVector(1,2,3,4))
///   root [9] anal->Vtx4PosMC.push_back(TLorentzVector(5,6,7,8))
///
/// etc...
///
///   root [10] t->Fill()
/// (Int_t)448
///
///   root [12] f->Write()
/// (Int_t)11238
///
/// And, the object can be read as follows
///
///   root [0] gSystem->Load("LArAnalysis_hxx.so");
///   root [1] anal = new LArAnalysis()
/// (class LArAnalysis*)0x1cb4850
///
///   root [2] f = TFile::Open("test.root")
/// (class TFile*)0x1e25ee0
///
///   root [3] LArTree->SetBranchAddress("LArAnalysis",&anal)
/// (const Int_t)0
///
///   root [5] LArTree->GetEvent()
/// (Int_t)448
///
///   root [6] anal->NVtxMC
/// (int)123
///
///   root [7] anal->Vtx4PosMC[0]->Print()
/// (x,y,z,t)=(1.000000,2.000000,3.000000,4.000000) (P,eta,phi,E)=(3.741657,1.103587,1.107149,4.000000)
///   root [10] anal->Vtx4PosMC[1]->Print()
/// (x,y,z,t)=(5.000000,6.000000,7.000000,8.000000) (P,eta,phi,E)=(10.488088,0.806083,0.876058,8.000000)
///
/// etc.
///


#include <TLorentzVector.h>
#include <TVector3.h>
#include "TObject.h"
///
/// The following directives are required for ROOT to understand how to address
/// STL vectors of vectors/strings/TVectors etc.
/// This is documented in the "Adding a class" chapter of the ROOT User's Guide
/// and an example implementation can be found in the hvector.C tutorial:
/// http://root.cern.ch/root/html/tutorials/tree/hvector.C.html
///

#ifndef __CINT__
template class std::vector<std::string>;
template class std::vector<std::vector<std::string> >;
template class std::vector<std::vector<bool> >;
template class std::vector<std::vector<int> >;
template class std::vector<std::vector<double> >;
template class std::vector<TVector3>;
template class std::vector<std::vector<TVector3> >;
template class std::vector<TLorentzVector>;
template class std::vector<std::vector<TLorentzVector> >;
#endif

///
/// The LArAnalysis class definintion.
/// This object will constitute a single TBranch.
///
class LArAnalysis : public TObject{

///
/// Abbreviations
///

typedef TLorentzVector TLVec;
typedef TVector3       T3Vec;

public:

  const int               kUnassigned;           ///< unassigned parameters initialised to this

  ///
  /// Event Bookeeping Information
  ///
  float                   EventEField;           ///< electric field strength
  int                     EventEventRun;         ///< run number
  int                     EventID;               ///< event ID
  int                     EventSubRun;           ///< sub run number
  float                   EventT0;               ///< t_0 of the event

  ///
  /// MC Truth Trajectories
  ///
  int                          NTrajMC;               ///< number of true trajectories in the event
  std::vector<int>             TrajNHitsMC;
  std::vector<int>             TrajChargeMC;          ///< true charge             [NTrajMC]
  std::vector<TLVec>           TrajEnd4MomMC;         ///< true momentum at end    [NTrajMC][px,py,pz,E]
  std::vector<TLVec>           TrajEnd4PosMC;         ///< true position at end    [NTrajMC][x,y,z,t]
  std::vector<int>             TrajIDMC;              ///< true ID                 [NTrajMC]
  std::vector<double>          TrajMassMC;            ///< true particle mass      [NTrajMC]
  std::vector<std::string>     TrajNameMC;            ///< true particle name      [NTrajMC]
  std::vector<int>             TrajPDGMC;             ///< true particle pdg code  [NTrajMC]
  std::vector<int>             TrajParentIDMC;        ///< true particle parent ID [NTrajMC]
  std::vector<TLVec>           TrajStart4MomMC;       ///< true momentum at start  [NTrajMC][px,py,pz,E]
  std::vector<TLVec>           TrajStart4PosMC;       ///< true position at start  [NTrajMC][x,y,z,t]
  std::vector<double>          TrajLengthMC;

  ///
  /// MC Truth Vertices
  ///
  int                               NVtxMC;                ///< number of true vertices in the event
  std::vector<TLVec>                Vtx4PosMC;             ///< true vertex position                       [NVtxMC][x,y,z,t]
  std::vector<int>                  VtxNPrimTrajMC;        ///< number of true trajectories at this vertex [NVtxMC]
  std::vector<std::vector<int> >    VtxPrimTrajIDsMC;      ///< true IDs of trajectories at this vertex    [NVtxMC][VtxNPrimTrajMC]
  std::vector<TLVec>                VtxTargetMomMC;        ///< true momentum of target                    [NVtxMC][px,py,pz,E]
  std::vector<int>                  VtxTargetPDGMC;        ///< pdg code of target                         [NVtxMC]
  std::vector<int>                  VtxReactionCodeMC;
  std::vector<TLVec>                VtxNuMomMC;
  std::vector<int>                  VtxNuPDGMC;
  
  ///
  /// Reconstructed Tracks
  ///
  int                          NTracks;               ///< number of reconstructed tracks
  std::vector<int>             TrackCharge;           ///< reconstructed charge            [NTracks]
  std::vector<TLVec>           TrackEnd4Pos;          ///< reconstructed position at end   [NTracks][x,y,z,t]
  std::vector<int>             TrackID;               ///< reconstructed ID                [NTracks]
  std::vector<TLVec>           TrackStart4Pos;        ///< reconstructed position at start [NTracks][x,y,z,t]
  std::vector<int>             TrackNHits;            ///< number of hits along this track [NTracks]
  std::vector<float>           TrackPitch;            ///< reconstructed track pitch       [NTracks]
  std::vector<double>          TrackLength;

  ///
  /// Reconstructed Clusters
  ///
  int                     NClusters;             ///< number of reconstructed clusters
  std::vector<int>        ClusterCharge;         ///< reconstructed charge                [NClusters]
  std::vector<int>        ClusterID;             ///< reconstructed ID                    [NClusters]
  std::vector<int>        ClusterNHits;          ///< number of hits along this track     [NClusters]
  std::vector<int>        ClusterTrackID;        ///< ID of track this cluster belongs to [NClusters]

  ///
  /// Reconstructed Hits
  ///
  int                     NHits;                      ///< number of reconstructed hits
  std::vector<int>             HitChannel;            ///< reconstructed hit channel         [NHits]
  std::vector<float>           HitCharge;             ///< reconstructed hit charge          [NHits]
  std::vector<int>             HitClusterID;          ///< ID of cluster this hit belongs to [NHits]
  std::vector<float>           HitPeakT;              ///< reconstructed hit peak time       [NHits]
  std::vector<int>             HitPlane;              ///< reconstructed hit plane           [NHits]
  std::vector<int>             HitTrackID;            ///< ID of track this hit belongs to   [NHits]
  std::vector<int>             HitWire;               ///< reconstructed hit wire            [NHits]
  std::vector<int>             HitTPC;                ///< reconstructed hit TPC            [NHits]
  std::vector<T3Vec>           Hit3Pos;               ///< reconstructed 3D-matched position [NHits]
  std::vector<bool>            HitIsMatched;          ///< whether hit has a valid 3D position [NHits]

  std::vector<TVectorD>        EigenValues;
  std::vector<TMatrixD>        EigenVectors;
  std::vector<TMatrixD>        Covariance;
  ///
  /// Reset parameter values and vector containers
  ///
  virtual void Clear();

  ///
  /// Default constructor/destructor
  ///
  LArAnalysis();
  virtual ~LArAnalysis();
  ClassDef(LArAnalysis,1);
};

#endif
