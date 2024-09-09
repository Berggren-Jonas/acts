// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
//#include "Acts/Core/Navigation/SurfaceCandidatesUpdaters.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
//#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
//#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"

#include <ctime>
#include <fstream>

/*
#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
*/
using namespace Acts;
using namespace Acts::Experimental;
//using namespace ActsExamples;
using namespace Acts::UnitLiterals;

namespace Acts::Test {

GeometryContext tContext;
MagneticFieldContext mfContext;
// AlgorithmContext alContext;

// A straw selector for the SurfaceCollector
struct StrawSelector {
  /// Call operator
  /// @param sf The input surface to be checked
  bool operator()(const Surface& sf) const {
    return (sf.type() == Surface::Straw);
  }
};

using ActionListType = ActionList<SurfaceCollector<StrawSelector>>;
using AbortListType = AbortList<Acts::EndOfWorldReached>;

// function that configures the mockup geometry and build the detector
//  for the navigation tests
std::shared_ptr<Acts::Experimental::DetectorVolume> buildMockup(
    MockupSectorBuilder::Config& mCfg) {
  auto mockup_chamberConfig_inner = MockupSectorBuilder::ChamberConfig();
  auto mockup_chamberConfig_middle = MockupSectorBuilder::ChamberConfig();
  auto mockup_chamberConfig_outer = MockupSectorBuilder::ChamberConfig();
}

/*
  //   mockup_config.gdmlPath =
  //       "
  //       ../../../../acts/Examples/Detectors/MuonSpectrometerMockupDetector/"
  //       "MuonChamber.gdml";
  //   mockup_config.NumberOfSectors = 8;
  //   mockup_config.binning = false;

  mockup_chamberConfig_inner.name = "Inner_Detector_Chamber";
  mockup_chamberConfig_inner.SensitiveNames = {"Inner_Skin"};
  mockup_chamberConfig_inner.PassiveNames = {"xx"};

  mockup_chamberConfig_middle.name = "Middle_Detector_Chamber";
  mockup_chamberConfig_middle.SensitiveNames = {"Middle_Skin"};
  mockup_chamberConfig_middle.PassiveNames = {"xx"};

  mockup_chamberConfig_outer.name = "Outer_Detector_Chamber";
  mockup_chamberConfig_outer.SensitiveNames = {"Outer_Skin"};
  mockup_chamberConfig_outer.PassiveNames = {"xx"};

  MockupSectorBuilder mockup_builder(mCfg);

  GeometryContext gctx;
  auto detectorVolume_inner_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_inner);

  auto detectorVolume_middle_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_middle);

  auto detectorVolume_outer_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_outer);

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detector_volumes = {detectorVolume_inner_chamber,
                          detectorVolume_middle_chamber,
                          detectorVolume_outer_chamber};

  auto detectorVolume_sector = mockup_builder.buildSector(detector_volumes);

  return detectorVolume_sector;
}

BOOST_AUTO_TEST_CASE(MuonMockupSpectrometer_TryAllNavigation) {
  // Build the detector that is going to be used for the propagation -the mockup
  // muon geometry

  // BUild three stations along z
  auto mockup_config = MockupSectorBuilder::Config();

  mockup_config.gdmlPath =
      " ../../../../acts/Examples/Detectors/MuonSpectrometerMockupDetector/"
      "MuonChamber.gdml";
  mockup_config.NumberOfSectors = 10;
  mockup_config.binning = true;
  // mockup_config.zOffset = 5000.;

  // auto stationpos = buildMockup(mockup_config);

  // mockup_config.zOffset=0;
  auto station0 = buildMockup(mockup_config);

  // mockup_config.zOffset = -5000.;
  // auto stationneg = buildMockup(mockup_config);

  auto bounds = std::make_unique<Acts::CylinderVolumeBounds>(
      0, station0->volumeBounds().values()[1],
      station0->volumeBounds().values()[2]);

  auto worldVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), tContext,
      "World_Detector_Volume",
      Transform3(Transform3::Identity() *
                 Acts::AngleAxis3(M_PI / 2, Acts::Vector3(0., 0., 1))),
      std::move(bounds), std::vector<std::shared_ptr<Acts::Surface>>{},
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>{
          station0},
      Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  worldVolume->assignGeometryId(Acts::GeometryIdentifier{}.setVolume(250));

  MockupSectorBuilder::drawSector(worldVolume, "sector_test.obj");

  auto detector_sector = Acts::Experimental::Detector::makeShared(
      "Detector", {worldVolume}, Acts::Experimental::tryRootVolumes());

  // Set the navigator for the propagator
  DetectorNavigator::Config navCfg;
  navCfg.detector = detector_sector.get();
  auto navigator = DetectorNavigator(
      navCfg, Acts::getDefaultLogger("DetectorNavigator",
                                     Acts::Logging::Level::WARNING));

  // Set the stepper for the propagator with a magnetic field
  auto bField = std::make_shared<ConstantBField>(
      Vector3{0 * Acts::UnitConstants::T, 0, 0});
  auto stepper = EigenStepper<>(bField);

  // Set the propagator options and the propagator
  auto options =
      PropagatorOptions<ActionListType, AbortListType>(tContext, mfContext);
  auto propagator =
      Propagator<Acts::EigenStepper<>, Acts::Experimental::DetectorNavigator>(
          stepper, navigator,
          Acts::getDefaultLogger("Propagator", Acts::Logging::Level::WARNING));

  options.pathLimit = 100_m;
  //  options.maxSteps = 20000;
  // options.maxStepSize = 1000_mm;

  // Test for different values of pT
  std::vector<float> pTValue = {1.};
  // Set the Random Engine

  std::cout << worldVolume->volumeBounds().values()[1] << std::endl;
  std::cout << worldVolume->volumeBounds().values()[2] << std::endl;
  std::cin.ignore();
  std::cout << std::atan(worldVolume->volumeBounds().values()[1] /
                         worldVolume->volumeBounds().values()[2])
            << std::endl;
  auto rMax = worldVolume->volumeBounds().values()[1];
  auto hlengthZ = worldVolume->volumeBounds().values()[2];
  float theta = std::acos(hlengthZ / rMax);

  double writtenEta{-1.f}, writtenPhi{-1.f}, writtenMom{-1.f}, cpuTime{-1.f},
      realTime{-1.f};
  unsigned int hittedSurf{0};

  std::unique_ptr<TFile> out_file{
      std::make_unique<TFile>("GridNavigation.root", "RECREATE")};
  std::unique_ptr<TTree> out_tree{
      std::make_unique<TTree>("PropagationTest", "Blub")};
  out_tree->Branch("particleEta", &writtenEta);
  out_tree->Branch("particlePhi", &writtenPhi);
  out_tree->Branch("partileMom", &writtenMom);
  out_tree->Branch("cpuTime", &cpuTime);
  out_tree->Branch("realTime", &realTime);
  out_tree->Branch("collectedHits", &hittedSurf);

  // Configure the particle generator for every pT value
  for (std::size_t i = 0; i < pTValue.size(); i++) {
    // std::unordered_map<string, string>
    // The file to store the event parameters
    std::ofstream outfile;
    std::ofstream position;
    outfile.open("TryAllEvents_pT_" + std::to_string(pTValue[i]) + ".txt");
    position.open("TryAllHitPoints_pT" + std::to_string(pTValue[i]) + ".txt");

    WhiteBoard eventStore(
        getDefaultLogger("Event_Store", Acts::Logging::Level::WARNING),
        {{"Particles_pT", std::to_string(pTValue[i])}});

    ParametricParticleGenerator::Config pCfg;
    pCfg.thetaMin = theta;
    pCfg.thetaMax = M_PI - theta;
    // std::cout<<"theta min="<<theta<<std::endl;
    // std::cout<<"theta max="<<M_PI - theta<<std::endl;
    // std::cin.ignore();

    // pCfg.thetaMin = 91*M_PI/180.;
    // pCfg.thetaMax = 91*M_PI/180.;
    // pCfg.thetaMin = 0.5*M_PI;
    // pCfg.thetaMax = 0.5*M_PI;

    pCfg.etaUniform = true;
    pCfg.phiMin = 0;
    pCfg.phiMax = 2 * M_PI;
    pCfg.pMin = pTValue[i] * Acts::UnitConstants::GeV;
    pCfg.pMax = pTValue[i] * Acts::UnitConstants::GeV;
    pCfg.pTransverse = true;
    pCfg.numParticles = 10000;
    int njobs = 100;

    for (int nj = 0; nj < njobs; nj++) {
      ParametricParticleGenerator pgenerator{pCfg};

      auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
          ActsExamples::RandomNumbers::Config{nj});

      AlgorithmContext alContext(0, i, eventStore);
      RandomEngine randomEng = rnd->spawnGenerator(alContext);

      auto particles = pgenerator(randomEng);
      // std::cout<<"particles size="<<particles.size()<<std::endl;
      // std::cin.ignore();

      // simulate the tracks using the parameters from the particle generator

      for (auto ip : particles) {
        Vector4 pos = ip.fourPosition();
        Vector3 mom = ip.momentum();
        ActsScalar pT = ip.transverseMomentum();
        auto eta = -std::log(std::tan(ip.theta() / 2));

        // std::cin.ignore();
        ActsScalar qOverp = ip.qOverP();
        // std::cout<<qOverp<<std::endl;
        // std::cin.ignore();
        ParticleHypothesis phypothesis = ip.hypothesis();
        CurvilinearTrackParameters start(pos, ip.phi(), ip.theta(), qOverp,
                                         std::nullopt, phypothesis);
        //outfile<<eta<<"\t";
        //outfile<<ip.phi()<<"\t";

        TStopwatch watch{};
        watch.Start();

        const auto& presult = propagator.propagate(start, options).value();

        watch.Stop();
        realTime = (watch.RealTime() * 1000.);
        cpuTime = (watch.CpuTime() * 1000.);

        auto& cSurfaces =
            presult.get<SurfaceCollector<StrawSelector>::result_type>();
        hittedSurf = cSurfaces.collected.size();
        writtenEta = eta;
        writtenPhi = ip.phi();
        writtenMom = pT;
        out_tree->Fill();

        // std::cout<<"pt="<<pT<<std::endl;
        // std::cout<<"theta="<<ip.theta()<<std::endl;
        // std::cout<<"phi="<<ip.phi()<<std::endl;
        // std::cout<<"collected
        // surfaces="<<cSurfaces.collected.size()<<std::endl;

        // outfile<<cSurfaces.collected.size()<<"\t";

        // for(auto& sc : cSurfaces.collected){
        //    auto hittedsurf = sc.surface;
        //    position<<ip.phi()<<"\t"<<sc.position.x()<<"\t"<<sc.position.y()<<"\t"<<sc.position.z()<<"\n";
        //  }

        // outfile<< 1000.0*(c_end - c_start)/CLOCKS_PER_SEC<<"\n";
        // std::cin.ignore();
      }
    }

    // outfile.close();
    // position.close();
  }
  out_file->WriteObject(out_tree.get(), out_tree->GetName());
  out_tree.reset();

  // define start parameters
  //   Vector4 pos(0, 4750., 0, 0);
  //   Vector3 mom(0, 10, 0);
  //   CurvilinearTrackParameters start(pos, mom, +1_e / mom.norm(),
  //   std::nullopt,
  //                                    Acts::ParticleHypothesis::pion());
  //   const auto& presult = propagator.propagate(start, options).value();

  //   auto& cSurfaces =
  //   presult.get<SurfaceCollector<StrawSelector>::result_type>();

  // we expect 20 straw surfaces collected (one from each tube layer of each
  // chamber)
  // BOOST_CHECK_EQUAL(cSurfaces.collected.size(), 20);
}

BOOST_AUTO_TEST_CASE(MuonMockupSpectrometer_GridNavigation) {}

*/
}  // namespace Acts::Test
