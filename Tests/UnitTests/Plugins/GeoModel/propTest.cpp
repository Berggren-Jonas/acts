// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

// clang-formal off
#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
// clang-formal on
#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
// #include
// "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"

#include <ctime>
#include <fstream>

#include <TFile.h>
#include <TStopwatch.h>
#include <TTree.h>
BOOST_AUTO_TEST_SUITE(GeoModelPlugin)

Acts::GeometryContext gContext;
Acts::MagneticFieldContext mfContext;
struct StrawSelector {
  /// Call operator
  /// @param sf The input surface to be checked
  bool operator()(const Acts::Surface& sf) const {
    return (sf.type() == Acts::Surface::Straw);
  }
};

using ActionListType = Acts::ActionList<Acts::SurfaceCollector<StrawSelector>>;
using AbortListType = Acts::AbortList<Acts::EndOfWorldReached>;

BOOST_AUTO_TEST_CASE(proTest) {
  // Setting parameters for conversion
  Acts::GeoModelDetectorObjectFactory::Config factoryCfg;
  factoryCfg.materialList = {"Aluminium"};
  factoryCfg.nameList = {"MDT", "Tube"};
  factoryCfg.convertBox = {"MDT"};
  Acts::GeoModelDetectorObjectFactory::Options factoryOpt;
  factoryOpt.queries = {"Muon"};
  Acts::GeoModelDetectorObjectFactory::Cache cache;
  Acts::GeoModelTree tree = Acts::GeoModelReader::readFromDb(
      "/home/cberggre/ATLAS-R3-MUONTEST_v3.db");

  // converting
  auto factory = Acts::GeoModelDetectorObjectFactory(factoryCfg);
  factory.construct(cache, gContext, tree, factoryOpt);
  auto sensSurfaces = cache.sensitiveSurfaces;
  auto boxes = cache.boundingBoxes;
  std::vector<std::shared_ptr<Acts::Surface>> surfaces(sensSurfaces.size());
  std::transform(
      sensSurfaces.begin(), sensSurfaces.end(), surfaces.begin(),
      [](const std::tuple<std::shared_ptr<Acts::GeoModelDetectorElement>,
                          std::shared_ptr<Acts::Surface>>& t) {
        return std::get<1>(t);
      });
  std::vector<int> ids(16, 1);
  for (int j = 0; j < boxes.size(); j++) {
    // Setting the layer idenifier according to phi sector and +- eta
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream stream(boxes[j]->name());
    while (std::getline(stream, token, '_')) {
      tokens.push_back(token);
    }
    bool negEta = boxes[j]->name().find("-") != std::string::npos;
    int phiSector;
    if (negEta) {
      phiSector = 8 + std::stoi(tokens[tokens.size() - 2]);
    } else {
      phiSector = std::stoi(tokens[tokens.size() - 2]);
    }
    boxes[j]->assignGeometryId(
        Acts::GeometryIdentifier{}.setLayer(phiSector).setVolume(
            ids[phiSector - 1]));
    // setting the identifier for the tubes
    auto tubes = boxes[j]->surfacePtrs();
    for (int k = 0; k < tubes.size(); k++) {
      tubes[k]->assignGeometryId(Acts::GeometryIdentifier{}
                                     .setLayer(phiSector)
                                     .setVolume(ids[phiSector - 1])
                                     .setSensitive(k + 1));
    }
    ids[phiSector - 1] = ids[phiSector - 1] + 1;
  }
  // for (int i=0;i<surfaces.size();i++){
  // surfaces[i]->assignGeometryId(Acts::GeometryIdentifier{}.setLayer(chId.first).setVolume(chId.second).setSensitive(++surfId));
  //}

  // constructing word volume
  auto bounds = std::make_unique<Acts::CylinderVolumeBounds>(0, 15000, 25000);
  auto worldVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), gContext,
      "World_Detector_Volume",
      Acts::Transform3(Acts::Transform3::Identity() *
                       Acts::AngleAxis3(M_PI / 2, Acts::Vector3(0., 0., 1))),
      std::move(bounds), std::vector<std::shared_ptr<Acts::Surface>>{}, boxes,
      Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());
  worldVolume->assignGeometryId(
      Acts::GeometryIdentifier{}.setVolume(surfaces.size() + boxes.size()));

  // geometry of the world volume
  auto rMax = worldVolume->volumeBounds().values()[1];
  auto hlengthZ = worldVolume->volumeBounds().values()[2];
  float theta = std::acos(hlengthZ / rMax);
  std::vector<float> pTValue = {1.};

  // iterate over pt values
  // generate particles
  ActsExamples::ParametricParticleGenerator::Config pCfg;
  pCfg.thetaMin = theta;
  pCfg.thetaMax = M_PI - theta;
  pCfg.etaUniform = true;
  pCfg.phiMin = 0;
  pCfg.phiMax = 2 * M_PI;
  pCfg.pMin = Acts::UnitConstants::GeV;
  pCfg.pMax = Acts::UnitConstants::GeV;
  pCfg.pTransverse = true;
  pCfg.numParticles = 10000;
  ActsExamples::WhiteBoard eventStore(
      Acts::getDefaultLogger("Event_Store", Acts::Logging::Level::WARNING),
      {{"Particles_pT", std::to_string(1)}});
  ActsExamples::ParametricParticleGenerator pgenerator{pCfg};
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
      ActsExamples::RandomNumbers::Config{static_cast<uint64_t>(0)});
  ActsExamples::AlgorithmContext alContext(0, 0, eventStore);
  ActsExamples::RandomEngine randomEng = rnd->spawnGenerator(alContext);
  auto particles = std::get<1>(pgenerator(randomEng));
  for (auto ip : particles) {
    Acts::Vector4 pos = ip.fourPosition();
    Acts::Vector3 mom = ip.momentum();
    Acts::ActsScalar pT = ip.transverseMomentum();
    auto eta = -std::log(std::tan(ip.theta() / 2));
    Acts::ActsScalar qOverp = ip.qOverP();
    Acts::ParticleHypothesis phypothesis = ip.hypothesis();
    std::cout << ip.phi() << " " << ip.theta() << std::endl;
    Acts::CurvilinearTrackParameters start({0, 0, 0, 0}, ip.phi(), ip.theta(),
                                           qOverp, std::nullopt, phypothesis);
    using Propagator = Acts::Propagator<Acts::StraightLineStepper,
                                        Acts::Experimental::DetectorNavigator>;
    // using PropagatorOptions = Propagator::Options<>;
    using PropagatorOptions =
        Propagator::Options<ActionListType, AbortListType>;
    // Set the stepper for the propagator with a magnetic field
    auto stepper = Acts::StraightLineStepper();
    Acts::Experimental::DetectorNavigator::Config navCfg;
    auto detector_sector = Acts::Experimental::Detector::makeShared(
        "Detector", {worldVolume}, Acts::Experimental::tryRootVolumes());
    navCfg.detector = detector_sector.get();
    auto navigator = Acts::Experimental::DetectorNavigator(
        navCfg, Acts::getDefaultLogger("DetectorNavigator",
                                       Acts::Logging::Level::WARNING));
    Propagator propagator(stepper, navigator);

    PropagatorOptions options(gContext, {});
    options.direction = Acts::Direction::Backward;
    const auto& presult = propagator.propagate(start, options).value();
    auto& cSurfaces =
        presult.get<Acts::SurfaceCollector<StrawSelector>::result_type>();
  }
}
BOOST_AUTO_TEST_SUITE_END()
