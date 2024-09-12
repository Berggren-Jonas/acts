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
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
//#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"

BOOST_AUTO_TEST_SUITE(GeoModelPlugin)

Acts::GeometryContext gContext;



BOOST_AUTO_TEST_CASE(proTest) {
  auto bounds = std::make_unique<Acts::CylinderVolumeBounds>(
      0, 15000, 25000);
  auto worldVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), gContext,
      "World_Detector_Volume",
      Acts::Transform3(Acts::Transform3::Identity() *
                 Acts::AngleAxis3(M_PI / 2, Acts::Vector3(0., 0., 1))),
      std::move(bounds), std::vector<std::shared_ptr<Acts::Surface>>{},
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>{},
      Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());
  //worldVolume->assignGeometryId(Acts::GeometryIdentifier{}.setVolume(250));
  auto rMax = worldVolume->volumeBounds().values()[1];
  auto hlengthZ = worldVolume->volumeBounds().values()[2];
  float theta = std::acos(hlengthZ / rMax);
  std::vector<float> pTValue = {1.};
  for (std::size_t i = 0; i < pTValue.size(); i++) {
    ActsExamples::ParametricParticleGenerator::Config pCfg;
    pCfg.thetaMin = theta;
    pCfg.thetaMax = M_PI - theta;
    ActsExamples::WhiteBoard eventStore(
        Acts::getDefaultLogger("Event_Store", Acts::Logging::Level::WARNING),
        {{"Particles_pT", std::to_string(pTValue[i])}});
    int njobs = 100;
    for (int nj = 0; nj < njobs; nj++) {
      ActsExamples::ParametricParticleGenerator pgenerator{pCfg};
      auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
          ActsExamples::RandomNumbers::Config{static_cast<uint64_t>(nj)});
      ActsExamples::AlgorithmContext alContext(0, i, eventStore);
      ActsExamples::RandomEngine randomEng = rnd->spawnGenerator(alContext);
      auto particles = pgenerator(randomEng);
      /*
      for (auto ip : particles) {
      }
      */
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()
