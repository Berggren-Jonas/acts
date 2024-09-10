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

BOOST_AUTO_TEST_SUITE(GeoModelPlugin)




BOOST_AUTO_TEST_CASE(proTest) {
  // create pars for conversion
  Acts::GeoModelDetectorObjectFactory::Config gmConfig;
  Acts::GeoModelDetectorObjectFactory::Options gmOptions;
  gmConfig.convertBox = {"Box"};
  Acts::GeometryContext gContext;
  Acts::MagneticFieldContext mfContext;
  Acts::GeoModelDetectorObjectFactory::Cache gmCache;

  // create factory instance
  Acts::GeoModelDetectorObjectFactory factory(gmConfig);
  auto gmTree = Acts::GeoModelReader::readFromDb("~/ATLAS-R3-MUONTEST_v3.db");
  factory.construct(gmCache, gContext, gmTree, gmOptions);
  const auto gmBoxes = gmCache.boundingBoxes;
  const auto gmSensSurfaces = gmCache.sensitiveSurfaces;

  std::vector<std::shared_ptr<Acts::Surface>> gmSurfaces(gmSensSurfaces.size());
  std::transform(gmSensSurfaces.begin(), gmSensSurfaces.end(), gmSurfaces.begin(),
                 [](const std::tuple<std::shared_ptr<Acts::GeoModelDetectorElement>,
                                     std::shared_ptr<Acts::Surface>>& t) {
                   return std::get<1>(t);
                 });

  auto bounds = std::make_shared<Acts::CylinderVolumeBounds>(0, 15000, 25000);
  auto worldVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), gContext,
      "World_Detector_Volume",
      Acts::Transform3(Acts::Transform3::Identity() *
                 Acts::AngleAxis3(M_PI / 2, Acts::Vector3(0., 0., 1))),
      bounds, gmSurfaces, gmBoxes,
      Acts::Experimental::tryAllSubVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  worldVolume->assignGeometryId(Acts::GeometryIdentifier{}.setVolume(250));//TODO what is that?

  auto detector_sector = Acts::Experimental::Detector::makeShared(
      "Detector", {worldVolume}, Acts::Experimental::tryRootVolumes());

  // Set the navigator for the propagator
  Acts::Experimental::DetectorNavigator::Config navCfg;
  navCfg.detector = detector_sector.get();
  auto navigator = Acts::Experimental::DetectorNavigator(
      navCfg, Acts::getDefaultLogger("DetectorNavigator",
                                     Acts::Logging::Level::WARNING));
  // Set the stepper for the propagator with a magnetic field
  auto bField = std::make_shared<Acts::ConstantBField>(
      Acts::Vector3{0 * Acts::UnitConstants::T, 0, 0});
  auto stepper = Acts::EigenStepper<>(bField);

  // Set the propagator options and the propagator
  auto options =
      Acts::PropagatorOptions<Acts::ActionList, Acts::AbortList>(gContext, mfContext);
  auto propagator =
      Acts::Propagator<Acts::EigenStepper<>, Acts::Experimental::DetectorNavigator>(
          stepper, navigator,
          Acts::getDefaultLogger("Propagator", Acts::Logging::Level::WARNING));
}

BOOST_AUTO_TEST_SUITE_END()
