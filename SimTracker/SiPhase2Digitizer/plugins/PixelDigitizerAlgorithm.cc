#include <iostream>
#include <cmath>

#include "SimTracker/SiPhase2Digitizer/plugins/PixelDigitizerAlgorithm.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineSimService.h"

#include "CondFormats/SiPixelObjects/interface/GlobalPixel.h"
#include "CondFormats/DataRecord/interface/SiPixelQualityRcd.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleSimRcd.h"

// Geometry
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "CLHEP/Random/RandGaussQ.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

using namespace edm;
using namespace sipixelobjects;

void PixelDigitizerAlgorithm::init(const edm::EventSetup& es) {
  if (use_ineff_from_db_)  // load gain calibration service fromdb...
    theSiPixelGainCalibrationService_->setESObjects(es);

  if (use_deadmodule_DB_)
    es.get<SiPixelQualityRcd>().get(SiPixelBadModule_);

  if (use_LorentzAngle_DB_)  // Get Lorentz angle from DB record
    es.get<SiPixelLorentzAngleSimRcd>().get(SiPixelLorentzAngle_);

  // gets the map and geometry from the DB (to kill ROCs)
  es.get<SiPixelFedCablingMapRcd>().get(fedCablingMap_);
  es.get<TrackerDigiGeometryRecord>().get(geom_);
}

PixelDigitizerAlgorithm::PixelDigitizerAlgorithm(const edm::ParameterSet& conf)
    : Phase2TrackerDigitizerAlgorithm(conf.getParameter<ParameterSet>("AlgorithmCommon"),
                                      conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")),
      odd_row_interchannelCoupling_next_row_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
                                                 .getParameter<double>("Odd_row_interchannelCoupling_next_row")),
      even_row_interchannelCoupling_next_row_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
                                                  .getParameter<double>("Even_row_interchannelCoupling_next_row")),
      odd_column_interchannelCoupling_next_column_(
          conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
              .getParameter<double>("Odd_column_interchannelCoupling_next_column")),
      even_column_interchannelCoupling_next_column_(
          conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
              .getParameter<double>("Even_column_interchannelCoupling_next_column")),
      timewalk_model(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<std::string>("TimewalkModelDataFile")) {
  pixelFlag_ = true;
  LogInfo("PixelDigitizerAlgorithm") << "Algorithm constructed "
                                     << "Configuration parameters:"
                                     << "Threshold/Gain = "
                                     << "threshold in electron Endcap = " << theThresholdInE_Endcap_
                                     << "threshold in electron Barrel = " << theThresholdInE_Barrel_ << " "
                                     << theElectronPerADC_ << " " << theAdcFullScale_ << " The delta cut-off is set to "
                                     << tMax_ << " pix-inefficiency " << addPixelInefficiency_;
}
PixelDigitizerAlgorithm::~PixelDigitizerAlgorithm() { LogDebug("PixelDigitizerAlgorithm") << "Algorithm deleted"; }
void PixelDigitizerAlgorithm::accumulateSimHits(std::vector<PSimHit>::const_iterator inputBegin,
                                                std::vector<PSimHit>::const_iterator inputEnd,
                                                const size_t inputBeginGlobalIndex,
                                                const uint32_t tofBin,
                                                const Phase2TrackerGeomDetUnit* pixdet,
                                                const GlobalVector& bfield) {
  // produce SignalPoint's for all SimHit's in detector
  // Loop over hits
  uint32_t detId = pixdet->geographicalId().rawId();
  size_t simHitGlobalIndex = inputBeginGlobalIndex;  // This needs to be stored to create the digi-sim link later

  // find the relevant hits
  std::vector<PSimHit> matchedSimHits;
  std::copy_if(inputBegin, inputEnd, std::back_inserter(matchedSimHits), [detId](auto const& hit) -> bool {
    return hit.detUnitId() == detId;
  });
  // loop over a much reduced set of SimHits
  for (auto& hit : matchedSimHits) {
    LogDebug("PixelDigitizerAlgorithm") << hit.particleType() << " " << hit.pabs() << " " << hit.energyLoss() << " "
                                        << hit.tof() << " " << hit.trackId() << " " << hit.processType() << " "
                                        << hit.detUnitId() << hit.entryPoint() << " " << hit.exitPoint();

    std::vector<DigitizerUtility::EnergyDepositUnit> ionization_points;
    std::vector<DigitizerUtility::SignalPoint> collection_points;

    // apply correction to tof
    double time = hit.tof() - pixdet->surface().toGlobal((hit).localPosition()).mag() / 30.;
    hit.setTof(time);

    // check if the hit arrived durring this bunch crossing
    if (hit.tof() >= theTofLowerCut_ && hit.tof() <= theTofUpperCut_) {
      primary_ionization(hit, ionization_points);  // fills ionization_points

      // transforms ionization_points -> collection_points
      drift(hit, pixdet, bfield, ionization_points, collection_points);

      // compute induced signal on readout elements and add to _signal
      // hit needed only for SimHit<-->Digi link
      induce_signal(hit, simHitGlobalIndex, tofBin, pixdet, collection_points);
    }
    ++simHitGlobalIndex;
  }
}

// ======================================================================
//
//  Add  Cross-talk contribution
//
// ======================================================================
void PixelDigitizerAlgorithm::add_cross_talk(const Phase2TrackerGeomDetUnit* pixdet) {
  if (!pixelFlag_)
    return;

  const Phase2TrackerTopology* topol = &pixdet->specificTopology();

  // cross-talk calculation valid for the case of 25x100 pixels
  const float pitch_first = 0.0025;
  const float pitch_second = 0.0100;

  // 0.5 um tolerance when comparing the pitch to accommodate the small changes in different TK geometrie (temporary fix)
  const double pitch_tolerance(0.0005);

  if (std::abs(topol->pitch().first - pitch_first) > pitch_tolerance ||
      std::abs(topol->pitch().second - pitch_second) > pitch_tolerance)
    return;

  uint32_t detID = pixdet->geographicalId().rawId();
  signal_map_type& theSignal = _signal[detID];
  signal_map_type signalNew;

  int numRows = topol->nrows();
  int numColumns = topol->ncolumns();

  for (auto& s : theSignal) {
    float signalInElectrons = s.second.ampl();  // signal in electrons

    auto hitChan = PixelDigi::channelToPixel(s.first);

    float signalInElectrons_odd_row_Xtalk_next_row = signalInElectrons * odd_row_interchannelCoupling_next_row_;
    float signalInElectrons_even_row_Xtalk_next_row = signalInElectrons * even_row_interchannelCoupling_next_row_;
    float signalInElectrons_odd_column_Xtalk_next_column =
        signalInElectrons * odd_column_interchannelCoupling_next_column_;
    float signalInElectrons_even_column_Xtalk_next_column =
        signalInElectrons * even_column_interchannelCoupling_next_column_;

    // subtract the charge which will be shared
    s.second.set(signalInElectrons - signalInElectrons_odd_row_Xtalk_next_row -
                 signalInElectrons_even_row_Xtalk_next_row - signalInElectrons_odd_column_Xtalk_next_column -
                 signalInElectrons_even_column_Xtalk_next_column);

    if (hitChan.first != 0) {
      auto XtalkPrev = std::make_pair(hitChan.first - 1, hitChan.second);
      int chanXtalkPrev = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second);
      if (hitChan.first % 2 == 1)
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_even_row_Xtalk_next_row, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_row_Xtalk_next_row, nullptr, -1.0));
    }
    if (hitChan.first < numRows - 1) {
      auto XtalkNext = std::make_pair(hitChan.first + 1, hitChan.second);
      int chanXtalkNext = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkNext.first, XtalkNext.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkNext.first, XtalkNext.second);
      if (hitChan.first % 2 == 1)
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_row_Xtalk_next_row, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_even_row_Xtalk_next_row, nullptr, -1.0));
    }

    if (hitChan.second != 0) {
      auto XtalkPrev = std::make_pair(hitChan.first, hitChan.second - 1);
      int chanXtalkPrev = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second);
      if (hitChan.second % 2 == 1)
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_even_column_Xtalk_next_column, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkPrev,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_column_Xtalk_next_column, nullptr, -1.0));
    }
    if (hitChan.second < numColumns - 1) {
      auto XtalkNext = std::make_pair(hitChan.first, hitChan.second + 1);
      int chanXtalkNext = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkNext.first, XtalkNext.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkNext.first, XtalkNext.second);
      if (hitChan.second % 2 == 1)
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_odd_column_Xtalk_next_column, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkNext,
                          DigitizerUtility::Amplitude(signalInElectrons_even_column_Xtalk_next_column, nullptr, -1.0));
    }
  }
  for (auto const& l : signalNew) {
    int chan = l.first;
    auto iter = theSignal.find(chan);
    if (iter != theSignal.end()) {
      iter->second += l.second.ampl();
    } else {
      theSignal.emplace(chan, DigitizerUtility::Amplitude(l.second.ampl(), nullptr, -1.0));
    }
  }
}

void PixelDigitizerAlgorithm::digitize(const Phase2TrackerGeomDetUnit* pixdet,
                                               std::map<int, DigitizerUtility::DigiSimInfo>& digi_map,
                                               const TrackerTopology* tTopo) {
  uint32_t detID = pixdet->geographicalId().rawId();
  auto it = _signal.find(detID);
  if (it == _signal.end())
    return;

  const signal_map_type& theSignal = _signal[detID];

  uint32_t Sub_detid = DetId(detID).subdetId();

  float theThresholdInE = 0.;
  float theHIPThresholdInE = 0.;
  // Define Threshold
  if (Sub_detid == PixelSubdetector::PixelBarrel || Sub_detid == StripSubdetector::TOB) {  // Barrel modules
    theThresholdInE = addThresholdSmearing_ ? smearedThreshold_Barrel_->fire()             // gaussian smearing
                                            : theThresholdInE_Barrel_;                     // no smearing
    theHIPThresholdInE = theHIPThresholdInE_Barrel_;
  } else {                                                                      // Forward disks modules
    theThresholdInE = addThresholdSmearing_ ? smearedThreshold_Endcap_->fire()  // gaussian smearing
                                            : theThresholdInE_Endcap_;          // no smearing
    theHIPThresholdInE = theHIPThresholdInE_Endcap_;
  }

  //  if (addNoise) add_noise(pixdet, theThresholdInE/theNoiseInElectrons_);  // generate noise
  if (addNoise_)
    add_noise(pixdet);  // generate noise
  if (addXtalk_)
    add_cross_talk(pixdet);
  if (addNoisyPixels_)
    add_noisy_cells(pixdet, theHIPThresholdInE / theElectronPerADC_);

  // Do only if needed
  if (addPixelInefficiency_ && !theSignal.empty()) {
    if (use_ineff_from_db_)
      pixel_inefficiency_db(detID);
    else
      pixel_inefficiency(subdetEfficiencies_, pixdet, tTopo);
  }
  if (use_module_killing_) {
    if (use_deadmodule_DB_)  // remove dead modules using DB
      module_killing_DB(pixdet);
    else  // remove dead modules using the list in cfg file
      module_killing_conf(detID);
  }

  // Digitize if the signal is greater than threshold
  for (auto const& s : theSignal) {
    const DigitizerUtility::Amplitude& sig_data = s.second;
    float signalInElectrons = sig_data.ampl();

    const auto& info_list = sig_data.simInfoList();
    const auto it = std::max_element(info_list.begin(), info_list.end());
    const DigitizerUtility::SimHitInfo* hit_info = it->second.get();
    if (hit_info) {
      double time = hit_info->time() + timewalk_model(signalInElectrons, theThresholdInE);
      if (time < theTofLowerCut_ || time > theTofUpperCut_)
        continue;
    }

    if (signalInElectrons >= theThresholdInE) {  // check threshold
      DigitizerUtility::DigiSimInfo info;
      info.sig_tot = convertSignalToAdc(detID, signalInElectrons, theThresholdInE);  // adc
      info.ot_bit = signalInElectrons > theHIPThresholdInE ? true : false;
      if (makeDigiSimLinks_) {
        for (auto const& l : sig_data.simInfoList()) {
          float charge_frac = l.first / signalInElectrons;
          if (l.first > -5.0)
            info.simInfoList.push_back({charge_frac, l.second.get()});
        }
      }
      digi_map.insert({s.first, info});
    }
  }
}

PixelDigitizerAlgorithm::TimewalkModel::TimewalkModel(const std::string& file_path) {
  try {
    std::ifstream file(file_path);
    parse_csv_line(file, input_charge);
    parse_csv_line(file, threshold);
    parse_csv_line(file, delay);
  }
  catch (std::exception& e) {
    throw cms::Exception("Configuration") << "Timewalk model data file (" << file_path << ") error: " << e.what();
  }
  if (delay.size() != input_charge.size() * threshold.size())
    throw cms::Exception("Configuration") << "Timewalk model data file (" << file_path << ") error: series have incompatible size.";
}

double PixelDigitizerAlgorithm::TimewalkModel::operator()(double q_in, double q_threshold) const {
  auto index_x = find_closest_index(input_charge, q_in);
  auto index_y = find_closest_index(threshold, q_threshold);
  return delay[index_x * threshold.size() + index_y];
}

template <class Stream>
void PixelDigitizerAlgorithm::TimewalkModel::parse_csv_line(Stream& stream, std::vector<double>& vec) {
  std::string line;
  std::getline(stream, line);
  std::istringstream ss(line);
  std::string value;
  while (std::getline(ss, value, ',')) {
    vec.push_back(std::stod(value));
  }
}

std::size_t PixelDigitizerAlgorithm::TimewalkModel::find_closest_index(const std::vector<double>& vec, double value) const {
    auto it = std::lower_bound(vec.begin(), vec.end(), value);
  
    if (it == vec.begin()) return 0;
    else if (it == vec.end()) return vec.size() - 1;
    else {
      auto it_upper = it;
      auto it_lower = --it;
      
      auto closest = (value - *it_lower > *it_upper - value) ? it_upper : it_lower;
      return std::distance(vec.begin(), closest);
    }
}
