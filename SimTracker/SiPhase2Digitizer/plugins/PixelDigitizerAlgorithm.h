#ifndef _SimTracker_SiPhase2Digitizer_PixelDigitizerAlgorithm_h
#define _SimTracker_SiPhase2Digitizer_PixelDigitizerAlgorithm_h

#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerAlgorithm.h"

class PixelDigitizerAlgorithm : public Phase2TrackerDigitizerAlgorithm {
private:
  class TimewalkModel {
  public:
    // construct the model from a file.
    // The file must contain 3 lines:
    // 1st line: comma-sepparated values of input charge
    // 2nd line: comma-sepparated values of threshold
    // 3rd line: comma-sepparated values of delay, one for each (input charge, threshold) pair
    TimewalkModel(const std::string& filename);

    // returns delay for given input charge and threshold
    double operator()(double q_in, double q_threshold) const;

  private:
    template <class Stream>
    void parse_csv_line(Stream& stream, std::vector<double>& vec);

    std::size_t find_closest_index(const std::vector<double>& vec, double value) const;

    std::vector<double> input_charge;
    std::vector<double> threshold;
    std::vector<double> delay;
  };

public:
  PixelDigitizerAlgorithm(const edm::ParameterSet& conf);
  ~PixelDigitizerAlgorithm() override;

  // initialization that cannot be done in the constructor
  void init(const edm::EventSetup& es) override;

  // void initializeEvent();
  // run the algorithm to digitize a single det
  void accumulateSimHits(const std::vector<PSimHit>::const_iterator inputBegin,
                         const std::vector<PSimHit>::const_iterator inputEnd,
                         const size_t inputBeginGlobalIndex,
                         const uint32_t tofBin,
                         const Phase2TrackerGeomDetUnit* pixdet,
                         const GlobalVector& bfield) override;
  void add_cross_talk(const Phase2TrackerGeomDetUnit* pixdet) override;

  void digitize(const Phase2TrackerGeomDetUnit* pixdet,
                std::map<int, DigitizerUtility::DigiSimInfo>& digi_map,
                const TrackerTopology* tTopo) override;

  // Addition four xtalk-related parameters to PixelDigitizerAlgorithm specific parameters initialized in Phase2TrackerDigitizerAlgorithm
  const double odd_row_interchannelCoupling_next_row_;
  const double even_row_interchannelCoupling_next_row_;
  const double odd_column_interchannelCoupling_next_column_;
  const double even_column_interchannelCoupling_next_column_;

  const TimewalkModel timewalk_model;
};
#endif
