import FWCore.ParameterSet.Config as cms

phase2TrackerDigitizer = cms.PSet(
# For the Digitizer
    accumulatorType = cms.string("Phase2TrackerDigitizer"),
    hitsProducer = cms.string('g4SimHits'),
    ROUList = cms.vstring(
        'TrackerHitsPixelBarrelLowTof', 
        'TrackerHitsPixelBarrelHighTof', 
        'TrackerHitsPixelEndcapLowTof', 
        'TrackerHitsPixelEndcapHighTof'),
    GeometryType = cms.string('idealForDigi'),
    isOTreadoutAnalog = cms.bool(False),#set this to true if you want analog readout for OT
# Common for Algos
    premixStage1 = cms.bool(False),
    AlgorithmCommon = cms.PSet(
      DeltaProductionCut = cms.double(0.03),
      makeDigiSimLinks = cms.untracked.bool(True),
    ),
# Specific parameters
#Pixel Digitizer Algorithm
    PixelDigitizerAlgorithm = cms.PSet(
      ElectronPerAdc = cms.double(600.0),
      ReadoutNoiseInElec = cms.double(0.0),
      ThresholdInElectrons_Barrel = cms.double(1200.0),
      ThresholdInElectrons_Endcap = cms.double(1200.0),
      AddThresholdSmearing = cms.bool(False),
      ThresholdSmearing_Barrel = cms.double(0.0),
      ThresholdSmearing_Endcap = cms.double(0.0),
      HIPThresholdInElectrons_Barrel = cms.double(1.0e10), # very high value to avoid Over threshold bit
      HIPThresholdInElectrons_Endcap = cms.double(1.0e10), # very high value to avoid Over threshold bit
      NoiseInElectrons = cms.double(0.0),
      Phase2ReadoutMode = cms.int32(-1), # Flag to decide Readout Mode :Digital(0) or Analog (linear TDR (-1), dual slope with slope parameters (+1,+2,+3,+4) with threshold subtraction
      AdcFullScale = cms.int32(15),
      TofUpperCut = cms.double(12.5),
      TofLowerCut = cms.double(-12.5),
      AddNoisyPixels = cms.bool(False),
      Alpha2Order = cms.bool(True),			#D.B.: second order effect, does not switch off magnetic field as described
      AddNoise = cms.bool(False),
      AddXTalk = cms.bool(False),			#D.B.
      InterstripCoupling = cms.double(0.0),	#D.B. # No need to be used in PixelDigitizerAlgorithm
      Odd_row_interchannelCoupling_next_row = cms.double(0.20),
      Even_row_interchannelCoupling_next_row = cms.double(0.0),
      Odd_column_interchannelCoupling_next_column = cms.double(0.0),
      Even_column_interchannelCoupling_next_column = cms.double(0.0),
      SigmaZero = cms.double(0.00037),  		#D.B.: 3.7um spread for 300um-thick sensor, renormalized in digitizerAlgo
      SigmaCoeff = cms.double(0),  		#S.D: setting SigmaCoeff=0 for IT-pixel
      ClusterWidth = cms.double(3),		#D.B.: this is used as number of sigmas for charge collection (3=+-3sigmas)
      LorentzAngle_DB = cms.bool(True),			
      TanLorentzAnglePerTesla_Endcap = cms.double(0.106),
      TanLorentzAnglePerTesla_Barrel = cms.double(0.106),
      KillModules = cms.bool(False),
      DeadModules_DB = cms.bool(False),
      DeadModules = cms.VPSet(),
      AddInefficiency = cms.bool(False),
      Inefficiency_DB = cms.bool(False),
      EfficiencyFactors_Barrel = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999 ),
      EfficiencyFactors_Endcap = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 
      0.999, 0.999 ),#Efficiencies kept as Side2Disk1,Side1Disk1 and so on
      CellsToKill = cms.VPSet(),
      
      # new parameters
      doValidation = cms.bool(True),
      ValidationOutputFile = cms.string('validationHistograms.root'),
      # RD53A
      TimewalkModel = cms.PSet(
        ThresholdValues = cms.vdouble(1000, 1200),
        Curves = cms.VPSet(
            cms.PSet(
                charge = cms.vdouble(1.000E3, 1.050E3, 1.100E3, 1.150E3, 1.200E3, 1.250E3, 1.300E3, 1.400E3, 1.500E3, 1.600E3, 1.700E3, 1.800E3, 1.900E3, 2.000E3, 4.000E3, 6.000E3, 8.000E3, 10.00E3, 15.00E3, 20.00E3, 25.00E3),
                delay = cms.vdouble(43.51, 34.410000000000004, 30.45, 27.73, 25.650000000000002, 23.96, 22.55, 20.29, 18.54, 17.11, 15.92, 14.92, 14.040000000000001, 13.270000000000001, 6.51, 4.234, 3.0300000000000002, 2.272, 1.188, 0.601, 0.2351)
            ),
            cms.PSet(
                charge = cms.vdouble(1.200E3, 1.250E3, 1.300E3, 1.350E3, 1.400E3, 1.450E3, 1.500E3, 1.600E3, 1.700E3, 1.800E3, 1.900E3, 2.000E3, 4.000E3, 6.000E3, 8.000E3, 10.00E3, 15.00E3, 20.00E3, 25.00E3),
                delay = cms.vdouble(42.78, 33.589999999999996, 29.87, 27.32, 25.37, 23.779999999999998, 22.44, 20.29, 18.6, 17.220000000000002, 16.07, 15.069999999999999, 7.032, 4.516, 3.214, 2.404, 1.2590000000000001, 0.6393, 0.2535)
            )
        )
      )
      # RD53B
    #   TimewalkModel = cms.PSet(
    #     ThresholdValues = cms.vdouble(1000, 1200),
    #     Curves = cms.VPSet(
    #         cms.PSet(
    #             charge = cms.vdouble(1.000E3, 1.050E3, 1.100E3, 1.200E3, 1.300E3, 1.400E3, 1.500E3, 1.600E3, 1.700E3, 1.800E3, 1.900E3, 2.000E3, 4.000E3, 6.000E3, 8.000E3, 10.00E3, 15.00E3, 20.00E3, 25.00E3),
    #             delay = cms.vdouble(28.0, 22.67, 20.18, 17.14, 15.17, 13.729999999999999, 12.59, 11.67, 10.9, 10.24, 9.662, 9.155, 4.598, 2.995, 2.122, 1.564, 0.7648, 0.3422, 0.08872000000000001)
    #         ),
    #         cms.PSet(
    #             charge = cms.vdouble(1.200E3, 1.250E3, 1.300E3, 1.400E3, 1.500E3, 1.600E3, 1.700E3, 1.800E3, 1.900E3, 2.000E3, 4.000E3, 6.000E3, 8.000E3, 10.00E3, 15.00E3, 20.00E3, 25.00E3),
    #             delay = cms.vdouble(28.45, 22.82, 20.41, 17.47, 15.55, 14.12, 13.01, 12.09, 11.31, 10.639999999999999, 5.1160000000000005, 3.317, 2.359, 1.754, 0.8986, 0.45, 0.18209999999999998)
    #         )
    #     )
    #   )
    ),
#Pixel in PS Module
    PSPDigitizerAlgorithm = cms.PSet(
      ElectronPerAdc = cms.double(135.0),
      ReadoutNoiseInElec = cms.double(200.0),#D.B.:Fill readout noise, including all readout chain, relevant for smearing
      ThresholdInElectrons_Barrel = cms.double(6300.), #(0.4 MIP = 0.4 * 16000 e)
      ThresholdInElectrons_Endcap = cms.double(6300.), #(0.4 MIP = 0.4 * 16000 e) 
      AddThresholdSmearing = cms.bool(True),
      ThresholdSmearing_Barrel = cms.double(630.0),
      ThresholdSmearing_Endcap = cms.double(630.0),
      HIPThresholdInElectrons_Barrel = cms.double(1.0e10), # very high value to avoid Over threshold bit
      HIPThresholdInElectrons_Endcap = cms.double(1.0e10), # very high value to avoid Over threshold bit
      NoiseInElectrons = cms.double(200),	         # 30% of the readout noise (should be changed in future)
      Phase2ReadoutMode = cms.int32(0), # Flag to decide Readout Mode :Digital(0) or Analog (linear TDR (-1), dual slope with slope parameters (+1,+2,+3,+4) with threshold subtraction
      AdcFullScale = cms.int32(255),
      TofUpperCut = cms.double(12.5),
      TofLowerCut = cms.double(-12.5),
      AddNoisyPixels = cms.bool(True),
      Alpha2Order = cms.bool(True),			#D.B.: second order effect, does not switch off magnetic field as described
      AddNoise = cms.bool(True),
      AddXTalk = cms.bool(True),			#D.B.
      InterstripCoupling = cms.double(0.05),	#D.B.
      SigmaZero = cms.double(0.00037),  		#D.B.: 3.7um spread for 300um-thick sensor, renormalized in digitizerAlgo
      SigmaCoeff = cms.double(1.80),  		#D.B.: to be confirmed with simulations in CMSSW_6.X
      ClusterWidth = cms.double(3),		#D.B.: this is used as number of sigmas for charge collection (3=+-3sigmas)
      LorentzAngle_DB = cms.bool(False),			
      TanLorentzAnglePerTesla_Endcap = cms.double(0.07),
      TanLorentzAnglePerTesla_Barrel = cms.double(0.07),
      KillModules = cms.bool(False),
      DeadModules_DB = cms.bool(False),
      DeadModules = cms.VPSet(),
      AddInefficiency = cms.bool(False),
      Inefficiency_DB = cms.bool(False),				
      EfficiencyFactors_Barrel = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999 ),
      EfficiencyFactors_Endcap = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 
      0.999, 0.999 ),#Efficiencies kept as Side2Disk1,Side1Disk1 and so on
      CellsToKill = cms.VPSet()
    ),
#Strip in PS module
    PSSDigitizerAlgorithm = cms.PSet(
      ElectronPerAdc = cms.double(135.0),
#D.B.:the noise should be a function of strip capacitance, roughly: ReadoutNoiseInElec=500+(64*Cdet[pF]) ~= 500+(64*1.5[cm])
      ReadoutNoiseInElec = cms.double(700.0),#D.B.:Fill readout noise, including all readout chain, relevant for smearing
      ThresholdInElectrons_Barrel = cms.double(6300.), #(0.4 MIP = 0.4 * 16000 e)
      ThresholdInElectrons_Endcap = cms.double(6300.), #(0.4 MIP = 0.4 * 16000 e)
      AddThresholdSmearing = cms.bool(True),
      ThresholdSmearing_Barrel = cms.double(630.0),
      ThresholdSmearing_Endcap = cms.double(630.0),
      HIPThresholdInElectrons_Barrel = cms.double(21000.), # 1.4 MIP considered as HIP
      HIPThresholdInElectrons_Endcap = cms.double(21000.), # 1.4 MIP considered as HIP 
      NoiseInElectrons = cms.double(700),	         # 30% of the readout noise (should be changed in future)
      Phase2ReadoutMode = cms.int32(0), # Flag to decide Readout Mode :Digital(0) or Analog (linear TDR (-1), dual slope with slope parameters (+1,+2,+3,+4) with threshold subtraction
      AdcFullScale = cms.int32(255),
      TofUpperCut = cms.double(12.5),
      TofLowerCut = cms.double(-12.5),
      AddNoisyPixels = cms.bool(True),
      Alpha2Order = cms.bool(True),			#D.B.: second order effect, does not switch off magnetic field as described
      AddNoise = cms.bool(True),
      AddXTalk = cms.bool(True),			#D.B.
      InterstripCoupling = cms.double(0.05),	#D.B.
      SigmaZero = cms.double(0.00037),  		#D.B.: 3.7um spread for 300um-thick sensor, renormalized in digitizerAlgo
      SigmaCoeff = cms.double(1.80),  		#D.B.: to be confirmed with simulations in CMSSW_6.X
      ClusterWidth = cms.double(3),		#D.B.: this is used as number of sigmas for charge collection (3=+-3sigmas)
      LorentzAngle_DB = cms.bool(False),			
      TanLorentzAnglePerTesla_Endcap = cms.double(0.07),
      TanLorentzAnglePerTesla_Barrel = cms.double(0.07),
      KillModules = cms.bool(False),
      DeadModules_DB = cms.bool(False),
      DeadModules = cms.VPSet(),
      AddInefficiency = cms.bool(False),
      Inefficiency_DB = cms.bool(False),				
      EfficiencyFactors_Barrel = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999 ),
      EfficiencyFactors_Endcap = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 
      0.999, 0.999 ),#Efficiencies kept as Side2Disk1,Side1Disk1 and so on
      CellsToKill = cms.VPSet()
    ),
#Two Strip Module
    SSDigitizerAlgorithm = cms.PSet(
      ElectronPerAdc = cms.double(135.0),
#D.B.:the noise should be a function of strip capacitance, roughly: ReadoutNoiseInElec=500+(64*Cdet[pF]) ~= 500+(64*1.5[cm])
      ReadoutNoiseInElec = cms.double(1000.0),#D.B.:Fill readout noise, including all readout chain, relevant for smearing
      ThresholdInElectrons_Barrel = cms.double(5800.), #D.B.: this should correspond to a threshold of 530mV    
      ThresholdInElectrons_Endcap = cms.double(5800.),
      AddThresholdSmearing = cms.bool(True),
      ThresholdSmearing_Barrel = cms.double(580.0),#D.B.: changed (~5mV peakToPeak --> 1.76mV rms) (was 210.0)
      ThresholdSmearing_Endcap = cms.double(580.0),#D.B.: changed (~5mV peakToPeak --> 1.76mV rms) (was 245.0)
      HIPThresholdInElectrons_Barrel = cms.double(1.0e10), # very high value to avoid Over threshold bit
      HIPThresholdInElectrons_Endcap = cms.double(1.0e10), # very high value to avoid Over threshold bit
      NoiseInElectrons = cms.double(1000),	         # 30% of the readout noise (should be changed in future)
      Phase2ReadoutMode = cms.int32(0), # Flag to decide Readout Mode :Digital(0) or Analog (linear TDR (-1), dual slope with slope parameters (+1,+2,+3,+4) with threshold subtraction
      AdcFullScale = cms.int32(255),
      TofUpperCut = cms.double(12.5),
      TofLowerCut = cms.double(-12.5),
      AddNoisyPixels = cms.bool(True),
      Alpha2Order = cms.bool(True),			#D.B.: second order effect, does not switch off magnetic field as described
      AddNoise = cms.bool(True),
      AddXTalk = cms.bool(True),			#D.B.
      InterstripCoupling = cms.double(0.05),	#D.B.
      SigmaZero = cms.double(0.00037),  		#D.B.: 3.7um spread for 300um-thick sensor, renormalized in digitizerAlgo
      SigmaCoeff = cms.double(1.80),  		#D.B.: to be confirmed with simulations in CMSSW_6.X
      ClusterWidth = cms.double(3),		#D.B.: this is used as number of sigmas for charge collection (3=+-3sigmas)
      LorentzAngle_DB = cms.bool(False),			
      TanLorentzAnglePerTesla_Endcap = cms.double(0.07),
      TanLorentzAnglePerTesla_Barrel = cms.double(0.07),
      KillModules = cms.bool(False),
      DeadModules_DB = cms.bool(False),
      DeadModules = cms.VPSet(),
      AddInefficiency = cms.bool(False),
      Inefficiency_DB = cms.bool(False),				
      EfficiencyFactors_Barrel = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999 ),
      EfficiencyFactors_Endcap = cms.vdouble(0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 
      0.999, 0.999 ),#Efficiencies kept as Side2Disk1,Side1Disk1 and so on
      CellsToKill = cms.VPSet(),
      HitDetectionMode = cms.int32(0),  # (0/1/2/3/4 => SquareWindow/SampledMode/LatchedMode/SampledOrLachedMode/HIPFindingMode)
      PulseShapeParameters = cms.vdouble(-3.0, 16.043703, 99.999857, 40.571650, 2.0, 1.2459094),
        CBCDeadTime = cms.double(0.0) # (2.7 ns deadtime in latched mode)
    )
)

# For premixing stage1
# - add noise as by default
# - do not add noisy pixels (to be done in stage2)
# - do not apply inefficiency (to be done in stage2)
# - disable threshold smearing
#
# For inner pixel
# - extend the dynamic range of ADCs
#
# For outer tracker
# - force analog readout to get the ADCs
#
# NOTE: It is currently assumed that all sub-digitizers have the same ElectronPerAdc.
from Configuration.ProcessModifiers.premix_stage1_cff import premix_stage1
_premixStage1ModifyDict = dict(
    premixStage1 = True,
    PixelDigitizerAlgorithm = dict(
        AddNoisyPixels = False,
        AddInefficiency = False,
        AddThresholdSmearing = False,
        ElectronPerAdc = phase2TrackerDigitizer.PSPDigitizerAlgorithm.ElectronPerAdc.value(),
        AdcFullScale = phase2TrackerDigitizer.PSPDigitizerAlgorithm.AdcFullScale.value(),
    ),
    PSPDigitizerAlgorithm = dict(
        AddNoisyPixels = False,
        AddInefficiency = False,
        AddThresholdSmearing = False,
    ),
    PSSDigitizerAlgorithm = dict(
        AddNoisyPixels = False,
        AddInefficiency = False,
        AddThresholdSmearing = False,
    ),
    SSDigitizerAlgorithm = dict(
        AddNoisyPixels = False,
        AddInefficiency = False,
        AddThresholdSmearing = False,
    ),
)
premix_stage1.toModify(phase2TrackerDigitizer, **_premixStage1ModifyDict)
