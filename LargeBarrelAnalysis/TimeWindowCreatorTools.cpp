/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file TimeWindowCreatorTools.cpp
 */

#include "TimeWindowCreatorTools.h"
#include "UniversalFileLoader.h"

using namespace std;

/**
 * Sorting method for Signal Channels by time value
 */
 void TimeWindowCreatorTools::sortByValue(vector<JPetSigCh>& input)
 {
   std::sort(input.begin(), input.end(),
     [] (JPetSigCh sigCh1, JPetSigCh sigCh2) {
       return sigCh1.getValue() < sigCh2.getValue();
     }
   );
 }

/**
 * Building all Signal Chnnels from one TDC
 */
vector<JPetSigCh> TimeWindowCreatorTools::buildSigChs(
  TDCChannel* tdcChannel, const JPetTOMBChannel& tombChannel,
  map<unsigned int, vector<double>>& timeCalibrationMap,
  map<unsigned int, vector<double>>& thresholdsMap,
  double maxTime, double minTime, bool setTHRValuesFromChannels,
  JPetStatistics& stats, bool saveHistos
){
  vector<JPetSigCh> allTDCSigChs;
  // Loop over all entries on leading edge in current TOMBChannel and create SigCh
  for (int j = 0; j < tdcChannel->GetLeadHitsNum(); j++) {
    auto leadTime = tdcChannel->GetLeadTime(j);
    if (leadTime > maxTime || leadTime < minTime ) { continue; }
    auto leadSigCh = generateSigCh(
      leadTime, tombChannel, timeCalibrationMap, thresholdsMap,
      JPetSigCh::Leading, setTHRValuesFromChannels
    );
    allTDCSigChs.push_back(leadSigCh);
    if (saveHistos){
      stats.getHisto1D(Form("pm_occupation_thr%d", tombChannel.getLocalChannelNumber()))
        ->Fill(tombChannel.getPM().getID());
    }
  }
  // Loop over all entries on trailing edge in current TOMBChannel and create SigCh
  for (int j = 0; j < tdcChannel->GetTrailHitsNum(); j++) {
    auto trailTime = tdcChannel->GetTrailTime(j);
    if (trailTime > maxTime || trailTime < minTime ) { continue; }
    auto trailSigCh = generateSigCh(
      trailTime, tombChannel, timeCalibrationMap, thresholdsMap,
      JPetSigCh::Trailing, setTHRValuesFromChannels
    );
    allTDCSigChs.push_back(trailSigCh);
    if (saveHistos){
      stats.getHisto1D(Form("pm_occupation_thr%d", tombChannel.getLocalChannelNumber()))
        ->Fill(tombChannel.getPM().getID());
    }
  }
  return allTDCSigChs;
}

/**
 * Method for investigation of repetated Leading edges - setting RecoFlag for each SigCh
 */
void TimeWindowCreatorTools::flagSigChs(
  vector<JPetSigCh>& inputSigChs, JPetStatistics& stats, bool saveHistos
) {
  for(unsigned int i=0; i<inputSigChs.size(); i++) {
    auto& sigCh1 = inputSigChs.at(i);
    if(sigCh1.getType() == JPetSigCh::Trailing) {
      sigCh1.setRecoFlag(JPetSigCh::Good);
      if(saveHistos){
        stats.getHisto1D("good_vs_bad_sigch")->Fill(1);
      }
      continue;
    }
    for(unsigned int j=i+1; j<inputSigChs.size(); j++) {
      auto& sigCh2 = inputSigChs.at(j);
      if(sigCh1.getPM().getID() != sigCh2.getPM().getID()) { continue; }
      if(sigCh1.getThresholdNumber() != sigCh2.getThresholdNumber()) { continue; }
      if(sigCh2.getType() == JPetSigCh::Leading) {
        sigCh1.setRecoFlag(JPetSigCh::Corrupted);
        if(saveHistos){
          stats.getHisto1D("good_vs_bad_sigch")->Fill(2);
          stats.getHisto1D("LL_per_PM")->Fill(sigCh1.getPM().getID());
          stats.getHisto1D("LL_per_THR")->Fill(sigCh1.getThresholdNumber());
          stats.getHisto1D("LL_time_diff")->Fill(sigCh2.getValue()-sigCh1.getValue());
        }
        break;
      } else if(sigCh2.getType()==JPetSigCh::Trailing) {
        sigCh1.setRecoFlag(JPetSigCh::Good);
        if(saveHistos){
          stats.getHisto1D("LT_time_diff")->Fill(sigCh2.getValue()-sigCh1.getValue());
        }
        break;
      }
    }
  }
}

/**
* Sets up Signal Channel fields
*/
JPetSigCh TimeWindowCreatorTools::generateSigCh(
  double tdcChannelTime, const JPetTOMBChannel& channel,
  map<unsigned int, vector<double>>& timeCalibrationMap,
  map<unsigned int, vector<double>>& thresholdsMap,
  JPetSigCh::EdgeType edge, bool setTHRValuesFromChannels
) {
  JPetSigCh sigCh;
  sigCh.setValue(1000.*(tdcChannelTime
    + UniversalFileLoader::getConfigurationParameter(timeCalibrationMap, channel.getChannel())
  ));
  sigCh.setType(edge);
  sigCh.setTOMBChannel(channel);
  sigCh.setPM(channel.getPM());
  sigCh.setFEB(channel.getFEB());
  sigCh.setTRB(channel.getTRB());
  sigCh.setDAQch(channel.getChannel());
  sigCh.setThresholdNumber(channel.getLocalChannelNumber());
  if(setTHRValuesFromChannels) {
    sigCh.setThreshold(channel.getThreshold());
  } else {
    sigCh.setThreshold(
      UniversalFileLoader::getConfigurationParameter(thresholdsMap, channel.getChannel())
    );
  }
  return sigCh;
}
