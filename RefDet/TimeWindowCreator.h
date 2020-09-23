/**
 *  @copyright Copyright 2019 The J-PET Framework Authors. All rights reserved.
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
 *  @file TimeWindowCreator.h
 */

#ifndef TIMEWINDOWCREATOR_H
#define TIMEWINDOWCREATOR_H

#include <JPetTimeWindow/JPetTimeWindow.h>
#include <JPetUserTask/JPetUserTask.h>
#include <JPetChannel/JPetChannel.h>
#include <JPetSigCh/JPetSigCh.h>
#include <map>
#include <set>

class JPetWriter;

#ifdef __CINT__
#define override
#endif

/**
 * @brief User Task: translate Unpacker EventIII data to JPetTimeWindow
 *
 * Task translates data from Unpacker file fomrat - EventIII to JPetTimeWindow.
 * Parameters for start and end time can be specified in user options, default are
 * provided. Moreover time calibration and threshold values injection can be
 * performed, if ASCII files of standard format were provided. In case of errors,
 * creation of Time Windows continues without this additional information.
 */
class TimeWindowCreator: public JPetUserTask
{
public:
	TimeWindowCreator(const char* name);
	virtual ~TimeWindowCreator();
	virtual bool init() override;
	virtual bool exec() override;
	virtual bool terminate() override;

protected:
	void saveSigChs(const std::vector<JPetSigCh>& sigChVec);
	void initialiseHistograms();
	void fillChannelHistos(
	  const JPetChannel& channel, JPetSigCh::EdgeType edge
	);
	const std::string kSaveControlHistosParamKey = "Save_Control_Histograms_bool";
	const std::string kMaxTimeParamKey = "TimeWindowCreator_MaxTime_float";
	const std::string kMinTimeParamKey = "TimeWindowCreator_MinTime_float";
	long long int fCurrEventNumber = 0;
	bool fSaveControlHistos = true;
	double fMinTime = -1.e6;
	double fMaxTime = 0.;
	static JPetSigCh generateSigCh(
		double time, const JPetChannel& channel, JPetSigCh::EdgeType edge
	);
 	bool readMappingForWLS();


  //      scin           side                  pos_in_matrix  thr     
  std::map<int,
	   std::map<JPetPM::Side,
		    std::map<int,
			     std::map<int,
				      std::map<JPetSigCh::EdgeType,
					       std::vector<float>> > > > > times;

  std::vector<std::pair<double, double>> fTagDtMeans;

  std::map<int, int> fTDCtoWLSIDside0;
  std::map<int, int> fTDCtoWLSIDside1;
  std::vector< std::map< int, int > > fMappingForWLS; 
};

#endif /* !TIMEWINDOWCREATOR_H */
