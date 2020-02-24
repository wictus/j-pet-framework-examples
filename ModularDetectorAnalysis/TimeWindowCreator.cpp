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
 *  @file TimeWindowCreator.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include "TimeWindowCreator.h"
#include <EventIII.h>
#include <tuple>
#include <algorithm>

using namespace jpet_options_tools;
using namespace std;

TimeWindowCreator::TimeWindowCreator(const char* name): JPetUserTask(name) {}

TimeWindowCreator::~TimeWindowCreator() {}

bool TimeWindowCreator::init()
{
  INFO("TimeSlot Creation Started");
  fOutputEvents = new JPetTimeWindow("JPetSigCh");
  
  // Reading values from the user options if available
  // Min allowed signal time
  if (isOptionSet(fParams.getOptions(), kMinTimeParamKey)) {
    fMinTime = getOptionAsFloat(fParams.getOptions(), kMinTimeParamKey);
  } else {
    WARNING(
      Form("No value of the %s parameter provided by the user. Using default value of %lf.",
        kMinTimeParamKey.c_str(), fMinTime
      )
    );
  }
  // Max allowed signal time
  if (isOptionSet(fParams.getOptions(), kMaxTimeParamKey)) {
    fMaxTime = getOptionAsFloat(fParams.getOptions(), kMaxTimeParamKey);
  } else {
    WARNING(
      Form("No value of the %s parameter provided by the user. Using default value of %lf.",
        kMaxTimeParamKey.c_str(), fMaxTime
      )
    );
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }

  /************************************************************************/
  /* My histos                                                            */
  /************************************************************************/
  getStatistics().createHistogram(new TH1F("leads_trails_per_thr",
					   "Leads - trails per threshold of pm",
					   9,-4.5, 4.5));
  getStatistics().createHistogram(new TH2F("leads_vs_trails_per_thr",
					   "Leads vs trails per threshold of pm"
					   ";Leads;Trails",
					   5, -0.5, 4.5,
					   5, -0.5, 4.5
					   ));
  getStatistics().createHistogram(new TH1F("subseq_leads_t_diff",
					   "Time difference between subsequent leading edges"
					   ";#Delta t [ns]",
					   1000000, -1000000, 1000000
					   ));
  getStatistics().createHistogram(new TH1F("subseq_trails_t_diff",
					   "Time difference between subsequent leading edges"
					   ";#Delta t [ns]",
					   1000000, -1000000, 1000000
					   ));
  getStatistics().createHistogram(new TH1F("neighbors_pm_t_diff_lead",
					   "Time difference between leads from neighboring PMs"
					   ";#Delta t [ns]",
					   1000000, -45000, 45000
					   ));
  getStatistics().createHistogram(new TH1F("neighbors_pm_t_diff_trail",
					   "Time difference between trails from neighboring PMs"
					   ";#Delta t [ns]",
					   1000000, -45000, 45000
					   ));
  
  getStatistics().createHistogram(new TH1F("tots_thr1",
					   "TOT threshold 1"
					   ";TOT [ns]",
					   10000, 0., 150.
					   ));
  getStatistics().createHistogram(new TH1F("tots_thr2",
					   "TOT threshold 2"
					   ";TOT [ns]",
					   10000, 0., 150.
					   ));

  getStatistics().createHistogram(new TH2F("coinc_channels",
					   "Coinciding channels"
					   ";Channel 1;Channel 2",
					   220, 2090.5, 2310.5,
					   220, 2090.5, 2310.5
					   ));

  for(int i=0;i<13;++i){
    getStatistics().createHistogram(new TH1F(Form("dt_scin_%d", i+1),
					     Form("Time diff scin %d; dt [ns]", i+1),
					     100, -15., 15.
					     ));
  }

  // checking coincidences within a single FTAB
  getStatistics().createHistogram(new TH1F("all_ftab_coinc_a",
					   "Time coincidences in the whole FTAB, side A"
					   ";#Delta t [ns]",
					   10000, 0., 100.
					   ));  

  getStatistics().createHistogram(new TH2F("scin_coinc_ftab_a",
					   "Scin-Scin coincidences, FTAB A"
					   ";2nd scin;1st scin",
					   13, 0.5, 13.5,
					   13, 0.5, 13.5
					   ));

  getStatistics().createHistogram(new TH2F("pm_coinc_ftab_a",
					   "PM-PM (in matrix) coincidences, FTAB A"
					   ";2nd PM pos in matrix;1st pm pos in matrix",
					   4, 0.5, 4.5,
					   4, 0.5, 4.5
					   ));
  
  getStatistics().createHistogram(new TH2F("thr_coinc_ftab_a",
					   "Thr-Thr coincidences, FTAB A"
					   ";2nd hit on threshold;1st hit on threshold",
					   2, 0.5, 2.5,
					   2, 0.5, 2.5
					   ));

  getStatistics().createHistogram(new TH1F("channel_after_masking",
					   "Channels surviving masking"
					   ";#Delta t [ns]",
					   230, 2090.5, 2320.5
					   ));

  // PM-PM coincidences in a single matrix, each possible combination
  for(int side=1; side <= 2; ++side){
    for(int scin=1; scin <= 13; ++scin){
      for(int thr=1; thr<=2; ++thr){
	for(int i=1; i<=4; ++i){
	  for(int j=i+1; j<=4; ++j){

	    getStatistics().createHistogram(new TH1F(Form("dt_PmPm_side_%d_scin_%d_pm_%d_pm_%d_thr%d",
							  side, scin, i, j, thr), ";#Delta t [ns]",
						     200, -10., 10.
						     )
					    );
	  }
	}
      }
    }
  }
  
  // Control histograms
  if (fSaveControlHistos) { initialiseHistograms(); }
  return true;
}

bool TimeWindowCreator::exec()
{
  if (auto event = dynamic_cast<EventIII* const> (fEvent)) {
    int kTDCChannels = event->GetTotalNTDCChannels();
    if (fSaveControlHistos){
      getStatistics().getHisto1D("sig_ch_per_time_slot")->Fill(kTDCChannels);
    }
    // Loop over all TDC channels in file
    auto tdcChannels = event->GetTDCChannelsArray();

    vector<JPetSigCh> mySigChs;

    for (int i = 0; i < kTDCChannels; ++i) {
      auto tdcChannel = dynamic_cast<TDCChannel* const> (tdcChannels->At(i));
      auto channelNumber = tdcChannel->GetChannel();
      
      //      if(channelNumber == 2199 || channelNumber == 2304 || channelNumber == 2297) continue;

      for (int j = 0; j < kTDCChannels; ++j) {
	auto tdcChannel2 = dynamic_cast<TDCChannel* const> (tdcChannels->At(j));
	auto channelNumber2 = tdcChannel2->GetChannel();

	if(channelNumber != channelNumber2)
	  getStatistics().getHisto2D("coinc_channels")->Fill(channelNumber, channelNumber2);
      }

      // Skip trigger signals - every 105th
      //      if (channelNumber % 105 == 0) continue;

      // Check if channel exists in database from loaded local file
      if (getParamBank().getChannels().count(channelNumber) == 0) {
        if (fSaveControlHistos){
          getStatistics().getHisto1D("channel_missing")->Fill(channelNumber);
        }
        WARNING(
          Form("DAQ Channel %d appears in data but does not exist in the detector setup.", channelNumber)
        );
        continue;
      }

      // Get channel for corresponding number
      auto& channel = getParamBank().getChannel(channelNumber);

      vector<JPetSigCh> allSigChs;

      // Loop over all entries on leading edge in current TDCChannel and create SigCh
      for (int j = 0; j < tdcChannel->GetLeadHitsNum(); j++) {
        auto leadTime = tdcChannel->GetLeadTime(j);
        if (leadTime > fMaxTime || leadTime < fMinTime ) {
	  continue; }
		  
        auto leadSigCh = generateSigCh(
				       leadTime, channel, JPetSigCh::Leading
        );
        allSigChs.push_back(leadSigCh);
	mySigChs.push_back(leadSigCh);
        if (fSaveControlHistos){
          fillChannelHistos(channel, JPetSigCh::Leading);
        }
      }

      // Loop over all entries on trailing edge in current TDCChannel and create SigCh
      for (int j = 0; j < tdcChannel->GetTrailHitsNum(); j++) {
        auto trailTime = tdcChannel->GetTrailTime(j);
        if (trailTime > fMaxTime || trailTime < fMinTime ) { continue; }

        auto trailSigCh = generateSigCh(
          trailTime, channel, JPetSigCh::Trailing
        );
        allSigChs.push_back(trailSigCh);
	mySigChs.push_back(trailSigCh);
        if (fSaveControlHistos){
          fillChannelHistos(channel, JPetSigCh::Trailing);
        }
      }
      
      // Save result
      saveSigChs(allSigChs);
    }


    /********************************************************************/
    /* My study                                                         */
    /********************************************************************/
    // time, scin, pm, thr
    vector<tuple<double, int, int, int >> side_a;
    vector<tuple<double, int, int, int >> side_b;

    // scin - vector of <time, type, pm>
    map<int, vector<tuple<double, JPetSigCh::EdgeType, int>>> times_by_scin_a;
    map<int, vector<tuple<double, JPetSigCh::EdgeType, int>>> times_by_scin_b;
	
    for(auto& sc: mySigChs){

      getStatistics().getHisto1D("channel_after_masking")->Fill(sc.getChannel().getID());
      
      int scin = sc.getChannel().getPM().getScin().getID();
      JPetPM::Side side = sc.getChannel().getPM().getSide();
      //      int scin = 207;

      int pm = sc.getChannel().getPM().getMatrixPosition();
      
      int thr = sc.getChannel().getThresholdNumber();
      JPetSigCh::EdgeType type = sc.getType();
      double time = sc.getTime();

      if(side == JPetPM::SideA && type==JPetSigCh::Leading){
	side_a.push_back(tie(time, scin, pm, thr));
      }else{
	side_b.push_back(tie(time, scin, pm, thr));
      }

      if(side == JPetPM::SideA){
	times_by_scin_a[scin].push_back(tie(time, type, pm));
      }else{
	times_by_scin_b[scin].push_back(tie(time, type, pm));
      }
            
      times[scin][side][pm][thr][type].push_back(time);
    }

    /*
    // try clustering of TDC hits
    // scin, pm - time of hit with good TOT
    map<pair<int, int>, vector<double>> hits_by_scin_a;
    map<pair<int, int>, vector<double>> hits_by_scin_b;
    
    for(auto& scin_times : times_by_scin_a){
      sort(scin_times.second.begin(), scin_times.second.end(), [](const tuple<double, JPetSigCh::EdgeType, int>& A, const tuple<double, JPetSigCh::EdgeType, int>& B){
								 return get<0>(A) < get<0>(B);
							       });

      int leads = 0;
      int trails = 0;
      for(auto& hit: scin_times.second){
	if(get<1>(hit) == JPetSigCh::Leading){
	  leads++;
	}else{
	  trails++;
	}
      }
      //      std::cout << "Side A Scin " << scin_times.first << " L: " << leads << " T: " << trails << std::endl;

      // try to match lead-trail
      for(auto& hit_l: scin_times.second){
	if(get<1>(hit_l) == JPetSigCh::Leading){

	  for(auto& hit_t: scin_times.second){
	    if(get<1>(hit_t) == JPetSigCh::Trailing){

	      double dt = (get<0>(hit_t) - get<0>(hit_l))/1000.;
	      if( get<2>(hit_t) == get<2>(hit_l) && dt > 0. && dt < 100. ){
		hits_by_scin_a[make_pair(scin_times.first, get<2>(hit_t))].push_back(get<0>(hit_l));
	      }
	      
	    }
	  }

	}      
      }


    }

    
    for(auto& scin_times : times_by_scin_b){
      sort(scin_times.second.begin(), scin_times.second.end(), [](const tuple<double, JPetSigCh::EdgeType, int>& A, const tuple<double, JPetSigCh::EdgeType, int>& B){
								 return get<0>(A) < get<0>(B);
							       });

      int leads = 0;
      int trails = 0;
      for(auto& hit: scin_times.second){
	if(get<1>(hit) == JPetSigCh::Leading){
	  leads++;
	}else{
	  trails++;
	}
      }
      //      std::cout << "Side B Scin " << scin_times.first << " L: " << leads << " T: " << trails << std::endl;

      // try to match lead-trail
      for(auto& hit_l: scin_times.second){
	if(get<1>(hit_l) == JPetSigCh::Leading){

	  for(auto& hit_t: scin_times.second){
	    if(get<1>(hit_t) == JPetSigCh::Trailing){

	      double dt = (get<0>(hit_t) - get<0>(hit_l))/1000.;
	      if( get<2>(hit_t) == get<2>(hit_l) && dt > 0. && dt < 100. ){
		hits_by_scin_b[make_pair(scin_times.first, get<2>(hit_t))].push_back(get<0>(hit_l));
	      }
	      
	    }
	  }

	}      
      }
      
    }

    // make clusters out of hits
    for(auto& spt: hits_by_scin_a){
      
    }
    */

    
    // look for light sharing or crosstalk
    sort(side_a.begin(), side_a.end(), [](const tuple<double, int, int, int >& A, const tuple<double, int, int, int >& B){
					 return get<0>(A) < get<0>(B);
				       });

    sort(side_b.begin(), side_b.end(), [](const tuple<double, int, int, int >& A, const tuple<double, int, int, int >& B){
					 return get<0>(A) < get<0>(B);
				       });

    double prev_time = 1.e20;
    int prev_scin = -1;
    int prev_pm = -1;
    int prev_thr = -1;

    for(auto& hit: side_a){
      double t = get<0>(hit);
      int scin = get<1>(hit) - 200;
      int pm = get<2>(hit);
      int thr = get<3>(hit);


      //      std::cout << scin << " : " << pm << " : " << thr << std::endl;
      
      getStatistics().getHisto1D("all_ftab_coinc_a")->Fill((t-prev_time)/1000.);
      if ( fabs(t-prev_time) < 1000.0 ){
	// check coincidences
	getStatistics().getHisto2D("scin_coinc_ftab_a")->Fill(scin, prev_scin);
	if(scin==prev_scin){
	  getStatistics().getHisto2D("pm_coinc_ftab_a")->Fill(pm, prev_pm);
	}
	if(scin!=prev_scin || pm != prev_pm){
	  getStatistics().getHisto2D("thr_coinc_ftab_a")->Fill(thr, prev_thr);
	}
      }

      prev_time = t;
      prev_scin = scin;
      prev_pm = pm; 
      prev_thr = thr;
      
    }

    

    
    // fill histos
    for(auto& kscin: times){ // scintillators
      for(auto& kside: kscin.second){// sides
	
	// check time differences of neighboring PMs
	for(auto& kpm1: kside.second){// pms 1
	  for(auto& kpm2: kside.second){// pms 2

	    if( kpm1.first >= kpm2.first ) continue;

	      //	    if(kpm1.first != 3 || kpm2.first !=4)

	    for(int thr=1;thr<=2;++thr){

	      if( kpm1.second[thr][JPetSigCh::Leading].size() == 1 && kpm2.second[thr][JPetSigCh::Leading].size() == 1 ){

		double tdiff = kpm1.second.at(thr)[JPetSigCh::Leading].front() - kpm2.second.at(thr)[JPetSigCh::Leading].front();
		TH1F * hh = getStatistics().getHisto1D(Form("dt_PmPm_side_%d_scin_%d_pm_%d_pm_%d_thr%d",
							  kside.first+1, kscin.first-200, kpm1.first, kpm2.first, thr));

		if(hh) hh->Fill(tdiff/1000.);
	      }
	    }
	    // leading
	    if( kpm1.second[1][JPetSigCh::Leading].size() == 1 && kpm2.second[1][JPetSigCh::Leading].size() == 1 ){
	      double tdiff = kpm1.second[1][JPetSigCh::Leading].front() - kpm2.second[1][JPetSigCh::Leading].front();
	      getStatistics().getHisto1D("neighbors_pm_t_diff_lead")->Fill(tdiff/1000.);
	    }
	    // trailing
	    if( kpm1.second[1][JPetSigCh::Trailing].size() == 1 && kpm2.second[1][JPetSigCh::Trailing].size() == 1 ){
	      double tdiff = kpm1.second[1][JPetSigCh::Trailing].front() - kpm2.second[1][JPetSigCh::Trailing].front();
	      getStatistics().getHisto1D("neighbors_pm_t_diff_trail")->Fill(tdiff/1000.);
	    }
	  }
	}
	
	for(auto& kpm: kside.second){// pms
	  for(auto& kthr: kpm.second){
	    getStatistics().getHisto1D("leads_trails_per_thr")->Fill((int)kthr.second[JPetSigCh::Leading].size() - (int)kthr.second[JPetSigCh::Trailing].size());
	    getStatistics().getHisto2D("leads_vs_trails_per_thr")->Fill(kthr.second[JPetSigCh::Leading].size(), kthr.second[JPetSigCh::Trailing].size());
	    
	    // check times
	    for(int k=1; k < kthr.second[JPetSigCh::Leading].size();++k){
	      double tdiff = kthr.second[JPetSigCh::Leading].at(k) - kthr.second[JPetSigCh::Leading].at(k-1);
	      getStatistics().getHisto1D("subseq_leads_t_diff")->Fill(tdiff/1000.);
	    }

	    for(int k=1; k < kthr.second[JPetSigCh::Trailing].size();++k){
	      double tdiff = kthr.second[JPetSigCh::Trailing].at(k) - kthr.second[JPetSigCh::Trailing].at(k-1);
	      getStatistics().getHisto1D("subseq_trails_t_diff")->Fill(tdiff/1000.);
	    }
	    
	    // check TOTs
	    if(kthr.second[JPetSigCh::Leading].size() == 1 &&
	       kthr.second[JPetSigCh::Trailing].size() == 1){
	      double tot = kthr.second[JPetSigCh::Trailing].front() -
		kthr.second[JPetSigCh::Leading].front();
	      getStatistics().getHisto1D(Form("tots_thr%d", kthr.first))->Fill(tot/1000.);
	    }
	    
	    
	  }
	}
      }	
    }


    // try to test the time resolution
    for(auto& kscin: times){
      auto& sideA = kscin.second[JPetPM::SideA];
      auto& sideB = kscin.second[JPetPM::SideB];
      
      for(int pos=1;pos<=4;++pos){
	auto& times_a = sideA[pos][2][JPetSigCh::Leading];
	auto& times_b = sideB[pos][2][JPetSigCh::Leading];
	
	//	std::cout << "Pos: " << pos << " A: " << times_a.size() << " B: " << times_b.size() << std::endl;
	
	if(times_a.size() == 1 && times_b.size() == 1){
	  //	  std::cout << kscin.first << std::endl;
	  //	  std::cout << (times_a.front() - times_b.front()) / 1000. << std::endl;
	  getStatistics().getHisto1D(Form("dt_scin_%d", kscin.first-200))->Fill((times_a.front() - times_b.front()) / 1000.);
	}
	
      }
      
    }

    
    times.clear();
    /********************************************************************/
    /* End of my study                                                  */
    /********************************************************************/
    

    fCurrEventNumber++;
  } else { return false; }
  return true;
}

/**
* Sets up Signal Channel fields
*/
void TimeWindowCreator::fillChannelHistos(
  const JPetChannel& channel, JPetSigCh::EdgeType edge
) {
  getStatistics().getHisto1D("channel_occ")->Fill(channel.getID());
  getStatistics().getHisto1D("channel_thrnum")->Fill(channel.getThresholdNumber());
  getStatistics().getHisto1D("matrix_occ")->Fill(channel.getPM().getMatrixPosition());
  getStatistics().getHisto1D("pm_occ")->Fill(channel.getPM().getID());
  getStatistics().getHisto1D("scin_occ")->Fill(channel.getPM().getScin().getID());
  getStatistics().getHisto1D("slot_occ")->Fill(channel.getPM().getScin().getID());

  if(edge == JPetSigCh::Leading){
    getStatistics().getHisto1D("channel_occ_leads")->Fill(channel.getID());
    getStatistics().getHisto1D("channel_thrnum_leads")->Fill(channel.getThresholdNumber());
    getStatistics().getHisto1D("matrix_occ_leads")->Fill(channel.getPM().getMatrixPosition());
    getStatistics().getHisto1D("pm_occ_leads")->Fill(channel.getPM().getID());
    getStatistics().getHisto1D("scin_occ_leads")->Fill(channel.getPM().getScin().getID());
    getStatistics().getHisto1D("slot_occ_leads")->Fill(channel.getPM().getScin().getID());
  } else if(edge == JPetSigCh::Trailing){
    getStatistics().getHisto1D("channel_occ_trails")->Fill(channel.getID());
    getStatistics().getHisto1D("channel_thrnum_trails")->Fill(channel.getThresholdNumber());
    getStatistics().getHisto1D("matrix_occ_trails")->Fill(channel.getPM().getMatrixPosition());
    getStatistics().getHisto1D("pm_occ_trails")->Fill(channel.getPM().getID());
    getStatistics().getHisto1D("scin_occ_trails")->Fill(channel.getPM().getScin().getID());
    getStatistics().getHisto1D("slot_occ_trails")->Fill(channel.getPM().getScin().getID());
  }

  if(channel.getPM().getSide() == JPetPM::SideA){
    getStatistics().getHisto1D("pm_occ_sides")->Fill(1);
  }else if(channel.getPM().getSide() == JPetPM::SideB){
    getStatistics().getHisto1D("pm_occ_sides")->Fill(2);
  }
}

/**
* Sets up Signal Channel fields
*/
JPetSigCh TimeWindowCreator::generateSigCh(
  double time, const JPetChannel& channel, JPetSigCh::EdgeType edge
) {
  JPetSigCh sigCh;
  sigCh.setTime(1000.*time);
  sigCh.setType(edge);
  sigCh.setChannel(channel);
  sigCh.setRecoFlag(JPetSigCh::Good);
  return sigCh;
}

bool TimeWindowCreator::terminate()
{
  INFO("TimeSlot Creation Ended");
  return true;
}

void TimeWindowCreator::saveSigChs(const vector<JPetSigCh>& sigChVec)
{
  //  for (auto & sigCh : sigChVec) { fOutputEvents->add<JPetSigCh>(sigCh); }
}

void TimeWindowCreator::initialiseHistograms(){

  getStatistics().createHistogram(
    new TH1F("sig_ch_per_time_slot", "Signal Channels Per Time Slot", 50, -0.5, 50.5)
  );
  getStatistics().getHisto1D("sig_ch_per_time_slot")
  ->GetXaxis()->SetTitle("Signal Channels in Time Slot");
  getStatistics().getHisto1D("sig_ch_per_time_slot")
  ->GetYaxis()->SetTitle("Number of Time Slots");

  // Channels
  getStatistics().createHistogram(
    new TH1F("channel_occ", "Channels occupation", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_occ")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_occ")->GetYaxis()->SetTitle("Number of SigCh");

  // 2204 2309
  // high occ: 2198, 2199, 2204, 2205, 2296, 2303, 2304, 2305, 2308, 2309
  getStatistics().createHistogram(
    new TH1F("channel_missing", "Channels missing in configuration", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_missing")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_missing")->GetYaxis()->SetTitle("Counts");

  getStatistics().createHistogram(
    new TH1F("channel_thrnum", "Channels threshold numbers", 4, 0.5, 4.5)
  );
  getStatistics().getHisto1D("channel_thrnum")->GetXaxis()->SetTitle("Channel Threshold Number");
  getStatistics().getHisto1D("channel_thrnum")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_thrnum_leads", "Channels threshold numbers LEADS", 4, 0.5, 4.5)
  );
  getStatistics().getHisto1D("channel_thrnum_leads")->GetXaxis()->SetTitle("Channel Threshold Number");
  getStatistics().getHisto1D("channel_thrnum_leads")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_thrnum_trails", "Channels threshold numbers TRAILS", 4, 0.5, 4.5)
  );
  getStatistics().getHisto1D("channel_thrnum_trails")->GetXaxis()->SetTitle("Channel Threshold Number");
  getStatistics().getHisto1D("channel_thrnum_trails")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_occ_leads", "Channels occupation - Leading channels", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_occ_leads")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_occ_leads")->GetYaxis()->SetTitle("Number of Lading SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_occ_trails", "Channels occupation - Trailing channels", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_occ_trails")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");

  // SiPMs
  getStatistics().createHistogram(
    new TH1F("pm_occ", "PMs occupation", 111, 399.5, 510.5)
  );
  getStatistics().getHisto1D("pm_occ")->GetXaxis()->SetTitle("PM ID");
  getStatistics().getHisto1D("pm_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("pm_occ_leads", "PMs occupation LEADS", 111, 399.5, 510.5)
  );
  getStatistics().getHisto1D("pm_occ_leads")->GetXaxis()->SetTitle("PM ID");
  getStatistics().getHisto1D("pm_occ_leads")->GetYaxis()->SetTitle("Number of Leading SigCh");

  getStatistics().createHistogram(
    new TH1F("pm_occ_trails", "PMs occupation TRAILS", 111, 399.5, 510.5)
  );
  getStatistics().getHisto1D("pm_occ_trails")->GetXaxis()->SetTitle("PM ID");
  getStatistics().getHisto1D("pm_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");

  getStatistics().createHistogram(
    new TH1F("pm_occ_sides", "PMs occupation of sides", 2, 0.5, 2.5)
  );
  getStatistics().getHisto1D("pm_occ_sides")->GetXaxis()->SetBinLabel(1, "SIDE A");
  getStatistics().getHisto1D("pm_occ_sides")->GetXaxis()->SetBinLabel(2, "SIDE B");
  getStatistics().getHisto1D("pm_occ_sides")->GetYaxis()->SetTitle("Number of SigCh");

  // Matrix position
  getStatistics().createHistogram(
    new TH1F("matrix_occ", "Position in matrix in PMs occupation", 5, -0.5, 4.5)
  );
  getStatistics().getHisto1D("matrix_occ")->GetXaxis()->SetTitle("Matrix position");
  getStatistics().getHisto1D("matrix_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("matrix_occ_leads", "Position in matrix in PMs occupation LEADS", 5, -0.5, 4.5)
  );
  getStatistics().getHisto1D("matrix_occ_leads")->GetXaxis()->SetTitle("Matrix position");
  getStatistics().getHisto1D("matrix_occ_leads")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("matrix_occ_trails", "Position in matrix in PMs occupation TRAILS", 5, -0.5, 4.5)
  );
  getStatistics().getHisto1D("matrix_occ_trails")->GetXaxis()->SetTitle("Matrix position");
  getStatistics().getHisto1D("matrix_occ_trails")->GetYaxis()->SetTitle("Number of SigCh");

  // Scins
  getStatistics().createHistogram(
    new TH1F("scin_occ", "Scins occupation", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("scin_occ")->GetXaxis()->SetTitle("SCIN ID");
  getStatistics().getHisto1D("scin_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("scin_occ_leads", "Scins occupation LEADS",16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("scin_occ_leads")->GetXaxis()->SetTitle("SCIN ID");
  getStatistics().getHisto1D("scin_occ_leads")->GetYaxis()->SetTitle("Number of Leading SigCh");

  getStatistics().createHistogram(
    new TH1F("scin_occ_trails", "Scins occupation TRAILS", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("scin_occ_trails")->GetXaxis()->SetTitle("SCIN ID");
  getStatistics().getHisto1D("scin_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");

  // Slots
  getStatistics().createHistogram(
    new TH1F("slot_occ", "Slots occupation", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("slot_occ")->GetXaxis()->SetTitle("SLOT ID");
  getStatistics().getHisto1D("slot_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("slot_occ_leads", "Slots occupation LEADS", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("slot_occ_leads")->GetXaxis()->SetTitle("SLOT ID");
  getStatistics().getHisto1D("slot_occ_leads")->GetYaxis()->SetTitle("Number of Leading SigCh");

  getStatistics().createHistogram(
    new TH1F("slot_occ_trails", "Slots occupation TRAILS", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("slot_occ_trails")->GetXaxis()->SetTitle("SLOT ID");
  getStatistics().getHisto1D("slot_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");
}
