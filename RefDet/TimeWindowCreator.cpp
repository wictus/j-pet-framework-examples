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
#include <array>
#include <cmath>

using namespace jpet_options_tools;
using namespace std;

TimeWindowCreator::TimeWindowCreator(const char* name): JPetUserTask(name) {}

TimeWindowCreator::~TimeWindowCreator() {}

double fixTime(double dt){
  int offset = round(dt / 2.7027);
  return dt - offset * 2.7027;
}

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
					   10000, -1000, 1000
					   ));
  getStatistics().createHistogram(new TH1F("subseq_trails_t_diff",
					   "Time difference between subsequent leading edges"
					   ";#Delta t [ns]",
					   1000000, -1000000, 1000000
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

  getStatistics().createHistogram(new TH1F("ref_tot_coin",
					   "TOT on RefDet in coinc."
					   ";TOT [ns]",
					   10000, 0., 150.
					   ));

  getStatistics().createHistogram(new TH1F("ref_tot",
					   "TOT on RefDet"
					   ";TOT [ns]",
					   10000, 0., 150.
					   ));

  getStatistics().createHistogram(new TH1F("ref_n_hits",
					   "No ref. det. hits in time window",
                                           10, -0.5, 9.5
					   ));

  getStatistics().createHistogram(
                                  new TH1F("h_pmt_no", "Observed PM", 4, 0.5, 4.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_strip_no", "Observed Scin", 13, 0.5, 13.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_thr_tdiff", "Time difference between A and B thresholds;#Delta t_{AB} [ns]", 1000, -10., 10.)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("h_leads_size_a", "Number of thr A leading edges in TW", 10, -0.5, 9.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_trails_size_a", "Number of thr A trailing edges in TW", 10, -0.5, 9.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_leads_size_b", "Number of thr B leading edges in TW", 10, -0.5, 9.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_trails_size_b", "Number of thr B trailing edges in TW", 10, -0.5, 9.5)
                                  );
  

  getStatistics().createHistogram(
                                  new TH1F("h_tot_a", "TOT thr A;TOT [ns]", 4000, -200., 200.)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_tot_b", "TOT thr B;TOT [ns]", 4000, -200., 200.)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("h_lead_trail_a", "Leads - trails, thr A;N_{LEAD} - N_{TRAIL}", 5, -2.5, 2.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_lead_trail_b", "Leads - trails, thr B;N_{LEAD} - N_{TRAIL}", 5, -2.5, 2.5)
                                  );


  getStatistics().createHistogram(new TH1F("cluster_size", "Cluster size (no. PMs)", 5, -0.5, 4.5));
  
  // time differences for different time calculation variants
  for(int scin=1; scin <= 13; ++scin){
                                                 
    for(int pm=1; pm<=4; ++pm){
      getStatistics().createHistogram(new TH1F(Form("dt_1pm_a_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               200, -20., 20.
                                               )
                                      );
      
      getStatistics().createHistogram(new TH1F(Form("dt_1pm_b_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               200, -20., 20.
                                               )
                                      );
      
      getStatistics().createHistogram(new TH1F(Form("dt_tag_1pm_a_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               200, -20., 20.
                                               )
                                            );
      
      getStatistics().createHistogram(new TH1F(Form("dt_tag_1pm_b_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               200, -20., 20.
                                               )
                                      );
      

      // histos of TOT for each channel separately
      getStatistics().createHistogram(new TH1F(Form("tot_scin_%d_pm_%d_thr_A",
                                                    scin, pm), "TOT;TOT [ns]",
                                               260, 0., 130.
                                               )
                                      );

      getStatistics().createHistogram(new TH1F(Form("tot_scin_%d_pm_%d_thr_B",
                                                    scin, pm), "TOT;TOT [ns]",
                                               260, 0., 130.
                                               )
                                      );

      getStatistics().createHistogram(new TH1F(Form("tot_coinc_scin_%d_pm_%d_thr_A",
                                                    scin, pm), "TOT;TOT [ns]",
                                               260, 0., 130.
                                               )
                                      );

      getStatistics().createHistogram(new TH1F(Form("tot_coinc_scin_%d_pm_%d_thr_B",
                                                    scin, pm), "TOT;TOT [ns]",
                                               260, 0., 130.
                                               )
                                      );
      
    }
  }

  getStatistics().createHistogram(new TH1F("Inter-module Dt", "Inter-module Dt;#Delta t [ns]",
                                           100000, -50000, 50000));
  
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
      
      //      Skip trigger signals - every 105th
      if (channelNumber % 105 == 104) continue;

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

      // Loop over all entries on leading edge in current TDCChannel and create SigCh
      for (int j = 0; j < tdcChannel->GetLeadHitsNum(); j++) {
        auto leadTime = tdcChannel->GetLeadTime(j);
        if (leadTime > fMaxTime || leadTime < fMinTime ) {
	  continue; }
		  
        auto leadSigCh = generateSigCh(
				       leadTime, channel, JPetSigCh::Leading
        );

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

	mySigChs.push_back(trailSigCh);
        if (fSaveControlHistos){
          fillChannelHistos(channel, JPetSigCh::Trailing);
        }
      }
      
    }


    /********************************************************************/
    /* My study                                                         */
    /********************************************************************/
    //
    using times = vector<double>;
    using timesByEdge = array<times, 2>;
    using timesByThreshold = array<timesByEdge, 3>;
    using timesByPM = array<timesByThreshold,5>;
    using timesByScin = array<timesByPM,28>;
    using timesBySide = array<timesByScin,2>;

    timesBySide all_sigchs;
    	
    for(auto& sc: mySigChs){
      
      int scin = sc.getChannel().getPM().getScin().getID() - 200;
      JPetPM::Side side = sc.getChannel().getPM().getSide();

      int pm = sc.getChannel().getPM().getMatrixPosition();

      int thr = sc.getChannel().getThresholdNumber();
      JPetSigCh::EdgeType type = sc.getType();
      double time = sc.getTime();

      getStatistics().getHisto1D("h_pmt_no")->Fill(pm);
      getStatistics().getHisto1D("h_strip_no")->Fill(scin);

      int side_no = 0; // 0 - Side A, 1 - Side B
      if(side == JPetPM::SideA){
        side_no = 0;
      }else{
        side_no = 1;
      }

      int edge = 0; // 0 - Leading, 1 - Trailing
      if(type==JPetSigCh::Leading){
        edge = 0;
      }else{
        edge = 1;
      }

      all_sigchs[side_no][scin][pm][thr][edge].push_back( time );  
    }

    /**************************************************************************/
    /* Actual clustering                                                      */
    /**************************************************************************/
    using hitsByPM = array<times,5>;
    using hitsByScin = array<hitsByPM,28>;
    using hitsBySide = array<hitsByScin,2>;    

    hitsBySide hits;
    
    for(int side=0;side<=1;++side){
      for(int scin=1;scin<=13;++scin){

        for(int pm = 1; pm <= 4; ++pm){
        
	  times thr_a_leads = all_sigchs[side][scin][pm][1][0]; 
	  times thr_b_leads = all_sigchs[side][scin][pm][2][0]; 

	  times thr_a_trails = all_sigchs[side][scin][pm][1][1]; 
	  times thr_b_trails = all_sigchs[side][scin][pm][2][1]; 



          if(scin <= 13){
            if( thr_a_leads.size() != 1 ||
                thr_b_leads.size() != 1 ||
                thr_a_trails.size() != 1 ||
                thr_b_trails.size() != 1
                ){
              continue;
            }
          }
          
          double tot_a = thr_a_trails.front() - thr_a_leads.front();
          double tot_b = thr_b_trails.front() - thr_b_leads.front();

          if(scin <= 13){
            getStatistics().getHisto1D(Form("tot_scin_%d_pm_%d_thr_A", scin, pm))->Fill(tot_a/1000.);
            getStatistics().getHisto1D(Form("tot_scin_%d_pm_%d_thr_B", scin, pm))->Fill(tot_b/1000.);
          }

          // time coincidence within 6 ns
          if( fabs(thr_a_leads.front() - thr_b_leads.front()) < 6000. ){
            
            hits[side][scin][pm].push_back(thr_a_leads.front() / 1000. );
            hits[side][scin][pm].push_back(thr_b_leads.front() / 1000.);
            
          }

	}
      }
    }
 
    bool was_tag = false;
    std::vector<std::pair<double, double>> tag_times;
    
    // look for tags on the reference detector
    times thr_a_leads = all_sigchs[0][20][1][1][0]; 
    times thr_a_trails = all_sigchs[0][20][1][1][1]; 
    times thr_b_leads = all_sigchs[0][20][1][2][0]; 
    getStatistics().getHisto1D("ref_n_hits")->Fill(thr_a_leads.size());

    if( thr_a_leads.size() == 0 ){ // no chance for coincidence
      return true;
    }

    for(int k=0;k<thr_a_leads.size() && k<thr_b_leads.size(); ++k){
      double tot = 1.e6;
      if(thr_a_leads.size() == thr_a_trails.size()){
        tot = (thr_a_trails.at(k) - thr_a_leads.at(k)) / 1000.;
      }
      getStatistics().getHisto1D("ref_tot")->Fill(tot);
      tag_times.push_back(std::make_pair(thr_a_leads.at(k) / 1000., tot));
    }
        
    // fill single-channel time diffs in case there was a tag on the "far" side
    for(int scin=1;scin<=13;++scin){
      for(int pm = 1; pm <= 4; ++pm){
        
        if(hits[0][scin][pm].size() != 2 || hits[1][scin][pm].size() != 2){
          continue;
        }

        double t_a = 0.5*(hits[0][scin][pm].at(0) + hits[1][scin][pm].at(0));
        double dt = hits[0][scin][pm].at(0) - hits[1][scin][pm].at(0);       
        
        getStatistics().getHisto1D(Form("dt_1pm_a_scin_%d_pm_%d",
                                        scin, pm))->Fill(dt);
        
        for(auto& tt: tag_times){

          double t = tt.first;
          double tot = tt.second;
          
          getStatistics().getHisto1D("Inter-module Dt")->Fill(t_a - t);

          if( fabs(t_a - t + 30.0) < 10.0 ){ // we had a coincidence with a tag
            
            // fill Dt in coincidence            

            getStatistics().getHisto1D(Form("dt_tag_1pm_a_scin_%d_pm_%d",
                                          scin, pm))->Fill(dt);        
            dt = hits[0][scin][pm].at(1) - hits[1][scin][pm].at(1);
            getStatistics().getHisto1D(Form("dt_tag_1pm_b_scin_%d_pm_%d",
                                          scin, pm))->Fill(dt);;

            getStatistics().getHisto1D("ref_tot_coin")->Fill(tot);
            
            // fill TOT in coincidence
            for(int side=0;side<=1;++side){
              times thr_a_leads = all_sigchs[side][scin][pm][1][0]; 
              times thr_b_leads = all_sigchs[side][scin][pm][2][0]; 
              
              times thr_a_trails = all_sigchs[side][scin][pm][1][1]; 
              times thr_b_trails = all_sigchs[side][scin][pm][2][1]; 
            
              if( !(thr_a_leads.size() != 1 ||
                    thr_b_leads.size() != 1 ||
                    thr_a_trails.size() != 1 ||
                    thr_b_trails.size() != 1) ){
                double tot_a = thr_a_trails.front() - thr_a_leads.front();
                double tot_b = thr_b_trails.front() - thr_b_leads.front();              
              
                getStatistics().getHisto1D(Form("tot_coinc_scin_%d_pm_%d_thr_A", scin, pm))->Fill(tot_a/1000.);
                getStatistics().getHisto1D(Form("tot_coinc_scin_%d_pm_%d_thr_B", scin, pm))->Fill(tot_b/1000.);
              }
            }     
          }
        }
      }
      
      
    }
    
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
