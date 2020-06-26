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

  getStatistics().createHistogram(new TH2F("n_good_tots", "Signals with good TOT;thr A;thr B",
                                           5, -0.5, 4.5, 5, -0.5, 4.5));
  
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
  for(int scin=1; scin <= 13; scin+=12){                                                
    for(int pm=1; pm<=4; ++pm){
      getStatistics().createHistogram(new TH1F(Form("dt_1pm_a_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               400, -20., 20.
                                               )
                                      );

      getStatistics().createHistogram(new TH1F(Form("dt_totcut_1pm_a_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               400, -20., 20.
                                               )
                                      );
      
      getStatistics().createHistogram(new TH1F(Form("dt_1pm_b_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               400, -20., 20.
                                               )
                                      );
      
      getStatistics().createHistogram(new TH1F(Form("dt_coinc_1pm_a_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               400, -20., 20.
                                               )
                                            );
      
      getStatistics().createHistogram(new TH1F(Form("dt_coinc_1pm_b_scin_%d_pm_%d",
                                                    scin, pm), "Dt;#Delta t [ns]",
                                               400, -20., 20.
                                               )
                                      );
      

      // histos of TOT for each channel separately
      getStatistics().createHistogram(new TH1F(Form("tot_side_0_scin_%d_pm_%d_thr_A",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );
      getStatistics().createHistogram(new TH1F(Form("tot_side_1_scin_%d_pm_%d_thr_A",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );

      getStatistics().createHistogram(new TH1F(Form("tot_side_0_scin_%d_pm_%d_thr_B",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );
      
      getStatistics().createHistogram(new TH1F(Form("tot_side_1_scin_%d_pm_%d_thr_B",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );

      getStatistics().createHistogram(new TH1F(Form("tot_coinc_side_0_scin_%d_pm_%d_thr_A",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );
      getStatistics().createHistogram(new TH1F(Form("tot_coinc_side_1_scin_%d_pm_%d_thr_A",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );

      getStatistics().createHistogram(new TH1F(Form("tot_coinc_side_0_scin_%d_pm_%d_thr_B",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );
      
      getStatistics().createHistogram(new TH1F(Form("tot_coinc_side_1_scin_%d_pm_%d_thr_B",
                                                    scin, pm), "TOT;TOT [ns]",
                                               100, 0., 200.
                                               )
                                      );

      
    }
  }

  getStatistics().createHistogram(new TH1F("Inter-module Dt", "Inter-module Dt;#Delta t [ns]",
                                           100000, -50000, 50000));

  getStatistics().createHistogram(new TH1F("inter-thr tdiff", "Inter-thredhold Dt;#Delta t [ns]",
                                           30000, -15, 15));

  for(int pm=1;pm<5;++pm){
    getStatistics().createHistogram(new TH2F(Form("dt_vs_dt_pm_%d", pm),
                                             Form("Dt scin 1 vs Dt scin 13, pm %d", pm),
                                             400, -20., 20.,
                                             400, -20., 20.
                                             ));
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
    using timesAndTOTs = vector<pair<double,double>>; // (lead time, TOT)
    using hitsByPM = array<timesAndTOTs,5>;
    using hitsByScin = array<hitsByPM,28>;
    using hitsBySide = array<hitsByScin,2>;    

    hitsBySide hits;

    // define TOT cuts
    array<array<array<array<double,2>,5>,28>,2> tot_cuts;
    tot_cuts[0][13][1][0] = 80;
    tot_cuts[0][13][1][1] = 70;
    tot_cuts[0][13][2][0] = 90;
    tot_cuts[0][13][2][1] = 90;
    tot_cuts[0][13][3][0] = 60;
    tot_cuts[0][13][3][1] = 60;
    tot_cuts[0][13][4][0] = 70;
    tot_cuts[0][13][4][1] = 60;

    tot_cuts[1][13][1][0] = 100;
    tot_cuts[1][13][1][1] = 90;
    tot_cuts[1][13][2][0] = 100;
    tot_cuts[1][13][2][1] = 110;
    tot_cuts[1][13][3][0] = 110;
    tot_cuts[1][13][3][1] = 130;
    tot_cuts[1][13][4][0] = 100;
    tot_cuts[1][13][4][1] = 120;

    tot_cuts[0][1][1][0] = 70;
    tot_cuts[0][1][1][1] = 80;
    tot_cuts[0][1][2][0] = 80;
    tot_cuts[0][1][2][1] = 80;
    tot_cuts[0][1][3][0] = 90;
    tot_cuts[0][1][3][1] = 70;
    tot_cuts[0][1][4][0] = 90;
    tot_cuts[0][1][4][1] = 80;

    tot_cuts[1][1][1][0] = 100;
    tot_cuts[1][1][1][1] = 100;
    tot_cuts[1][1][2][0] = 60;
    tot_cuts[1][1][2][1] = 60;
    tot_cuts[1][1][3][0] = 40;
    tot_cuts[1][1][3][1] = 40;
    tot_cuts[1][1][4][0] = 40;
    tot_cuts[1][1][4][1] = 40;

    
    for(int side=0;side<=1;++side){
      for(int scin=1;scin<=13;scin+=12){

        for(int pm = 1; pm <= 4; ++pm){
        
	  times thr_a_leads = all_sigchs[side][scin][pm][1][0]; 
	  times thr_b_leads = all_sigchs[side][scin][pm][2][0]; 

	  times thr_a_trails = all_sigchs[side][scin][pm][1][1]; 
	  times thr_b_trails = all_sigchs[side][scin][pm][2][1]; 

          // check reasonable TOT on A
          std::vector<std::pair<double, double>> sigs_a;
          for(double t_a_l : thr_a_leads){
            for(double t_a_t : thr_a_trails){
              double dt = (t_a_t - t_a_l) / 1000.;
              if( dt > 0. && dt < 200. ){ // reasonable single-threshold TOT
                sigs_a.push_back(std::make_pair(t_a_l / 1000., dt));
              }
              getStatistics().getHisto1D(Form("tot_side_%d_scin_%d_pm_%d_thr_A", side, scin, pm))->Fill(dt);
            }
          }
       
          // check reasonable TOT on B
          std::vector<std::pair<double, double>> sigs_b;
          for(double t_b_l : thr_b_leads){
            for(double t_b_t : thr_b_trails){
              double dt = (t_b_t - t_b_l) / 1000.;
              if( dt > 0. && dt < 200. ){ // reasonable single-threshold TOT
                sigs_b.push_back(std::make_pair(t_b_l / 1000., dt));
              }
              getStatistics().getHisto1D(Form("tot_side_%d_scin_%d_pm_%d_thr_B", side, scin, pm))->Fill(dt);
            }
          }
                    
          getStatistics().getHisto2D("n_good_tots")->Fill(sigs_a.size(), sigs_b.size());

          for(auto& tta : sigs_a){
            for(auto& ttb : sigs_b){
              double dt = tta.first - ttb.first;
              getStatistics().getHisto1D("inter-thr tdiff")->Fill(dt);
              if(fabs(dt) < 4.0){
                hits[side][scin][pm].push_back(tta);
              }
            } 
          }
         
	}
      }
    }

    /**********************************************************************/
    /* Look for coincidences between sides on a single strip              */
    /**********************************************************************/
    using allHitInfo=vector<array<double,5>>; // (t_hit, t_left, t_tight, TOT_left, TOT_right)
    array<allHitInfo,5> scin_1_hits;
    array<allHitInfo,5> scin_13_hits;
    
    for(int scin=1;scin<=13;scin+=12){
      for(int pm = 1; pm <= 4; ++pm){
        for(auto& tt_left: hits[0][scin][pm]){
          for(auto& tt_right: hits[1][scin][pm]){
          
            double dt = tt_right.first - tt_left.first;
            getStatistics().getHisto1D(Form("dt_1pm_a_scin_%d_pm_%d",
                                            scin, pm))->Fill(dt);          

            // fill dt after TOT cut
            double tot_left = tt_left.second;
            double tot_right = tt_right.second;
            if( tot_left > tot_cuts[0][scin][pm][0] / 2. && tot_left < tot_cuts[0][scin][pm][0] &&
                tot_right > tot_cuts[1][scin][pm][0] / 2. && tot_right < tot_cuts[1][scin][pm][0]
                ){
              getStatistics().getHisto1D(Form("dt_totcut_1pm_a_scin_%d_pm_%d",
                                              scin, pm))->Fill(dt);          
            }
            
            double t = 0.5*(tt_right.first + tt_left.first); // hit time

            if( fabs(dt) < 20.0 ){
              if(scin==1){
                scin_1_hits[pm].push_back({t, tt_left.first, tt_right.first, tt_left.second, tt_right.second});
              }
              if(scin==13){
                scin_13_hits[pm].push_back({t, tt_left.first, tt_right.first, tt_left.second, tt_right.second});
              }
            }
            
            
          }
        }
      }
    }

    // fix for missing signals from side B of scin 1 pm 1
    for(auto& tt: scin_1_hits[2]){
      scin_1_hits[1].push_back({tt[0], 0., 0., 0., 0.});
    }
      
    /*********************************************************************/
    /* Look for coincidences between the strips                          */
    /*********************************************************************/

    for(int pm_1 = 1; pm_1 <= 4; ++pm_1){
      for(int pm_13 = 1; pm_13 <= 4; ++pm_13){
        for(auto& tt_1: scin_1_hits[pm_1]){
          for(auto& tt_13: scin_13_hits[pm_13]){

            double dt = tt_1[0] - tt_13[0];

            getStatistics().getHisto1D("Inter-module Dt")->Fill(dt);

            if( fabs(dt) < 8.0 ){ // inter-strip coincidence

              // fill strip time diferences in coincidence
              if(tt_1[4] > fTOTmedianCuts[1][JPetPM::SideB][pm_1] &&
                 tt_1[3] > fTOTmedianCuts[1][JPetPM::SideA][pm_1]
                 ){
                getStatistics().getHisto1D(Form("dt_coinc_1pm_a_scin_%d_pm_%d",
                                                1, pm_1))->Fill(tt_1[2] - tt_1[1]);
              }
              if(tt_13[4] > fTOTmedianCuts[13][JPetPM::SideB][pm_13] &&
                 tt_13[3] > fTOTmedianCuts[13][JPetPM::SideA][pm_13]
                 ){
                getStatistics().getHisto1D(Form("dt_coinc_1pm_a_scin_%d_pm_%d",
                                                13, pm_13))->Fill(tt_13[2] - tt_13[1]);          
              }

              // fill Dt vs Dt
              if(pm_1 == pm_13){
                getStatistics().getHisto2D(Form("dt_vs_dt_pm_%d", pm_1))->Fill(tt_1[2] - tt_1[1], tt_13[2] - tt_13[1]);
              }
              
              // fill TOT in coincidence
              getStatistics().getHisto1D(Form("tot_coinc_side_1_scin_%d_pm_%d_thr_A",
                                              1, pm_1))->Fill(tt_1[4]);
              getStatistics().getHisto1D(Form("tot_coinc_side_0_scin_%d_pm_%d_thr_A",
                                              1, pm_1))->Fill(tt_1[3]);
              getStatistics().getHisto1D(Form("tot_coinc_side_1_scin_%d_pm_%d_thr_A",
                                              13, pm_13))->Fill(tt_13[4]);
              getStatistics().getHisto1D(Form("tot_coinc_side_0_scin_%d_pm_%d_thr_A",
                                              13, pm_13))->Fill(tt_13[3]);
                            
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

  fTOTmedianCuts[1][JPetPM::SideA][1] = 55.0;
  fTOTmedianCuts[1][JPetPM::SideA][2] = 55.0;
  fTOTmedianCuts[1][JPetPM::SideA][3] = 55.0;
  fTOTmedianCuts[1][JPetPM::SideA][4] = 55.0;

  fTOTmedianCuts[13][JPetPM::SideA][1] = 60.0;
  fTOTmedianCuts[13][JPetPM::SideA][2] = 65.0;
  fTOTmedianCuts[13][JPetPM::SideA][3] = 50.0;
  fTOTmedianCuts[13][JPetPM::SideA][4] = 55.0;

  fTOTmedianCuts[1][JPetPM::SideB][1] = 20.0;
  fTOTmedianCuts[1][JPetPM::SideB][2] = 20.0;
  fTOTmedianCuts[1][JPetPM::SideB][3] = 20.0;
  fTOTmedianCuts[1][JPetPM::SideB][4] = 30.0;

  fTOTmedianCuts[13][JPetPM::SideB][1] = 65.0;
  fTOTmedianCuts[13][JPetPM::SideB][2] = 65.0;
  fTOTmedianCuts[13][JPetPM::SideB][3] = 75.0;
  fTOTmedianCuts[13][JPetPM::SideB][4] = 75.0;
  
  
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
