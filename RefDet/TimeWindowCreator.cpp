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
					   10000, -5, 5
					   ));
  getStatistics().createHistogram(new TH1F("subseq_trails_t_diff",
					   "Time difference between subsequent trailing edges"
					   ";#Delta t [ns]",
					   10000, -5, 5
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
                                  new TH1F("h_strip_no", "Observed Scin", 52, 0.5, 52.5)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("h_wls_pm_no", "Responding WLS PM", 64, 0.5, 64.5)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("h_wls_scin_no_side0", "Responding WLS Scin, side0", 13, 0.5, 13.5)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("h_wls_scin_no_side1", "Responding WLS Scin, side1", 13, 0.5, 13.5)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("h_wls_pm_pos", "Weighted PM position", 1000, 0., 64.)
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
  
  // time differences for strips
  for(int scin=1; scin <= 4*13; scin++){

    if (scin >= 14 && scin <= 26) {
      continue;
    }
    
    getStatistics().createHistogram(new TH1F(Form("dt_scin_%d", scin), "Dt;#Delta t [ns]",
                                             400, -20., 20.
                                             )
                                    );           
  }

  // histos of TOT for each channel separately
  for (int scin = 1; scin <= 4 * 13; scin += 1) {
    for (int pm = 1; pm <= 4; ++pm) {
      getStatistics().createHistogram(
          new TH1F(Form("tot_side_0_scin_%d_pm_%d_thr_A", scin, pm),
                   "TOT;TOT [ns]", 100, 0., 200.));
      getStatistics().createHistogram(
          new TH1F(Form("tot_side_1_scin_%d_pm_%d_thr_A", scin, pm),
                   "TOT;TOT [ns]", 100, 0., 200.));
    }
  }
  
  getStatistics().createHistogram(new TH2F(
                                           "tag_vs_module_dt",
                                           "Dt on tagging scin vs on module scin;"
                                           "#Delta t_{MOD} [ns];#Delta t_{TAG} [ns]",
                                           400, -20., 20., 400, -20., 20.)
                                  );

  getStatistics().createHistogram(new TH2F(
                                           "module_l1_vs_l2_dt",
                                           "Dt on module scin in L1 vs in L2;"
                                           "#Delta t_{L2} [ns];#Delta t_{L1} [ns]",
                                           400, -20., 20., 400, -20., 20.)
                                  );
  
  
  getStatistics().createHistogram(new TH1F("Inter-module Dt", "Inter-module Dt;#Delta t [ns]",
                                           100000, -50000, 50000));

  getStatistics().createHistogram(new TH1F("inter-thr tdiff", "Inter-thredhold Dt;#Delta t [ns]",
                                           30000, -15, 15));

  getStatistics().createHistogram(new TH2F("mult_lead_trail", "TDC signals on trails vs leads; N. lead; N. trail", 9, -0.5, 8.5, 9, -0.5, 8.5));

  getStatistics().createHistogram(new TH2F("2lead_1trail", ";T-L1 [ns];T-L2 [ns]", 1000, -100, 100, 1000, -100, 100));
  getStatistics().createHistogram(new TH2F("2trail_1lead", ";T1-L [ns];T2-L [ns]", 1000, -100, 100, 1000, -100, 100));
  
  
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
    using timesByScin = array<timesByPM,60>;
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
    using hitsByScin = array<hitsByPM,60>;
    using hitsBySide = array<hitsByScin,2>;    

    hitsBySide hits;
    
    for(int side=0;side<=1;++side){
      for(int scin=1;scin<=52;scin+=1){

        for(int pm = 1; pm <= 4; ++pm){
        
	  times thr_a_leads = all_sigchs[side][scin][pm][1][0]; 
	  times thr_b_leads = all_sigchs[side][scin][pm][2][0]; 

	  times thr_a_trails = all_sigchs[side][scin][pm][1][1]; 
	  times thr_b_trails = all_sigchs[side][scin][pm][2][1];

          // fill auxilliary histos to understand the issue
          getStatistics().getHisto2D("mult_lead_trail")->Fill(thr_a_leads.size(), thr_a_trails.size());

          if (thr_a_leads.size() == 2 && thr_a_trails.size() == 1) {
            getStatistics().getHisto2D("2lead_1trail")->Fill(thr_a_trails.front() - thr_a_leads.front(),
                                                             thr_a_trails.front() - thr_a_leads.back());
          }
          if (thr_a_leads.size() == 1 && thr_a_trails.size() == 2) {
            getStatistics().getHisto2D("2trail_1lead")->Fill(thr_a_trails.front() - thr_a_leads.front(),
                                                             thr_a_trails.back() - thr_a_leads.front());
          }

          for (int k=1; k < thr_a_leads.size(); ++k) {
            getStatistics().getHisto1D("subseq_leads_t_diff")->Fill((thr_a_leads[k] - thr_a_leads[k-1])/1000.);
          }
          // artificially fill negative subseq_tdiff to count single-signal cases
          if ( thr_a_leads.size() == 1) {
            getStatistics().getHisto1D("subseq_leads_t_diff")->Fill(-1000.);
          }
          
          for (int k=1; k < thr_a_trails.size(); ++k) {
            getStatistics().getHisto1D("subseq_trails_t_diff")->Fill((thr_a_trails[k] - thr_a_trails[k-1])/1000.);
          }
          // artificially fill negative subseq_tdiff to count single-signal cases
          if ( thr_a_trails.size() == 1) {
            getStatistics().getHisto1D("subseq_trails_t_diff")->Fill(-1000.);
          }

          std::vector<times*> lists = {&thr_a_leads, &thr_b_leads, &thr_a_trails, &thr_b_trails};
          
          std::for_each(lists.begin(), lists.end(),
                        [](std::vector<double>* vec) {
                          if (vec->size() > 1) {
                            vec->erase(vec->begin()+1, vec->end());
                          }
                        });
          
          // check reasonable TOT on A
          std::vector<std::pair<double, double>> sigs_a;
          for(double t_a_l : thr_a_leads){
            for(double t_a_t : thr_a_trails){
              double dt = (t_a_t - t_a_l) / 1000.;
              if( dt > 0. && dt < 200. ){ // reasonable single-threshold TOT
                sigs_a.push_back(std::make_pair(t_a_l / 1000.
                                                , dt));
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
    array<allHitInfo,60> scin_hits;

    for (int scin = 1; scin <= 52; scin += 1) {

      if (scin >= 14 && scin <= 26) {
        continue;
      }

      for (int pm_left = 1; pm_left <= 4; ++pm_left) {
        for (int pm_right = 1; pm_right <= 4; ++pm_right) {
          for (auto &tt_left : hits[0][scin][pm_left]) {
            for (auto &tt_right : hits[1][scin][pm_right]) {

              double dt = tt_right.first - tt_left.first;

              double tot_left = tt_left.second;
              double tot_right = tt_right.second;

              double t = 0.5 * (tt_right.first + tt_left.first); // hit time
              if (fabs(dt) < 20.0) {

                getStatistics().getHisto1D(Form("dt_scin_%d", scin))->Fill(dt);

                scin_hits[scin].push_back({t, tt_left.first, tt_right.first,
                    tt_left.second, tt_right.second});
                
              }
            }
          }
        }
      }
    }
    
    /*********************************************************************/
    /* Look for coincidences between the strips                          */
    /*********************************************************************/

     for (int tag_scin = 1; tag_scin <= 13; ++tag_scin) {
       for (int mod_scin = 27; mod_scin <= 52 ; ++mod_scin) {

         
         for (auto &h1 : scin_hits[mod_scin]) {
           bool coinc_found = false;
           for (auto &h2 : scin_hits[tag_scin]) {

             // soft collimation based on Dt from the tagging scins
             /*
             if (fabs((h2[2] - h2[1]) - fTagDtMeans[tag_scin - 1].first) >
                 0.5 * fTagDtMeans[tag_scin - 1].second) {
               continue;
             }
             */

             if (fabs(h1[0] - h2[0]) < 20.0) {

               // fill Dt plots
               getStatistics().getHisto2D("tag_vs_module_dt")->Fill(h1[2] - h1[1], h2[2] - h2[1]);
               
               // tagger - module coincidenceg in the scintillators!
               coinc_found = true;
               
               // study WLS response coincident with that

               std::vector<std::pair<int, double>> wls_pm_coinc;

               for(int side=0; side<=1; ++side){
                 for (int wls_scin = 14; wls_scin <= 26; ++wls_scin) {
                   for (int pm = 1; pm <= 4; ++pm) {

                     for (auto &wls_sig : hits[side][wls_scin][pm]) {

                       if (fabs(wls_sig.first - h1[0]) < 20.0) {

                         int scin_no = wls_scin - 13;
                         int pm_no = pm;
                         if (side == 0) {
                           pm_no += (scin_no - (scin_no<=9 ? 4 : 5))*4;
                         } else {
                           pm_no += 32 + (13-scin_no)*4;
                         }

                         getStatistics().getHisto1D("h_wls_pm_no")->Fill(pm_no);
                         getStatistics().getHisto1D(Form("h_wls_scin_no_side%d", side))->Fill(scin_no);

                         wls_pm_coinc.push_back({pm_no, wls_sig.second});
                       }
                       
                     }
                   }
                 }
               }

               /*
               bool was_data = false;
               for (auto &wls_pm : wls_pm_coinc) {
                 if(was_data) std::cout << " | ";
                 std::cout << wls_pm.first << " : " << wls_pm.second;
                 was_data = true;
               }
               if (was_data) {
                 std::cout << std::endl;
               }
               */

               // calculate weighted position
               double pos = 0;
               double tot_sum = 0;
               for (auto &wls_pm : wls_pm_coinc) {
                 pos += wls_pm.first * wls_pm.second;
                 tot_sum += wls_pm.second;
               }
               pos /= tot_sum;               

               getStatistics().getHisto1D("h_wls_pm_pos")->Fill(pos);
               
               break;
             }
           }           
           if(coinc_found) break;
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

  // fTagDtMeans =   {
  //   {-2.79883, 2.10773}, {-3.17619, 1.94407}, {-2.26928, 1.35067},
  //   {-3.47174, 1.92962}, {-5.5183, 1.44181},  {-3.56402, 2.07303},
  //   {-3.65391, 1.28607}, {-3.4557, 2.33178},  {-2.10269, 1.43848},
  //   {-3.99234, 1.84888}, {-4.76713, 1.15378}, {-4.27078, 1.74875},
  //   {-4.21927, 1.63704}
  // };
    
  getStatistics()
      .createHistogram(new TH1F("sig_ch_per_time_slot",
                                "Signal Channels Per Time Slot", 50, -0.5,
                                50.5));
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
