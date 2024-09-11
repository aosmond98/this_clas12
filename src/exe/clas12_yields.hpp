
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"

// Helper function to check if a string contains a substring (case-sensitive)
bool contains(const std::string& str, const std::string& substr) {
  return str.find(substr) != std::string::npos;
}

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id, const std::string& output_filename) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  // Determine the beam energy from the environment variable
  float beam_energy = 24.0;
  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Determine the data processing type (gen, rec, or exp) based on output filename
  bool is_gen_data = contains(output_filename, "gen");
  bool is_rec_data = contains(output_filename, "rec");
  bool is_exp_data = contains(output_filename, "exp");

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = is_gen_data || is_rec_data ? std::make_shared<Branches12>(_chain, true) : std::make_shared<Branches12>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  size_t total_twopion_events = 0;

  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);
    
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
      
      // ----- Process Generated Data -----
      if (is_gen_data) {
        // ----- Generated reaction class -----
        auto mc_event = std::make_shared<MCReaction>(data, beam_energy);
        for (int part = 1; part < data->mc_npart(); part++) {
          if (data->mc_pid(part) == PIP) {
            mc_event->SetMCPip(part);
          } else if (data->mc_pid(part) == PROTON) {
            mc_event->SetMCProton(part);
          } else if (data->mc_pid(part) == PIM) {
            mc_event->SetMCPim(part);
          }
        }

          // ----- Generated data output -----
        csv_data output;
        output.event = current_event;
        output.w_mc = mc_event->W_mc();
        output.q2_mc = mc_event->Q2_mc();
        output.weight_gen = mc_event->weight();

        _sync->write(output);
        total++;  // Increment for all events when processing gen data
      }

      // ----- Process Reconstructed Data -----
      else if (is_rec_data) {
        int statusPim = -9999;
        int statusPip = -9999;
        int statusProt = -9999;
        float vertex_hadron[3][3];

        // Make cuts
        if (data->mc_npart() < 1) continue;
        auto dt = std::make_shared<Delta_T>(data);
        auto cuts = std::make_shared<rga_Cuts>(data);
        if (!cuts->ElectronCuts()) continue;
        
        total++;  // Increment only if the event is processed with rec cuts

        // ----- Reconstructed reaction class -----
        auto event = std::make_shared<Reaction>(data, beam_energy);
        for (int part = 1; part < data->gpart(); part++) {
          dt->dt_calc(part);
          if (cuts->IsProton(part)) {
            event->SetProton(part);
            statusProt = abs(data->status(part));
          } else if (cuts->IsPip(part)) {
            event->SetPip(part);
            statusPip = abs(data->status(part));
          } else if (cuts->IsPim(part)) {
            event->SetPim(part);
            statusPim = abs(data->status(part));
          } else {
            event->SetOther(part);
          }
        }
        
        // if (event->TwoPion_missingPim()) {
        // if (event->TwoPion_missingPip()) {
        // if (event->TwoPion_missingProt()) {
        if (event->TwoPion_exclusive()) {
          if (event->W() > 1.25 && event->W() < 2.55 && event->Q2() > 1.5 && event->Q2() < 30.0 && event->weight() > 0.0) {
            total_twopion_events++;

            // ----- Reconstructed data output -----

            csv_data output;

            output.event = current_event;
            output.w = event->W();
            output.q2 = event->Q2();
            output.weight_rec = event->weight();

            _sync->write(output);
          }
        }
      }
    }
    
    std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
    std::cout << " total no of events = " << total << std::endl;
    std::cout << " total no of twopion events = " << total_twopion_events << std::endl;

    return num_of_events;
  }

#endif
