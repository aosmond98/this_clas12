#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"

// Helper function to check if a string contains a substring (case-sensitive)
static bool contains(const std::string& str, const std::string& substr) {
  return str.find(substr) != std::string::npos;
}

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram>& _hists, int thread_id, const std::string& output_filename) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  // Determine the beam energy from the environment variable or set default value
  float beam_energy = 10.6;  // Default value
  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));



  // Determine the data processing type (rec or exp) based on output filename
  bool is_rec_data = contains(output_filename, "rec");
  bool is_exp_data = contains(output_filename, "exp");

  // Determine the topology based on output filename
  bool is_topology_excl = contains(output_filename, "excl");
  bool is_topology_mProt = contains(output_filename, "mProt");
  bool is_topology_mPip = contains(output_filename, "mPip");
  bool is_topology_mPim = contains(output_filename, "mPim");


  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = is_rec_data 
                ? std::make_shared<Branches12>(_chain, true)  // For rec
                : std::make_shared<Branches12>(_chain);       // For exp

  // Total number of events "Processed"
  size_t total = 0;
  size_t total_twopion_events = 0;

  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);

    // Print progress for the 0th thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    // Skip events for rec data with no particles
    if (is_rec_data && data->mc_npart() < 1) continue;

    // If we pass electron cuts the event is processed
    total++; // Increment processed events

    // Make cuts
    auto dt = std::make_shared<Delta_T>(data);
    auto cuts = std::make_shared<rga_Cuts>(data);
    if (!cuts->ElectronCuts()) continue;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data, beam_energy, is_rec_data ? "rec" : "exp");

    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);
      _hists->Fill_MomVsBeta(data, part, event);
      _hists->Fill_deltat_pi(data, dt, part, event);
      _hists->Fill_deltat_prot(data, dt, part, event);
      _hists->Fill_MomVsMM2(data, part, event); // test moving these down?

      // Check particle IDs and fill the reaction class
      if (cuts->IsProton(part)) {
        event->SetProton(part);
      } else if (cuts->IsPip(part)) {
        event->SetPip(part);
      } else if (cuts->IsPim(part)) {
        event->SetPim(part);
      } else {
        event->SetOther(part);
      }
    }

    double q2_min_analysis = -1.0, q2_max_analysis = 30.0;
    double w_min_analysis = 1.0, w_max_analysis = 2.5;

    // Dynamically set W and Q2 limits based on BEAM_E
    if (getenv("BEAM_E") != NULL) {
      double beam_energy = atof(getenv("BEAM_E"));
      if (beam_energy < 3) {
        q2_max_analysis = 1.0;
        w_max_analysis = 3.5;
        w_min_analysis = 0.9;
      } else if (beam_energy < 11) {
        q2_min_analysis = 0.0;
        q2_max_analysis = 12.0;
        w_min_analysis = 1.0;
        w_max_analysis = 2.5;
      } else if (beam_energy < 24) {
        q2_min_analysis = 2.0;
        q2_max_analysis = 30.0;
        w_min_analysis = 1.0;
        w_max_analysis = 2.5;
      }
    }

    // Update the condition to use the dynamically set limits
    if ((is_topology_excl && event->TwoPion_exclusive()) ||
        (is_topology_mProt && event->TwoPion_missingProt()) ||
        (is_topology_mPip && event->TwoPion_missingPip()) ||
        (is_topology_mPim && event->TwoPion_missingPim())) {
      if (event->W() > w_min_analysis && event->W() < w_max_analysis && 
          event->Q2() > q2_min_analysis && event->Q2() < q2_max_analysis && 
          event->weight() > 0.0) {
        _hists->Fill_WvsQ2(event, data);
        _hists->Fill_MM2(event, data);
        _hists->Fill_MM2withbins(event);
        total_twopion_events++;
      }
    }
  }

  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  std::cout << " total no of events = " << total << std::endl;
  std::cout << " total no of twopion events = " << total_twopion_events << std::endl;

  return num_of_events;
}

#endif