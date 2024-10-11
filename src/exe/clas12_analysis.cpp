#include <future>
#include <thread>
#include "clas12_analysis.hpp"

int main(int argc, char** argv) {
  // Need this to make sure root doesn't break
  ROOT::EnableThreadSafety();

  // Set the number of threads
  int NUM_THREADS = 4;
  if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
  if (NUM_THREADS > argc - 2) NUM_THREADS = 1;  // Ensure thread count doesn't exceed input files

  // Parse beam energy from environment variable (default to 10.6)
  float beam_energy = 10.6;
  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Make a vector of vectors of strings the size of the number of threads
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);

  // Get the output file name
  std::string outfilename;

  // Check if there are enough arguments
  if (argc >= 2) {
    // First argument is the output file
    outfilename = argv[1];
    // All other files are split evenly by the number of threads
    for (int i = 2; i < argc; i++) {
      infilenames[i % NUM_THREADS].push_back(argv[i]);
    }
  } else {
    std::cerr << "Not enough arguments! Usage: ./program outputfile inputfiles...\n";
    return 1;
  }

  // Make a set of threads (Futures are special threads which return a value)
  std::future<size_t> threads[NUM_THREADS];

  // Define events to be used to get Hz later
  size_t events = 0;

  // Make your histograms object as a shared pointer that all the threads will share
  auto hists = std::make_shared<Histogram>(outfilename);

  // Lambda function to run files on each thread
  auto run_files = [&hists, &outfilename](std::vector<std::string> inputs, int thread_id) {
    // Make a new chain to process for this thread
    auto chain = std::make_shared<TChain>("clas12");
    // Add every file to the chain
    for (auto in : inputs) chain->Add(in.c_str());

    // Run the function over each thread, passing the output filename to determine rec/exp
    return run<Cuts>(std::move(chain), hists, thread_id, outfilename);
  };

  // Start timer
  auto start = std::chrono::high_resolution_clock::now();

  // For each thread, run the function asynchronously
  for (size_t i = 0; i < NUM_THREADS; i++) {
    threads[i] = std::async(run_files, infilenames.at(i), i);
  }

  // Gather results from each thread
  for (size_t i = 0; i < NUM_THREADS; i++) {
    events += threads[i].get();  // Get the number of events processed by this thread
  }

  // Timer and Hz calculator functions that print at the end
  std::cout.imbue(std::locale(""));  // Puts commas in the output
  std::chrono::duration<double> elapsed_full = std::chrono::high_resolution_clock::now() - start;
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;

  return 0;
}




// #include <future>
// #include <thread>
// #include "clas12_analysis.hpp"

// int main(int argc, char** argv) {
//   // Need this to make sure root doesn't break
//   ROOT::EnableThreadSafety();

//   int NUM_THREADS = 4;
//   if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
//   if (NUM_THREADS > argc - NUM_THREADS) NUM_THREADS = 1;

//   // Make a vector of vectors of strings the size of the number of threads
//   std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
//   // Get the output file name
//   std::string outfilename;

//   if (argc >= 2) {
//     // First argument is the output file
//     outfilename = argv[1];
//     // All other files are split evently by the under of threads
//     for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
//   } else {
//     return 1;
//   }

//   // Make a set of threads (Futures are special threads which return a value)
//   std::future<size_t> threads[NUM_THREADS];

//   // Define events to be used to get Hz later
//   size_t events = 0;

//   // Make your histograms object as a shared pointer that all the threads will have
//   auto hists = std::make_shared<Histogram>(outfilename);

//   auto run_files = [&hists](std::vector<std::string> inputs, int thread_id) {
//     // Called once for each thread
//     // Make a new chain to process for this thread
//     auto chain = std::make_shared<TChain>("clas12");
//     // Add every file to the chain
//     for (auto in : inputs) chain->Add(in.c_str());
//     // Run the function over each thread
//     return run<Cuts>(std::move(chain), hists, thread_id);
//   };

//   // Start timer
//   auto start = std::chrono::high_resolution_clock::now();
//   // For each thread
//   for (size_t i = 0; i < NUM_THREADS; i++) {
//     // Set the thread to run a task A-Syncroisly
//     // The function we run is the first argument (run_files)
//     // The functions areruments are all the remaining arguments
//     threads[i] = std::async(run_files, infilenames.at(i), i);
//   }

//   // For each thread
//   for (size_t i = 0; i < NUM_THREADS; i++) {
//     // Get the information from the thread in this case how many events each thread actually computed
//     events += threads[i].get();
//   }

//   // Timer and Hz calculator functions that print at the end
//   std::cout.imbue(std::locale(""));  // Puts commas in
//   std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
//   std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
//   std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
//   return 0;
// }

// // #include <future>
// // #include <thread>
// // #include "clas12_analysis.hpp"

// // int main(int argc, char** argv) {
// //   // Need this to make sure root doesn't break
// //   ROOT::EnableThreadSafety();

// //   int NUM_THREADS = 4;
// //   if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
// //   if (NUM_THREADS > argc - NUM_THREADS) NUM_THREADS = 1;

// //   // Parse the beam energy from the command line
// //   float beam_energy = 24.0;  // Default value
// //   if (argc > 1) {
// //     beam_energy = std::stof(argv[1]);
// //   }

// //   // Make a vector of vectors of strings the size of the number of threads
// //   std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
// //   // Get the output file name
// //   std::string outfilename;

// //   if (argc >= 3) {
// //     // Second argument is the output file
// //     outfilename = argv[2];
// //     // All other files are split evenly by the number of threads
// //     for (int i = 3; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
// //   } else {
// //     return 1;
// //   }

// //   // Make a set of threads (Futures are special threads which return a value)
// //   std::future<size_t> threads[NUM_THREADS];

// //   // Define events to be used to get Hz later
// //   size_t events = 0;

// //   // Make your histograms object as a shared pointer that all the threads will have
// //   auto hists = std::make_shared<Histogram>(outfilename);

// //   auto run_files = [&hists, beam_energy](std::vector<std::string> inputs, int thread_id) {
// //     // Called once for each thread
// //     // Make a new chain to process for this thread
// //     auto chain = std::make_shared<TChain>("clas12");
// //     // Add every file to the chain
// //     for (auto in : inputs) chain->Add(in.c_str());
// //     // Run the function over each thread
// //     return run<Cuts>(std::move(chain), hists, thread_id, beam_energy);
// //   };

// //   // Start timer
// //   auto start = std::chrono::high_resolution_clock::now();
// //   // For each thread
// //   for (size_t i = 0; i < NUM_THREADS; i++) {
// //     // Set the thread to run a task A-Synchronously
// //     // The function we run is the first argument (run_files)
// //     // The functions arguments are all the remaining arguments
// //     threads[i] = std::async(run_files, infilenames.at(i), i);
// //   }

// //   // For each thread
// //   for (size_t i = 0; i < NUM_THREADS; i++) {
// //     // Get the information from the thread in this case how many events each thread actually computed
// //     events += threads[i].get();
// //   }

// //   // Timer and Hz calculator functions that print at the end
// //   std::cout.imbue(std::locale(""));  // Puts commas in
// //   std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
// //   std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
// //   std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
// //   return 0;
// // }
