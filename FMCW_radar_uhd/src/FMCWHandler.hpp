#ifndef FMCWHANDLER
#define FMCWHANDLER


//c standard library
    #include <iostream>
    #include <cstdlib>
    #include <string>
    #include <complex>
    #include <csignal>
    #include <thread>

//Radar Class
    #include "RADAR.hpp"
//JSON class
    #include <nlohmann/json.hpp>

using RADAR_namespace::RADAR;
using json = nlohmann::json;

namespace FMCWHandler_namespace {
    
    template<typename data_type>
    class FMCWHandler {
        private:
        
            bool radar_enabled;

            json radar_config;
            json fmcw_config;

            RADAR<data_type> Radar;

            //to support multiple runs
            bool multiple_runs;
            bool evaluate_spoofing_enabled;
            bool evaluate_parameter_estimation_enabled;
            std::string radar_tx_files_folder_path;
            size_t num_runs;

        public:
            /**
             * @brief Construct a new FMCWHandler object
             * 
             * @param fmcw_config_obj a JSON config object for the FMCW simulation
             * @param radar_config_obj a JSON config object for the victim
             * @param run (default false) on true, runs the experiment
             */
            FMCWHandler(json fmcw_config_obj, json radar_config_obj, bool run = false)
                :fmcw_config(fmcw_config_obj),
                radar_config(radar_config_obj),
                Radar() //use the default constructor
                {
                    
                    if(check_config())
                    {
                        get_enabled_status();
                        get_multiple_runs_configuration();
                        init_FMCW_devices();
                        if (run)
                        {
                            run_FMCW();
                        } 
                    }
            }
            
            /**
             * @brief Check the json config files to make sure all necessary parameters are included
             * 
             * @return true - JSON is all good and has required elements
             * @return false - JSON is missing certain fields
             */
            bool check_config(){
                bool config_good = true;

                //Radar enabled
                if(fmcw_config["Radar_enabled"].is_null()){
                    std::cerr << "FMCWHandler::check_config: Radar_enabled not in JSON" <<std::endl;
                    config_good = false;
                }

                //performing multiple runs
                if(fmcw_config["multiple_runs"]["perform_multiple_runs"].is_null()){
                    std::cerr << "FMCWHandler::check_config: perform_multiple_runs not in JSON" <<std::endl;
                    config_good = false;
                }else if (fmcw_config["multiple_runs"]["perform_multiple_runs"].get<bool>())
                {
                    //the number of runs to perform
                    if(fmcw_config["multiple_runs"]["num_runs"].is_null()){
                        std::cerr << "FMCWHandler::check_config: num_runs not in JSON" <<std::endl;
                        config_good = false;
                    }
                    //the number of runs to perform
                    if(fmcw_config["multiple_runs"]["evaluate_spoofing"].is_null()){
                        std::cerr << "FMCWHandler::check_config: evaluate_spoofing not in JSON" <<std::endl;
                        config_good = false;
                    }
                    //the number of runs to perform
                    if(fmcw_config["multiple_runs"]["evaluate_parameter_estimation"].is_null()){
                        std::cerr << "FMCWHandler::check_config: evaluate_parameter_estimation not in JSON" <<std::endl;
                        config_good = false;
                    }
                }
                return config_good;
            }

            /**
             * @brief Get the enabled status from the JSON
             * 
             */
            void get_enabled_status(void){
                radar_enabled = fmcw_config["Radar_enabled"].get<bool>();
            }

            /**
             * @brief Setup whether or not to perform multiple runs. If multiple runs are requested,
             * also get the number of runs
             * 
             */
            void get_multiple_runs_configuration(void){
                multiple_runs = fmcw_config["multiple_runs"]["perform_multiple_runs"].get<bool>();
                if(multiple_runs)
                {
                    num_runs = fmcw_config["multiple_runs"]["num_runs"].get<size_t>();
                    evaluate_spoofing_enabled = fmcw_config["multiple_runs"]["evaluate_spoofing"].get<bool>();
                    evaluate_parameter_estimation_enabled = fmcw_config["multiple_runs"]["evaluate_parameter_estimation"].get<bool>();
                }
            }

            /**
             * @brief Initialize the victim based on its enable status
             * 
             */
            void init_FMCW_devices(void){
                if (radar_enabled)
                {
                    //initialize the victim device
                    Radar.init(radar_config);
                }
            }
            
            /**
             * @brief Run the FMCW simulation with the vicitm in a separate thread
             * 
             * @param run_number (defaults to 0) the run number if multiple runs are being performed
             */
            void run_FMCW_trial(size_t run_number = 0){

                if(multiple_runs)
                {
                    std::cout << std::endl << std::endl << "FMCWHandler::run_FMCW_trial: running trial: " << run_number << std::endl;
                }
                else
                {
                    std::cout << "FMCWHandler::run_FMCW_trial: running trial" << std::endl;
                }

                if (radar_enabled)
                {
                    //initialize victim
                    Radar.initialize_radar(multiple_runs,run_number,evaluate_spoofing_enabled,evaluate_parameter_estimation_enabled);

                    //run the victim
                    Radar.run_RADAR();
                }
                else
                {
                    std::cout << "FMCWHandler::run_FMCW_trial: victim not enabled, no run performed" << std::endl;
                }

                if(multiple_runs)
                {
                    std::cout << "FMCWHandler::run_FMCW_trial: trial " << run_number << " complete" << std::endl;
                }
                else
                {
                    std::cout << "FMCWHandler::run_FMCW_trial: trial complete" << std::endl;
                }
            }

            /**
             * @brief Runs the FMCW implementation (runs multiple trials if multiple_trials is true)
             * 
             */
            void run_FMCW(void){
                
                if (multiple_runs)
                {
                    for (size_t i = 1; i <= num_runs; i++)
                    {
                        run_FMCW_trial(i);
                    }
                    
                }
                else //only a single run requested
                {
                    run_FMCW_trial();
                }
                
                
            }
    };
}

#endif