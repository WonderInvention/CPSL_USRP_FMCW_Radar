# CPSL Radar FMCW Code

The following repository was created to implement a simple FMCW radar on the USRP with offline processing using a MATLAB script

## Install Required Dependencies

To run the code, the following dependencies must also be installed:
1. nlohmann_json (at least 3.10.5)
2. uhd (at least 4.1.0)
3. c++ (at least C++ 11)

Instructions to install each package are available here:

### 1. Install nlohmann_json
1. clone the nlohmann_json repository
```
git clone https://github.com/nlohmann/json.git
```

2. create a build directory and navigate to it
```
mkdir build
cd build
```

3. run cmake with the path to the source directory
```
cmake ../json
```

4. run make to build the library
```
make
```

5. install the library
```
sudo make install
```

### 2. Install uhd from source

I installed uhd from source. This can be accomplished by the instructions on the [UHD Install Guide](https://files.ettus.com/manual/page_build_guide.html)


## Building C++ Code
To build the c++ code used to interact with the USRP radar, perform the following steps:

1. Clone the CPSL_USRP_FMCW_Radar git repository
```
git clone https://github.com/davidmhunt/CPSL_USRP_FMCW_Radar.git
```

2. Navigate to the FMCW_radar_uhd folder
```
cd CPSL_USRP_FMCW_Radar/FMCW_radar_uhd
```

3. make the build folder
```
mkdir build
```

4. build the project
```
cmake ..
make
```

## Running Experiments on the USRP

Running experiments can be accomplished using the following steps:
1. Start experiment in MATLAB
    1. (optional) create a config.json file for the radar
    2. Start running experiment in MATLAB
2. Run experiment on USRP
    1. (optional) update config.json file for radar
    2. Rebuild c++ code for updated experiment
    3. Run experiment
3. Process streamed data from USRP in MATLAB

### 1. Start Experiment in MATLAB

#### Loading an FMCW radar configuration
1. Open up the [DIAG_USRP_Run_Radar.mlx](MATLAB/DIAG_USRP_Run_Radar.mlx) matlab live script file.
2. In the first section ("Initialize the FMCW Configuration), update the config_folder_path variable with the path to the CPSL_USRP_FMCW_Radar/MATLAB/config_files folder on your machine.The [config_files](MATLAB/config_files/) folder contains several configurations that I have used in the past, however you can also create your own configuration as well. To do so, complete the following steps:
    1. 