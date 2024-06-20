/*
BSD 3-Clause License

Copyright (c) Alliance for Sustainable Energy, LLC. See also https://github.com/NREL/ssc/blob/develop/LICENSE
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Flat plate collector subcomponent model
#include "core.h"
//#include "tckernel.h"


#include "common.h"
#include "flat_plate_solar_collector.h"
#include <ctime>
#include <algorithm>
#include <iterator>

// signed/unsigned mismatch
#pragma warning (disable : 4388)

static var_info _cm_vtab_flat_plate_subcomponent[] = {

    // Weather Reader
    { SSC_INPUT,        SSC_STRING,      "file_name",                 "Local weather file with path",                                                     "none",         "",               "weather",        "?",                       "LOCAL_FILE",            "" },
    { SSC_INPUT,        SSC_TABLE,       "solar_resource_data",       "Weather resource data in memory",                                                  "",             "",               "weather",        "?",                       "",                      "" },

    // Case Parameters
    { SSC_INPUT,        SSC_NUMBER,      "start_step",                "Hour of year",                                                                     "-",            "",               "",               "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "m_dot",                     "Fluid mass flow rate into flat plate collector",                                   "kg/s",         "",               "",               "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_in",                      "Fluid temperature into flat plate collector",                                      "C",            "",               "",               "*",                       "",                      "" },
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         
    // Flat Plate Collectors                                                                                                                                                                                                             
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_tested_frta",    "FRta from certification testing",                                                  "none",         "",               "solar_field",    "?=0.733",                 "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_tested_frul",    "FRUL from certification testing",                                                  "W/m2-K",       "",               "solar_field",    "?=3.41",                  "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_tested_iam",     "IAM from certification testing",                                                   "none",         "",               "solar_field",    "?=-0.06",                 "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_tested_area_coll", "Area of collector used for certification testing",                               "m2",           "",               "solar_field",    "?=3.73",                  "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_tested_m_dot",   "Mass flow used during certification testing",                                      "kg/s",         "",               "solar_field",    "?=0.0498",                "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_tested_heat_capacity", "Heat capacity of fluid used for certification testing",                      "kJ/kg-K",      "",               "solar_field",    "?=4.182",                 "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_azimuth",        "Azimuth of flat plate collectors, clockwise from North",                           "deg",          "",               "solar_field",    "?=180",                   "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_tilt",           "Tilt of flat plate collectors",                                                    "deg",          "",               "solar_field",    "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plate_albedo",         "Ground reflectance factor",                                                        "0..1",         "",               "SWH",            "*",                       "",                      "" },

    { SSC_INPUT,        SSC_NUMBER,      "flat_plates_in_series",     "Number of flat plate collectors in series",                                        "",             "",               "trough_field",   "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "flat_plates_in_parallel",   "Number of flat plate collectors in parallel",                                      "",             "",               "trough_field",   "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "Fluid_FPC",                 "Flat plate array HTF fluid ID number",                                             "none",         "",               "solar_field",    "*",                       "",                      "" },
    { SSC_INPUT,        SSC_MATRIX,      "fl_props_FPC",              "User defined field fluid property data",                                           "-",            "",               "solar_field",    "*",                       "",                      "" },


    // Piping Parameters
    { SSC_INPUT,        SSC_NUMBER,      "pipe_diam",                 "Pipe diameter",                                                                    "m",            "",               "solar_field",    "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "pipe_k",                    "Pipe insulation thermal conductivity",                                             "W/m2 K",       "",               "solar_field",    "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "pipe_insul",                "Pipe insulation thickness",                                                        "m",            "",               "solar_field",    "*",                       "",                      "" },
    { SSC_INPUT,        SSC_NUMBER,      "pipe_length",               "Piping total length",                                                              "m",            "",               "solar_field",    "*",                       "",                      "" },



    // Outputs
    { SSC_OUTPUT,       SSC_NUMBER,      "T_out",                     "Temperature leaving solar collector array (and piping)",                           "C",            "",               "solar_field",    "*",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "Q_gain",                    "Heat gained in solar collector array",                                             "MWt",          "",               "solar_field",    "*",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "Q_loss",                    "Heat lost in solar collector array (and piping)",                                  "MWt",          "",               "solar_field",    "*",                       "",                      "" },




var_info_invalid };



class cm_flat_plate_subcomponent : public compute_module
{
public:

    cm_flat_plate_subcomponent()
    {
        add_var_info(_cm_vtab_flat_plate_subcomponent);
    }

    void exec()
    {
        // Uncomment following 2 lines to write cmod inputs to LK script
        //FILE* fp = fopen("flat_plate_cmod_to_lk.lk", "w");
        //write_cmod_to_lk_script(fp, m_vartab);

        // ********************************
        // ********************************
        // Weather reader
        // ********************************
        // ********************************
        C_csp_weatherreader weather_reader;
        C_csp_solver::S_sim_setup sim_setup;
        int n_steps_fixed;
        int steps_per_hour;
        {
            if (is_assigned("file_name")) {
                weather_reader.m_weather_data_provider = make_shared<weatherfile>(as_string("file_name"));
                if (weather_reader.m_weather_data_provider->has_message()) log(weather_reader.m_weather_data_provider->message(), SSC_WARNING);
            }
            if (is_assigned("solar_resource_data")) {
                weather_reader.m_weather_data_provider = make_shared<weatherdata>(lookup("solar_resource_data"));
                if (weather_reader.m_weather_data_provider->has_message()) log(weather_reader.m_weather_data_provider->message(), SSC_WARNING);
            }

            weather_reader.m_filename = as_string("file_name");
            weather_reader.m_trackmode = 0;
            //weather_reader.m_tilt = 0.0;
            //weather_reader.m_azimuth = 0.0;
            // Initialize to get weather file info
            weather_reader.init();
            if (weather_reader.has_error()) throw exec_error("trough_physical", weather_reader.get_error());

            // Set steps per hour
            double nhourssim = 8760.0;                                  //[hr] Number of hours to simulate

            sim_setup.m_sim_time_start = 0.0;                           //[s] starting first hour of year
            sim_setup.m_sim_time_end = nhourssim * 3600.;                 //[s] full year simulation

            steps_per_hour = 1;                                     //[-]

            int n_wf_records = (int)weather_reader.m_weather_data_provider->nrecords();
            steps_per_hour = n_wf_records / 8760;                       //[-]

            n_steps_fixed = steps_per_hour * 8760;                    //[-]
            sim_setup.m_report_step = 3600.0 / (double)steps_per_hour;  //[s]
        }

        // Parse HTF Properties
        HTFProperties htf_props;
        {
            int fluid_id = as_integer("Fluid_FPC");
            util::matrix_t<double> fl_props = as_matrix("fl_props_FPC");
            string error_msg;

            if (fluid_id != HTFProperties::User_defined)
            {
                if (!htf_props.SetFluid(fluid_id))
                {
                    throw(C_csp_exception("Field HTF code is not recognized", "Trough Collector Solver"));
                }
            }
            else if (fluid_id == HTFProperties::User_defined)
            {
                int n_rows = (int)fl_props.nrows();
                int n_cols = (int)fl_props.ncols();
                if (n_rows > 2 && n_cols == 7)
                {
                    if (!htf_props.SetUserDefinedFluid(fl_props))
                    {
                        
                        error_msg = util::format(htf_props.UserFluidErrMessage(), n_rows, n_cols);
                        throw(C_csp_exception(error_msg, "Flat Plate Collector Solver"));
                    }
                }
                else
                {
                    error_msg = util::format("The user defined field HTF table must contain at least 3 rows and exactly 7 columns. The current table contains %d row(s) and %d column(s)", n_rows, n_cols);
                    throw(C_csp_exception(error_msg, "Flat Plate Collector Solver"));
                }
            }
            else
            {
                throw(C_csp_exception("HTF code is not recognized", "Flat Plate Solver"));
            }
        }

        // Initialize flat plate collectors
        CollectorTestSpecifications collector_test_specifications;
        collector_test_specifications.FRta = as_double("flat_plate_tested_frta");
        collector_test_specifications.FRUL = as_double("flat_plate_tested_frul");
        collector_test_specifications.iam = as_double("flat_plate_tested_iam");
        collector_test_specifications.area_coll = as_double("flat_plate_tested_area_coll");             // [m2]
        collector_test_specifications.m_dot = as_double("flat_plate_tested_m_dot");                     // [kg/s]   
        collector_test_specifications.heat_capacity = as_double("flat_plate_tested_heat_capacity");     // [kJ/kg-K]

        CollectorLocation collector_location;
        collector_location.latitude = weather_reader.ms_solved_params.m_lat;		//[deg]
        collector_location.longitude = weather_reader.ms_solved_params.m_lon;	    //[deg]
        collector_location.timezone = weather_reader.ms_solved_params.m_tz;	    	//[hr]

        CollectorOrientation collector_orientation;
        collector_orientation.azimuth = as_double("flat_plate_azimuth");
        collector_orientation.tilt = as_double("flat_plate_tilt");

        ArrayDimensions array_dimensions;
        array_dimensions.num_in_series = as_integer("flat_plates_in_series");
        array_dimensions.num_in_parallel = as_integer("flat_plates_in_parallel");

        Pipe inlet_pipe(as_double("pipe_diam"), as_double("pipe_k"), as_double("pipe_insul"), as_double("pipe_length") / 2.0);
        Pipe outlet_pipe(inlet_pipe);

        FlatPlateArray flat_plate_array_ = FlatPlateArray(collector_test_specifications, collector_location,
            collector_orientation, array_dimensions,
            inlet_pipe, outlet_pipe);

        // Define Time Step
        double time_step = 3600;          // [s]
        double start_time_step = as_double("start_step");   // start time (number of steps from janurary 1)

        C_csp_solver_sim_info sim_info;
        sim_info.ms_ts.m_step = time_step;       // [s] Size of timestep

        // Read weather for time step
        weather_reader.read_time_step(start_time_step, sim_info);
        C_csp_weatherreader::S_outputs weather = weather_reader.ms_outputs;

        // Flat plate array
        tm datetime;
        datetime.tm_year = weather.m_year - 1900;  // years since 1900
        datetime.tm_mon = weather.m_month - 1;     // months since Jan. (Jan. = 0)
        datetime.tm_mday = weather.m_day;
        datetime.tm_hour = weather.m_hour;
        datetime.tm_min = weather.m_minute;
        datetime.tm_sec = 0.;
        ExternalConditions case_conditions;
        case_conditions.weather.ambient_temp = weather.m_tdry;
        case_conditions.weather.dni = weather.m_beam;
        case_conditions.weather.dhi = weather.m_diffuse;
        case_conditions.weather.ghi = weather.m_global;
        case_conditions.weather.wind_speed = weather.m_wspd;
        case_conditions.weather.wind_direction = weather.m_wdir;

        case_conditions.inlet_fluid_flow.m_dot = as_double("m_dot");    // [kg/s]
        case_conditions.inlet_fluid_flow.temp = as_double("T_in");      // [C]
        case_conditions.inlet_fluid_flow.fluid = htf_props;
        case_conditions.albedo = as_double("flat_plate_albedo");        // []

        case_conditions.inlet_fluid_flow.specific_heat = case_conditions.inlet_fluid_flow.fluid.Cp(case_conditions.inlet_fluid_flow.temp+273.15);

        // Simulate
        HeatAndTempInOut sim_results = flat_plate_array_.HeatFlowsAndOutletTemp(datetime, case_conditions);

        // Assign Outputs
        assign("T_out", sim_results.T_out);             // [C]
        assign("Q_gain", sim_results.Q_gain / 1.E3);    // [MWt]
        assign("Q_loss", sim_results.Q_loss / 1.E3);    // [MWt]

        int x = 0;
    }

    template <typename T>
    void set_vector(const std::string& name, const vector<T> vec)
    {
        int size = vec.size();
        ssc_number_t* alloc_vals = allocate(name, size);
        for (int i = 0; i < size; i++)
            alloc_vals[i] = vec[i];    // []
    }

};

DEFINE_MODULE_ENTRY(flat_plate_subcomponent, "Flat plate subcomponent", 1)
