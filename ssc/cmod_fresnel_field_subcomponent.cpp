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

// Trough CSP - physical model
#include "core.h"
//#include "tckernel.h"

// for adjustment factors
#include "common.h"

#include "csp_solver_core.h"
#include "csp_solver_fresnel_collector_receiver.h"
#include <ctime>
#include <algorithm>
#include <iterator>

// signed/unsigned mismatch
#pragma warning (disable : 4388)

static var_info _cm_vtab_fresnel_field_subcomponent[] = {

    /* VARTYPE          DATATYPE         NAME                         LABEL                                                                               UNITS           META              GROUP             REQUIRED_IF                CONSTRAINTS         UI_HINTS*/

    { SSC_INPUT,        SSC_NUMBER,      "sim_type",                  "1 (default): timeseries, 2: design only",                                          "",             "",               "System Control", "?=1",                    "",                     ""},

    // Weather Reader
    { SSC_INPUT,        SSC_STRING,      "file_name",                 "Local weather file with path",                                                     "none",         "",               "weather",        "*",                       "LOCAL_FILE",          "" },


    // Case Parameters
    { SSC_INPUT,        SSC_NUMBER,      "time_step",                 "Length of time step",                                                              "s",            "",               "",               "*",                      "",                     ""},
    { SSC_INPUT,        SSC_NUMBER,      "start_step",                "Number of time steps (from beginning of year) to start",                           "-",            "",               "",               "*",                      "",                     ""},
    { SSC_INPUT,        SSC_NUMBER,      "T_htf_in",                  "Temperature of HTF into field",                                                    "C",            "",               "",               "*",                      "",                     ""},
    { SSC_INPUT,        SSC_NUMBER,      "field_mode",                "Field operation mode (0,1,2,3:OFF,OFFNOSTARTUP,STARTUP,ON)",                       "",             "",               "",               "*",                      "",                     ""},
    { SSC_INPUT,        SSC_NUMBER,      "defocus_field_control",     "Field defocus control (0-1)",                                                      "-",            "",               "",               "field_mode=3",           "",                     "" },


    // System Design

    { SSC_INPUT,        SSC_NUMBER,      "solar_mult_in",             "Solar multiple Input",                                                             "",             "",               "System_Design",  "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "total_Ap_in",               "Field aperture Input",                                                             "m3",           "",               "System_Design",  "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "solar_mult_or_Ap",          "Design using specified solar mult or field aperture",                              "m3",           "",               "System_Design",  "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_loop_in_des",             "Design loop inlet temperature",                                                    "C",            "",               "System_Design",  "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_loop_out",                "Target loop outlet temperature",                                                   "C",            "",               "System_Design",  "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "I_bn_des",                  "Solar irradiation at design",                                                      "W/m2",         "",               "System_Design",  "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "tshours",                   "Equivalent full-load thermal storage hours",                                       "hr",           "",               "System_Design",  "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "q_pb_design",               "Design heat input to power block",                                                 "MWt",          "",               "System_Design",  "*",                      "",                     "" },


    // Solar Field

    { SSC_INPUT,        SSC_NUMBER,      "nMod",                      "Number of collector modules in a loop",                                            "",             "",               "Solar_Field",    "*",                      "INTEGER",              "" },
    { SSC_INPUT,        SSC_NUMBER,      "eta_pump",                  "HTF pump efficiency",                                                              "",             "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "HDR_rough",                 "Header pipe roughness",                                                            "m",            "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "theta_stow",                "stow angle",                                                                       "deg",          "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "theta_dep",                 "deploy angle",                                                                     "deg",          "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "FieldConfig",               "Number of subfield headers",                                                       "",             "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "m_dot_htfmin",              "Minimum loop HTF flow rate",                                                       "kg/s",         "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "m_dot_htfmax",              "Maximum loop HTF flow rate",                                                       "kg/s",         "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "Fluid",                     "Field HTF fluid number",                                                           "",             "",               "Solar_Field",    "*",                      "INTEGER",              "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_fp",                      "Freeze protection temperature (heat trace activation temperature)",                "C",            "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "V_hdr_max",                 "Maximum HTF velocity in the header at design",                                     "m/s",          "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "V_hdr_min",                 "Minimum HTF velocity in the header at design",                                     "m/s",          "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "Pipe_hl_coef",              "Loss coefficient from the header - runner pipe - and non-HCE piping",              "W/m2-K",       "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "mc_bal_hot",                "The heat capacity of the balance of plant on the hot side",                        "kWht/K-MWt",   "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "mc_bal_cold",               "The heat capacity of the balance of plant on the cold side",                       "kWht/K-MWt",   "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "mc_bal_sca",                "Non-HTF heat capacity associated with each SCA - per meter basis",                 "Wht/K-m",      "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "water_per_wash",            "Water usage per wash",                                                             "L/m2_aper",    "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "washes_per_year",           "Mirror washing frequency",                                                         "none",         "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "rec_htf_vol",               "Volume of HTF in a single collector unit per unit aperture area",                  "L/m2-ap",      "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_amb_sf_des",              "Ambient design-point temperature for the solar field",                             "C",            "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "V_wind_des",                "Design-point wind velocity",                                                       "m/s",          "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_MATRIX,      "field_fl_props",            "Fluid property data",                                                              "",             "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "SCA_drives_elec",           "Tracking power in Watts per SCA drive",                                            "W/module",     "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "land_mult",                 "Non-solar field land area multiplier",                                             "-",            "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_startup",                 "Power block startup temperature",                                                  "C",            "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "rec_su_delay",              "Fixed startup delay time for the receiver",                                        "hr",           "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "rec_qf_delay",              "Energy-based receiver startup delay (fraction of rated thermal power)",            "-",            "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "p_start",                   "Collector startup energy, per SCA",                                                "kWe-hr",       "",               "Solar_Field",    "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "L_rnr_pb",                  "Length of runner pipe in power block",                                             "m",            "",               "Solar_Field",    "*",                      "",                     "" },

    // Collector and Receiver

    { SSC_INPUT,        SSC_NUMBER,      "ColAz",                     "Collector azimuth angle",                                                          "deg",          "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "opt_model",                 "The optical model",                                                                "",             "",               "Col_Rec",        "*",                      "INTEGER",              "" },
    { SSC_INPUT,        SSC_NUMBER,      "A_aperture",                "Reflective aperture area of the collector",                                        "m2",           "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "reflectivity",              "Solar-weighted mirror reflectivity value",                                         "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "TrackingError",             "Tracking error derate",                                                            "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "GeomEffects",               "Geometry effects derate",                                                          "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "Dirt_mirror",               "User-defined dirt on mirror derate",                                               "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "Error",                     "User-defined general optical error derate",                                        "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "L_mod",                     "The length of the collector module",                                               "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "IAM_T_coefs",               "Incidence angle modifier coefficients - transversal plane",                        "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "IAM_L_coefs",               "Incidence angle modifier coefficients - longitudinal plane",                       "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_MATRIX,      "OpticalTable",              "Values of the optical efficiency table",                                           "",             "",               "Col_Rec",        "*",                      "",                     "" },
                                     
    { SSC_INPUT,        SSC_NUMBER,      "rec_model",                 "Receiver model type (1=Polynomial ; 2=Evac tube)",                                 "",             "",               "Col_Rec",        "*",                      "INTEGER",              "" },
    { SSC_INPUT,        SSC_ARRAY,       "HCE_FieldFrac",             "The fraction of the field occupied by this HCE type",                              "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "D_abs_in",                  "The inner absorber tube diameter",                                                 "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "D_abs_out",                 "The outer absorber tube diameter",                                                 "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "D_glass_in",                "The inner glass envelope diameter",                                                "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "D_glass_out",               "The outer glass envelope diameter",                                                "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "D_plug",                    "The diameter of the absorber flow plug (optional)",                                "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "Flow_type",                 "The flow type through the absorber",                                               "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "Rough",                     "Roughness of the internal surface",                                                "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "alpha_env",                 "Envelope absorptance",                                                             "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_MATRIX,      "epsilon_abs_1",             "Absorber emittance - HCE variation 1",                                             "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_MATRIX,      "epsilon_abs_2",             "Absorber emittance - HCE variation 2",                                             "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_MATRIX,      "epsilon_abs_3",             "Absorber emittance - HCE variation 3",                                             "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_MATRIX,      "epsilon_abs_4",             "Absorber emittance - HCE variation 4",                                             "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "alpha_abs",                 "Absorber absorptance",                                                             "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "Tau_envelope",              "Envelope transmittance",                                                           "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "epsilon_glass",             "Glass envelope emissivity",                                                        "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "GlazingIntactIn",           "The glazing intact flag",                                                          "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "P_a",                       "Annulus gas pressure",                                                             "torr",         "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "AnnulusGas",                "Annulus gas type (1=air; 26=Ar; 27=H2)",                                           "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "AbsorberMaterial",          "Absorber material type",                                                           "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "Shadowing",                 "Receiver bellows shadowing loss factor",                                           "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "dirt_env",                  "Loss due to dirt on the receiver envelope",                                        "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "Design_loss",               "Receiver heat loss at design",                                                     "W/m",          "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "L_mod_spacing",             "Piping distance between sequential modules in a loop",                             "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "L_crossover",               "Length of crossover piping in a loop",                                             "m",            "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "HL_T_coefs",                "HTF temperature-dependent heat loss coefficients",                                 "W/m-K",        "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "HL_w_coefs",                "Wind-speed-dependent heat loss coefficients",                                      "W/m-(m/s)",    "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "DP_nominal",                "Pressure drop across a single collector assembly at design",                       "bar",          "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_ARRAY,       "DP_coefs",                  "Pressure drop mass flow based part-load curve",                                    "",             "",               "Col_Rec",        "*",                      "",                     "" },
    { SSC_INPUT,        SSC_NUMBER,      "nRecVar",                   "Number of receiver variations",                                                    "",             "",               "Col_Rec",        "?=4",                    "INTEGER",              "" },

        // Optional Component Initialization (state at start of first timestep)
    // Trough field
    { SSC_INPUT,        SSC_NUMBER,      "rec_op_mode_initial",       "Initial receiver operating mode 0: off, 1: startup, 2: on",                        "-",            "",               "System Control", "",                       "",                    "" },
    { SSC_INPUT,        SSC_NUMBER,      "defocus_initial",           "Initial receiver defocus",                                                         "-",            "",               "System Control", "",                       "",                    "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_in_loop_initial",         "Initial loop inlet, cold header and cold runner fluid temperature",                "C",            "",               "System Control", "",                       "",                    "" },
    { SSC_INPUT,        SSC_NUMBER,      "T_out_loop_initial",        "Initial loop outlet, hot header and hot runner fluid temperature",                 "C",            "",               "System Control", "",                       "",                    "" },
    { SSC_INPUT,        SSC_ARRAY,       "T_out_scas_initial",        "Initial SCA outlet temperatures",                                                  "C",            "",               "System Control", "",                       "",                    "" },


    // OUTPUTS
        // Design Point Outputs

         // Solar Field
    { SSC_OUTPUT,       SSC_NUMBER,      "A_loop",                           "Aperture of a single loop",                                            "m2",           "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "loop_opt_eff",                     "Loop optical efficiency at design",                                    "",             "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "loop_therm_eff",                   "Loop thermal efficiency at design",                                    "",             "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "loop_eff",                         "Total loop conversion efficiency at design",                           "",             "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "sm1_aperture",                     "Total required aperture, SM=1",                                        "m2",           "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "sm1_nLoops",                       "Required number of loops, SM=1",                                       "",             "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "total_tracking_power",             "Design tracking power",                                                "MW",           "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "A_field",                          "Total field aperture",                                                 "m2",           "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "q_field_des_actual",               "Design-point thermal power from the solar field limited by mass flow", "MW",           "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "q_field_des_ideal",                "Design-point thermal power from the solar field with no limit",        "MW",           "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "field_area",                       "Solar field area",                                                     "acres",        "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "total_land_area",                  "Total land area",                                                      "acres",        "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "field_htf_min_temp",               "Minimum field htf temp",                                               "C",            "",         "Power Cycle",                    "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "field_htf_max_temp",               "Maximum field htf temp",                                               "C",            "",         "Power Cycle",                    "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "mdot_field_des",                   "Field design HTF mass flow rate",                                      "kg/s",         "",         "Receiver",                       "*",                                                                "",              "" },
                                         
    { SSC_OUTPUT,       SSC_NUMBER,      "dP_field_des_SS",                  "Steady State Field design total pressure drop",                        "bar",          "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "Q_field_des_SS",                   "Steady State Field design thermal power",                              "MWt",           "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "T_field_out_des_SS",               "Steady State Field design outlet temperature",                         "C",            "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "m_dot_des_SS",                     "Steady State Field mass flow rate",                                    "kg/s",         "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "m_dot_loop_des_SS",                "Steady State Loop mass flow rate",                                     "kg/s",         "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "V_hdr_min_des_SS",                 "Steady State min header velocity",                                     "m/s",          "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "V_hdr_max_des_SS",                 "Steady State max header velocity",                                     "m/s",          "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "eta_optical_des_SS",               "Steady State optical efficiency",                                      "",             "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "therm_eff_des_SS",                 "Steady State field optical efficiency",                                "",             "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "eff_des_SS",                       "Steady State field total efficiency",                                  "",             "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "W_dot_pump_des_SS",                "Steady State field pumping power",                                           "MWe",             "",          "Receiver",                       "*",                                                                "",              "" },
                                         
                                         
    { SSC_OUTPUT,       SSC_NUMBER,      "T_loop_out_des_SS",                "Steady State loop design outlet temperature",                          "C",            "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "Q_loop_des_SS",                    "Steady State loop design thermal power",                               "MWt",          "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "therm_eff_loop_des_SS",            "Steady State loop optical efficiency",                                 "",             "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "eff_loop_des_SS",                  "Steady State loop total efficiency",                                   "",             "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "W_dot_pump_des_SS",                "Steady State field pumping power",                                     "MWe",          "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "Q_loss_receiver_des_SS",           "Steady State field heat loss from receiver",                           "MWt",          "",          "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "Q_loss_hdr_rnr_des_SS",            "Steady State field heat loss from headers and runners",                "MWt",          "",          "Receiver",                       "*",                                                                "",              "" },
                                         
                                         
    // Collector and Receiver            
    { SSC_OUTPUT,       SSC_NUMBER,      "DP_pressure_loss",                 "Total loop pressure loss at design",                                   "bar",          "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "avg_dt_des",                       "Average field temp difference at design",                              "C",            "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "hl_des",                           "Heat loss at design",                                                  "W/m",          "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "opt_derate",                       "Receiver optical derate",                                              "",             "",         "Receiver",                       "*",                                                                "",              "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "opt_normal",                       "Collector optical loss at normal incidence",                           "",             "",         "Receiver",                       "*",                                                                "",              "" },

    // Hourly

        // Solar Field
    { SSC_OUTPUT,       SSC_NUMBER,      "EqOpteff",                         "Field optical efficiency before defocus",                              "",             "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "SCAs_def",                         "Field fraction of focused SCAs",                                       "",             "",         "solar_field",    "sim_type=1",                       "",                      "" },
                        
    { SSC_OUTPUT,       SSC_NUMBER,      "q_inc_sf_tot",                     "Field thermal power incident",                                         "MWt",          "",         "solar_field",    "*",                                "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "q_dot_rec_inc",                    "Receiver thermal power incident",                                      "MWt",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "q_dot_rec_thermal_loss",           "Receiver thermal losses",                                              "MWt",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "q_dot_rec_abs",                    "Receiver thermal power absorbed",                                      "MWt",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "rec_thermal_eff",                  "Receiver thermal efficiency",                                          "",             "",         "solar_field",    "sim_type=1",                       "",                      "" },
                        
    { SSC_OUTPUT,       SSC_NUMBER,      "q_dot_piping_loss",                "Field piping thermal losses",                                          "MWt",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "e_dot_field_int_energy",           "Field change in material/htf internal energy",                         "MWt",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "q_dot_htf_sf_out",                 "Field thermal power leaving in HTF",                                   "MWt",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "q_dot_freeze_prot",                "Field freeze protection required",                                     "MWt",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
                        
    { SSC_OUTPUT,       SSC_NUMBER,      "m_dot_loop",                       "Receiver mass flow rate",                                              "kg/s",         "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "m_dot_field_recirc",               "Field total mass flow recirculated",                                   "kg/s",         "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "m_dot_field_delivered",            "Field total mass flow delivered",                                      "kg/s",         "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "T_field_cold_in",                  "Field timestep-averaged inlet temperature",                            "C",            "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "T_rec_cold_in",                    "Loop timestep-averaged inlet temperature",                             "C",            "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "T_rec_hot_out",                    "Loop timestep-averaged outlet temperature",                            "C",            "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "T_field_hot_out",                  "Field timestep-averaged outlet temperature",                           "C",            "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "deltaP_field",                     "Field pressure drop",                                                  "bar",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
                        
    { SSC_OUTPUT,       SSC_NUMBER,      "W_dot_sca_track",                  "Field collector tracking power",                                       "MWe",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "W_dot_field_pump",                 "Field htf pumping power",                                              "MWe",          "",         "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "recirculating",                    "Field recirculating (bypass valve open)",                              "-",            "",         "solar_field",    "sim_type=1",                       "",                      "" },

    { SSC_OUTPUT,       SSC_NUMBER,      "rec_op_mode_final",         "Final receiver operating mode (0,1,2,3:OFF,OFFNOSTARTUP,STARTUP,ON)",              "-",            "",               "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "defocus_final",             "Defocus final",                                                                    "-",            "",               "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "T_in_loop_final",           "Final loop inlet, cold header and cold runner fluid temperature",                  "C",            "",               "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "T_out_loop_final",          "Final loop outlet, hot header and hot runner fluid temperature",                   "C",            "",               "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_ARRAY,       "T_out_scas_last_final",     "Final SCA outlet temperatures",                                                    "C",            "",               "solar_field",    "sim_type=1",                       "",                      "" },
    { SSC_OUTPUT,       SSC_NUMBER,      "time_required_su",          "Time required to startup (only used for startup mode)",                            "s",            "",               "solar_field",    "sim_type=1",                       "",                      "" },



var_info_invalid };

class cm_fresnel_field_subcomponent : public compute_module
{
public:

    cm_fresnel_field_subcomponent()
    {
        add_var_info(_cm_vtab_fresnel_field_subcomponent);
    }

    void exec()
    {
        // Uncomment following 2 lines to write cmod inputs to LK script
        //FILE* fp = fopen("fresnel_iph_cmod_to_lk.lk", "w");
        //write_cmod_to_lk_script(fp, m_vartab);

        // Common Parameters
        int sim_type = as_number("sim_type");
        double q_dot_pc_des = as_double("q_pb_design");         //[MWt] HEAT SINK design thermal power

        // *****************************************************
        // System Design Parameters


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
            weather_reader.m_tilt = 0.0;
            weather_reader.m_azimuth = 0.0;
            // Initialize to get weather file info
            weather_reader.init();
            if (weather_reader.has_error()) throw exec_error("fresnel_physical", weather_reader.get_error());

            // Set up ssc output arrays
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


        // ********************************
        // ********************************
        // Solar field, fresnel
        // ********************************
        // ********************************
        C_csp_fresnel_collector_receiver c_fresnel;
        {
            // Inputs
            {
                c_fresnel.m_solar_mult_or_Ap = as_integer("solar_mult_or_Ap");
                c_fresnel.m_solar_mult_in = as_double("solar_mult_in");
                c_fresnel.m_total_Ap_in = as_double("total_Ap_in");

                c_fresnel.m_nMod = as_integer("nMod");
                c_fresnel.m_nRecVar = as_integer("nRecVar");

                c_fresnel.m_eta_pump = as_number("eta_pump");
                c_fresnel.m_HDR_rough = as_number("HDR_rough");
                c_fresnel.m_theta_stow = as_number("theta_stow");
                c_fresnel.m_theta_dep = as_number("theta_dep");
                c_fresnel.m_FieldConfig = as_integer("FieldConfig");
                c_fresnel.m_T_startup = as_number("T_startup");

                // Set P_ref = q_pb_design and eta_ref = 1 
                c_fresnel.m_P_ref = as_double("q_pb_design") * 1e6;
                c_fresnel.m_eta_ref = 1;

                c_fresnel.m_m_dot_htfmin = as_number("m_dot_htfmin");
                c_fresnel.m_m_dot_htfmax = as_number("m_dot_htfmax");
                c_fresnel.m_T_loop_in_des = as_number("T_loop_in_des");

                c_fresnel.m_T_loop_out_des = as_number("T_loop_out");
                c_fresnel.m_Fluid = as_integer("Fluid");

                c_fresnel.m_field_fl_props = as_matrix("field_fl_props");
                c_fresnel.m_T_fp = as_number("T_fp");
                c_fresnel.m_I_bn_des = as_number("I_bn_des");
                c_fresnel.m_V_hdr_max = as_number("V_hdr_max");
                c_fresnel.m_V_hdr_min = as_number("V_hdr_min");
                c_fresnel.m_Pipe_hl_coef = as_number("Pipe_hl_coef");
                c_fresnel.m_SCA_drives_elec = as_number("SCA_drives_elec");
                c_fresnel.m_ColAz = as_number("ColAz");

                c_fresnel.m_mc_bal_hot = as_number("mc_bal_hot");
                c_fresnel.m_mc_bal_cold = as_number("mc_bal_cold");
                c_fresnel.m_mc_bal_sca = as_number("mc_bal_sca");

                c_fresnel.m_opt_model = as_integer("opt_model");

                c_fresnel.m_A_aperture = as_number("A_aperture");
                c_fresnel.m_reflectivity = as_number("reflectivity");
                c_fresnel.m_TrackingError = as_number("TrackingError");
                c_fresnel.m_GeomEffects = as_number("GeomEffects");
                c_fresnel.m_Dirt_mirror = as_number("Dirt_mirror");
                c_fresnel.m_Error = as_number("Error");
                c_fresnel.m_L_mod = as_number("L_mod");

                c_fresnel.m_IAM_T_coefs = as_vector_double("IAM_T_coefs");
                c_fresnel.m_IAM_L_coefs = as_vector_double("IAM_L_coefs");
                c_fresnel.m_OpticalTable = as_matrix("OpticalTable");
                c_fresnel.m_rec_model = as_integer("rec_model");

                c_fresnel.m_HCE_FieldFrac = as_vector_double("HCE_FieldFrac");
                c_fresnel.m_D_abs_in = as_vector_double("D_abs_in");
                c_fresnel.m_D_abs_out = as_vector_double("D_abs_out");
                c_fresnel.m_D_glass_in = as_vector_double("D_glass_in");
                c_fresnel.m_D_glass_out = as_vector_double("D_glass_out");
                c_fresnel.m_D_plug = as_vector_double("D_plug");
                c_fresnel.m_Flow_type = as_vector_double("Flow_type");
                c_fresnel.m_Rough = as_vector_double("Rough");
                c_fresnel.m_alpha_env = as_vector_double("alpha_env");

                c_fresnel.m_epsilon_abs_1 = as_matrix_transpose("epsilon_abs_1");
                c_fresnel.m_epsilon_abs_2 = as_matrix_transpose("epsilon_abs_2");
                c_fresnel.m_epsilon_abs_3 = as_matrix_transpose("epsilon_abs_3");
                c_fresnel.m_epsilon_abs_4 = as_matrix_transpose("epsilon_abs_4");

                c_fresnel.m_alpha_abs = as_vector_double("alpha_abs");
                c_fresnel.m_Tau_envelope = as_vector_double("Tau_envelope");
                c_fresnel.m_epsilon_glass = as_vector_double("epsilon_glass");
                c_fresnel.m_GlazingIntact = as_vector_bool("GlazingIntactIn");

                c_fresnel.m_P_a = as_vector_double("P_a");

                c_fresnel.m_AnnulusGas = as_vector_double("AnnulusGas");
                c_fresnel.m_AbsorberMaterial = as_vector_double("AbsorberMaterial");
                c_fresnel.m_Shadowing = as_vector_double("Shadowing");
                c_fresnel.m_dirt_env = as_vector_double("dirt_env");
                c_fresnel.m_Design_loss = as_vector_double("Design_loss");

                c_fresnel.m_L_mod_spacing = as_number("L_mod_spacing");
                c_fresnel.m_L_crossover = as_number("L_crossover");
                c_fresnel.m_HL_T_coefs = as_vector_double("HL_T_coefs");
                c_fresnel.m_HL_w_coefs = as_vector_double("HL_w_coefs");

                c_fresnel.m_DP_nominal = as_number("DP_nominal");
                c_fresnel.m_DP_coefs = as_vector_double("DP_coefs");
                c_fresnel.m_rec_htf_vol = as_number("rec_htf_vol");

                c_fresnel.m_L_rnr_pb = as_number("L_rnr_pb"); // No power block line length
                c_fresnel.m_rec_su_delay = as_number("rec_su_delay");
                c_fresnel.m_rec_qf_delay = as_number("rec_qf_delay");
                c_fresnel.m_p_start = as_number("p_start");

                c_fresnel.m_V_wind_des = as_number("V_wind_des");
                c_fresnel.m_T_amb_sf_des = as_number("T_amb_sf_des");

                // Check initialization variables
                if (is_assigned("rec_op_mode_initial")) {
                    c_fresnel.m_operating_mode_initial = (C_csp_collector_receiver::E_csp_cr_modes)as_integer("rec_op_mode_initial");
                }
                if (is_assigned("defocus_initial")) {
                    c_fresnel.m_defocus_initial = as_integer("defocus_initial");
                }
                if (is_assigned("T_in_loop_initial")) {
                    c_fresnel.m_T_in_loop_initial = as_double("T_in_loop_initial") + 273.15; // [K] Convert to Kelvin
                }
                if (is_assigned("T_out_loop_initial")) {
                    c_fresnel.m_T_out_loop_initial = as_double("T_out_loop_initial") + 273.15; // [K] Convert to Kelvin
                }
                if (is_assigned("T_out_scas_initial")) {
                    size_t n_T_out_scas_last_initial = -1;
                    ssc_number_t* T_out_scas_last_initial = as_array("T_out_scas_initial", &n_T_out_scas_last_initial);
                    std::copy(T_out_scas_last_initial, T_out_scas_last_initial + n_T_out_scas_last_initial, back_inserter(c_fresnel.m_T_out_scas_last_initial));
                    for (double& temp : c_fresnel.m_T_out_scas_last_initial)
                        temp += 273.15; // [K] Convert to Kelvin
                }
            }

            // Calculate solar multiple (needed for other component constructors)
            // Need latitude from weather reader
            weather_reader.init();
            c_fresnel.design_solar_mult(weather_reader.ms_solved_params.m_lat);
        }

        // Initialize fresnel Field
        C_csp_collector_receiver::S_csp_cr_init_inputs init_inputs;
        init_inputs.m_latitude = weather_reader.ms_solved_params.m_lat;		//[deg]
        init_inputs.m_longitude = weather_reader.ms_solved_params.m_lon;	//[deg]
        init_inputs.m_tz = weather_reader.ms_solved_params.m_tz;	    	//[hr]
        init_inputs.m_shift = weather_reader.ms_solved_params.m_shift;		//[deg]
        init_inputs.m_elev = weather_reader.ms_solved_params.m_elev;		//[m]
        C_csp_collector_receiver::S_csp_cr_solved_params cr_solved_params;

        c_fresnel.init(init_inputs, cr_solved_params);

        // Output Design Point Calculations
        {
            // System Design Calcs
            //double eta_ref = as_double("eta_ref");                          //[-]
            //double W_dot_cycle_des = as_double("P_ref");                    //[MWe]
            double tshours = as_double("tshours");                          //[-]
            double solar_mult_des = c_fresnel.m_solar_mult;
            double q_pb_design = as_double("q_pb_design");

            double q_dot_pc_des = q_pb_design;               //[MWt]
            double Q_tes = q_dot_pc_des * tshours;                          //[MWt-hr]

            double mdot_field_des = c_fresnel.m_m_dot_design;          // [kg/s]

            double avg_T_des = (c_fresnel.m_T_loop_in_des + c_fresnel.m_T_loop_out_des) / 2.0;

            // Solar Field

            double W_dot_col_tracking_des = c_fresnel.get_tracking_power();                 // [MWe]
            double A_loop = c_fresnel.m_A_loop;                                             // [m2]
            double loop_opt_eff = c_fresnel.m_loop_opt_eff;
            double loop_therm_eff = c_fresnel.m_loop_therm_eff;
            double loop_eff = c_fresnel.m_loop_eff;
            double sm1_aperture = c_fresnel.m_Ap_sm1;                                       // [m2]
            double sm1_nLoops = c_fresnel.m_nLoops_sm1;
            double total_tracking_power = c_fresnel.m_W_dot_sca_tracking_nom;               // [MW]
            double A_field = c_fresnel.m_Ap_tot;                                            // [m2]
            double q_field_des = c_fresnel.m_q_design_actual / 1e6;                         // [MW]
            double q_field_des_ideal = c_fresnel.m_q_design_ideal / 1e6;                    // [MW]

            double field_area = A_field / 4046.85642;                                       // [acres] (convert m2 to acre)
            double land_mult = as_double("land_mult");
            double total_land_area = field_area * land_mult;                                // [acres]

            double field_htf_min_temp = c_fresnel.m_htfProps.min_temp() - 273.15;           // [C]
            double field_htf_max_temp = c_fresnel.m_htfProps.max_temp() - 273.15;           // [C]

            // steady state results
            double dP_field_des_SS = c_fresnel.m_dP_des_SS;                                 // [bar]
            double Q_field_des_SS = c_fresnel.m_Q_field_des_SS / 1e6;                       // [MW]
            double T_field_out_des_SS = c_fresnel.m_T_field_out_des_SS;                     // [C]
            double m_dot_des_SS = c_fresnel.m_m_dot_des_SS;                                 // [kg/s]
            double m_dot_loop_des_SS = c_fresnel.m_m_dot_loop_des_SS;                       // [kg/s]
            double V_hdr_min_des_SS = c_fresnel.m_V_hdr_min_des_SS;                         // [m/s]
            double V_hdr_max_des_SS = c_fresnel.m_V_hdr_max_des_SS;                         // [m/s]
            double eta_optical_des_SS = c_fresnel.m_eta_optical_des_SS;
            double therm_eff_des_SS = c_fresnel.m_therm_eff_des_SS;
            double eff_des_SS = c_fresnel.m_eff_des_SS;
            double W_dot_pump_des_SS = c_fresnel.m_W_dot_pump_des_SS;                       // [MWe]

            double T_loop_out_des_SS = c_fresnel.m_T_loop_out_des_SS;                       // [C]
            double Q_loop_des_SS = c_fresnel.m_Q_loop_des_SS / 1e6;                         // [MW]
            double therm_eff_loop_des_SS = c_fresnel.m_therm_eff_loop_des_SS;
            double eff_loop_des_SS = c_fresnel.m_eff_loop_des_SS;
            double Q_loss_receiver_des_SS = c_fresnel.m_Q_loss_receiver_des_SS;             // [MWt]
            double Q_loss_hdr_rnr_des_SS = c_fresnel.m_Q_loss_hdr_rnr_des_SS;               // [MWt]

            // Assign
            {
                assign("A_loop", A_loop);
                assign("loop_opt_eff", loop_opt_eff);
                assign("loop_therm_eff", loop_therm_eff);
                assign("loop_eff", loop_eff);
                assign("sm1_aperture", sm1_aperture);
                assign("sm1_nLoops", sm1_nLoops);
                assign("total_tracking_power", total_tracking_power);
                assign("A_field", A_field);
                assign("q_field_des_actual", q_field_des);
                assign("q_field_des_ideal", q_field_des_ideal);
                assign("field_area", field_area);
                assign("total_land_area", total_land_area);
                assign("field_htf_min_temp", field_htf_min_temp);
                assign("field_htf_max_temp", field_htf_max_temp);
                assign("mdot_field_des", mdot_field_des);

                assign("dP_field_des_SS", dP_field_des_SS);
                assign("Q_field_des_SS", Q_field_des_SS);
                assign("T_field_out_des_SS", T_field_out_des_SS);
                assign("m_dot_des_SS", m_dot_des_SS);
                assign("m_dot_loop_des_SS", m_dot_loop_des_SS);
                assign("V_hdr_min_des_SS", V_hdr_min_des_SS);
                assign("V_hdr_max_des_SS", V_hdr_max_des_SS);
                assign("eta_optical_des_SS", eta_optical_des_SS);
                assign("therm_eff_des_SS", therm_eff_des_SS);
                assign("eff_des_SS", eff_des_SS);
                assign("W_dot_pump_des_SS", W_dot_pump_des_SS);

                assign("T_loop_out_des_SS", T_loop_out_des_SS);
                assign("Q_loop_des_SS", Q_loop_des_SS);
                assign("therm_eff_loop_des_SS", therm_eff_loop_des_SS);
                assign("eff_loop_des_SS", eff_loop_des_SS);

                assign("Q_loss_receiver_des_SS", Q_loss_receiver_des_SS);
                assign("Q_loss_hdr_rnr_des_SS", Q_loss_hdr_rnr_des_SS);
            }

            // Collector and Receiver
            double DP_pressure_loss = c_fresnel.m_nMod * c_fresnel.m_DP_nominal;            // [bar]
            double avg_dt_des = c_fresnel.m_dT_des;                                         // [C]
            double hl_des = c_fresnel.m_hl_des;                                             // [W/m]
            double opt_derate = c_fresnel.m_opt_derate;
            double opt_normal = c_fresnel.m_opt_normal;

            // Assign
            {
                assign("DP_pressure_loss", DP_pressure_loss);
                assign("avg_dt_des", avg_dt_des);
                assign("hl_des", hl_des);
                assign("opt_derate", opt_derate);
                assign("opt_normal", opt_normal);
            }

        }

        // Return if only called for design point
        if (sim_type != 1)
            return;

        // Simulate
        int field_mode = as_integer("field_mode");
        double T_htf_inlet = as_double("T_htf_in");    // [C]

        // Define Time Step
        double time_step = as_double("time_step");          // [s]
        double start_time_step = as_double("start_step");   // start time (number of steps from janurary 1)

        C_csp_solver_sim_info sim_info;
        sim_info.ms_ts.m_step = time_step;       // [s] Size of timestep

        // Read weather for time step
        weather_reader.read_time_step(start_time_step, sim_info);
        weather_reader.ms_outputs;

        // HTF State
        C_csp_solver_htf_1state htf_state;
        htf_state.m_temp = T_htf_inlet;             // [C] Inlet Temp

        C_csp_collector_receiver::S_csp_cr_out_solver cr_out_solver; // Output class

        // OFF
        if (field_mode == 0)
        {
            c_fresnel.off(weather_reader.ms_outputs, htf_state, cr_out_solver, sim_info);
            c_fresnel.converged();
        }
        // STARTUP
        if (field_mode == 2)
        {
            c_fresnel.startup(weather_reader.ms_outputs, htf_state, cr_out_solver, sim_info);
            c_fresnel.converged();
        }
        // ON
        if (field_mode == 3)
        {
            double q_dot_elec_to_CR_heat = 0;   // [MWt]
            double field_control = as_double("defocus_field_control");           // [-] Defocus control (1 is no defocus)

            c_fresnel.on(weather_reader.ms_outputs, htf_state, q_dot_elec_to_CR_heat, field_control,
                cr_out_solver, sim_info);
            c_fresnel.converged();  
        }
        

        // Allocate fresnel outputs
        {
            
            this->assign("EqOpteff", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_EQUIV_OPT_ETA_TOT));
            this->assign("SCAs_def", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_DEFOCUS));

            this->assign("q_inc_sf_tot", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_Q_DOT_INC_SF_TOT));
            this->assign("q_dot_rec_inc", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_Q_DOT_REC_INC));
            this->assign("q_dot_rec_thermal_loss", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_Q_DOT_REC_THERMAL_LOSS));
            this->assign("q_dot_rec_abs", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_Q_DOT_REC_ABS));
            this->assign("rec_thermal_eff", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_REC_THERMAL_EFF));

            this->assign("q_dot_piping_loss", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_Q_DOT_PIPING_LOSS));
            this->assign("e_dot_field_int_energy", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_E_DOT_INTERNAL_ENERGY));
            this->assign("q_dot_htf_sf_out", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_Q_DOT_HTF_OUT));
            this->assign("q_dot_freeze_prot", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_Q_DOT_FREEZE_PROT));

            this->assign("m_dot_loop", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_M_DOT_LOOP));
            this->assign("m_dot_field_recirc", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_M_DOT_FIELD_RECIRC));
            this->assign("m_dot_field_delivered", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_M_DOT_FIELD_DELIVERED));
            this->assign("T_field_cold_in", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_T_FIELD_COLD_IN));
            this->assign("T_rec_cold_in", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_T_REC_COLD_IN));
            this->assign("T_rec_hot_out", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_T_REC_HOT_OUT));
            this->assign("T_field_hot_out", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_T_FIELD_HOT_OUT));
            this->assign("deltaP_field", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_PRESSURE_DROP));

            this->assign("W_dot_sca_track", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_W_DOT_SCA_TRACK));
            this->assign("W_dot_field_pump", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_W_DOT_PUMP));
            this->assign("recirculating", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_IS_RECIRCULATING));

            this->assign("rec_op_mode_final", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_REC_OP_MODE_FINAL));
            this->assign("T_in_loop_final", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_T_IN_LOOP_FINAL) - 273.15);
            this->assign("T_out_loop_final", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_T_OUT_LOOP_FINAL) - 273.15);

            std::vector<double> T_out_scas_last_final_K = c_fresnel.get_scas_outlet_temps();   // [K]
            std::vector<double> T_out_scas_last_final_C;
            for (double temp : T_out_scas_last_final_K)
                T_out_scas_last_final_C.push_back(temp - 273.15);                             // [C]
            
            ssc_number_t* p_T_out_scas_last_final = allocate("T_out_scas_last_final", T_out_scas_last_final_C.size());
            std::copy(T_out_scas_last_final_C.begin(), T_out_scas_last_final_C.end(), p_T_out_scas_last_final);

            this->assign("time_required_su", cr_out_solver.m_time_required_su); //[s]
            this->assign("defocus_final", c_fresnel.mc_reported_outputs.value(C_csp_fresnel_collector_receiver::E_DEFOCUS));

        }

        // Output
        
        

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

DEFINE_MODULE_ENTRY(fresnel_field_subcomponent, "Physical fresnel field subcomponent", 1)
