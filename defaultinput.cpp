/*
 * CLASS
 * Copyright (c) 2010-2013 Meteorology and Air Quality section, Wageningen University and Research centre
 * Copyright (c) 2011-2013 Jordi Vila-Guerau de Arellano
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Bart van Stratum
 * Copyright (c) 2011-2013 Kees van den Dries
 *
 * This file is part of CLASS
 *
 * CLASS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * CLASS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with CLASS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "mainwindow.h"

void MainWindow::readdefaultinput()
{
  // temporary function to fill run form with default values

  // model test
  defaultinput.dt         = 60.;      // time step [s]
  defaultinput.runtime    = 43200.;   // total run time [s]
  defaultinput.sinperiod  = 43200.;   // period for sinusoidal heat fluxes [s]

  // mixed-layer input
  defaultinput.sw_ml      = true;     // mixed-layer model switch
  defaultinput.sw_ftcws   = false;    // compensate FT warming due to subsidence?
  defaultinput.sw_shearwe = false;    // Include shear effect entrainment
  defaultinput.h          = 200.;     // initial ABL height [m]
  defaultinput.Ps         = 101300.;  // surface pressure [Pa]
  defaultinput.omegas     = 0.;       // large scale vertical velocity [m s-1]
  defaultinput.fc         = 1.e-4;    // Coriolis parameter [m s-1]

  defaultinput.theta      = 288.;     // initial mixed-layer potential temperature [K]
  defaultinput.dtheta     = 1.;       // initial temperature jump at h [K]
  defaultinput.gammatheta = 0.006;    // free atmosphere potential temperature lapse rate [K m-1]
  defaultinput.advtheta   = 0.;       // advection of heat [K s-1]
  defaultinput.beta       = 0.2;      // entrainment ratio for virtual heat [-]
  defaultinput.wtheta     = 0.1;      // surface kinematic heat flux [K m s-1]
  defaultinput.sw_wtheta  = false;    // switch for sinusoidal wtheta

  defaultinput.q          = 0.008;    // initial mixed-layer specific humidity [kg kg-1]
  defaultinput.dq         = -0.001;   // initial specific humidity jump at h [kg kg-1]
  defaultinput.gammaq     = 0.;       // free atmosphere specific humidity lapse rate [kg kg-1 m-1]
  defaultinput.advq       = 0.;       // advection of moisture [kg kg-1 s-1]
  defaultinput.wq         = 0.0001;   // surface kinematic moisture flux [kg kg-1 m s-1]
  defaultinput.sw_wq      = false;    // switch for sinusoidal wq

  defaultinput.sw_wind    = false;     // prognostic wind switch
  defaultinput.u          = 6.;       // initial mixed-layer u-wind speed [m s-1]
  defaultinput.du         = 4.;       // initial u-wind jump at h [m s-1]
  defaultinput.gammau     = 0.;       // free atmosphere u-wind speed lapse rate [s-1]
  defaultinput.advu       = 0.;       // advection of u-wind [m s-2]

  defaultinput.v          = -4.0;     // initial mixed-layer u-wind speed [m s-1]
  defaultinput.dv         = 4.0;      // initial u-wind jump at h [m s-1]
  defaultinput.gammav     = 0.;       // free atmosphere v-wind speed lapse rate [s-1]
  defaultinput.advv       = 0.;       // advection of v-wind [m s-2]

  // BvS; a scalar, without the need for the chemistry scheme :)
  defaultinput.sca        = 0.;       // initial mixed-layer scalar [kg kg-1]
  defaultinput.dsca       = 0.;       // initial scalar jump at h [kg kg-1]
  defaultinput.gammasca   = 0.;       // free atmosphere scalar lapse rate [kg kg-1 m-1]
  defaultinput.advsca     = 0.;       // advection of scalar [kg kg-1 s-1]
  defaultinput.wsca       = 0.;       // surface kinematic scalar flux [kg kg-1 m s-1]

  defaultinput.CO2        = 422.;     // initial mixed-layer CO2 [ppm]
  defaultinput.dCO2       = -44.;     // initial CO2 jump at h [ppm]
  defaultinput.gammaCO2   = 0.;       // free atmosphere CO2 lapse rate [ppm]
  defaultinput.advCO2     = 0.;       // advection of CO2 [ppm]
  defaultinput.wCO2       = 0.;       // surface kinematic CO2 flux [ppm]

  // surface layer input
  defaultinput.sw_sl      = false;    // surface layer switch
  defaultinput.ustar      = 0.3;      // surface friction velocity [m s-1]
  defaultinput.z0m        = 0.05;     // roughness length for momentum [m]
  defaultinput.z0h        = 0.01;     // roughness length for scalars [m]
  
  // radiation parameters
  defaultinput.sw_rad     = false;    // radiation switch
  defaultinput.lat        = 51.97;    // latitude [deg]
  defaultinput.lon        = -4.93;    // longitude [deg]
  defaultinput.doy        = 268.;     // day of the year [-]
  defaultinput.tstart     = 6.8;      // time of the day [h UTC]
  defaultinput.cc         = 0.0;      // cloud cover fraction [-]
  defaultinput.Q          = 400.;      // net radiation [W m-2]
  
  // land surface parameters
  defaultinput.sw_ls      = false;    // land surface switch
  defaultinput.sw_jarvis  = true;     // Jarvis / A-Gs switch
  defaultinput.C3C4       = 0;        // C3 or C4 vegetation
  defaultinput.sw_sea     = false;    // land / sea switch
  defaultinput.wg         = 0.21;     // volumetric water content top soil layer [m3 m-3]
  defaultinput.w2         = 0.21;     // volumetric water content deeper soil layer [m3 m-3]
  defaultinput.cveg       = 0.9;      // vegetation fraction [-]
  defaultinput.Tsoil      = 285.;     // temperature top soil layer [K]
  defaultinput.T2         = 286.;     // temperature deeper soil layer [K]
  defaultinput.a          = 0.219;    // Clapp and Hornberger retention curve parameter a
  defaultinput.b          = 4.90;     // Clapp and Hornberger retention curve parameter b
  defaultinput.p          = 4.;       // Clapp and Hornberger retention curve parameter c
  defaultinput.CGsat      = 3.56e-6;  // saturated soil conductivity for heat
  
  defaultinput.wsat       = 0.472;    // saturated volumetric water content ECMWF config [-]
  defaultinput.wfc        = 0.323;    // volumetric water content field capacity [-]
  defaultinput.wwilt      = 0.171;    // volumetric water content wilting point [-]
  
  defaultinput.C1sat      = 0.132;
  defaultinput.C2ref      = 1.8;
  
  defaultinput.LAI        = 2.;       // leaf area index [-]
  defaultinput.gD         = 0.0;      // correction factor transpiration for VPD [-]
  defaultinput.rsmin      = 40.;      // minimum resistance transpiration [s m-1]
  defaultinput.rssoilmin  = 50.;      // minimun resistance soil evaporation [s m-1]
  defaultinput.alpha      = 0.25;     // surface albedo [-]
  
  defaultinput.Ts         = 290.;     // initial surface temperature [K]
  
  defaultinput.Wmax       = 0.0002;   // thickness of water layer on wet vegetation [m]
  defaultinput.Wl         = 0.0000;   // equivalent water layer depth for wet vegetation [m]
  
  defaultinput.Lambda     = 5.9;      // thermal diffusivity skin layer [-]

  // shallow-cumulus
  defaultinput.sw_cu      = true;    // shallow-cumulus switch
  defaultinput.sw_curad   = false;    // link ac -> cc -> radiation

  defaultinput.sw_plume   = true;    // calculate plume statistics
  defaultinput.sw_cin     = false;    // use CIN value from plume model for cumulus inhibition
  defaultinput.sw_ft_storage = false;
  // stratocumulus
  defaultinput.dFz        = 0.;       // Cloud top radiative divergence (stratocumulus)

  defaultinput.phi_cu     = 1.;
  defaultinput.wcld_fact  = 1.;
  defaultinput.hstore     = 3e3;
  //chemistry
  defaultinput.sw_chem    = false;
  defaultinput.sw_chem_constant  = true;
  defaultinput.sw_photo_constant = true;
  defaultinput.nsc        = 22;
  defaultinput.reactions  = defaultreactions; // CvH copy addresses, all runs will have therefore the same struct as reference

  defaultinput.rsize      = 28;
  defaultinput.csize      = 22;

  defaultinput.P_ref      = 100000;
  defaultinput.Tcbl_ref   = 298.;
  defaultinput.Tfc_ref    = 298.;
  defaultinput.qcbl_ref   = 0.010;
  defaultinput.qfc_ref    = 0.010;
  defaultinput.tod_ref    = 12;
  defaultinput.stocoef    = 0.;

  defaultinput.sc[0]      = 0.;
  defaultinput.dsc[0]     = 0.;
  defaultinput.gammasc[0] = 0.;
  defaultinput.advsc[0]   = 0.;
  defaultinput.wsc[0]     = 1.0;
  defaultinput.sw_wsc[0]  = 0;

  defaultinput.sc[1]      = 30.;
  defaultinput.dsc[1]     = 0.;
  defaultinput.gammasc[1] = 0.;
  defaultinput.advsc[1]   = 0.;
  defaultinput.wsc[1]     = 0.;
  defaultinput.sw_wsc[1]  = 0;

  defaultinput.sc[2]      = 0.;
  defaultinput.dsc[2]     = 0.;
  defaultinput.gammasc[2] = 0.;
  defaultinput.advsc[2]   = 0.;
  defaultinput.wsc[2]     = 0.;
  defaultinput.sw_wsc[2]  = 0;

  defaultinput.sc[3]      = 0.4;
  defaultinput.dsc[3]     = 0.;
  defaultinput.gammasc[3] = 0.;
  defaultinput.advsc[3]   = 0.;
  defaultinput.wsc[3]     = 0.;
  defaultinput.sw_wsc[3]  = 0;

  defaultinput.sc[4]      = 0.6;
  defaultinput.dsc[4]     = 0.;
  defaultinput.gammasc[4] = 0.;
  defaultinput.advsc[4]   = 0.;
  defaultinput.wsc[4]     = 0;
  defaultinput.sw_wsc[4]  = 0;

  defaultinput.sc[5]      = 1724.;
  defaultinput.dsc[5]     = 0.;
  defaultinput.gammasc[5] = 0.;
  defaultinput.advsc[5]   = 0.;
  defaultinput.wsc[5]     = 0.;
  defaultinput.sw_wsc[5]  = 0;

  defaultinput.sc[6]      = 0.;
  defaultinput.dsc[6]     = 0.;
  defaultinput.gammasc[6] = 0.;
  defaultinput.advsc[6]   = 0.;
  defaultinput.wsc[6]     = 0.;
  defaultinput.sw_wsc[6]  = 0;

  defaultinput.sc[7]      = 0.;
  defaultinput.dsc[7]     = 0.;
  defaultinput.gammasc[7] = 0.;
  defaultinput.advsc[7]   = 0.;
  defaultinput.wsc[7]     = 0.;
  defaultinput.sw_wsc[7]  = 0;

  defaultinput.sc[8]      = 0.;
  defaultinput.dsc[8]     = 0.;
  defaultinput.gammasc[8] = 0.;
  defaultinput.advsc[8]   = 0.;
  defaultinput.wsc[8]     = 0.;
  defaultinput.sw_wsc[8]  = 0;

  defaultinput.sc[9]      = 0.;
  defaultinput.dsc[9]     = 0.;
  defaultinput.gammasc[9] = 0.;
  defaultinput.advsc[9]   = 0.;
  defaultinput.wsc[9]     = .5;
  defaultinput.sw_wsc[9]  = 0;

  defaultinput.sc[10]      = 0.;
  defaultinput.dsc[10]     = 0.;
  defaultinput.gammasc[10] = 0.;
  defaultinput.advsc[10]   = 0.;
  defaultinput.wsc[10]     = 0.;
  defaultinput.sw_wsc[10]  = 0;

  defaultinput.sc[11]      = 0.;
  defaultinput.dsc[11]     = 0.;
  defaultinput.gammasc[11] = 0.;
  defaultinput.advsc[11]   = 0.;
  defaultinput.wsc[11]     = 0.;
  defaultinput.sw_wsc[11]  = 0;

  defaultinput.sc[12]      = 0.;
  defaultinput.dsc[12]     = 0.;
  defaultinput.gammasc[12] = 0.;
  defaultinput.advsc[12]   = 0.;
  defaultinput.wsc[12]     = 0.;
  defaultinput.sw_wsc[12]  = 0;

  defaultinput.sc[13]      = 100.;
  defaultinput.dsc[13]     = 0.;
  defaultinput.gammasc[13] = 0.;
  defaultinput.advsc[13]   = 0.;
  defaultinput.wsc[13]     = 0.;
  defaultinput.sw_wsc[13]  = 0;

  defaultinput.sc[14]      = 0.;
  defaultinput.dsc[14]     = 0.;
  defaultinput.gammasc[14] = 0.;
  defaultinput.advsc[14]   = 0.;
  defaultinput.wsc[14]     = 0.;
  defaultinput.sw_wsc[14]  = 0;

  defaultinput.sc[15]      = 0.;
  defaultinput.dsc[15]     = 0.;
  defaultinput.gammasc[15] = 0.;
  defaultinput.advsc[15]   = 0.;
  defaultinput.wsc[15]     = 0.;
  defaultinput.sw_wsc[15]  = 0;

  defaultinput.sc[16]      = 0.2e9;
  defaultinput.dsc[16]     = 0.;
  defaultinput.gammasc[16] = 0.;
  defaultinput.advsc[16]   = 0.;
  defaultinput.wsc[16]     = 0.;
  defaultinput.sw_wsc[16]  = 0;

  defaultinput.sc[17]      = 0.8e9;
  defaultinput.dsc[17]     = 0.;
  defaultinput.gammasc[17] = 0.;
  defaultinput.advsc[17]   = 0.;
  defaultinput.wsc[17]     = 0.;
  defaultinput.sw_wsc[17]  = 0;

  defaultinput.sc[18]      = 0.;
  defaultinput.dsc[18]     = 0.;
  defaultinput.gammasc[18] = 0.;
  defaultinput.advsc[18]   = 0.;
  defaultinput.wsc[18]     = 0.;
  defaultinput.sw_wsc[18]  = 0;

  defaultinput.sc[19]      = 0.;
  defaultinput.dsc[19]     = 0.;
  defaultinput.gammasc[19] = 0.;
  defaultinput.advsc[19]   = 0.;
  defaultinput.wsc[19]     = 0.;
  defaultinput.sw_wsc[19]  = 0;

  defaultinput.sc[20]      = 0.;
  defaultinput.dsc[20]     = 0.;
  defaultinput.gammasc[20] = 0.;
  defaultinput.advsc[20]   = 0.;
  defaultinput.wsc[20]     = 0.;
  defaultinput.sw_wsc[20]  = 0;

  defaultinput.sc[21]      = 0.;
  defaultinput.dsc[21]     = 0.;
  defaultinput.gammasc[21] = 0.;
  defaultinput.advsc[21]   = 0.;
  defaultinput.wsc[21]     = 0.;
  defaultinput.sw_wsc[21]  = 0;


  // TEMP FOR CHEMISTRY
  defaultreactions[0].rname           = "r01";
  defaultreactions[0].RadDep          = 1;
  defaultreactions[0].func1           = 2;
  defaultreactions[0].nr_chem_inp     = 1;
  defaultreactions[0].nr_chem_outp    = 2;
  defaultreactions[0].inp[0].coef     = 1.00;
  defaultreactions[0].inp[0].cname    = "O3";
  defaultreactions[0].inp[0].chem_nr  = 1;
  defaultreactions[0].inp[0].index    = 1;
  defaultreactions[0].outp[0].coef    = 1.00;
  defaultreactions[0].outp[0].cname   = "O1D";
  defaultreactions[0].outp[0].chem_nr = 2;
  defaultreactions[0].outp[0].index   = 2;
  defaultreactions[0].outp[1].coef    = 1.00;
  defaultreactions[0].outp[1].cname   = "O2";
  defaultreactions[0].outp[1].chem_nr = 16;
  defaultreactions[0].outp[1].index   = 16;
  defaultreactions[0].A               = +3.171e-4;
  defaultreactions[0].B               = -1.809;
  defaultreactions[0].C               = +1.200e+00;
  defaultreactions[0].D               = +1.000e+00;
  defaultreactions[0].E               = +1.000e+00;
  defaultreactions[0].F               = +1.000e+00;
  defaultreactions[0].G               = +1.000e+00;

  defaultreactions[1].rname           = "r02";
  defaultreactions[1].RadDep          = 0;
  defaultreactions[1].func1           = 2;
  defaultreactions[1].nr_chem_inp     = 2;
  defaultreactions[1].nr_chem_outp    = 1;
  defaultreactions[1].inp[0].coef     = 1.00;
  defaultreactions[1].inp[0].cname    = "O1D";
  defaultreactions[1].inp[0].chem_nr  = 2;
  defaultreactions[1].inp[0].index    = 2;
  defaultreactions[1].inp[1].coef     = 1.00;
  defaultreactions[1].inp[1].cname    = "H2O";
  defaultreactions[1].inp[1].chem_nr  = 14;
  defaultreactions[1].inp[1].index    = 14;
  defaultreactions[1].outp[0].coef    = 2.00;
  defaultreactions[1].outp[0].cname   = "OH";
  defaultreactions[1].outp[0].chem_nr = 11;
  defaultreactions[1].outp[0].index   = 11;
  defaultreactions[1].A               = +1.630e-10;
  defaultreactions[1].B               = +6.000e+01;
  defaultreactions[1].C               = +1.000e+00;
  defaultreactions[1].D               = +1.000e+00;
  defaultreactions[1].E               = +1.000e+00;
  defaultreactions[1].F               = +1.000e+00;
  defaultreactions[1].G               = +1.000e+00;

  defaultreactions[2].rname           = "r03";
  defaultreactions[2].RadDep          = 0;
  defaultreactions[2].func1           = 2;
  defaultreactions[2].nr_chem_inp     = 2;
  defaultreactions[2].nr_chem_outp    = 1;
  defaultreactions[2].inp[0].coef     = 1.00;
  defaultreactions[2].inp[0].cname    = "O1D";
  defaultreactions[2].inp[0].chem_nr  = 2;
  defaultreactions[2].inp[0].index    = 2;
  defaultreactions[2].inp[1].coef     = 1.00;
  defaultreactions[2].inp[1].cname    = "N2";
  defaultreactions[2].inp[1].chem_nr  = 17;
  defaultreactions[2].inp[1].index    = 17;
  defaultreactions[2].outp[0].coef    = 1.00;
  defaultreactions[2].outp[0].cname   = "O3";
  defaultreactions[2].outp[0].chem_nr = 1;
  defaultreactions[2].outp[0].index   = 1;
  defaultreactions[2].A               = +2.150e-11;
  defaultreactions[2].B               = +1.100e+02;
  defaultreactions[2].C               = +1.000e+00;
  defaultreactions[2].D               = +1.000e+00;
  defaultreactions[2].E               = +1.000e+00;
  defaultreactions[2].F               = +1.000e+00;
  defaultreactions[2].G               = +1.000e+00;

  defaultreactions[3].rname           = "r04";
  defaultreactions[3].RadDep          = 0;
  defaultreactions[3].func1           = 2;
  defaultreactions[3].nr_chem_inp     = 2;
  defaultreactions[3].nr_chem_outp    = 1;
  defaultreactions[3].inp[0].coef     = 1.00;
  defaultreactions[3].inp[0].cname    = "O1D";
  defaultreactions[3].inp[0].chem_nr  = 2;
  defaultreactions[3].inp[0].index    = 2;
  defaultreactions[3].inp[1].coef     = 1.00;
  defaultreactions[3].inp[1].cname    = "O2";
  defaultreactions[3].inp[1].chem_nr  = 16;
  defaultreactions[3].inp[1].index    = 16;
  defaultreactions[3].outp[0].coef    = 1.00;
  defaultreactions[3].outp[0].cname   = "O3";
  defaultreactions[3].outp[0].chem_nr = 1;
  defaultreactions[3].outp[0].index   = 1;
  defaultreactions[3].A               = +3.300e-11;
  defaultreactions[3].B               = +5.500e+01;
  defaultreactions[3].C               = +1.000e+00;
  defaultreactions[3].D               = +1.000e+00;
  defaultreactions[3].E               = +1.000e+00;
  defaultreactions[3].F               = +1.000e+00;
  defaultreactions[3].G               = +1.000e+00;

  defaultreactions[4].rname           = "r05";
  defaultreactions[4].RadDep          = 1;
  defaultreactions[4].func1           = 2;
  defaultreactions[4].nr_chem_inp     = 1;
  defaultreactions[4].nr_chem_outp    = 2;
  defaultreactions[4].inp[0].coef     = 1.00;
  defaultreactions[4].inp[0].cname    = "NO2";
  defaultreactions[4].inp[0].chem_nr  = 4;
  defaultreactions[4].inp[0].index    = 4;
  defaultreactions[4].outp[0].coef    = 1.00;
  defaultreactions[4].outp[0].cname   = "NO";
  defaultreactions[4].outp[0].chem_nr = 3;
  defaultreactions[4].outp[0].index   = 3;
  defaultreactions[4].outp[1].coef    = 1.00;
  defaultreactions[4].outp[1].cname   = "O3";
  defaultreactions[4].outp[1].chem_nr = 1;
  defaultreactions[4].outp[1].index   = 1;
  defaultreactions[4].A               = +1.930e-02;
  defaultreactions[4].B               = -5.750e-01;
  defaultreactions[4].C               = +1.000e+00;
  defaultreactions[4].D               = +1.000e+00;
  defaultreactions[4].E               = +1.000e+00;
  defaultreactions[4].F               = +1.000e+00;
  defaultreactions[4].G               = +1.000e+00;

  defaultreactions[5].rname           = "r06";
  defaultreactions[5].RadDep          = 1;
  defaultreactions[5].func1           = 2;
  defaultreactions[5].nr_chem_inp     = 1;
  defaultreactions[5].nr_chem_outp    = 1;
  defaultreactions[5].inp[0].coef     = 1.00;
  defaultreactions[5].inp[0].cname    = "CH2O";
  defaultreactions[5].inp[0].chem_nr  = 6;
  defaultreactions[5].inp[0].index    = 6;
  defaultreactions[5].outp[0].coef    = 1.00;
  defaultreactions[5].outp[0].cname   = "HO2";
  defaultreactions[5].outp[0].chem_nr = 12;
  defaultreactions[5].outp[0].index   = 12;
  defaultreactions[5].A               = +2.16e-4;
  defaultreactions[5].B               = -7.78e-01;
  defaultreactions[5].C               = +1.000e+00;
  defaultreactions[5].D               = +1.000e+00;
  defaultreactions[5].E               = +1.000e+00;
  defaultreactions[5].F               = +1.000e+00;
  defaultreactions[5].G               = +1.000e+00;

  defaultreactions[6].rname           = "r07";
  defaultreactions[6].RadDep          = 0;
  defaultreactions[6].func1           = 1;
  defaultreactions[6].nr_chem_inp     = 2;
  defaultreactions[6].nr_chem_outp    = 1;
  defaultreactions[6].inp[0].coef     = 1.00;
  defaultreactions[6].inp[0].cname    = "OH";
  defaultreactions[6].inp[0].chem_nr  = 11;
  defaultreactions[6].inp[0].index    = 11;
  defaultreactions[6].inp[1].coef     = 1.00;
  defaultreactions[6].inp[1].cname    = "CO";
  defaultreactions[6].inp[1].chem_nr  = 13;
  defaultreactions[6].inp[1].index    = 13;
  defaultreactions[6].outp[0].coef    = 1.00;
  defaultreactions[6].outp[0].cname   = "HO2";
  defaultreactions[6].outp[0].chem_nr = 12;
  defaultreactions[6].outp[0].index   = 12;
  defaultreactions[6].A               = +2.400e-13;
  defaultreactions[6].B               = +1.000e+00;
  defaultreactions[6].C               = +1.000e+00;
  defaultreactions[6].D               = +1.000e+00;
  defaultreactions[6].E               = +1.000e+00;
  defaultreactions[6].F               = +1.000e+00;
  defaultreactions[6].G               = +1.000e+00;

  defaultreactions[7].rname           = "r08";
  defaultreactions[7].RadDep          = 0;
  defaultreactions[7].func1           = 2;
  defaultreactions[7].nr_chem_inp     = 2;
  defaultreactions[7].nr_chem_outp    = 1;
  defaultreactions[7].inp[0].coef     = 1.00;
  defaultreactions[7].inp[0].cname    = "OH";
  defaultreactions[7].inp[0].chem_nr  = 11;
  defaultreactions[7].inp[0].index    = 11;
  defaultreactions[7].inp[1].coef     = 1.00;
  defaultreactions[7].inp[1].cname    = "CH4";
  defaultreactions[7].inp[1].chem_nr  = 5;
  defaultreactions[7].inp[1].index    = 5;
  defaultreactions[7].outp[0].coef    = 1.00;
  defaultreactions[7].outp[0].cname   = "CH3O2";
  defaultreactions[7].outp[0].chem_nr = 7;
  defaultreactions[7].outp[0].index   = 7;
  defaultreactions[7].A               = +2.450e-12;
  defaultreactions[7].B               = -1.775e+03;
  defaultreactions[7].C               = +1.000e+00;
  defaultreactions[7].D               = +1.000e+00;
  defaultreactions[7].E               = +1.000e+00;
  defaultreactions[7].F               = +1.000e+00;
  defaultreactions[7].G               = +1.000e+00;

  defaultreactions[8].rname           = "r09";
  defaultreactions[8].func1           = 1;
  defaultreactions[8].nr_chem_inp     = 2;
  defaultreactions[8].nr_chem_outp    = 1;
  defaultreactions[8].inp[0].coef     = 1.00;
  defaultreactions[8].inp[0].cname    = "OH";
  defaultreactions[8].inp[0].chem_nr  = 11;
  defaultreactions[8].inp[0].index    = 11;
  defaultreactions[8].inp[1].coef     = 1.00;
  defaultreactions[8].inp[1].cname    = "ISO";
  defaultreactions[8].inp[1].chem_nr  = 9;
  defaultreactions[8].inp[1].index    = 9;
  defaultreactions[8].outp[0].coef    = 1.00;
  defaultreactions[8].outp[0].cname   = "RO2";
  defaultreactions[8].outp[0].chem_nr = 10;
  defaultreactions[8].outp[0].index   = 10;
  defaultreactions[8].A               = +1.000e-10;
  defaultreactions[8].B               = +1.000e+00;
  defaultreactions[8].C               = +1.000e+00;
  defaultreactions[8].D               = +1.000e+00;
  defaultreactions[8].E               = +1.000e+00;
  defaultreactions[8].F               = +1.000e+00;
  defaultreactions[8].G               = +1.000e+00;

  defaultreactions[9].rname           = "r10";
  defaultreactions[9].RadDep          = 0;
  defaultreactions[9].func1           = 1;
  defaultreactions[9].nr_chem_inp     = 2;
  defaultreactions[9].nr_chem_outp    = 2;
  defaultreactions[9].inp[0].coef     = 1.00;
  defaultreactions[9].inp[0].cname    = "OH";
  defaultreactions[9].inp[0].chem_nr  = 11;
  defaultreactions[9].inp[0].index    = 11;
  defaultreactions[9].inp[1].coef     = 1.00;
  defaultreactions[9].inp[1].cname    = "MVK";
  defaultreactions[9].inp[1].chem_nr  = 8;
  defaultreactions[9].inp[1].index    = 8;
  defaultreactions[9].outp[0].coef    = 1.00;
  defaultreactions[9].outp[0].cname   = "HO2";
  defaultreactions[9].outp[0].chem_nr = 12;
  defaultreactions[9].outp[0].index   = 12;
  defaultreactions[9].outp[1].coef    = 1.00;
  defaultreactions[9].outp[1].cname   = "CH2O";
  defaultreactions[9].outp[1].chem_nr = 6;
  defaultreactions[9].outp[1].index   = 6;
  defaultreactions[9].A               = +2.400e-11;
  defaultreactions[9].B               = +1.000e+00;
  defaultreactions[9].C               = +1.000e+00;
  defaultreactions[9].D               = +1.000e+00;
  defaultreactions[9].E               = +1.000e+00;
  defaultreactions[9].F               = +1.000e+00;
  defaultreactions[9].G               = +1.000e+00;

  defaultreactions[10].rname           = "r11";
  defaultreactions[10].RadDep          = 0;
  defaultreactions[10].func1           = 2;
  defaultreactions[10].nr_chem_inp     = 2;
  defaultreactions[10].nr_chem_outp    = 2;
  defaultreactions[10].inp[0].coef     = 1.00;
  defaultreactions[10].inp[0].cname    = "OH";
  defaultreactions[10].inp[0].chem_nr  = 11;
  defaultreactions[10].inp[0].index    = 11;
  defaultreactions[10].inp[1].coef     = 1.00;
  defaultreactions[10].inp[1].cname    = "HO2";
  defaultreactions[10].inp[1].chem_nr  = 12;
  defaultreactions[10].inp[1].index    = 12;
  defaultreactions[10].outp[0].coef    = 1.00;
  defaultreactions[10].outp[0].cname   = "H2O";
  defaultreactions[10].outp[0].chem_nr = 14;
  defaultreactions[10].outp[0].index   = 14;
  defaultreactions[10].outp[1].coef    = 1.00;
  defaultreactions[10].outp[1].cname   = "O2";
  defaultreactions[10].outp[1].chem_nr = 16;
  defaultreactions[10].outp[1].index   = 16;
  defaultreactions[10].A               = +4.800e-11;
  defaultreactions[10].B               = +2.500e+02;
  defaultreactions[10].C               = +1.000e+00;
  defaultreactions[10].D               = +1.000e+00;
  defaultreactions[10].E               = +1.000e+00;
  defaultreactions[10].F               = +1.000e+00;
  defaultreactions[10].G               = +1.000e+00;

  defaultreactions[11].rname           = "r12";
  defaultreactions[11].RadDep          = 0;
  defaultreactions[11].func1           = 2;
  defaultreactions[11].nr_chem_inp     = 2;
  defaultreactions[11].nr_chem_outp    = 2;
  defaultreactions[11].inp[0].coef     = 1.00;
  defaultreactions[11].inp[0].cname    = "OH";
  defaultreactions[11].inp[0].chem_nr  = 11;
  defaultreactions[11].inp[0].index    = 11;
  defaultreactions[11].inp[1].coef     = 1.00;
  defaultreactions[11].inp[1].cname    = "H2O2";
  defaultreactions[11].inp[1].chem_nr  = 19;
  defaultreactions[11].inp[1].index    = 19;
  defaultreactions[11].outp[0].coef    = 1.00;
  defaultreactions[11].outp[0].cname   = "H2O";
  defaultreactions[11].outp[0].chem_nr = 14;
  defaultreactions[11].outp[0].index   = 14;
  defaultreactions[11].outp[1].coef    = 1.00;
  defaultreactions[11].outp[1].cname   = "HO2";
  defaultreactions[11].outp[1].chem_nr = 12;
  defaultreactions[11].outp[1].index   = 12;
  defaultreactions[11].A               = +2.900e-12;
  defaultreactions[11].B               = -1.600e+02;
  defaultreactions[11].C               = +1.000e+00;
  defaultreactions[11].D               = +1.000e+00;
  defaultreactions[11].E               = +1.000e+00;
  defaultreactions[11].F               = +1.000e+00;
  defaultreactions[11].G               = +1.000e+00;

  defaultreactions[12].rname           = "r13";
  defaultreactions[12].RadDep          = 0;
  defaultreactions[12].func1           = 2;
  defaultreactions[12].nr_chem_inp     = 2;
  defaultreactions[12].nr_chem_outp    = 2;
  defaultreactions[12].inp[0].coef     = 1.00;
  defaultreactions[12].inp[0].cname    = "HO2";
  defaultreactions[12].inp[0].chem_nr  = 12;
  defaultreactions[12].inp[0].index    = 12;
  defaultreactions[12].inp[1].coef     = 1.00;
  defaultreactions[12].inp[1].cname    = "NO";
  defaultreactions[12].inp[1].chem_nr  = 3;
  defaultreactions[12].inp[1].index    = 3;
  defaultreactions[12].outp[0].coef    = 1.00;
  defaultreactions[12].outp[0].cname   = "OH";
  defaultreactions[12].outp[0].chem_nr = 11;
  defaultreactions[12].outp[0].index   = 11;
  defaultreactions[12].outp[1].coef    = 1.00;
  defaultreactions[12].outp[1].cname   = "NO2";
  defaultreactions[12].outp[1].chem_nr = 4;
  defaultreactions[12].outp[1].index   = 4;
  defaultreactions[12].A               = +3.500e-12;
  defaultreactions[12].B               = +2.500e+02;
  defaultreactions[12].C               = +1.000e+00;
  defaultreactions[12].D               = +1.000e+00;
  defaultreactions[12].E               = +1.000e+00;
  defaultreactions[12].F               = +1.000e+00;
  defaultreactions[12].G               = +1.000e+00;

  defaultreactions[13].rname           = "r14";
  defaultreactions[13].RadDep          = 0;
  defaultreactions[13].func1           = 2;
  defaultreactions[13].nr_chem_inp     = 2;
  defaultreactions[13].nr_chem_outp    = 3;
  defaultreactions[13].inp[0].coef     = 1.00;
  defaultreactions[13].inp[0].cname    = "CH3O2";
  defaultreactions[13].inp[0].chem_nr  = 7;
  defaultreactions[13].inp[0].index    = 7;
  defaultreactions[13].inp[1].coef     = 1.00;
  defaultreactions[13].inp[1].cname    = "NO";
  defaultreactions[13].inp[1].chem_nr  = 3;
  defaultreactions[13].inp[1].index    = 3;
  defaultreactions[13].outp[0].coef    = 1.00;
  defaultreactions[13].outp[0].cname   = "HO2";
  defaultreactions[13].outp[0].chem_nr = 12;
  defaultreactions[13].outp[0].index   = 12;
  defaultreactions[13].outp[1].coef    = 1.00;
  defaultreactions[13].outp[1].cname   = "NO2";
  defaultreactions[13].outp[1].chem_nr = 4;
  defaultreactions[13].outp[1].index   = 4;
  defaultreactions[13].outp[2].coef    = 1.00;
  defaultreactions[13].outp[2].cname   = "CH2O";
  defaultreactions[13].outp[2].chem_nr = 6;
  defaultreactions[13].outp[2].index   = 6;
  defaultreactions[13].A               = +2.800e-12;
  defaultreactions[13].B               = +3.000e+02;
  defaultreactions[13].C               = +1.000e+00;
  defaultreactions[13].D               = +1.000e+00;
  defaultreactions[13].E               = +1.000e+00;
  defaultreactions[13].F               = +1.000e+00;
  defaultreactions[13].G               = +1.000e+00;

  defaultreactions[14].rname           = "r15";
  defaultreactions[14].RadDep          = 0;
  defaultreactions[14].func1           = 1;
  defaultreactions[14].nr_chem_inp     = 2;
  defaultreactions[14].nr_chem_outp    = 4;
  defaultreactions[14].inp[0].coef     = 1.00;
  defaultreactions[14].inp[0].cname    = "RO2";
  defaultreactions[14].inp[0].chem_nr  = 10;
  defaultreactions[14].inp[0].index    = 10;
  defaultreactions[14].inp[1].coef     = 1.00;
  defaultreactions[14].inp[1].cname    = "NO";
  defaultreactions[14].inp[1].chem_nr  = 3;
  defaultreactions[14].inp[1].index    = 3;
  defaultreactions[14].outp[0].coef    = 1.00;
  defaultreactions[14].outp[0].cname   = "HO2";
  defaultreactions[14].outp[0].chem_nr = 12;
  defaultreactions[14].outp[0].index   = 12;
  defaultreactions[14].outp[1].coef    = 1.00;
  defaultreactions[14].outp[1].cname   = "NO2";
  defaultreactions[14].outp[1].chem_nr = 4;
  defaultreactions[14].outp[1].index   = 4;
  defaultreactions[14].outp[2].coef    = 1.00;
  defaultreactions[14].outp[2].cname   = "MVK";
  defaultreactions[14].outp[2].chem_nr = 8;
  defaultreactions[14].outp[2].index   = 8;
  defaultreactions[14].outp[3].coef    = 1.00;
  defaultreactions[14].outp[3].cname   = "CH2O";
  defaultreactions[14].outp[3].chem_nr = 6;
  defaultreactions[14].outp[3].index   = 6;
  defaultreactions[14].A               = +1.000e-11;
  defaultreactions[14].B               = +1.000e+00;
  defaultreactions[14].C               = +1.000e+00;
  defaultreactions[14].D               = +1.000e+00;
  defaultreactions[14].E               = +1.000e+00;
  defaultreactions[14].F               = +1.000e+00;
  defaultreactions[14].G               = +1.000e+00;

  defaultreactions[15].rname           = "r16";
  defaultreactions[15].RadDep          = 0;
  defaultreactions[15].func1           = 2;
  defaultreactions[15].nr_chem_inp     = 2;
  defaultreactions[15].nr_chem_outp    = 2;
  defaultreactions[15].inp[0].coef     = 1.00;
  defaultreactions[15].inp[0].cname    = "OH";
  defaultreactions[15].inp[0].chem_nr  = 11;
  defaultreactions[15].inp[0].index    = 11;
  defaultreactions[15].inp[1].coef     = 1.00;
  defaultreactions[15].inp[1].cname    = "CH2O";
  defaultreactions[15].inp[1].chem_nr  = 6;
  defaultreactions[15].inp[1].index    = 6;
  defaultreactions[15].outp[0].coef    = 1.00;
  defaultreactions[15].outp[0].cname   = "HO2";
  defaultreactions[15].outp[0].chem_nr = 12;
  defaultreactions[15].outp[0].index   = 12;
  defaultreactions[15].outp[1].coef    = 1.00;
  defaultreactions[15].outp[1].cname   = "CO";
  defaultreactions[15].outp[1].chem_nr = 13;
  defaultreactions[15].outp[1].index   = 13;
  defaultreactions[15].A               = +5.500e-12;
  defaultreactions[15].B               = +1.250e+02;
  defaultreactions[15].C               = +1.000e+00;
  defaultreactions[15].D               = +1.000e+00;
  defaultreactions[15].E               = +1.000e+00;
  defaultreactions[15].F               = +1.000e+00;
  defaultreactions[15].G               = +1.000e+00;

  defaultreactions[16].rname           = "r17";
  defaultreactions[16].RadDep          = 0;
  defaultreactions[16].func1           = 6;
  defaultreactions[16].nr_chem_inp     = 1;
  defaultreactions[16].nr_chem_outp    = 1;
  defaultreactions[16].inp[0].coef     = 2.00;
  defaultreactions[16].inp[0].cname    = "HO2";
  defaultreactions[16].inp[0].chem_nr  = 12;
  defaultreactions[16].inp[0].index    = 12;
  defaultreactions[16].outp[0].coef    = 1.00;
  defaultreactions[16].outp[0].cname   = "H2O2";
  defaultreactions[16].outp[0].chem_nr = 19;
  defaultreactions[16].outp[0].index   = 19;
  defaultreactions[16].A               = +2.200e-13;
  defaultreactions[16].B               = +6.000e+02;
  defaultreactions[16].C               = +1.900e-33;
  defaultreactions[16].D               = +9.800e+02;
  defaultreactions[16].E               = +1.400e-21;
  defaultreactions[16].F               = +2.200e+03;
  defaultreactions[16].G               = +1.000e+00;

  defaultreactions[17].rname           = "r18";
  defaultreactions[17].RadDep          = 0;
  defaultreactions[17].func1           = 2;
  defaultreactions[17].nr_chem_inp     = 2;
  defaultreactions[17].nr_chem_outp    = 1;
  defaultreactions[17].inp[0].coef     = 1.00;
  defaultreactions[17].inp[0].cname    = "CH3O2";
  defaultreactions[17].inp[0].chem_nr  = 7;
  defaultreactions[17].inp[0].index    = 7;
  defaultreactions[17].inp[1].coef     = 1.00;
  defaultreactions[17].inp[1].cname    = "HO2";
  defaultreactions[17].inp[1].chem_nr  = 12;
  defaultreactions[17].inp[1].index    = 12;
  defaultreactions[17].outp[0].coef    = 1.00;
  defaultreactions[17].outp[0].cname   = "Product";
  defaultreactions[17].outp[0].chem_nr = 15;
  defaultreactions[17].outp[0].index   = 15;
  defaultreactions[17].A               = +4.100e-13;
  defaultreactions[17].B               = +7.500e+02;
  defaultreactions[17].C               = +1.000e+00;
  defaultreactions[17].D               = +1.000e+00;
  defaultreactions[17].E               = +1.000e+00;
  defaultreactions[17].F               = +1.000e+00;
  defaultreactions[17].G               = +1.000e+00;

  defaultreactions[18].rname           = "r19";
  defaultreactions[18].RadDep          = 0;
  defaultreactions[18].func1           = 1;
  defaultreactions[18].nr_chem_inp     = 2;
  defaultreactions[18].nr_chem_outp    = 2;
  defaultreactions[18].inp[0].coef     = 1.00;
  defaultreactions[18].inp[0].cname    = "RO2";
  defaultreactions[18].inp[0].chem_nr  = 10;
  defaultreactions[18].inp[0].index    = 10;
  defaultreactions[18].inp[1].coef     = 1.00;
  defaultreactions[18].inp[1].cname    = "HO2";
  defaultreactions[18].inp[1].chem_nr  = 12;
  defaultreactions[18].inp[1].index    = 12;
  defaultreactions[18].outp[0].coef    = 0.00;
  defaultreactions[18].outp[0].cname   = "OH";
  defaultreactions[18].outp[0].chem_nr = 11;
  defaultreactions[18].outp[0].index   = 11;
  defaultreactions[18].outp[1].coef    = 1.00;
  defaultreactions[18].outp[1].cname   = "Product";
  defaultreactions[18].outp[1].chem_nr = 15;
  defaultreactions[18].outp[1].index   = 15;
  defaultreactions[18].A               = +1.500e-11;
  defaultreactions[18].B               = +7.500e+02;
  defaultreactions[18].C               = +1.000e+00;
  defaultreactions[18].D               = +1.000e+00;
  defaultreactions[18].E               = +1.000e+00;
  defaultreactions[18].F               = +1.000e+00;
  defaultreactions[18].G               = +1.000e+00;

  defaultreactions[19].rname           = "r20";
  defaultreactions[19].RadDep          = 0;
  defaultreactions[19].func1           = 2;
  defaultreactions[19].nr_chem_inp     = 2;
  defaultreactions[19].nr_chem_outp    = 1;
  defaultreactions[19].inp[0].coef     = 1.00;
  defaultreactions[19].inp[0].cname    = "OH";
  defaultreactions[19].inp[0].chem_nr  = 11;
  defaultreactions[19].inp[0].index    = 11;
  defaultreactions[19].inp[1].coef     = 1.00;
  defaultreactions[19].inp[1].cname    = "NO2";
  defaultreactions[19].inp[1].chem_nr  = 4;
  defaultreactions[19].inp[1].index    = 4;
  defaultreactions[19].outp[0].coef    = 1.00;
  defaultreactions[19].outp[0].cname   = "HNO3";
  defaultreactions[19].outp[0].chem_nr = 18;
  defaultreactions[19].outp[0].index   = 18;
  defaultreactions[19].A               = +3.500e-12;
  defaultreactions[19].B               = +3.400e+02;
  defaultreactions[19].C               = +1.000e+00;
  defaultreactions[19].D               = +1.000e+00;
  defaultreactions[19].E               = +1.000e+00;
  defaultreactions[19].F               = +1.000e+00;
  defaultreactions[19].G               = +1.000e+00;

  defaultreactions[20].rname           = "r21";
  defaultreactions[20].RadDep          = 0;
  defaultreactions[20].func1           = 2;
  defaultreactions[20].nr_chem_inp     = 2;
  defaultreactions[20].nr_chem_outp    = 1;
  defaultreactions[20].inp[0].coef     = 1.00;
  defaultreactions[20].inp[0].cname    = "NO";
  defaultreactions[20].inp[0].chem_nr  = 3;
  defaultreactions[20].inp[0].index    = 3;
  defaultreactions[20].inp[1].coef     = 1.00;
  defaultreactions[20].inp[1].cname    = "O3";
  defaultreactions[20].inp[1].chem_nr  = 1;
  defaultreactions[20].inp[1].index    = 1;
  defaultreactions[20].outp[0].coef    = 1.00;
  defaultreactions[20].outp[0].cname   = "NO2";
  defaultreactions[20].outp[0].chem_nr = 4;
  defaultreactions[20].outp[0].index   = 4;
  defaultreactions[20].A               = +3.000e-12;
  defaultreactions[20].B               = -1.500e+03;
  defaultreactions[20].C               = +1.000e+00;
  defaultreactions[20].D               = +1.000e+00;
  defaultreactions[20].E               = +1.000e+00;
  defaultreactions[20].F               = +1.000e+00;
  defaultreactions[20].G               = +1.000e+00;

  defaultreactions[21].rname           = "r22";
  defaultreactions[21].RadDep          = 0;
  defaultreactions[21].func1           = 2;
  defaultreactions[21].nr_chem_inp     = 2;
  defaultreactions[21].nr_chem_outp    = 1;
  defaultreactions[21].inp[0].coef     = 1.00;
  defaultreactions[21].inp[0].cname    = "NO";
  defaultreactions[21].inp[0].chem_nr  = 3;
  defaultreactions[21].inp[0].index    = 3;
  defaultreactions[21].inp[1].coef     = 1.00;
  defaultreactions[21].inp[1].cname    = "NO3";
  defaultreactions[21].inp[1].chem_nr  = 20;
  defaultreactions[21].inp[1].index    = 20;
  defaultreactions[21].outp[0].coef    = 2.00;
  defaultreactions[21].outp[0].cname   = "NO2";
  defaultreactions[21].outp[0].chem_nr = 4;
  defaultreactions[21].outp[0].index   = 4;
  defaultreactions[21].A               = +1.800e-11;
  defaultreactions[21].B               = +1.100e+02;
  defaultreactions[21].C               = +1.000e+00;
  defaultreactions[21].D               = +1.000e+00;
  defaultreactions[21].E               = +1.000e+00;
  defaultreactions[21].F               = +1.000e+00;
  defaultreactions[21].G               = +1.000e+00;

  defaultreactions[22].rname           = "r23";
  defaultreactions[22].RadDep          = 0;
  defaultreactions[22].func1           = 2;
  defaultreactions[22].nr_chem_inp     = 2;
  defaultreactions[22].nr_chem_outp    = 2;
  defaultreactions[22].inp[0].coef     = 1.00;
  defaultreactions[22].inp[0].cname    = "NO2";
  defaultreactions[22].inp[0].chem_nr  = 4;
  defaultreactions[22].inp[0].index    = 4;
  defaultreactions[22].inp[1].coef     = 1.00;
  defaultreactions[22].inp[1].cname    = "O3";
  defaultreactions[22].inp[1].chem_nr  = 1;
  defaultreactions[22].inp[1].index    = 1;
  defaultreactions[22].outp[0].coef    = 1.00;
  defaultreactions[22].outp[0].cname   = "NO3";
  defaultreactions[22].outp[0].chem_nr = 20;
  defaultreactions[22].outp[0].index   = 20;
  defaultreactions[22].outp[1].coef    = 1.00;
  defaultreactions[22].outp[1].cname   = "O2";
  defaultreactions[22].outp[1].chem_nr = 16;
  defaultreactions[22].outp[1].index   = 16;
  defaultreactions[22].A               = +1.400e-13;
  defaultreactions[22].B               = -2.470e+03;
  defaultreactions[22].C               = +1.000e+00;
  defaultreactions[22].D               = +1.000e+00;
  defaultreactions[22].E               = +1.000e+00;
  defaultreactions[22].F               = +1.000e+00;
  defaultreactions[22].G               = +1.000e+00;

  defaultreactions[23].rname           = "r24";
  defaultreactions[23].RadDep          = 0;
  defaultreactions[23].func1           = 4;
  defaultreactions[23].nr_chem_inp     = 2;
  defaultreactions[23].nr_chem_outp    = 1;
  defaultreactions[23].inp[0].coef     = 1.00;
  defaultreactions[23].inp[0].cname    = "NO2";
  defaultreactions[23].inp[0].chem_nr  = 4;
  defaultreactions[23].inp[0].index    = 4;
  defaultreactions[23].inp[1].coef     = 1.00;
  defaultreactions[23].inp[1].cname    = "NO3";
  defaultreactions[23].inp[1].chem_nr  = 20;
  defaultreactions[23].inp[1].index    = 20;
  defaultreactions[23].outp[0].coef    = 1.00;
  defaultreactions[23].outp[0].cname   = "N2O5";
  defaultreactions[23].outp[0].chem_nr = 21;
  defaultreactions[23].outp[0].index   = 21;
  defaultreactions[23].A               = +3.600e-30;
  defaultreactions[23].B               = -4.100e+00;
  defaultreactions[23].C               = +0.000e+00;
  defaultreactions[23].D               = +1.900e-12;
  defaultreactions[23].E               = +2.000e-01;
  defaultreactions[23].F               = +0.000e+00;
  defaultreactions[23].G               = +3.500e-01;

  defaultreactions[24].rname           = "r25";
  defaultreactions[24].RadDep          = 0;
  defaultreactions[24].func1           = 5;
  defaultreactions[24].nr_chem_inp     = 1;
  defaultreactions[24].nr_chem_outp    = 2;
  defaultreactions[24].inp[0].coef     = 1.00;
  defaultreactions[24].inp[0].cname    = "N2O5";
  defaultreactions[24].inp[0].chem_nr  = 21;
  defaultreactions[24].inp[0].index    = 21;
  defaultreactions[24].outp[0].coef    = 1.00;
  defaultreactions[24].outp[0].cname   = "NO3";
  defaultreactions[24].outp[0].chem_nr = 20;
  defaultreactions[24].outp[0].index   = 20;
  defaultreactions[24].outp[1].coef    = 1.00;
  defaultreactions[24].outp[1].cname   = "NO2";
  defaultreactions[24].outp[1].chem_nr = 4;
  defaultreactions[24].outp[1].index   = 4;
  defaultreactions[24].A               = +1.300e-03;
  defaultreactions[24].B               = -3.500e+00;
  defaultreactions[24].C               = -1.100e+04;
  defaultreactions[24].D               = +9.700e+14;
  defaultreactions[24].E               = +1.000e-01;
  defaultreactions[24].F               = -1.108e+04;
  defaultreactions[24].G               = +3.500e-01;

  defaultreactions[25].rname           = "r26";
  defaultreactions[25].RadDep          = 0;
  defaultreactions[25].func1           = 2;
  defaultreactions[25].nr_chem_inp     = 2;
  defaultreactions[25].nr_chem_outp    = 1;
  defaultreactions[25].inp[0].coef     = 1.00;
  defaultreactions[25].inp[0].cname    = "N2O5";
  defaultreactions[25].inp[0].chem_nr  = 21;
  defaultreactions[25].inp[0].index    = 21;
  defaultreactions[25].inp[1].coef     = 1.00;
  defaultreactions[25].inp[1].cname    = "H2O";
  defaultreactions[25].inp[1].chem_nr  = 14;
  defaultreactions[25].inp[1].index    = 14;
  defaultreactions[25].outp[0].coef    = 2.00;
  defaultreactions[25].outp[0].cname   = "HNO3";
  defaultreactions[25].outp[0].chem_nr = 18;
  defaultreactions[25].outp[0].index   = 18;
  defaultreactions[25].A               = +2.500e-22;
  defaultreactions[25].B               = +0.000e+00;
  defaultreactions[25].C               = +1.000e+00;
  defaultreactions[25].D               = +1.000e+00;
  defaultreactions[25].E               = +1.000e+00;
  defaultreactions[25].F               = +1.000e+00;
  defaultreactions[25].G               = +1.000e+00;

  defaultreactions[26].rname           = "r27";
  defaultreactions[26].RadDep          = 0;
  defaultreactions[26].func1           = 7;
  defaultreactions[26].nr_chem_inp     = 2;
  defaultreactions[26].nr_chem_outp    = 2;
  defaultreactions[26].inp[0].coef     = 1.00;
  defaultreactions[26].inp[0].cname    = "N2O5";
  defaultreactions[26].inp[0].chem_nr  = 21;
  defaultreactions[26].inp[0].index    = 21;
  defaultreactions[26].inp[1].coef     = 2.00;
  defaultreactions[26].inp[1].cname    = "H2O";
  defaultreactions[26].inp[1].chem_nr  = 14;
  defaultreactions[26].inp[1].index    = 14;
  defaultreactions[26].outp[0].coef    = 2.00;
  defaultreactions[26].outp[0].cname   = "HNO3";
  defaultreactions[26].outp[0].chem_nr = 18;
  defaultreactions[26].outp[0].index   = 18;
  defaultreactions[26].outp[1].coef    = 1.00;
  defaultreactions[26].outp[1].cname   = "H2O";
  defaultreactions[26].outp[1].chem_nr = 14;
  defaultreactions[26].outp[1].index   = 14;
  defaultreactions[26].A               = +1.800e-39;
  defaultreactions[26].B               = +1.000e+00;
  defaultreactions[26].C               = +0.000e+00;
  defaultreactions[26].D               = +0.000e+00;
  defaultreactions[26].E               = +1.000e+00;
  defaultreactions[26].F               = +1.000e+00;
  defaultreactions[26].G               = +1.000e+00;

  defaultreactions[27].rname           = "r28";
  defaultreactions[27].RadDep          = 0;
  defaultreactions[27].func1           = 2;
  defaultreactions[27].nr_chem_inp     = 2;
  defaultreactions[27].nr_chem_outp    = 2;
  defaultreactions[27].inp[0].coef     = 1.00;
  defaultreactions[27].inp[0].cname    = "OH";
  defaultreactions[27].inp[0].chem_nr  = 11;
  defaultreactions[27].inp[0].index    = 11;
  defaultreactions[27].inp[1].coef     = 1.00;
  defaultreactions[27].inp[1].cname    = "O3";
  defaultreactions[27].inp[1].chem_nr  = 1;
  defaultreactions[27].inp[1].index    = 1;
  defaultreactions[27].outp[0].coef    = 1.00;
  defaultreactions[27].outp[0].cname   = "HO2";
  defaultreactions[27].outp[0].chem_nr = 12;
  defaultreactions[27].outp[0].index   = 12;
  defaultreactions[27].outp[1].coef    = 1.00;
  defaultreactions[27].outp[1].cname   = "O2";
  defaultreactions[27].outp[1].chem_nr = 16;
  defaultreactions[27].outp[1].index   = 16;
  defaultreactions[27].A               = +1.700e-12;
  defaultreactions[27].B               = -9.400e+02;
  defaultreactions[27].C               = +1.000e+00;
  defaultreactions[27].D               = +1.000e+00;
  defaultreactions[27].E               = +1.000e+00;
  defaultreactions[27].F               = +1.000e+00;
  defaultreactions[27].G               = +1.000e+00;

  // when adding reactions, increase rsize in defaultinput.cpp, in modelinput.cpp and in mainwindow.h
}
