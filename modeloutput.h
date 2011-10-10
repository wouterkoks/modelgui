#include <string>

class outputvar
{
public:
  double *data;
  std::string name;
  std::string unit;
  std::string description;
  std::string id;
};

class modeloutput
{
public:
  outputvar t;         // time [h]
  outputvar tutc;      // time in UTC [h]

  // mixed-layer variables
  outputvar h;         // CBL height [m]
  outputvar Ps;        // surface pressure [hPa]
  outputvar ws;        // large scale vertical velocity [m s-1]
  outputvar beta;      // entrainment ratio for virtual heat [-]
  outputvar lcl;       // lifting condensation level [m]
  outputvar we;        // entrainment velocity [m s-1]
  outputvar RH;        // Relative humidity at T=theta [-]
  outputvar RHtop;     // Relative humidity at mixed-layer top [-]

  // temperature
  outputvar theta;     // initial mixed-layer potential temperature [K]
  outputvar thetav;    // initial mixed-layer virtual potential temperature [K]
  outputvar dtheta;    // initial potential temperature jump at h [K]
  outputvar dthetav;   // initial virtual potential temperature jump at h [K]
  outputvar gammatheta;// free atmosphere potential temperature lapse rate [K m-1]
  outputvar sigmatheta;// mixed-layer top potential temperature std dev [K]
  outputvar advtheta;  // advection of heat [K s-1]
  outputvar wtheta;    // surface kinematic heat flux [K m s-1]
  outputvar wthetae;   // entrainment kinematic heat flux [K m s-1]
  outputvar wthetav;   // surface kinematic virtual heat flux [K m s-1]
  outputvar wthetaM;   // mass-flux kinematic heat flux [K m s-1]

  // moisture
  outputvar q;         // mixed-layer specific humidity [kg kg-1]
  //double *qsat;      // mixed-layer saturated specific humidity [kg kg-1]
  //double *e;         // mixed-layer vapor pressure [Pa]
  //double *esat;      // mixed-layer saturated vapor pressure [Pa]
  outputvar dq;        // initial specific humidity jump at h [kg kg-1]
  outputvar gammaq;    // free atmosphere specific humidity lapse rate [kg kg-1 m-1]
  outputvar sigmaq;    // mixed-layer top specific humidity std dev [kg kg-1]
  outputvar advq;      // advection of moisture [kg kg-1 s-1]
  outputvar wq;        // surface kinematic moisture flux [kg kg-1 m s-1]
  outputvar wqe;       // entrainment kinematic moisture flux [kg kg-1 m s-1]
  outputvar wqM;       // mass-flux kinematic moisture flux [kg kg-1 m s-1]

  // wind
  outputvar u;         // initial mixed-layer u-wind speed [m s-1]
  outputvar du;        // initial u-wind jump at h [m s-1]
  outputvar gammau;    // free atmosphere u-wind speed lapse rate [s-1]
  outputvar advu;      // advection of u-wind [m s-2]
  outputvar uw;        // surface u-momentum flux [m2 s-2]
  outputvar uwe;       // entrainment u-momentum flux [m2 s-2]

  outputvar v;         // initial mixed-layer u-wind speed [m s-1]
  outputvar dv;        // initial u-wind jump at h [m s-1]
  outputvar gammav;    // free atmosphere v-wind speed lapse rate [s-1]
  outputvar advv;      // advection of v-wind [m s-2]
  outputvar vw;        // surface v-momentum flux [m2 s-2]
  outputvar vwe;       // entrainment v-momentum flux [m2 s-2]

  // BvS; a scalar...
  outputvar sca;       // mixed-layer scalar [kg kg-1]
  outputvar dsca;      // initial scalar jump at h [kg kg-1]
  outputvar gammasca;  // free atmosphere scalar lapse rate [kg kg-1 m-1]
  outputvar advsca;    // advection of scalar [kg kg-1 s-1]
  outputvar wsca;      // surface kinematic scalar flux [kg kg-1 m s-1]
  outputvar wscae;     // entrainment kinematic scalar flux [kg kg-1 m s-1]
  outputvar wscaM;     // mass-flux kinematic scalar flux [kg kg-1 m s-1]
  outputvar sigmasca;  // mixed-layer top scalar std dev [kg kg-1]

  outputvar CO2;       // initial mixed-layer CO2 [ppm]
  outputvar dCO2;      // initial CO2 jump at h [ppm]
  outputvar gammaCO2;  // free atmosphere CO2 lapse rate [ppm]
  outputvar advCO2;    // advection of CO2 [ppm]
  outputvar wCO2;      // surface kinematic CO2 flux [ppm]
  outputvar wCO2A;     //
  outputvar wCO2R;     //
  outputvar wCO2e;     // entrainment kinematic CO2 flux [ppm]
  outputvar wCO2M;     // mass-flux kinematic CO2 flux [ppm]
  outputvar sigmaCO2;  // mixed-layer top CO2 std dev [kg kg-1]

  // surface-layer
  outputvar ustar;     // friction velocity [m s-1]
  outputvar L;         // Obukhov length [m]
  outputvar Rib;       // Bulk Richardson number [-]
  outputvar ra;        // aerodynamic resistance [m s-1]
  outputvar Cm;        // drag coefficient for momentum [-]
  outputvar Cs;        // drag coefficient for scalars [-]

  // radiation
  outputvar Swin;      // Incoming short wave radiation [W m-2]
  outputvar Swout;     // Outgoing short wave radiation [W m-2]
  outputvar Lwin;      // Incoming long wave radiation [W m-2]
  outputvar Lwout;     // Outgoing long wave radiation [W m-2]
  outputvar Q;         // Net radiation [W m-2]

  // land and soil
  outputvar wg;        // soil moisture [m3 m-3]
  outputvar Tsoil;     // soil temperature [K]
  outputvar Ts;        // surface temperature [K]
  outputvar Wl;        // liquid water on vegetation [m]
  outputvar rs;        // surface resistance [s m-1]

  outputvar H;         // sensible heat flux [W m-2]
  outputvar LE;        // latent heat flux [W m-2]
  outputvar G;         // ground heat flux [W m-2]

  // shallow-cumulus
  outputvar ac;        // cloud core fraction [-]
  outputvar M;         // mass-flux (/rho) [m s-1]

  // vertical profiles
  outputvar zprof;
  outputvar thetaprof;
  outputvar scaprof;
  outputvar wthetaprof;
  outputvar wthetavprof;
  outputvar wqprof;
  outputvar qprof;

  //chemistry
  outputvar *sc;        // mixed-layer specific humidity [kg kg-1]
  outputvar phi;        // photostationary state [-]
  outputvar k_r05;      // Photolysis rate of reaction r05

  modeloutput(int,int);
  void reset(int);
  void reload(int,int);
};
