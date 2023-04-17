#include <iostream>
#include <algorithm>
#include <iterator>

#include <string>
#include <cmath>
#include <vector>
#include <tuple>

#include "model.h"

using namespace std;

void model::initmlm() {
    h          = 400;
    theta      = 297;
    gammatheta = 3e-3;
    theta_ft0  = 298.6;
    q_ft0      = 0.01597982;
    q          = 15e-3;
    gammaq     = -3.5e-6;
    sigmaq2       = 5e-6;

    dz_ep      = 10;
    imax_ep    = 500;
    ent_corr_factor_ep = 0.7;
    wstar      = 40;

    tv_const = Rv / Rd - 1;
}




double model::epmodel::calc_sat(double temp, double pref) {
    double esat = 610.78 * exp(17.2694 * (temp - 273.16)/(temp - 35.86));
    double qsat = Rd / Rv * esat / (pref + (1 - Rd / Rv) * esat);
    return qsat;
}

void model::epmodel::initmodel(){
    Rd = mlm->Rd;
    Rv = mlm->Rv;
    Lv = mlm->Lv;
    cp = mlm->cp;
    g = mlm->g;
    tv_const = Rv / Rd - 1;
    Ps = mlm->Ps;

    z = 0;
    parcel.qt = mlm->q + 0.51 * sqrt(mlm->sigmaq2);  // parcel state at mixed-layer top
    cout << "qtp = "<< parcel.qt << endl;
    parcel.thl = mlm->theta;
    parcel.w = mlm->wstar;
    parcel.B = 0;
    parcel.cin = 0;

    pres = Ps;

    exner = 1;
    env.thl = mlm->theta;
    env.qt = mlm->q;

    theta_ft0 = mlm->theta_ft0;
    q_ft0 = mlm->q_ft0;
    gammatheta = mlm->gammatheta;
    gammaq = mlm->gammaq;
    dz = mlm->dz_ep;

    parcel_dry = false;
    inhibited = false;
    above_lcl = false;
    above_lfc = false;

    c1 = 0.5;  // vertical velocity constants from Simpson and Wiggert (1969), Jakob and Siebesma (2003)
    c2 = 0.3333;
}

void model::epmodel::runmodel(){
    initmodel();
    imax = mlm->imax_ep;
    for (int i = 0; i < imax; i++) {
        entrainment();  // effect of lateral entrainment rate on parcel thl, qt
        z += dz;
        get_env_stats(); 	// determine state of static environment at current height level
        calc_thermo();   // parcel thermodynamics
        cout << z << endl;


        if ((parcel.thv < env.thv) and (z > mlm->h)){  // calculate decrease in vertical velocity in CIN layer
            calc_w();
        }

        cont_bool = check_if_stop(); // check if the simulation can stop
        if (!cont_bool) {
            data_zsize = i;  // size of vertical profile data
            output_struct output_inst = get_output();
            break;
        }
        save_zstep();  // save data of current height level to create vertical profile
    }
}

void model::epmodel::get_env_stats() {
    if (z <= mlm->h) {
        env.thl = mlm->theta;
        env.qt  = mlm->q;
    }
    else {
        env.thl = theta_ft0 + mlm->gammatheta * z;
        env.qt  = q_ft0     + mlm->gammaq     * z;
    }
    calc_pres();
    env.temp = exner * env.thl;

    env.qsat = calc_sat(env.temp, pres);
    env.thv = env.thl * (1 + tv_const * env.qt);
}

void model::epmodel::calc_thermo() {
    // parcel thermodynamics
    double T_old = 0;  // overwritten later
    double T_i = parcel.thl * exner;
    int i = 0;
    double eps = 1e-8;
    while ((abs(T_i - T_old) > 1e-6) && (i < 100)) {  // hard-coded error margin
        T_old = T_i;
        double dfdT = (func(T_old + eps) - func(T_old - eps)) / (2 * eps);
        T_i -=  func(T_old) / dfdT;
        cout << T_i << endl;
        i++;
    }

    parcel.qsat = calc_sat(T_i, pres);
    parcel.ql = max(0., parcel.qt - parcel.qsat);
    parcel.thv = (T_i / exner) * (1 + tv_const * parcel.qt - (1 + tv_const) * parcel.ql);
}

void model::epmodel::calc_pres() {
    double Tv0 = theta_ft0 * (1 + tv_const * q_ft0);
    double gtv = (1 + tv_const * q_ft0) * (gammatheta - g / cp) + tv_const * theta_ft0 * gammaq;
    pres = Ps * pow(1 + gtv * z / Tv0, -g / (Rd * gtv));
    exner = pow(pres / Ps, Rd / cp);
}

double model::epmodel::func(double temp) {
    double qsat_f = calc_sat(temp, pres);
    double ql_f = max(0., parcel.qt - qsat_f);
    double function = temp - exner * parcel.thl - Lv * ql_f / cp;
    return function;
}

void model::epmodel::entrainment() {
    if (z < mlm->h) {
        ent = 0;
    } else if (above_lcl) {
        ent = mlm->ent_corr_factor_ep * 1.15 * env.qt / (env.qsat * ((z - z_lcl) + 300));  //parametrization of lateral entrainment rate by Lu et al. (2018)

        if (ent > 3e-3) {
            ent = 3e-3;  //artificial upper limit
        }
    } else {
        ent = 3e-3; // specify constant entrainment rate between h and LCL since the parametrization does not work in this region
    }
    parcel.qt  -=  ent * dz * (parcel.qt - env.qt);
    parcel.thl -=  ent * dz * (parcel.thl - env.thl);
}

void model::epmodel::calc_w() {
    parcel.B    = g * (parcel.thv - env.thv) / env.thv;
    parcel.w   += dz * (- c1 * ent * parcel.w + c2 * parcel.B / parcel.w);
    parcel.cin -= parcel.B * dz;
}

bool model::epmodel::check_if_stop() {
    // check if the simulation should stop
    if (!above_lcl) {
        if (parcel.ql > 0) {
            above_lcl = true;
            z_lcl = z;
        }
    }
    if (!above_lfc) {
        above_lfc = (above_lcl && parcel.thv > env.thv && z > mlm->h);
    }
    if (parcel.w < 0) {
        parcel.w = 0;
        inhibited = true;
    }
    cond_list = {inhibited, above_lfc}; // list of reasons to stop the loop
    cont_bool = find(begin(cond_list), end(cond_list), true) == end(cond_list);
    return cont_bool;
}

model::epmodel::output_struct model::epmodel::get_output() {
    // save resulting vertical velocity etc.
    output.w_lfc = parcel.w;
    output.cin = parcel.cin;
    output.thvp = parcel.thv;
    output.thve = env.thv;
    return output;
}

void model::epmodel::save_zstep() {
    // save vertical profiles
    out_z.w.push_back(parcel.w);
    out_z.cin.push_back(parcel.cin);
}


