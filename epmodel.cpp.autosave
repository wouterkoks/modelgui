#include <iostream>
#include <algorithm>
#include <iterator>

#include <string>
#include <cmath>
#include <vector>
#include <tuple>

#include "model.h"
#include "epmodel.h"

using namespace std;


double epmodel::calc_sat(double temp, double pref) {
    double esat = 610.78 * exp(17.2694 * (temp - 273.16)/(temp - 35.86));
    double qsat = Rd / Rv * esat / (pref + (1 - Rd / Rv) * esat);
    return qsat;
}


void epmodel::init(){
    Rd = 287;
    Rv = 461.5;
    Lv = 2.5e6;
    cp = 1005;
    g = 9.81;
    tv_const = Rv / Rd - 1;

    c1 = 0.5;  // vertical velocity constants from Simpson and Wiggert (1969), Jakob and Siebesma (2003)
    c2 = 0.33333;

    // initialize first height level
    z = 0;
    pres = input.Ps;
    env.thl = input.theta;
    env.qt = input.q;

    parcel.qt = input.q + input.phi_cu * sqrt(input.sigmaq2);
    parcel.thl = input.theta;
    parcel.w = input.wstar;
    parcel.cin = 0;
    parcel.B = 0;

    inhibited = false;
    above_lfc = false;
    above_lcl = false;
}

void epmodel::runmodel(){
    for (int i = 0; i < input.imax; i++) {
        entrainment();  // effect of lateral entrainment rate on parcel thl, qt
        z += input.dz;
        get_env_stats(); 	// determine state of static environment at current height level
        calc_thermo();   // parcel thermodynamics

        if ((parcel.thv < env.thv) and (z >= input.h)){  // calculate decrease in vertical velocity in CIN layer
            calc_w();
        }
        cont_bool = check_if_stop(); // check if the simulation can stop
        if (!cont_bool) {
            data_zsize = i;  // size of vertical profile data
            output = get_output();
            break;
        }
        save_zstep();  // save data of current height level to create vertical profile
    }
}

void epmodel::get_env_stats() {
    if (z < input.h) {
        env.thl = input.theta;
        env.qt  = input.q;
    }
    else {
        env.thl = input.theta_ft0 + input.gammatheta * z;
        env.qt  = input.q_ft0     + input.gammaq     * z;
        if (input.sw_ft_storage) {
            if (z < input.h + input.hstore) {
                env.thl += input.Stheta / input.hstore;
                env.qt  += input.Sq     / input.hstore;
            }
        }
    }
    calc_pres();
    env.temp = exner * env.thl;
    env.qsat = calc_sat(env.temp, pres);
    env.thv = env.thl * (1 + tv_const * env.qt);
}

void epmodel::calc_thermo() {
    // parcel thermodynamics
    double T_old = 0;  // overwritten later
    double T_i = parcel.thl * exner;
    int i = 0;
    double eps = 1e-8;
    while ((abs(T_i - T_old) > 1e-6) && (i < 100)) {  // hard-coded error margin
        T_old = T_i;
        double dfdT = (func(T_old + eps) - func(T_old - eps)) / (2 * eps);
        T_i -=  func(T_old) / dfdT;
        i++;
        if (i == 100) {
            cout << "ERROR: Calc_thermo loop did not converge!!" << endl;
        }
    }
    parcel.qsat = calc_sat(T_i, pres);
    parcel.ql = max(0., parcel.qt - parcel.qsat);
    parcel.thv = (T_i / exner) * (1 + tv_const * parcel.qt - (1 + tv_const) * parcel.ql);
}

double epmodel::func(double temp) {
    double qsat_f = calc_sat(temp, pres);
    double ql_f = max(0., parcel.qt - qsat_f);
    double function = temp - exner * parcel.thl - Lv * ql_f / cp;
    return function;
}

void epmodel::calc_pres() {
    double Tv0 = input.theta_ft0 * (1 + tv_const * input.q_ft0);
    double gtv = (1 + tv_const * input.q_ft0) * (input.gammatheta - g / cp) + tv_const * input.theta_ft0 * input.gammaq;
    pres = input.Ps * pow(1 + gtv * z / Tv0, -g / (Rd * gtv));
    exner = pow(pres / input.Ps, Rd / cp);
}

void epmodel::entrainment() {
    if (z < input.h) {
        ent = 0;
    } else if (above_lcl) {
        ent = input.ent_corr_factor * 1.15 * env.qt / (env.qsat * ((z - z_lcl) + 300));  //parametrization of lateral entrainment rate by Lu et al. (2018)
    } else {
        ent = 3e-3; // specify constant entrainment rate between h and LCL since the parametrization does not work in this region
    }
    parcel.qt  -=  ent * input.dz * (parcel.qt - env.qt);
    parcel.thl -=  ent * input.dz * (parcel.thl - env.thl);
}

void epmodel::calc_w() {
    // calculate vertical velocity based on constants from Simpson and Wiggert (1969), Jakob and Siebesma (2003)
    parcel.B    = g * (parcel.thv - env.thv) / env.thv;
    parcel.w   += input.dz * (- c1 * ent * parcel.w + c2 * parcel.B / parcel.w);
    parcel.cin -= parcel.B * input.dz;
}

bool epmodel::check_if_stop() {
    // check if the simulation should stop
    if (!above_lcl) {
        if (parcel.ql > 0) {
            above_lcl = true;
            z_lcl = z;
        }
    }
    above_lfc = (above_lcl && (parcel.thv > env.thv) && (z > input.h));

    if (parcel.w < 0) {
        parcel.w = 0;
        inhibited = true;
    }
    cond_list = {inhibited, above_lfc}; // list of reasons to stop the loop
    cont_bool = find(begin(cond_list), end(cond_list), true) == end(cond_list);  // boolean that determines whether simulation should continue
    return cont_bool;
}

epmodel::output_struct epmodel::get_output() {
    // save resulting vertical velocity etc.
    output_struct output;
    output.w_lfc = parcel.w;
    output.cin = parcel.cin;
    output.thvp = parcel.thv;
    output.thve = env.thv;
    output.inhibited = inhibited;
    return output;
}

void epmodel::save_zstep() {
    // save vertical profiles
    out_z.w.push_back(parcel.w);
    out_z.cin.push_back(parcel.cin);
}
