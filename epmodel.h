#ifndef EPMODEL_H
#define EPMODEL_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <string>


using namespace std;

class epmodel {
public:
    void runmodel();

    void init();
    struct input_struct {
        double Ps;
        double dz;
        double imax;
        double h;
        double theta;
        double q;
        double gammatheta;
        double gammaq;
        double theta_ft0;
        double q_ft0;
        double q2m;
        double wstar;
        double ent_corr_factor;
    } input;
    struct output_struct {
        double B;
        double w_lfc;
        double cin;
        double thvp;
        double thve;
    } output;
    output_struct get_output();
    struct out_z {
        std::vector<double> w;
        std::vector<double> cin;
        vector<double> thvp;
        vector<double> thve;
    } out_z;
    int data_zsize;

    void print_output();

private:
    double Rd;
    double Rv;
    double Lv;
    double cp;
    double g;
    double tv_const;
    double c1;
    double c2;

    void initmodel();
    void get_env_stats();
    void calc_thermo();
    void entrainment();
    void calc_w();
    void save_zstep();
    bool check_if_stop();
    double func(double T);
    void calc_pres();
    double calc_sat(double temp, double pref);

    double pres;
    double exner;
    double z;
    double dz;

    double ent;
    double z_lcl;

    bool cont_bool;
    bool parcel_dry;
    bool above_lcl;
    bool above_lfc;
    bool inhibited;


    struct parcel {
        double qt;
        double thl;
        double thv;
        double qsat;
        double ql;
        double w;
        double B;
        double cin;
        double dw; //remove
    } parcel;

    struct env {
        double qt;
        double thl;
        double thv;
        double qsat;
        double temp;
    } env;

    std::vector<bool> cond_list;
};


#endif
