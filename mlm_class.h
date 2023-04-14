#ifndef MLM_H
#define MLM_H

class mlm_class {
public:
    double theta;
    double gammatheta;
    double theta_ft0;
    double q_ft0;
    double gammaq;
    double q2_h;
    double q;
    double h;
    double dz_ep;
    double imax_ep;
    double ent_corr_factor_ep;
    double wstar;
    void initmlm();
};

#endif
