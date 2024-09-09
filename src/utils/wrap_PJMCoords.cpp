#include <cmath>

namespace py = pybind11;

typedef PyVector GCA;
typedef PyVector GCY;
typedef PyVector LSR;
typedef PyVector HCA;
typedef PyVector HGP;
typedef PyVector HEQ;
typedef PyVector bool6;
typedef PyVector vec6;
typedef PyMatrix mat33;

class OmniCoords {
  bool6 know;
  double Rsun, zsun, vcsun, Usun, Vsun, Wsun, epoch;
  mat33 EtoP;
  void SetTrans();
  void HEQfromHCA();
  void HGPfromHCA();
  void HCAfromLSR();
  void LSRfromGCA();
  void GCAfromGCY();
  void Backward(int);
  void HCAfromHEQ();
  void HCAfromHGP();
  void LSRfromHCA();
  void GCAfromLSR();
  void GCYfromGCA();
  void Forward(int);
  vec6 rv[6];

 public:
  OmniCoords();
  ~OmniCoords() {};
  
  void change_sol_pos(double, double);
  void change_vc(double);
  void change_vsol(double, double, double);
  void set_SBD10();
  void set_DB98();
  void change_epoch(double);

  void give_sol_pos(double&, double&);
  double give_Rsun() { return Rsun; }
  double give_zsun() { return zsun; }
  void give_vc(double&);
  double give_vcsun() { return vcsun; }
  void give_vsol(double&, double&, double&);
  void give_epoch(double&);

  vec6 give_HEQ() { return give(0); }
  vec6 give_HGP() { return give(1); }
  vec6 give_HCA() { return give(2); }
  vec6 give_LSR() { return give(3); }
  vec6 give_GCA() { return give(4); }
  vec6 give_GCY() { return give(5); }
  vec6 give_HEQ_units() { return give_units(0); }
  vec6 give_HGP_units() { return give_units(1); }
  vec6 give_HCA_units() { return give_units(2); }
  vec6 give_LSR_units() { return give_units(3); }
  vec6 give_GCA_units() { return give_units(4); }
  vec6 give_GCY_units() { return give_units(5); }
  vec6 give(int);
  vec6 give_units(int);
  void take_HEQ(vec6&&);
  void take_HGP(vec6&&);
  void take_HCA(vec6&&);
  void take_LSR(vec6&&);
  void take_GCA(vec6&&);
  void take_GCY(vec6&&);
  void take_HEQ_units(vec6&&);
  void take_HGP_units(vec6&&);
  void take_HCA_units(vec6&&);
  void take_LSR_units(vec6&&);
  void take_GCA_units(vec6&&);
  void take_GCY_units(vec6&&);

  vec6 HEQfromHGP(vec6&& sHGP) { take_HGP(std::move(sHGP)); return give_HEQ(); }
  vec6 HEQfromHCA(vec6&& sHCA) { take_HCA(std::move(sHCA)); return give_HEQ(); }
  vec6 HEQfromLSR(vec6&& sLSR) { take_LSR(std::move(sLSR)); return give_HEQ(); }
  vec6 HEQfromGCA(vec6&& sGCA) { take_GCA(std::move(sGCA)); return give_HEQ(); }
  vec6 HEQfromGCY(vec6&& sGCY) { take_GCY(std::move(sGCY)); return give_HEQ(); }

  vec6 HGPfromHEQ(vec6&& sHEQ) { take_HEQ(std::move(sHEQ)); return give_HGP(); }
  vec6 HGPfromHCA(vec6&& sHCA) { take_HCA(std::move(sHCA)); return give_HGP(); }
  vec6 HGPfromLSR(vec6&& sLSR) { take_LSR(std::move(sLSR)); return give_HGP(); }
  vec6 HGPfromGCA(vec6&& sGCA) { take_GCA(std::move(sGCA)); return give_HGP(); }
  vec6 HGPfromGCY(vec6&& sGCY) { take_GCY(std::move(sGCY)); return give_HGP(); }

  vec6 HCAfromHEQ(vec6&& sHEQ) { take_HEQ(std::move(sHEQ)); return give_HCA(); }
  vec6 HCAfromHGP(vec6&& sHGP) { take_HGP(std::move(sHGP)); return give_HCA(); }
  vec6 HCAfromLSR(vec6&& sLSR) { take_LSR(std::move(sLSR)); return give_HCA(); }
  vec6 HCAfromGCA(vec6&& sGCA) { take_GCA(std::move(sGCA)); return give_HCA(); }
  vec6 HCAfromGCY(vec6&& sGCY) { take_GCY(std::move(sGCY)); return give_HCA(); }

  vec6 LSRfromHEQ(vec6&& sHEQ) { take_HEQ(std::move(sHEQ)); return give_LSR(); }
  vec6 LSRfromHGP(vec6&& sHGP) { take_HGP(std::move(sHGP)); return give_LSR(); }
  vec6 LSRfromHCA(vec6&& sHCA) { take_HCA(std::move(sHCA)); return give_LSR(); }
  vec6 LSRfromGCA(vec6&& sGCA) { take_GCA(std::move(sGCA)); return give_LSR(); }
  vec6 LSRfromGCY(vec6&& sGCY) { take_GCY(std::move(sGCY)); return give_LSR(); }

  vec6 GCAfromHEQ(vec6&& sHEQ) { take_HEQ(std::move(sHEQ)); return give_GCA(); }
  vec6 GCAfromHGP(vec6&& sHGP) { take_HGP(std::move(sHGP)); return give_GCA(); }
  vec6 GCAfromHCA(vec6&& sHCA) { take_HCA(std::move(sHCA)); return give_GCA(); }
  vec6 GCAfromLSR(vec6&& sLSR) { take_LSR(std::move(sLSR)); return give_GCA(); }
  vec6 GCAfromGCY(vec6&& sGCY) { take_GCY(std::move(sGCY)); return give_GCA(); }

  vec6 GCYfromHEQ(vec6&& sHEQ) { take_HEQ(std::move(sHEQ)); return give_GCY(); }
  vec6 GCYfromHGP(vec6&& sHGP) { take_HGP(std::move(sHGP)); return give_GCY(); }
  vec6 GCYfromHCA(vec6&& sHCA) { take_HCA(std::move(sHCA)); return give_GCY(); }
  vec6 GCYfromLSR(vec6&& sLSR) { take_LSR(std::move(sLSR)); return give_GCY(); }
  vec6 GCYfromGCA(vec6&& sGCA) { take_GCA(std::move(sGCA)); return give_GCY(); }
};

inline OmniCoords::OmniCoords()
    : know(py::make_tuple(false, false, false, false, false, false)),
      EtoP(3, 3),
      rv{vec6(py::make_tuple(0, 0, 0, 0, 0, 0)),
         vec6(py::make_tuple(0, 0, 0, 0, 0, 0)),
         vec6(py::make_tuple(0, 0, 0, 0, 0, 0)),
         vec6(py::make_tuple(0, 0, 0, 0, 0, 0)),
         vec6(py::make_tuple(0, 0, 0, 0, 0, 0)),
         vec6(py::make_tuple(0, 0, 0, 0, 0, 0))} {
	epoch = 2000.;
	SetTrans();
	Rsun = Rsun;
	zsun = zsun;
	vcsun = vcsun;
	Usun = usun;
	Vsun = vsun;
	Wsun = wsun;

}

inline void OmniCoords::take_HEQ(vec6&& tHEQ) {
  rv[0] = std::move(tHEQ);
  know.__setitem__(0, py::cast(true));
}

inline void OmniCoords::take_HEQ_units(vec6&& tHEQ) {
  vec6 tmp = std::move(tHEQ);
  tmp.__setitem__(1, py::cast(tmp.__getitem__(1).cast<double>() * Units::degree));
  tmp.__setitem__(2, py::cast(tmp.__getitem__(2).cast<double>() * Units::degree));
  tmp.__setitem__(3, py::cast(tmp.__getitem__(3).cast<double>() * Units::kms));
  tmp.__setitem__(4, py::cast(tmp.__getitem__(4).cast<double>() * Units::masyr));
  tmp.__setitem__(5, py::cast(tmp.__getitem__(5).cast<double>() * Units::masyr));
  rv[0] = std::move(tmp);
  know.__setitem__(0, py::cast(true));
}

inline void OmniCoords::take_HGP(vec6&& tHGP) {
  rv[1] = std::move(tHGP);
  know.__setitem__(1, py::cast(true));
}

inline void OmniCoords::take_HGP_units(vec6&& tHGP) {
  vec6 tmp = std::move(tHGP);
  tmp.__setitem__(1, py::cast(tmp.__getitem__(1).cast<double>() * Units::degree));
  tmp.__setitem__(2, py::cast(tmp.__getitem__(2).cast<double>() * Units::degree));
  tmp.__setitem__(3, py::cast(tmp.__getitem__(3).cast<double>() * Units::kms));
  tmp.__setitem__(4, py::cast(tmp.__getitem__(4).cast<double>() * Units::masyr));
  tmp.__setitem__(5, py::cast(tmp.__getitem__(5).cast<double>() * Units::masyr));
  rv[1] = std::move(tmp);
  know.__setitem__(1, py::cast(true));
}

inline void OmniCoords::take_HCA(vec6&& tHCA) {
  rv[2] = std::move(tHCA);
  know.__setitem__(2, py::cast(true));
}

inline void OmniCoords::take_HCA_units(vec6&& tHCA) {
  vec6 tmp = std::move(tHCA);
  tmp.__setitem__(1, py::cast(tmp.__getitem__(1).cast<double>() * Units::degree));
  tmp.__setitem__(2, py::cast(tmp.__getitem__(2).cast<double>() * Units::degree));
  tmp.__setitem__(3, py::cast(tmp.__getitem__(3).cast<double>() * Units::kms));
  tmp.__setitem__(4, py::cast(tmp.__getitem__(4).cast<double>() * Units::masyr));
  tmp.__setitem__(5, py::cast(tmp.__getitem__(5).cast<double>() * Units::masyr));
  rv[2] = std::move(tmp);
  know.__setitem__(2, py::cast(true));
}

inline void OmniCoords::take_LSR(vec6&& tLSR) {
  rv[3] = std::move(tLSR);
  know.__setitem__(3, py::cast(true));
}

inline void OmniCoords::take_LSR_units(vec6&& tLSR) {
  vec6 tmp = std::move(tLSR);
  tmp.__setitem__(1, py::cast(tmp.__getitem__(1).cast<double>() * Units::degree));
  tmp.__setitem__(2, py::cast(tmp.__getitem__(2).cast<double>() * Units::degree));
  tmp.__setitem__(3, py::cast(tmp.__getitem__(3).cast<double>() * Units::kms));
  tmp.__setitem__(4, py::cast(tmp.__getitem__(4).cast<double>() * Units::masyr));
  tmp.__setitem__(5, py::cast(tmp.__getitem__(5).cast<double>() * Units::masyr));
  rv[3] = std::move(tmp);
  know.__setitem__(3, py::cast(true));
}

inline void OmniCoords::take_GCA(vec6&& tGCA) {
  rv[4] = std::move(tGCA);
  know.__setitem__(4, py::cast(true));
}

inline void OmniCoords::take_GCA_units(vec6&& tGCA) {
  vec6 tmp = std::move(tGCA);
  tmp.__setitem__(1, py::cast(tmp.__getitem__(1).cast<double>() * Units::degree));
  tmp.__setitem__(2, py::cast(tmp.__getitem__(2).cast<double>() * Units::degree));
  tmp.__setitem__(3, py::cast(tmp.__getitem__(3).cast<double>() * Units::kms));
  tmp.__setitem__(4, py::cast(tmp.__getitem__(4).cast<double>() * Units::masyr));
  tmp.__setitem__(5, py::cast(tmp.__getitem__(5).cast<double>() * Units::masyr));
  rv[4] = std::move(tmp);
  know.__setitem__(4, py::cast(true));
}

inline void OmniCoords::take_GCY(vec6&& tGCY) {
  rv[5] = std::move(tGCY);
  know.__setitem__(5, py::cast(true));
}

inline void OmniCoords::take_GCY_units(vec6&& tGCY) {
  vec6 tmp = std::move(tGCY);
  tmp.__setitem__(1, py::cast(tmp.__getitem__(1).cast<double>() * Units::degree));
  tmp.__setitem__(2, py::cast(tmp.__getitem__(2).cast<double>() * Units::degree));
  tmp.__setitem__(3, py::cast(tmp.__getitem__(3).cast<double>() * Units::kms));
  tmp.__setitem__(4, py::cast(tmp.__getitem__(4).cast<double>() * Units::masyr));
  tmp.__setitem__(5, py::cast(tmp.__getitem__(5).cast<double>() * Units::masyr));
  rv[5] = std::move(tmp);
  know.__setitem__(5, py::cast(true));
}

void OmniCoords::GCYfromGCA() {
    double s, c;
    rv[5].__setitem__(0, py::float_(hypot(rv[4].__getitem__(0).cast<double>(), rv[4].__getitem__(1).cast<double>())));
    rv[5].__setitem__(1, rv[4].__getitem__(2));
    if (rv[5].__getitem__(0).cast<double>()) {
        s = rv[4].__getitem__(1).cast<double>() / rv[5].__getitem__(0).cast<double>();
        c = rv[4].__getitem__(0).cast<double>() / rv[5].__getitem__(0).cast<double>();
        rv[5].__setitem__(2, py::float_((s > 0.0) ? acos(c) : 2 * M_PI - acos(c)));
    } else {
        s = 0.0;
        c = 1.0;
        rv[5].__setitem__(2, py::float_(0.0));
    }
    rv[5].__setitem__(3, py::float_(c * rv[4].__getitem__(3).cast<double>() + s * rv[4].__getitem__(4).cast<double>()));
    rv[5].__setitem__(5, py::float_(-s * rv[4].__getitem__(3).cast<double>() + c * rv[4].__getitem__(4).cast<double>()));
    rv[5].__setitem__(4, rv[4].__getitem__(5));
    know.__setitem__(5, py::bool_(true));
}

void OmniCoords::GCAfromGCY() {
    double s = sin(rv[5].__getitem__(2).cast<double>());
    double c = cos(rv[5].__getitem__(2).cast<double>());
    rv[4].__setitem__(0, py::float_(c * rv[5].__getitem__(0).cast<double>()));
    rv[4].__setitem__(1, py::float_(s * rv[5].__getitem__(0).cast<double>()));
    rv[4].__setitem__(2, rv[5].__getitem__(1));
    rv[4].__setitem__(3, py::float_(c * rv[5].__getitem__(3).cast<double>() - s * rv[5].__getitem__(5).cast<double>()));
    rv[4].__setitem__(4, py::float_(s * rv[5].__getitem__(3).cast<double>() + c * rv[5].__getitem__(5).cast<double>()));
    rv[4].__setitem__(5, rv[5].__getitem__(4));
    know.__setitem__(4, py::bool_(true));
}

void OmniCoords::LSRfromGCA() {
    static double z = 0.0, s = 0.0, c = 1.0;
    rv[3].__setitem__(0, py::float_(Rsun - rv[4].__getitem__(0).cast<double>()));   
    rv[3].__setitem__(1, py::float_(-rv[4].__getitem__(1).cast<double>()));        
    rv[3].__setitem__(2, py::float_(rv[4].__getitem__(2).cast<double>() - zsun));
    rv[3].__setitem__(3, py::float_(-rv[4].__getitem__(3).cast<double>()));
    rv[3].__setitem__(4, py::float_(vcsun - rv[4].__getitem__(4).cast<double>()));
    rv[3].__setitem__(5, rv[4].__getitem__(5));
    if (zsun) { 
        double t;
        if (z != zsun) {
            z = zsun;
            t = hypot(zsun, Rsun);
            s = zsun / t;
            c = Rsun / t;
        }
        t = rv[3].__getitem__(0).cast<double>();
        rv[3].__setitem__(0, py::float_(c * t - s * rv[3].__getitem__(2).cast<double>()));
        rv[3].__setitem__(2, py::float_(s * t + c * rv[3].__getitem__(2).cast<double>()));
        t = rv[3].__getitem__(3).cast<double>();
        rv[3].__setitem__(3, py::float_(c * t - s * rv[3].__getitem__(5).cast<double>()));
        rv[3].__setitem__(5, py::float_(s * t + c * rv[3].__getitem__(5).cast<double>()));
    }
    know.__setitem__(3, py::bool_(true));
}

void OmniCoords::GCAfromLSR() {
    static double z = 0.0, s = 0.0, c = 1.0;
    if (zsun) { 
        PyVector in = std::move(rv[3]);
        double t;
        if (z != zsun) {
            z = zsun;
            t = hypot(zsun, Rsun);
            s = zsun / t;
            c = Rsun / t;
        }
        t = in.__getitem__(0).cast<double>();
        in.__setitem__(0, py::float_(c * t + s * in.__getitem__(2).cast<double>()));
        in.__setitem__(2, py::float_(-s * t + c * in.__getitem__(2).cast<double>()));
        t = in.__getitem__(3).cast<double>();
        in.__setitem__(3, py::float_(c * t + s * in.__getitem__(5).cast<double>()));
        in.__setitem__(5, py::float_(-s * t + c * in.__getitem__(5).cast<double>()));
        rv[4].__setitem__(0, py::float_(Rsun - in.__getitem__(0).cast<double>()));
        rv[4].__setitem__(1, py::float_(-in.__getitem__(1).cast<double>()));
        rv[4].__setitem__(2, py::float_(in.__getitem__(2).cast<double>() + zsun));
        rv[4].__setitem__(3, py::float_(-in.__getitem__(3).cast<double>()));
        rv[4].__setitem__(4, py::float_(vcsun - in.__getitem__(4).cast<double>()));
        rv[4].__setitem__(5, in.__getitem__(5));
    } else {
        rv[4].__setitem__(0, py::float_(Rsun - rv[3].__getitem__(0).cast<double>()));
        rv[4].__setitem__(1, py::float_(-rv[3].__getitem__(1).cast<double>()));
        rv[4].__setitem__(2, rv[3].__getitem__(2));
        rv[4].__setitem__(3, py::float_(-rv[3].__getitem__(3).cast<double>()));
        rv[4].__setitem__(4, py::float_(vcsun - rv[3].__getitem__(4).cast<double>()));
        rv[4].__setitem__(5, rv[3].__getitem__(5));
    }
    know.__setitem__(4, py::bool_(true));
}

void OmniCoords::LSRfromHCA() {
    rv[3].__setitem__(0, rv[2].__getitem__(0));
    rv[3].__setitem__(1, rv[2].__getitem__(1));
    rv[3].__setitem__(2, rv[2].__getitem__(2));
    rv[3].__setitem__(3, py::float_(rv[2].__getitem__(3).cast<double>() + Usun));
    rv[3].__setitem__(4, py::float_(rv[2].__getitem__(4).cast<double>() + Vsun));
    rv[3].__setitem__(5, py::float_(rv[2].__getitem__(5).cast<double>() + Wsun));
    know.__setitem__(3, py::bool_(true));
}

void OmniCoords::HCAfromLSR() {
    rv[2].__setitem__(0, rv[3].__getitem__(0));
    rv[2].__setitem__(1, rv[3].__getitem__(1));
    rv[2].__setitem__(2, rv[3].__getitem__(2));
    rv[2].__setitem__(3, py::float_(rv[3].__getitem__(3).cast<double>() - Usun));
    rv[2].__setitem__(4, py::float_(rv[3].__getitem__(4).cast<double>() - Vsun));
    rv[2].__setitem__(5, py::float_(rv[3].__getitem__(5).cast<double>() - Wsun));
    know.__setitem__(2, py::bool_(true));
}

void OmniCoords::HGPfromHCA() {
    double R = hypot(rv[2].__getitem__(0).cast<double>(), rv[2].__getitem__(1).cast<double>());
    rv[1].__setitem__(0, py::float_(hypot(R, rv[2].__getitem__(2).cast<double>())));

    double cl = (R == 0.0) ? 1.0 : rv[2].__getitem__(0).cast<double>() / R;
    double sl = (R == 0.0) ? 0.0 : rv[2].__getitem__(1).cast<double>() / R;
    double cb = (rv[1].__getitem__(0).cast<double>() == 0.0) ? 1.0 : R / rv[1].__getitem__(0).cast<double>();
    double sb = (rv[1].__getitem__(0).cast<double>() == 0.0) ? 0.0 : rv[2].__getitem__(2).cast<double>() / rv[1].__getitem__(0).cast<double>();
    double temp = cl * rv[2].__getitem__(3).cast<double>() + sl * rv[2].__getitem__(4).cast<double>();

    rv[1].__setitem__(1, py::float_((sl < 0.0) ? TPi - acos(cl) : acos(cl)));
    rv[1].__setitem__(2, py::float_(asin(sb)));
    rv[1].__setitem__(3, py::float_(cb * temp + sb * rv[2].__getitem__(5).cast<double>()));
    rv[1].__setitem__(4, py::float_((rv[1].__getitem__(0).cast<double>() == 0.0) ? 0.0 : (cl * rv[2].__getitem__(4).cast<double>() - sl * rv[2].__getitem__(3).cast<double>()) / rv[1].__getitem__(0).cast<double>()));
    rv[1].__setitem__(5, py::float_((rv[1].__getitem__(0).cast<double>() == 0.0) ? 0.0 : (cb * rv[2].__getitem__(5).cast<double>() - sb * temp) / rv[1].__getitem__(0).cast<double>()));

    know.__setitem__(1, py::bool_(true));
}

void OmniCoords::HCAfromHGP() {
    double cl = cos(rv[1].__getitem__(1).cast<double>());
    double sl = sin(rv[1].__getitem__(1).cast<double>());
    double cb = cos(rv[1].__getitem__(2).cast<double>());
    double sb = sin(rv[1].__getitem__(2).cast<double>());
    double R = cb * rv[1].__getitem__(0).cast<double>();
    double vl = rv[1].__getitem__(0).cast<double>() * rv[1].__getitem__(4).cast<double>();
    double vb = rv[1].__getitem__(0).cast<double>() * rv[1].__getitem__(5).cast<double>();
    double temp = cb * rv[1].__getitem__(3).cast<double>() - sb * vb;

    rv[2].__setitem__(0, py::float_(cl * R));
    rv[2].__setitem__(1, py::float_(sl * R));
    rv[2].__setitem__(2, py::float_(sb * rv[1].__getitem__(0).cast<double>()));
    rv[2].__setitem__(3, py::float_(cl * temp - sl * vl));
    rv[2].__setitem__(4, py::float_(sl * temp + cl * vl));
    rv[2].__setitem__(5, py::float_(sb * rv[1].__getitem__(3).cast<double>() + cb * vb));

    know.__setitem__(2, py::bool_(true));
}
void OmniCoords::SetTrans() {
    const double A[6] = {0.01118086087, 6.770713945e-6, -6.738910167e-10,
                         1.463555541e-6, -1.667759063e-9, 8.720828496e-8},
                 B[6] = {0.01118086087, 6.770713945e-6, -6.738910167e-10,
                         5.307158404e-6, 3.199770295e-10, 8.825063437e-8},
                 C[6] = {0.009717173455, -4.136915141e-6, -1.052045688e-9,
                         -2.06845757e-6, -1.052045688e-9, -2.028121072e-7};

    if (epoch == 1991.25) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                py::tuple index = py::make_tuple(i, j);
                EtoP.__setitem__(index, py::float_(EtoPJ1991[i][j]));
                if (EtoPJ1991[i][j] == 0.0) {
                    std::cerr << "Something has gone wrong in OmniCoords.\n"
                              << "I've seen this before when OmniCoords declared as "
                              << "global variable\n";
                    exit(1);
                }
            }
        }
        return;
    }

    if (epoch == 2000.0) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                py::tuple index = py::make_tuple(i, j);
                EtoP.__setitem__(index, py::float_(EtoPJ2000[i][j]));
                if (EtoPJ2000[i][j] == 0.0) {
                    std::cerr << "Something has gone badly wrong in OmniCoords.\n"
                              << "I've seen this before when OmniCoords declared as "
                              << "global variable\n";
                    exit(1);
                }
            }
        }
        return;
    }

    double T = 0.01 * (epoch - 2000.0);
    double t = -T;
    double zt = (((A[5] * t + (A[4] * T + A[3])) * t + ((A[2] * T + A[1]) * T + A[0])) * t);
    double th = (((B[5] * t + (B[4] * T + B[3])) * t + ((B[2] * T + B[1]) * T + B[0])) * t);
    double Ze = (((C[5] * t + (C[4] * T + C[3])) * t + ((C[2] * T + C[1]) * T + C[0])) * t);
    double czt = cos(zt);
    double szt = sin(zt);
    double cth = cos(th);
    double sth = sin(th);
    double cZe = cos(Ze);
    double sZe = sin(Ze);

    double P[3][3] = {
        {czt * cth * cZe - szt * sZe, -czt * cth * sZe - szt * cZe, -czt * sth},
        {szt * cth * cZe + czt * sZe, -szt * cth * sZe + czt * cZe, -szt * sth},
        {szt * cZe, -sth * sZe, cth}
    };
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			EtoP.set_value(i, j, py::float_(0.0));
			for (int k = 0; k < 3; ++k) {
				double current_value = EtoP.getValueAt(i, j).cast<double>();
				double newValue = current_value + EtoPJ2000[i][k] * P[k][j];
				EtoP.set_value(i, j, py::float_(newValue));
			}
		}
	}

}
void OmniCoords::HEQfromHCA() {
    int i, j;
    py::list h = py::list(6); 

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            h[i] = h[i].cast<double>() + rv[2].__getitem__(j).cast<double>() * EtoP.getValueAt(i, j).cast<double>();
        	h[i + 3] = h[i + 3].cast<double>() + rv[2].__getitem__(j + 3).cast<double>() * EtoP.getValueAt(i, j).cast<double>();
		}
    }

    double R = hypot(h[0].cast<double>(), h[1].cast<double>());
    rv[0].__setitem__(0, py::float_(hypot(R, h[2].cast<double>())));

    double ca = (R == 0.0) ? 1.0 : h[0].cast<double>() / R;
    double sa = (R == 0.0) ? 0.0 : h[1].cast<double>() / R;
    double cd = (rv[0].__getitem__(0).cast<double>() == 0.0) ? 1.0 : R / rv[0].__getitem__(0).cast<double>();
    double sd = (rv[0].__getitem__(0).cast<double>() == 0.0) ? 0.0 : h[2].cast<double>() / rv[0].__getitem__(0).cast<double>();
    double temp = ca * h[3].cast<double>() + sa * h[4].cast<double>();

    rv[0].__setitem__(1, py::float_((sa < 0.0) ? TPi - acos(ca) : acos(ca)));
    rv[0].__setitem__(2, py::float_(asin(sd)));
    rv[0].__setitem__(3, py::float_(cd * temp + sd * h[5].cast<double>()));
    rv[0].__setitem__(4, py::float_((rv[0].__getitem__(0).cast<double>() == 0.0) ? 0.0 : (ca * h[4].cast<double>() - sa * h[3].cast<double>()) / rv[0].__getitem__(0).cast<double>()));
    rv[0].__setitem__(5, py::float_((rv[0].__getitem__(0).cast<double>() == 0.0) ? 0.0 : (cd * h[5].cast<double>() - sd * temp) / rv[0].__getitem__(0).cast<double>()));

    know.__setitem__(0, py::bool_(true));
}

void OmniCoords::HCAfromHEQ() {
    int i, j;
    py::list h = py::list(6); // Initialize a Python list with 6 elements

    double ca = cos(rv[0].__getitem__(1).cast<double>());
    double sa = sin(rv[0].__getitem__(1).cast<double>());
    double cd = cos(rv[0].__getitem__(2).cast<double>());
    double sd = sin(rv[0].__getitem__(2).cast<double>());
    double R = cd * rv[0].__getitem__(0).cast<double>();
    double va = rv[0].__getitem__(0).cast<double>() * rv[0].__getitem__(4).cast<double>();
    double vd = rv[0].__getitem__(0).cast<double>() * rv[0].__getitem__(5).cast<double>();
    double temp = cd * rv[0].__getitem__(3).cast<double>() - sd * vd;

    rv[2].__setitem__(0, py::float_(0.0));
    rv[2].__setitem__(1, py::float_(0.0));
    rv[2].__setitem__(2, py::float_(0.0));
    rv[2].__setitem__(3, py::float_(0.0));
    rv[2].__setitem__(4, py::float_(0.0));
    rv[2].__setitem__(5, py::float_(0.0));

    h[0] = py::float_(ca * R);
    h[1] = py::float_(sa * R);
    h[2] = py::float_(sd * rv[0].__getitem__(0).cast<double>());
    h[3] = py::float_(ca * temp - sa * va);
    h[4] = py::float_(sa * temp + ca * va);
    h[5] = py::float_(sd * rv[0].__getitem__(3).cast<double>() + cd * vd);

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
        	rv[2].__setitem__(i, py::float_(rv[2].__getitem__(i).cast<double>() + h[j].cast<double>() * EtoP.getValueAt(i, j).cast<double>()));
			rv[2].__setitem__(i + 3, py::float_(rv[2].__getitem__(i + 3).cast<double>() + h[j + 3].cast<double>() * EtoP.getValueAt(i, j).cast<double>()));
		}
    }

    know.__setitem__(2, py::bool_(true));
}

void init_pjmcoords(py::module_ &torus) {
	py::class_<OmniCoords>(torus, "OmniCoords")
        .def(py::init<>())
		.def(py::init<>())
        .def("change_sol_pos", &OmniCoords::change_sol_pos)  
        .def("change_vc", &OmniCoords::change_vc)            
        .def("change_vsol", &OmniCoords::change_vsol)        
        .def("set_SBD10", &OmniCoords::set_SBD10)            
        .def("set_DB98", &OmniCoords::set_DB98)              
        .def("change_epoch", &OmniCoords::change_epoch)	
		.def("give_sol_pos", [](OmniCoords &self) {
            double r, z;
            self.give_sol_pos(r, z);
            return py::make_tuple(r, z);
        })
        .def("give_Rsun", &OmniCoords::give_Rsun)
        .def("give_zsun", &OmniCoords::give_zsun)
        .def("give_vc", [](OmniCoords &self) {
            double vc;
            self.give_vc(vc);
            return vc;
        })
        .def("give_vcsun", &OmniCoords::give_vcsun)
        .def("give_vsol", [](OmniCoords &self) {
            double U, V, W;
            self.give_vsol(U, V, W);
            return py::make_tuple(U, V, W);
        })
        .def("give_epoch", [](OmniCoords &self) {
            double epoch;
            self.give_epoch(epoch);
            return epoch;
        })

		.def("give_HEQ", &OmniCoords::give_HEQ)
        .def("give_HGP", &OmniCoords::give_HGP)
        .def("give_HCA", &OmniCoords::give_HCA)
        .def("give_LSR", &OmniCoords::give_LSR)
        .def("give_GCA", &OmniCoords::give_GCA)
        .def("give_GCY", &OmniCoords::give_GCY)
        .def("give_HEQ_units", &OmniCoords::give_HEQ_units)
        .def("give_HGP_units", &OmniCoords::give_HGP_units)
        .def("give_HCA_units", &OmniCoords::give_HCA_units)
        .def("give_LSR_units", &OmniCoords::give_LSR_units)
        .def("give_GCA_units", &OmniCoords::give_GCA_units)
        .def("give_GCY_units", &OmniCoords::give_GCY_units)
		.def("give", &OmniCoords::give)                  
        .def("give_units", &OmniCoords::give_units)
		.def("take_HEQ", [](OmniCoords &self, vec6 &v) {  
            self.take_HEQ(std::move(v));
        })	
		.def("take_HGP", [](OmniCoords &self, vec6 &v) { 
			self.take_HGP(std::move(v)); 
		})
        .def("take_HCA", [](OmniCoords &self, vec6 &v) { 
			self.take_HCA(std::move(v)); 
		})
        .def("take_LSR", [](OmniCoords &self, vec6 &v) { 
			self.take_LSR(std::move(v)); 
		})
        .def("take_GCA", [](OmniCoords &self, vec6 &v) { 
			self.take_GCA(std::move(v)); 
		})
        .def("take_GCY", [](OmniCoords &self, vec6 &v) { 
			self.take_GCY(std::move(v)); 
		})
        .def("take_HEQ_units", [](OmniCoords &self, vec6 &v) { 
			self.take_HEQ_units(std::move(v)); 
		})
        .def("take_HGP_units", [](OmniCoords &self, vec6 &v) { 
			self.take_HGP_units(std::move(v)); 
		})
        .def("take_HCA_units", [](OmniCoords &self, vec6 &v) { 
			self.take_HCA_units(std::move(v)); 
		})
        .def("take_LSR_units", [](OmniCoords &self, vec6 &v) { 
			self.take_LSR_units(std::move(v)); 
		})
        .def("take_GCA_units", [](OmniCoords &self, vec6 &v) { 
			self.take_GCA_units(std::move(v)); 
		})
        .def("take_GCY_units", [](OmniCoords &self, vec6 &v) { 
			self.take_GCY_units(std::move(v)); 
		})
		.def("HEQfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromHCA(std::move(v));
        })		
		.def("HEQfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromHCA(std::move(v));
        })
        .def("HEQfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromLSR(std::move(v));
        })
        .def("HEQfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromGCA(std::move(v));
        })
        .def("HEQfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromGCY(std::move(v));
        })
        .def("HGPfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromHEQ(std::move(v));
        })
        .def("HGPfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromHCA(std::move(v));
        })
        .def("HGPfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromLSR(std::move(v));
        })
        .def("HGPfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromGCA(std::move(v));
        })
        .def("HGPfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromGCY(std::move(v));
        })
        .def("HCAfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromHEQ(std::move(v));
        })
        .def("HCAfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromHGP(std::move(v));
        })
        .def("HCAfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromLSR(std::move(v));
        })
        .def("HCAfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromGCA(std::move(v));
        })
        .def("HCAfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromGCY(std::move(v));
        })
        .def("LSRfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromHEQ(std::move(v));
        })
        .def("LSRfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromHGP(std::move(v));
        })
        .def("LSRfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromHCA(std::move(v));
        })
        .def("LSRfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromGCA(std::move(v));
        })
        .def("LSRfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromGCY(std::move(v));
        })
        .def("GCAfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromHEQ(std::move(v));
        })
        .def("GCAfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromHGP(std::move(v));
        })
        .def("GCAfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromHCA(std::move(v));
        })
        .def("GCAfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromLSR(std::move(v));
        })
        .def("GCAfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromGCY(std::move(v));
        })
        .def("GCYfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromHEQ(std::move(v));
        })
        .def("GCYfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromHGP(std::move(v));
        })
        .def("GCYfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromHCA(std::move(v));
        })
        .def("GCYfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromLSR(std::move(v));
        })
        .def("GCYfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromGCA(std::move(v));
        })
		.def("HEQfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromHCA(std::move(v));
        })
        .def("HEQfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromLSR(std::move(v));
        })
        .def("HEQfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromGCA(std::move(v));
        })
        .def("HEQfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.HEQfromGCY(std::move(v));
        })
        .def("HGPfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromHEQ(std::move(v));
        })
        .def("HGPfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromHCA(std::move(v));
        })
        .def("HGPfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromLSR(std::move(v));
        })
        .def("HGPfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromGCA(std::move(v));
        })
        .def("HGPfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.HGPfromGCY(std::move(v));
        })
        .def("HCAfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromHEQ(std::move(v));
        })
        .def("HCAfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromHGP(std::move(v));
        })
        .def("HCAfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromLSR(std::move(v));
        })
        .def("HCAfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromGCA(std::move(v));
        })
        .def("HCAfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.HCAfromGCY(std::move(v));
        })
        .def("LSRfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromHEQ(std::move(v));
        })
        .def("LSRfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromHGP(std::move(v));
        })
        .def("LSRfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromHCA(std::move(v));
        })
        .def("LSRfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromGCA(std::move(v));
        })
        .def("LSRfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.LSRfromGCY(std::move(v));
        })
        .def("GCAfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromHEQ(std::move(v));
        })
        .def("GCAfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromHGP(std::move(v));
        })
        .def("GCAfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromHCA(std::move(v));
        })
        .def("GCAfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromLSR(std::move(v));
        })
        .def("GCAfromGCY", [](OmniCoords &self, vec6 &v) {
            return self.GCAfromGCY(std::move(v));
        })
        .def("GCYfromHEQ", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromHEQ(std::move(v));
        })
        .def("GCYfromHGP", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromHGP(std::move(v));
        })
        .def("GCYfromHCA", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromHCA(std::move(v));
        })
        .def("GCYfromLSR", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromLSR(std::move(v));
        })
        .def("GCYfromGCA", [](OmniCoords &self, vec6 &v) {
            return self.GCYfromGCA(std::move(v));
        });
}

