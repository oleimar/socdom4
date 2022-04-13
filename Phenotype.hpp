#ifndef PHENOTYPE_HPP
#define PHENOTYPE_HPP

#include <string>
#include <array>
#include <ostream>
#include <istream>
#include <cmath>

// The EvoDom program runs evolutionary simulations
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//************************** struct LearnPars ******************************

// This struct stores an individual's current learning parameters, while
// learning from dominance interactions; the struct also stores the number of
// wins, losses and draws

template<typename T>
struct LearnPars {
    using val_type = T;
    using par_type = std::vector<T>;
    using ui_type = std::vector<unsigned>;
    LearnPars(unsigned lg_size,
              T w0 = 0, T th0 = 0) :
        w(lg_size, w0),
        th(lg_size, th0),
        wins(lg_size, 0),
        losses(lg_size, 0),
        draws(lg_size, 0) {}
    par_type w;
    par_type th;
    ui_type wins;
    ui_type losses;
    ui_type draws;
};


//************************* struct Phenotype ******************************

// Assumptions about GenType:
// types:
//   val_type
// member functions:
//   val_type Value()

// In addition to the genotypic trait values, consisting of w0, th0, g0, ga0,
// alphw, alphth, beta, v, gf, mf, cif, this class also stores the individual
// quality q, the 'perceived own fighting ability' y, the 'displayed fighting
// ability' z, the average per round payoff, the current 'Elo score', the
// current value of the group learning parameters (when saved this will be the
// value after the specified number of rounds of interaction during a
// generation), the number of rounds of interaction, individual number in
// group, the group number, the individual's sex, and whether it is alive.

template<typename GenType>
struct Phenotype {
// public:
    using flt = float;
    using v_type = std::vector<flt>;
    using lp_type = LearnPars<flt>;
    using gen_type = GenType;
    using val_type = typename gen_type::val_type;
    Phenotype(unsigned lg_size = 0,
              const gen_type& gt = gen_type()) :
        lp(lg_size) { Assign(gt); }
    void Assign(const gen_type& gt);
    // estimated value
    flt vhat(flt xi, unsigned j) const
    {
        return gf*lp.w[inum] + (1 - gf)*lp.w[j] + g0*xi;
    }
    // logit of probability of choosing A (fight)
    flt l(flt xi, unsigned j) const
    {
        return gf*lp.th[inum] + (1 - gf)*lp.th[j] + ga0*xi;
    }
    // probability of choosing A (fight)
    flt p(flt xi, unsigned j) const
    {
        flt h = gf*lp.th[inum] + (1 - gf)*lp.th[j] + ga0*xi;
        return 1/(1 + std::exp(-h));
    }
    bool Female() const { return female; }
    // public data members
    flt w0;      // parameter for estimated reward at start of generation
    flt th0;     // parameter for preference at start of generation
    flt g0;      // parameter for estimated reward at start of generation
    flt ga0;     // parameter for preference at start of generation
    flt alphw;   // learning rate
    flt alphth;  // learning rate
    flt beta;    // bystander learning rate factor
    flt v;       // perceived value of aggressive behaviour
    flt gf;      // generalisation factor
    flt mf;      // memory factor (between years)
    flt cif;     // strength of "competitive interference"
    flt dq;      // individual quality deviation (same for each age)
    flt q;       // individual quality (fighting ability)
    flt y;       // own perception of quality
    flt z;       // displayed quality
    flt dmg;     // damage: factor for mortality rates
    flt ares;    // acquired resources (current year)
    flt B0;      // multiplier for own decrease from own interference
    flt B1;      // multiplier for subordinate decrease from own interference
    flt D0;      // multiplier for 'personal' decrease of expected reproduction
    flt D1;      // multiplier for 'group' decrease of expected reproduction
    flt elor;    // version of Elo rating
    lp_type lp;         // learning parameters
    unsigned age;       // the age of the individual
    unsigned dage;      // the age of the individual in the year of death
    unsigned nInts;     // number of interactions experienced
    unsigned nRnds;     // number of rounds experienced
    unsigned nAA;       // number of AA interactions experienced
    unsigned ndom;      // number of group members dominated over
    unsigned nOffspr;   // lifetime number of offspring (so far)
    unsigned inum;      // individual number in local group
    unsigned yr1;       // the year the individual was yearling
    unsigned gnum;      // local group number
    int rank;           // dominance rank (1 is highest)
    bool female;
    bool alive;
};

template<typename GenType>
void Phenotype<GenType>::Assign(const gen_type& gt)
{
    val_type val = gt.Value();
    // assume val is a vector with components corresponding to the traits w0,
    // th0, g0, ga0, alphw, alphth, beta, v, gf, mf, cif
    w0 = val[0];
    th0 = val[1];
    g0 = val[2];
    ga0 = val[3];
    alphw = val[4];
    alphth = val[5];
    beta = val[6];
    v = val[7];
    gf = val[8];
    mf = val[9];
    cif = val[10];
    dq = 0;
    q = 0;
    y = 0;
    z = 0;
    dmg = 0;
    ares = 0;
    B0 = 1;
    B1 = 1;
    D0 = 1;
    D1 = 1;
    elor = 0;
    for (unsigned j = 0; j < lp.w.size(); ++j) {
        lp.w[j] = w0;
        lp.th[j] = th0;
        lp.wins[j] = 0;
        lp.losses[j] = 0;
        lp.draws[j] = 0;
    }
    age = 0;
    dage = 0;
    nInts = 0;
    nRnds = 0;
    nAA = 0;
    ndom = 0;
    nOffspr = 0;
    inum = 0;
    yr1 = 0;
    gnum = 0;
    rank = 0;
    female = true;
    alive = false;
}

#endif // PHENOTYPE_HPP
