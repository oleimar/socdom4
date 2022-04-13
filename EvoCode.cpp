#include "cpptoml.h"   // to read input parameters from TOML file
#include "EvoCode.hpp"
#include "hdf5code.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#ifdef PARA_RUN
#include <omp.h>
#endif

// The EvoDom program runs evolutionary simulations
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//************************** Read and ReadArr ****************************

// convenience functions to read from TOML input file

// this template function can be used for any type of single value
template<typename T>
void Get(std::shared_ptr<cpptoml::table> infile,
         T& value, const std::string& name)
{
    auto val = infile->get_as<T>(name);
    if (val) {
        value = *val;
    } else {
        std::cerr << "Read failed for identifier " << name << "\n";
    }
}

// this template function can be used for a vector or array (but there is no
// checking how many elements are read)
template<typename It>
void GetArr(std::shared_ptr<cpptoml::table> infile,
            It beg, const std::string& name)
{
    using valtype = typename std::iterator_traits<It>::value_type;
    auto vp = infile->get_array_of<valtype>(name);
    if (vp) {
        std::copy(vp->begin(), vp->end(), beg);
    } else {
        std::cerr << "Read failed for identifier " << name << "\n";
    }
}


//************************** class EvoInpData ****************************

EvoInpData::EvoInpData(const char* filename) :
      OK(false)
{
    auto idat = cpptoml::parse_file(filename);
    Get(idat, max_num_thrds, "max_num_thrds");
    Get(idat, num_loci, "num_loci");
    Get(idat, G, "G");
    Get(idat, M, "M");
    Get(idat, n1, "n1");
    Get(idat, np1, "np1");
    Get(idat, numI1, "numI1");
    Get(idat, irndx, "irndx");
    Get(idat, minrnds, "minrnds");
    Get(idat, domrnds, "domrnds");
    Get(idat, drawrnds, "drawrnds");
    Get(idat, numV1, "numV1");
    Get(idat, numyrs, "numyrs");
    Get(idat, nrci, "nrci");
    Get(idat, bystl, "bystl");
    Get(idat, degloc, "degloc");
    Get(idat, areslim, "areslim");
    Get(idat, V0, "V0");
    Get(idat, cdmg, "cdmg");
    Get(idat, c0, "c0");
    Get(idat, c1, "c1");
    Get(idat, mu1, "mu1");
    Get(idat, eta1, "eta1");
    Get(idat, b0, "b0");
    Get(idat, b1, "b1");
    Get(idat, phi, "phi");
    Get(idat, d0, "d0");
    Get(idat, d1, "d1");
    Get(idat, psi, "psi");
    Get(idat, a0, "a0");
    Get(idat, a1, "a1");
    Get(idat, sq, "sq");
    Get(idat, sqm, "sqm");
    Get(idat, sigma, "sigma");
    Get(idat, sigp, "sigp");
    Get(idat, pmax, "pmax");
    Get(idat, betelo, "betelo");
    Get(idat, pld, "pld");
    V1.resize(numV1);
    GetArr(idat, V1.begin(), "V1");
    Qm.resize(M);
    GetArr(idat, Qm.begin(), "Qm");
    mut_rate.resize(num_loci);
    GetArr(idat, mut_rate.begin(), "mut_rate");
    SD.resize(num_loci);
    GetArr(idat, SD.begin(), "SD");
    max_val.resize(num_loci);
    GetArr(idat, max_val.begin(), "max_val");
    min_val.resize(num_loci);
    GetArr(idat, min_val.begin(), "min_val");
    rho.resize(num_loci);
    GetArr(idat, rho.begin(), "rho");
    Get(idat, hist, "hist");
    Get(idat, read_from_file, "read_from_file");
    if (read_from_file) {
        Get(idat, h5InName, "h5InName");
    } else {
        all0.resize(num_loci);
        GetArr(idat, all0.begin(), "all0");
    }
    Get(idat, h5OutName, "h5OutName");
    if (hist) {
        Get(idat, h5HistName, "h5HistName");
    }
    InpName = std::string(filename);
    OK = true;
}


//****************************** Class Evo *****************************

Evo::Evo(const EvoInpData& eid) :
    id{eid},
    num_loci{id.num_loci},
    G{id.G},
    M{id.M},
    n1{id.n1},
    np1{id.np1},
    Nl{n1*M},
    Npl{np1*M},
    N1{n1*G},
    Np1{np1*G},
    N{Nl*G},
    Np{Npl*G},
    numI1{id.numI1},
    irndx{id.irndx},
    minrnds{id.minrnds},
    domrnds{id.domrnds},
    drawrnds{id.drawrnds},
    numV1{id.numV1},
    numyrs{id.numyrs},
    nrci{id.nrci},
    bystl{id.bystl},
    degloc{static_cast<flt>(id.degloc)},
    areslim{static_cast<flt>(id.areslim)},
    V0{static_cast<flt>(id.V0)},
    cdmg{static_cast<flt>(id.cdmg)},
    c0{static_cast<flt>(id.c0)},
    c1{static_cast<flt>(id.c1)},
    mu1{static_cast<flt>(id.mu1)},
    eta1{static_cast<flt>(id.eta1)},
    b0{static_cast<flt>(id.b0)},
    b1{static_cast<flt>(id.b1)},
    phi{static_cast<flt>(id.phi)},
    d0{static_cast<flt>(id.d0)},
    d1{static_cast<flt>(id.d1)},
    psi{static_cast<flt>(id.psi)},
    a0{static_cast<flt>(id.a0)},
    a1{static_cast<flt>(id.a1)},
    sq{static_cast<flt>(id.sq)},
    sqm{static_cast<flt>(id.sqm)},
    sigma{static_cast<flt>(id.sigma)},
    sigp{static_cast<flt>(id.sigp)},
    pmax{static_cast<flt>(id.pmax)},
    betelo{static_cast<flt>(id.betelo)},
    pld{static_cast<flt>(id.pld)},
    V1{id.V1.begin(), id.V1.end()},
    Qm{id.Qm.begin(), id.Qm.end()},
    hist{id.hist},
    num_thrds{1}
{
    // decide on number of threads for parallel processing
#ifdef PARA_RUN
    num_thrds = omp_get_max_threads();
    if (num_thrds > id.max_num_thrds) num_thrds = id.max_num_thrds;
    // check that there is at least one subpopulation per thread
    // if (num_thrds > nsp) num_thrds = nsp;
    std::cout << "Number of threads: "
              << num_thrds << '\n';
#endif
    // generate one seed for each thread
    sds.resize(num_thrds);
    std::random_device rd;
    for (unsigned i = 0; i < num_thrds; ++i) {
        // set up thread-local to be random number engines
        sds[i] = rd();
        vre.push_back(rand_eng(sds[i]));
    }
    for (unsigned i = 0; i < num_thrds; ++i) {
        // set up thread-local to be mutation records, with thread-local engine
        // and parameters controlling mutation, segregation and recombination
        rand_eng& eng = vre[i];
        mut_rec_type mr(eng, num_loci);
        for (unsigned l = 0; l < num_loci; ++l) {
            mr.mut_rate[l] = static_cast<flt>(id.mut_rate[l]);
            mr.SD[l] = static_cast<flt>(id.SD[l]);
            mr.max_val[l] = static_cast<flt>(id.max_val[l]);
            mr.min_val[l] = static_cast<flt>(id.min_val[l]);
            mr.rho[l] = static_cast<flt>(id.rho[l]);
        }
        vmr.push_back(mr);
    }

    int threadn = 0;

    // Note concerning thread safety: in order to avoid possible problems with
    // multiple threads, the std::vector containers f_pop, m_pop, and stat are
    // allocated once and for all here, and thread-local data are then copied
    // into position in these (thus avoiding potentially unsafe push_back and
    // insert).

    // Create N female "placeholder individuals" in population (based on the
    // constructor these are not alive). The convention is that individuals in
    // local group k (k = 0, ..., G-1) are found as alive phenotypes at
    // f_pop[j] with j = k*Nl, k*Nl+1, ..., k*Nl+Nl-1, i.e. j = k*Nl+i with i =
    // 0, ..., Nl-1. These values of i (i.e., i = 0, ..., Nl-1) correspond to
    // the original values of inum in each group gnum = k, assigned to
    // individuals when read_from_file is false. Further, when an individual
    // with age beyond the maximum is replaced by a yearling (in UpdateAgeStr),
    // the new individual 'inherits' the old inum.
    gam_type gam(num_loci);
    ind_type f_indi(Nl, gam);
    f_pop.resize(N, f_indi);
    // create Np male "placeholder genotypes" as male population
    gen_type gen(gam);
    m_pop.resize(Np, gen);
    // also create female and male offspring placeholders
    f_offspr.resize(N1, gen);
    m_offspr.resize(Np1, gen);
    // history stats
    if (hist) {
        stat.reserve(irndx);
    }

    // check if population data should be read from file
    if (id.read_from_file) {
        // Read_pop(id.InName);
        h5_read_pop(id.h5InName);
        UpdateAgeStr(0, vre[threadn]);
    } else {
        rand_eng& eng = vre[threadn];
        rand_norm nrq(0, sq);
        rand_norm nrqm(0, sqm);
        // construct all individuals with the same genotypes
        gam_type gam(num_loci); // starting gamete
        for (unsigned l = 0; l < num_loci; ++l) {
            gam.gamdat[l] = static_cast<flt>(id.all0[l]);
        }
        unsigned j = 0;
        unsigned jp = 0;
        for (unsigned k = 0; k < G; ++k) { // local groups
            for (unsigned m = 0; m < M; ++m) { // cohorts in groups
                // construct n1 females (genotype and phenotype) in each cohort
                for (unsigned i = 0; i < n1; ++i) { // females in cohort
                    ind_type ind(Nl, gam);
                    phen_type& ph = ind.phenotype;
                    ph.dq = nrq(eng);
                    ph.q = Qm[m] + ph.dq + nrqm(eng);
                    ph.y = a0*ph.q;
                    ph.z = a1*ph.q;
                    ph.gnum = k;         // set local group number
                    ph.age = m;          // set age
                    ph.dage = M;         // death age if survival over all ages
                    ph.inum = m*n1 + i;  // set individual number
                    ph.female = true;
                    ph.alive = true;
                    f_pop[j++] = ind;
                }
                // construct np1 male genotypes in each cohort
                for (unsigned ip = 0; ip < np1; ++ip) { // males in cohort
                    gen_type gen(gam);
                    m_pop[jp++] = gen;
                }
            }
        }
    }
}

void Evo::Run()
{
    Timer timer(std::cout);
    timer.Start();
    ProgressBar PrBar(std::cout, numyrs);
    // run through years
    // Time sequence within a year: (i)
    for (unsigned yr = 0; yr < numyrs; ++yr) {
        // use parallel for processing over the male-male interactions in local
        // groups, within a year
#pragma omp parallel for num_threads(num_thrds)
        for (unsigned k = 0; k < G; ++k) {
#ifdef PARA_RUN
            int threadn = omp_get_thread_num();
#else
            int threadn = 0;
#endif
            // thread-local random number engine
            rand_eng& eng = vre[threadn];
            // thread-local mutation record
            mut_rec_type& mr = vmr[threadn];
            // distribution needed
            rand_uni uni(0, 1);
            // thread-local container for local group individuals (females)
            vind_type lfpop;
            lfpop.reserve(Nl);
            for (unsigned i = 0; i < Nl; ++i) {
                unsigned j = k*Nl + i;
                lfpop.push_back(f_pop[j]);
            }
            // thread-local container for local group phenotypes
            vph_type lgph;
            lgph.reserve(Nl);
            // vector of indices in lfpop of individuals in container lgph
            ui_type lgi;
            lgi.reserve(Nl);
            for (unsigned i = 0; i < Nl; ++i) {
                phen_type& ph = lfpop[i].phenotype;
                if (ph.alive) { // only include alive individuals
                    lgph.push_back(ph);
                    lgi.push_back(i);
                }
            }
            // get history only for single threaded and final M years
            bool lhist = false;
            if (num_thrds == 1 && yr >= numyrs - M && hist) {
                lhist = true;
            }
            // form competition group
            ui_type cgi;
            for (unsigned i = 0; i < lgph.size(); ++i) {
                if (lgph[i].alive) { // this check is not really needed
                    cgi.push_back(i);
                }
            }
            unsigned g = cgi.size(); // competition group size
            if (g > 0) {
                // randomly shuffle indices (to avoid order effects)
                std::shuffle(cgi.begin(), cgi.end(), eng);
                vph_type cgph;
                cgph.reserve(g);
                for (unsigned i = 0; i < g; ++i) {
                    cgph.push_back(lgph[cgi[i]]);
                }
                // construct competition group
                cg_type cg(g, numI1,
                    irndx, minrnds, domrnds, drawrnds, yr,
                    cdmg, sigma, sigp, pmax, betelo, cgph, bystl, lhist);
                cg.Interact(eng);
                vph_type& memb = cg.Get_memb();


                const ui_type& ndom = cg.Get_ndom();
                // Apply mortality from damage (so far, all members are alive)
                int nalv = 0;      // count number alive
                ui_type ialv;      // indices in memb of alive individuals
                for (unsigned i = 0; i < g; ++i) {
                    phen_type& ph = memb[i];
                    if (uni(eng) >= std::exp(-c0*ph.dmg)) {
                        ph.alive = false;
                        ph.dage = ph.age;
                        ph.rank = 0;
                        ph.ndom = 0;
                    } else {
                        ialv.push_back(i);
                        ++nalv;
                    }
                }
                // assign resources; the basic assumption is that (neglecting
                // mortality), for a competition group there are amounts V1[k]
                // based on dominance positions k, and also the amount V0 each
                // receives irrespective of dominance; ndom[i] is the number of
                // individuals among the g that individual i has established
                // dominance over (thus, maximal value of ndom[i] is g-1);
                // however, with mortality dominance positions and resource
                // assignments should be distributed among alive individuals
                if (nalv > 1) {
                    ui_type idx(nalv, 0);
                    std::iota(idx.begin(), idx.end(), 0);
                    // sort the indices in idx based on ndom[ialv], from highest
                    // to lowest
                    std::stable_sort(idx.begin(), idx.end(),
                        [&ndom, &ialv](unsigned i1, unsigned i2)
                            { return ndom[ialv[i1]] > ndom[ialv[i2]]; });
                    // // NOTE; test printout of ndom and idx
                    // for (unsigned i = 0; i < g; ++i) {
                    //     phen_type& mi = memb[i];
                    //     std::cout << mi.inum + 1 << "  ";
                    // }
                    // std::cout << '\n';
                    // for (unsigned i = 0; i < g; ++i) {
                    //     std::cout << ndom[i] << "  ";
                    // }
                    // std::cout << '\n';
                    // for (unsigned i = 0; i < g; ++i) {
                    //     phen_type& mii = memb[idx[i]];
                    //     std::cout << mii.inum + 1 << "  ";
                    // }
                    // std::cout << '\n';
                    // std::cout << '\n';
                    // // NOTE: end test printout
                    if (ndom[idx[0]] == 0) {
                        // no dominance hierarchy, construct random positions
                        std::shuffle(idx.begin(), idx.end(), eng);
                    }
                    // apportion V1 values, from V1[numV1 -1] and downwards, in
                    // dominance order (NOTE: there is a possible issue with
                    // tied positions)
                    for (int i = 0; i < nalv; ++i) {
                        phen_type& mii = memb[ialv[idx[i]]];
                        flt v1 = V1[numV1 - 1 - i];
                        mii.ares = v1;
                        mii.rank = i + 1;
                    }
                    // adjust acquired resources for interference by dominants
                    // in subordinate resource acquisition
                    for (int i = 0; i < nalv; ++i) {
                        phen_type& mii = memb[ialv[idx[i]]];
                        mii.B0 = std::exp(-b0*mii.cif);
                        mii.B1 = 1 - phi + phi*std::exp(-b1*mii.cif);
                        if (i + 1 < nalv) { // there are subordinates
                            mii.ares *= mii.B0; // decrease of own ares
                            int jhi = std::min(nalv - i, nrci + 1);
                            for (int j = 1; j < jhi; ++j) {
                                // reduction of subordinate ares
                                phen_type& mij = memb[ialv[idx[i + j]]];
                                mij.ares *= mii.B1;
                            }
                        }
                    }
                } else if (nalv > 0) { // only one alive group member
                    flt v1 = V1[numV1 - 1];
                    memb[ialv[0]].ares = v1;
                    memb[ialv[0]].rank = 1;
                }
                // adjust for decrease in expected reproduction from damage and
                // death
                flt D1 = 1; // 'collective' decrease
                for (unsigned i = 0; i < g; ++i) {
                    const phen_type& ph = memb[i];
                    if (ph.alive) {
                        D1 *= 1 - psi + psi*std::exp(-d1*ph.dmg);
                    } else {
                        // effect of absence of individual
                        D1 *= 1 - psi;
                    }
                }
                for (unsigned i = 0; i < g; ++i) {
                    phen_type& ph = memb[i];
                    if (ph.alive) {
                        ph.D0 = std::exp(-d0*ph.dmg);
                        ph.D1 = D1;
                        ph.ares *= ph.D0*ph.D1;
                    }
                }
                if (lhist) {
                    // append history
                    const vs_type& st = cg.Get_stat();
                    stat.insert(stat.end(), st.begin(), st.end());
                }
                // put individuals in the mating group back into lgph
                for (unsigned i = 0; i < g; ++i) {
                    lgph[cgi[i]] = memb[i];
                }
            } else {
                // NOTE: all individuals are dead, and to avoid failure of
                // reproduction the current year, allocate expected reproduction
                // to all lowest age individuals in lgph
                for (unsigned i = 0; i < lgph.size(); ++i) {
                    if (lgph[i].age == 0) {
                        lgph[i].ares = 1;
                    }
                }
            }

            // put female phenotypes from lgph into the correct local group
            // individuals
            for (unsigned i = 0; i < lgph.size(); ++i) {
                lfpop[lgi[i]].phenotype = lgph[i];
            }
            // copy females from lfpop back to f_pop
            for (unsigned i = 0; i < Nl; ++i) {
                unsigned j = k*Nl + i;
                f_pop[j] = lfpop[i];
            }
        } // end of parallel for (over local groups)

        // reproduction in entire population
        SelectReproduce(vmr[0]);
        // damage healing and survival from mating season and reproduction to
        // the next year; run through f_pop; check survival
        rand_uni uni(0, 1);
        for (unsigned i = 0; i < N; ++i) {
            phen_type& ph = f_pop[i].phenotype;
            if (ph.alive) {
                // first do healing
                ph.dmg *= 1 - eta1;
                // then do survival
                if (uni(vre[0]) >= std::exp(-mu1 - c1*ph.dmg)) {
                    ph.alive = false;
                    ph.dage = ph.age;
                }
            }
        }
        // if not final year, update age structure
        if (yr < numyrs - 1) {
            UpdateAgeStr(yr + 1, vre[0]);
        }
        // all set to start next year
        ++PrBar;
    }

    PrBar.Final();
    timer.Stop();
    timer.Display();
    h5_write_pop(id.h5OutName);
    if (hist) {
        Set_stat_tm();
        h5_write_hist(id.h5HistName);
    }
}

void Evo::SelectReproduce(mut_rec_type& mr)
{
    // get discrete distribution with female expected reproduction as weights
    v_type wei(N);
    flt tot_wei = 0;
    for (unsigned i = 0; i < N; ++i) {
        const phen_type& ph = f_pop[i].phenotype;
        // set weights for ares below a limit to zero
        wei[i] = (ph.ares < areslim) ? 0 : ph.ares;
        tot_wei += wei[i];
    }
    rand_uni uni(0, 1);
    // get distributions for 'outside option' reproduction
    v_type wei0(N);
    flt tot_wei0 = 0;
    for (unsigned i = 0; i < N; ++i) {
        const phen_type& ph = f_pop[i].phenotype;
        wei0[i] = (ph.alive == 1) ? ph.D0 : 0;
        tot_wei0 += wei0[i];
    }
    if (tot_wei > 0) {
        if (degloc > 0) {
            // adjust weights in wei
            for (unsigned k = 0; k < G; ++k) { // local groups
                flt tot_weil = 0;
                for (unsigned i = 0; i < Nl; ++i) {
                    tot_weil += wei[k*Nl + i];
                }
                if (tot_weil > 0) {
                    flt mult = degloc*tot_wei/(G*tot_weil) + (1 - degloc);
                    for (unsigned i = 0; i < Nl; ++i) {
                        wei[k*Nl + i] *= mult;
                    }
                }
                // for tot_weil == 0, let group reproduction fail (not really
                // purely local competition for degloc == 1)
            }
        }
        rand_discr dsf(wei.begin(), wei.end());
        rand_discr dsf0(wei0.begin(), wei0.end());
        // get uniform integer distribution for male genotypes
        rand_ui dsm(0, Np - 1);
        // get N1 female offspring
        if (pld > 0) {
            // get female offspring for each next-year group
            for (unsigned k = 0; k < G; ++k) { // local groups
                // get acquired resource weights for group
                v_type weil(Nl);
                flt tot_weil = 0;
                for (unsigned i = 0; i < Nl; ++i) {
                    weil[i] = wei[k*Nl + i];
                    tot_weil += weil[i];
                }
                rand_discr dsfl(weil.begin(), weil.end());
                // get resource-weight adjusted probability of female offspring
                // being locally derived (does not matter is prloc > 1)
                flt prloc = pld*G*tot_weil/tot_wei;
                // get n1 female offspring
                for (unsigned j = 0; j < n1; ++j) {
                    // find mother of individual to be constructed
                    unsigned imat = 0;
                    if (uni(mr.eng) < prloc) { // local mother
                        imat = k*Nl + dsfl(mr.eng);
                    } else {
                        if (uni(mr.eng) < V0) { // check for 'outside option'
                            imat = dsf0(mr.eng);
                        } else {
                            imat = dsf(mr.eng);
                        }
                    }
                    ind_type& matind = f_pop[imat];
                    // find father genotype for individual to be constructed
                    unsigned ipat = dsm(mr.eng);
                    gen_type& patind = m_pop[ipat];
                    // construct offspring genotype
                    gen_type offspr(matind.GetGamete(mr),
                                    patind.GetGamete(mr));
                    // copy new female genotype to f_offspr
                    f_offspr[k*n1 + j] = offspr;
                    // update mother's record of reproductive success
                    matind.phenotype.nOffspr += 1;
                }
            }
        } else { // global dispersal
            for (unsigned j = 0; j < N1; ++j) {
                // find mother for individual to be constructed
                unsigned imat = 0;
                if (uni(mr.eng) < V0) { // check for 'outside option'
                    imat = dsf0(mr.eng);
                } else {
                    imat = dsf(mr.eng);
                }
                ind_type& matind = f_pop[imat];
                // find father genotype for individual to be constructed
                unsigned ipat = dsm(mr.eng);
                gen_type& patind = m_pop[ipat];
                // construct offspring genotype
                gen_type offspr(matind.GetGamete(mr), patind.GetGamete(mr));
                // copy new female genotype to f_offspr
                f_offspr[j] = offspr;
                // update mother's record of reproductive success
                matind.phenotype.nOffspr += 1;
            }
        }
        // get Np1 male offspring
        for (unsigned j = 0; j < Np1; ++j) {
            // find mother for individual to be constructed
            unsigned imat = 0;
            if (uni(mr.eng) < V0) { // check for 'outside option'
                imat = dsf0(mr.eng);
            } else {
                imat = dsf(mr.eng);
            }
            ind_type& matind = f_pop[imat];
            // find father genotype for individual to be constructed
            unsigned ipat = dsm(mr.eng);
            gen_type& patind = m_pop[ipat];
            // construct offspring genotype
            gen_type offspr(matind.GetGamete(mr), patind.GetGamete(mr));
            // copy new male genotype to m_offspr
            m_offspr[j] = offspr;
            // update mother's record of reproductive success
            matind.phenotype.nOffspr += 1;
        }
    } else {
        std::cout << "NOTE: reproduction failed\n";
    }
}

void Evo::Forget(phen_type& ph)
{
    // move learning parameters of the phenotype towards 'naive' state
    lp_type& lp = ph.lp;
    // deviations of w and th from initial value are
    // multiplied by memory factor
    for (unsigned j = 0; j < lp.w.size(); ++j) {
        // lp.w[j] = ph.w0 + ph.mf*(lp.w[j] - ph.w0);
        // lp.th[j] = ph.th0 + ph.mf*(lp.th[j] - ph.th0);
        lp.w[j] *= ph.mf;
        lp.th[j] *= ph.mf;
    }
}

void Evo::UpdateAgeStr(unsigned yr1, rand_eng& eng)
{
    rand_uni uni(0, 1);
    rand_norm nrq(0, sq);
    rand_norm nrqm(0, sqm);
    // first increase the age of each individual in f_pop, and put those
    // individuals who reach beyond maximum age as dead;  for those below
    // maximum age, adjust the fighting ability, and put ares and rank for all
    // to zero
    for (unsigned i = 0; i < N; ++i) {
        phen_type& ph = f_pop[i].phenotype;
        ph.age += 1;
        if (ph.age >= M) {
            ph.alive = false;
        } else {
            // update fighting ability and apply forgetting for alive
            // individuals
            if (ph.alive) {
                ph.q = Qm[ph.age] + ph.dq + nrqm(eng);
                ph.y = a0*ph.q;
                ph.z = a1*ph.q;
                Forget(ph);
            }
        }
        ph.ares = 0;
        ph.rank = 0;
        ph.ndom = 0;
    }
    // then copy female offspring genotypes from f_offspr into the positions of
    // individuals who reach beyond maximum age and update the learning
    // parameters of other local group members
    unsigned j = 0;
    for (unsigned i = 0; i < N; ++i) {
        phen_type& ph = f_pop[i].phenotype;
        if (ph.age >= M) {
            gen_type& gen = f_offspr[j];
            f_pop[i].genotype = gen;
            // keep individual and local group numbers
            unsigned inum = ph.inum;
            unsigned gnum = ph.gnum;
            ph.Assign(gen);
            // count number of replacements
            ++j;
            // set aspects of phenotype not set correctly by Assign
            ph.dq = nrq(eng);
            ph.q = Qm[ph.age] + ph.dq + nrqm(eng); // ph.age is zero
            ph.y = a0*ph.q;
            ph.z = a1*ph.q;
            ph.dage = M;
            ph.inum = inum;
            ph.yr1 = yr1;
            ph.gnum = gnum;
            ph.alive = true;
            // update relevant learning parameters of older local group females
            for (unsigned k = gnum*Nl; k < (gnum + 1)*Nl; ++k) {
                phen_type& phl = f_pop[k].phenotype;
                if (phl.age > 0 && phl.age < M && phl.inum != inum) {
                    lp_type& lp = phl.lp;
                    lp.w[inum] = phl.w0;
                    lp.th[inum] = phl.th0;
                    lp.wins[inum] = 0;
                    lp.losses[inum] = 0;
                    lp.draws[inum] = 0;
                }
            }
        }
    }
    // NOTE: test printout
    // if (j != N1) {
    //     std::cout << "j = " << j << " problem UpdateAgeStr\n";
    // }
    // then copy (move) males in m_pop to correspond to increased age
    for (unsigned i = Np - 1; i > Npl - 1; --i) {
        m_pop[i] = m_pop[i - Npl];
    }
    // then copy male offspring into the positions corresponding to the
    // youngest age
    for (unsigned i = 0; i < Np1; ++i) {
        m_pop[i] = m_offspr[i];
    }
}

void Evo::Set_stat_tm()
{
    // compute and assign tminyr and tm fields of history in stat
    unsigned hlen = stat.size();
    if (hlen > 0) {
        unsigned kh = 0; // counter for stat elements
        // run through history
        while (kh < hlen) {
            stat_type& st = stat[kh];
            unsigned yr = st.yr;
            unsigned cnum = st.cnum;
            unsigned g = st.g;
            // number of interactions per pair
            unsigned nI = numI1;
            // number of interactions in current competition group
            unsigned nInt = nI*g*(g - 1)/2;
            flt dInt = 1.0/nInt;
            flt tminyr = cnum*dInt;
            st.tminyr = tminyr;
            st.tm = st.yr + tminyr;
            ++kh;
        }
    }
}

void Evo::h5_read_pop(const std::string& infilename)
{
    // read data and put in pop
    h5R h5(infilename);
    // female gametes
    std::vector<v_type> fgams(N, v_type(num_loci));
    // read maternal gametes of females
    h5.read_flt_arr("FlMatGam", fgams);
    for (unsigned i = 0; i < N; ++i) {
        gam_type& gam = f_pop[i].genotype.mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = fgams[i][l];
        }
    }
    // read paternal gametes of females
    h5.read_flt_arr("FlPatGam", fgams);
    for (unsigned i = 0; i < N; ++i) {
        gam_type& gam = f_pop[i].genotype.pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = fgams[i][l];
        }
    }
    // male gametes
    std::vector<v_type> mgams(Np, v_type(num_loci));
    // read maternal gametes of males
    h5.read_flt_arr("MlMatGam", mgams);
    for (unsigned i = 0; i < Np; ++i) {
        gam_type& gam = m_pop[i].mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = mgams[i][l];
        }
    }
    // read paternal gametes of males
    h5.read_flt_arr("MlPatGam", mgams);
    for (unsigned i = 0; i < Np; ++i) {
        gam_type& gam = m_pop[i].pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = mgams[i][l];
        }
    }
    // female offspring gametes
    std::vector<v_type> fogams(N1, v_type(num_loci));
    // read maternal gametes of female offspring
    h5.read_flt_arr("FloMatGam", fogams);
    for (unsigned i = 0; i < N1; ++i) {
        gam_type& gam = f_offspr[i].mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = fogams[i][l];
        }
    }
    // read paternal gametes of female offspring
    h5.read_flt_arr("FloPatGam", fogams);
    for (unsigned i = 0; i < N1; ++i) {
        gam_type& gam = f_offspr[i].pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = fogams[i][l];
        }
    }
    // male offspring gametes
    std::vector<v_type> mogams(Np1, v_type(num_loci));
    // read maternal gametes of male offspring
    h5.read_flt_arr("MloMatGam", mogams);
    for (unsigned i = 0; i < Np1; ++i) {
        gam_type& gam = m_offspr[i].mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = mogams[i][l];
        }
    }
    // read paternal gametes of male offspring
    h5.read_flt_arr("MloPatGam", mogams);
    for (unsigned i = 0; i < Np1; ++i) {
        gam_type& gam = m_offspr[i].pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = mogams[i][l];
        }
    }
    v_type fval(N);
    // w0
    h5.read_flt("w0", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.w0 = fval[i];
    }
    // th0
    h5.read_flt("th0", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.th0 = fval[i];
    }
    // g0
    h5.read_flt("g0", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.g0 = fval[i];
    }
    // ga0
    h5.read_flt("ga0", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.ga0 = fval[i];
    }
    // alphw
    h5.read_flt("alphw", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.alphw = fval[i];
    }
    // alphth
    h5.read_flt("alphth", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.alphth = fval[i];
    }
    // beta
    h5.read_flt("beta", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.beta = fval[i];
    }
    // v
    h5.read_flt("v", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.v = fval[i];
    }
    // gf
    h5.read_flt("gf", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.gf = fval[i];
    }
    // mf
    h5.read_flt("mf", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.mf = fval[i];
    }
    // cif
    h5.read_flt("cif", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.cif = fval[i];
    }
    // dq
    h5.read_flt("dq", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.dq = fval[i];
    }
    // q
    h5.read_flt("q", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.q = fval[i];
    }
    // y
    h5.read_flt("y", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.y = fval[i];
    }
    // z
    h5.read_flt("z", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.z = fval[i];
    }
    // dmg
    h5.read_flt("dmg", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.dmg = fval[i];
    }
    // ares
    h5.read_flt("ares", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.ares = fval[i];
    }
    // B0
    h5.read_flt("B0", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.B0 = fval[i];
    }
    // B1
    h5.read_flt("B1", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.B1 = fval[i];
    }
    // D0
    h5.read_flt("D0", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.D0 = fval[i];
    }
    // D1
    h5.read_flt("D1", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.D1 = fval[i];
    }
    // elor
    h5.read_flt("elor", fval);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.elor = fval[i];
    }
    // std::vector to hold unsigned int member
    ui_type uival(N);
    // age
    h5.read_uint("age", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.age = uival[i];
    }
    // dage
    h5.read_uint("dage", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.dage = uival[i];
    }
    // nInts
    h5.read_uint("nInts", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.nInts = uival[i];
    }
    // nRnds
    h5.read_uint("nRnds", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.nRnds = uival[i];
    }
    // nAA
    h5.read_uint("nAA", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.nAA = uival[i];
    }
    // ndom
    h5.read_uint("ndom", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.ndom = uival[i];
    }
    // nOffspr
    h5.read_uint("nOffspr", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.nOffspr = uival[i];
    }
    // inum
    h5.read_uint("inum", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.inum = uival[i];
    }
    // yr1
    h5.read_uint("yr1", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.yr1 = uival[i];
    }
    // gnum
    h5.read_uint("gnum", uival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.gnum = uival[i];
    }
    // std::vector to hold int (or bool) member
    std::vector<int> ival(N);
    // rank
    h5.read_int("rank", ival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.rank = ival[i];
    }
    // female
    h5.read_int("female", ival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.female = ival[i];
    }
    // alive
    h5.read_int("alive", ival);
    for (unsigned i = 0; i < N; ++i) {
        f_pop[i].phenotype.alive = ival[i];
    }
    // read learning parameters w and th
    std::vector<v_type> pars(N, v_type(Nl));
    // read w parameters
    h5.read_flt_arr("w", pars);
    for (unsigned i = 0; i < N; ++i) {
        v_type& par = f_pop[i].phenotype.lp.w;
        for (unsigned j = 0; j < Nl; ++j) {
            par[j] = pars[i][j];
        }
    }
    // read th parameters
    h5.read_flt_arr("th", pars);
    for (unsigned i = 0; i < N; ++i) {
        v_type& par = f_pop[i].phenotype.lp.th;
        for (unsigned j = 0; j < Nl; ++j) {
            par[j] = pars[i][j];
        }
    }
    // read parameters wins, losses, and draws
    std::vector<ui_type> uipars(N, ui_type(Nl));
    // read wins parameters
    h5.read_uint_arr("wins", uipars);
    for (unsigned i = 0; i < N; ++i) {
        ui_type& uipar = f_pop[i].phenotype.lp.wins;
        for (unsigned j = 0; j < Nl; ++j) {
            uipar[j] = uipars[i][j];
        }
    }
    // read losses parameters
    h5.read_uint_arr("losses", uipars);
    for (unsigned i = 0; i < N; ++i) {
        ui_type& uipar = f_pop[i].phenotype.lp.losses;
        for (unsigned j = 0; j < Nl; ++j) {
            uipar[j] = uipars[i][j];
        }
    }
    // read draws parameters
    h5.read_uint_arr("draws", uipars);
    for (unsigned i = 0; i < N; ++i) {
        ui_type& uipar = f_pop[i].phenotype.lp.draws;
        for (unsigned j = 0; j < Nl; ++j) {
            uipar[j] = uipars[i][j];
        }
    }
}

void Evo::h5_write_pop(const std::string& outfilename) const
{
    h5W h5(outfilename);
    // female gametes
    std::vector<v_type> fgams(N, v_type(num_loci));
    // write maternal gametes of females
    for (unsigned i = 0; i < N; ++i) {
        const gam_type& gam = f_pop[i].genotype.mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            fgams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("FlMatGam", fgams);
    // write paternal gametes of females
    for (unsigned i = 0; i < N; ++i) {
        const gam_type& gam = f_pop[i].genotype.pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            fgams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("FlPatGam", fgams);
    // male gametes
    std::vector<v_type> mgams(Np, v_type(num_loci));
    // write maternal gametes of males
    for (unsigned i = 0; i < Np; ++i) {
        const gam_type& gam = m_pop[i].mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            mgams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("MlMatGam", mgams);
    // write paternal gametes of males
    for (unsigned i = 0; i < Np; ++i) {
        const gam_type& gam = m_pop[i].pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            mgams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("MlPatGam", mgams);
    // female offspring gametes
    std::vector<v_type> fogams(N1, v_type(num_loci));
    // write maternal gametes of female offspring
    for (unsigned i = 0; i < N1; ++i) {
        const gam_type& gam = f_offspr[i].mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            fogams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("FloMatGam", fogams);
    // write paternal gametes of female offspring
    for (unsigned i = 0; i < N1; ++i) {
        const gam_type& gam = f_offspr[i].pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            fogams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("FloPatGam", fogams);
    // male offspring gametes
    std::vector<v_type> mogams(Np1, v_type(num_loci));
    // write maternal gametes of male offspring
    for (unsigned i = 0; i < Np1; ++i) {
        const gam_type& gam = m_offspr[i].mat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            mogams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("MloMatGam", mogams);
    // write paternal gametes of male offspring
    for (unsigned i = 0; i < Np1; ++i) {
        const gam_type& gam = m_offspr[i].pat_gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            mogams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("MloPatGam", mogams);
    // write members of phenotypes
    // std::vector to hold flt member
    v_type fval(N);
    // w0
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.w0; });
    h5.write_flt("w0", fval);
    // th0
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.th0; });
    h5.write_flt("th0", fval);
    // g0
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.g0; });
    h5.write_flt("g0", fval);
    // ga0
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.ga0; });
    h5.write_flt("ga0", fval);
    // alphw
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.alphw; });
    h5.write_flt("alphw", fval);
    // alphth
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.alphth; });
    h5.write_flt("alphth", fval);
    // beta
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.beta; });
    h5.write_flt("beta", fval);
    // v
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.v; });
    h5.write_flt("v", fval);
    // gf
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.gf; });
    h5.write_flt("gf", fval);
    // mf
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.mf; });
    h5.write_flt("mf", fval);
    // cif
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.cif; });
    h5.write_flt("cif", fval);
    // dq
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.dq; });
    h5.write_flt("dq", fval);
    // q
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.q; });
    h5.write_flt("q", fval);
    // y
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.y; });
    h5.write_flt("y", fval);
    // z
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.z; });
    h5.write_flt("z", fval);
    // dmg
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.dmg; });
    h5.write_flt("dmg", fval);
    // ares
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.ares; });
    h5.write_flt("ares", fval);
    // B0
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.B0; });
    h5.write_flt("B0", fval);
    // B1
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.B1; });
    h5.write_flt("B1", fval);
    // D0
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.D0; });
    h5.write_flt("D0", fval);
    // D1
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.D1; });
    h5.write_flt("D1", fval);
    // elor
    std::transform(f_pop.begin(), f_pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.elor; });
    h5.write_flt("elor", fval);
    // std::vector to hold unsigned int member
    ui_type uival(N);
    // age
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.age; });
    h5.write_uint("age", uival);
    // dage
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.dage; });
    h5.write_uint("dage", uival);
    // nInts
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.nInts; });
    h5.write_uint("nInts", uival);
    // nRnds
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.nRnds; });
    h5.write_uint("nRnds", uival);
    // nAA
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.nAA; });
    h5.write_uint("nAA", uival);
    // ndom
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.ndom; });
    h5.write_uint("ndom", uival);
    // nOffspr
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.nOffspr; });
    h5.write_uint("nOffspr", uival);
    // inum
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.inum; });
    h5.write_uint("inum", uival);
    // yr1
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.yr1; });
    h5.write_uint("yr1", uival);
    // gnum
    std::transform(f_pop.begin(), f_pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.gnum; });
    h5.write_uint("gnum", uival);
    // std::vector to hold int (or bool) member
    std::vector<int> ival(N);
    // rank
    std::transform(f_pop.begin(), f_pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.rank; });
    h5.write_int("rank", ival);
    // female
    std::transform(f_pop.begin(), f_pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.female; });
    h5.write_int("female", ival);
    // alive
    std::transform(f_pop.begin(), f_pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.alive; });
    h5.write_int("alive", ival);
    // write learning parameters w and th
    std::vector<v_type> pars(N, v_type(Nl));
    // write w parameters
    for (unsigned i = 0; i < N; ++i) {
        const v_type& par = f_pop[i].phenotype.lp.w;
        for (unsigned j = 0; j < Nl; ++j) {
            pars[i][j] = par[j];
        }
    }
    h5.write_flt_arr("w", pars);
    // write th parameters
    for (unsigned i = 0; i < N; ++i) {
        const v_type& par = f_pop[i].phenotype.lp.th;
        for (unsigned j = 0; j < Nl; ++j) {
            pars[i][j] = par[j];
        }
    }
    h5.write_flt_arr("th", pars);
    // write parameters wins, losses, and draws
    std::vector<ui_type> uipars(N, ui_type(Nl));
    // write wins parameters
    for (unsigned i = 0; i < N; ++i) {
        const ui_type& uipar = f_pop[i].phenotype.lp.wins;
        for (unsigned j = 0; j < Nl; ++j) {
            uipars[i][j] = uipar[j];
        }
    }
    h5.write_uint_arr("wins", uipars);
    // write losses parameters
    for (unsigned i = 0; i < N; ++i) {
        const ui_type& uipar = f_pop[i].phenotype.lp.losses;
        for (unsigned j = 0; j < Nl; ++j) {
            uipars[i][j] = uipar[j];
        }
    }
    h5.write_uint_arr("losses", uipars);
    // write draws parameters
    for (unsigned i = 0; i < N; ++i) {
        const ui_type& uipar = f_pop[i].phenotype.lp.draws;
        for (unsigned j = 0; j < Nl; ++j) {
            uipars[i][j] = uipar[j];
        }
    }
    h5.write_uint_arr("draws", uipars);
}

void Evo::h5_write_hist(const std::string& histfilename) const
{
    h5W h5(histfilename);
    unsigned hlen = stat.size();
    // std::vector to hold unsigned int member
    ui_type uival(hlen);
    // gnum
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.gnum; });
    h5.write_uint("gnum", uival);
    // yr
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.yr; });
    h5.write_uint("yr", uival);
    // g
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.g; });
    h5.write_uint("g", uival);
    // cnum
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.cnum; });
    h5.write_uint("cnum", uival);
    // i
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.i; });
    h5.write_uint("i", uival);
    // j
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.j; });
    h5.write_uint("j", uival);
    // yr1i
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.yr1i; });
    h5.write_uint("yr1i", uival);
    // yr1j
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.yr1j; });
    h5.write_uint("yr1j", uival);
    // agei
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.agei; });
    h5.write_uint("agei", uival);
    // agej
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.agej; });
    h5.write_uint("agej", uival);
    // ui
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.ui; });
    h5.write_uint("ui", uival);
    // uj
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.uj; });
    h5.write_uint("uj", uival);
    // irnds
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.irnds; });
    h5.write_uint("irnds", uival);
    // nAA
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.nAA; });
    h5.write_uint("nAA", uival);
    // nAAi
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.nAAi; });
    h5.write_uint("nAAi", uival);
    // nAAj
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.nAAj; });
    h5.write_uint("nAAj", uival);
    // ndomi
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.ndomi; });
    h5.write_uint("ndomi", uival);
    // ndomj
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.ndomj; });
    h5.write_uint("ndomj", uival);
    // std::vector to hold int member
    std::vector<int> ival(hlen);
    // winij
    std::transform(stat.begin(), stat.end(), ival.begin(),
                   [](const stat_type& st) -> int
                   { return st.winij; });
    h5.write_int("winij", ival);
    // std::vector to hold flt member
    v_type fval(hlen);
    // qi
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.qi; });
    h5.write_flt("qi", fval);
    // qj
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.qj; });
    h5.write_flt("qj", fval);
    // dmgi
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.dmgi; });
    h5.write_flt("dmgi", fval);
    // dmgj
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.dmgj; });
    h5.write_flt("dmgj", fval);
    // lij
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.lij; });
    h5.write_flt("lij", fval);
    // lji
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.lji; });
    h5.write_flt("lji", fval);
    // pij
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.pij; });
    h5.write_flt("pij", fval);
    // pji
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.pji; });
    h5.write_flt("pji", fval);
    // wii
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.wii; });
    h5.write_flt("wii", fval);
    // wij
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.wij; });
    h5.write_flt("wij", fval);
    // wjj
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.wjj; });
    h5.write_flt("wjj", fval);
    // wji
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.wji; });
    h5.write_flt("wji", fval);
    // thii
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.thii; });
    h5.write_flt("thii", fval);
    // thij
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.thij; });
    h5.write_flt("thij", fval);
    // thjj
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.thjj; });
    h5.write_flt("thjj", fval);
    // thji
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.thji; });
    h5.write_flt("thji", fval);
    // elori
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.elori; });
    h5.write_flt("elori", fval);
    // elorj
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.elorj; });
    h5.write_flt("elorj", fval);
    // tminyr
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.tminyr; });
    h5.write_flt("tminyr", fval);
    // tm
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.tm; });
    h5.write_flt("tm", fval);
}
