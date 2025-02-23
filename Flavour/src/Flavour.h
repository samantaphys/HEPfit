/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOUR_H
#define FLAVOUR_H

class StandardModel;
class HeffDF2;
class HeffDF1_Plepnu;
class HeffDB1;
class HeffDC1;
class HeffDS1;
class MVll;
class MPll;
class MVgamma;
class MVlnu;
class MPlnu;
#include "QCD.h"
#include <boost/tuple/tuple.hpp>
#include <memory>

/**
 * @class Flavour
 * @ingroup Flavour
 * @brief The parent class in Flavour for calculating all the Wilson coefficients for various Flavor Violating processes.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The Flavour class aggregates the Wilson coefficients for each of the processes generated by the model by calling the Hamiltonians.
 */

class Flavour {
public:

    /**
     * @brief The constructor.
     * @param[in] SM_i a reference to an object of the class StandardModel
     */
    Flavour(const StandardModel& SM_i);

    /**
     * @brief The member that returns an object of the class HeffDF2.
     * @return returns the Hamiltonian for the \f$ \Delta F = 2 \f$ processes
     *
     */
    HeffDF2& getHDF2() const;
    
    
    /**
     * @brief The member that returns an object of the class HeffDF1_Plepnu.
     * @return returns the Hamiltonian for the \f$ \Delta F = 1 \f$ processes pion leptonic decay
     *
     */
    HeffDF1_Plepnu& getHDF1_Plepnu() const;

    /**
     * @brief The member that returns an object of the class HeffDS1.
     * @return returns the Hamiltonian for the \f$ \Delta S = 1 \f$ processes.
     *
     */
    HeffDS1& getHDS1() const;
    
    
    /**
     * @brief The member that returns an object of the class HeffDC1.
     * @return returns the Hamiltonian for the \f$ \Delta C = 1 \f$ processes.
     *
     */
    HeffDC1& getHDC1() const;

    /**
     * @brief The member that returns an object of the class HeffDB1.
     * @return returns the Hamiltonian for the \f$ \Delta B = 1 \f$ processes.
     *
     */
    HeffDB1& getHDB1() const;

    /**
     * @brief Computes the Wilson coefficient for the process \f$ B_d \to \mu \mu \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the Wilson coefficients for the process \f$ B_d \to \mu \mu \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBd(double mu, schemes scheme = NDR) const;

    gslpp::vector<gslpp::complex>** ComputeCoeffDS1pnunu() const;

    gslpp::vector<gslpp::complex>** ComputeCoeffDS1mumu() const;

    /**
     * @brief Computes the Wilson coefficient for the process \f$ B_s \to \mu \mu \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the Wilson coefficients for the process \f$ B_s \to \mu \mu \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBs(double mu, schemes scheme = NDR, bool SM = false) const;

    gslpp::vector<gslpp::complex>** ComputeCoeffdd(double mu, schemes scheme = NDR) const;

    gslpp::vector<gslpp::complex>** ComputeCoeffK(double mu, schemes scheme = NDR) const;

    gslpp::vector<gslpp::complex>** ComputeCoeffmK(double mu, schemes scheme = NDR) const;

    gslpp::vector<gslpp::complex>** ComputeCoeffDS1PPv(double mu, schemes scheme = NDR) const;
    gslpp::vector<gslpp::complex>** ComputeCoeffDS1PPz(double muc, schemes scheme = NDR) const;
    
    /**
     * @brief Computes the Wilson coefficient for the process \f$ B_s \to \mu \mu \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the Wilson coefficients for the process \f$ B_s \to \mu \mu \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffsmumu(double mu, schemes scheme = NDR) const;

    /**
     * @brief Computes the Wilson coefficient for the process \f$ B_d \to \mu \mu \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the Wilson coefficients for the process \f$ B_d \to \mu \mu \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffdmumu(double mu, schemes scheme = NDR) const;

    
    gslpp::vector<gslpp::complex>** ComputeCoeffcleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) const;
    
    
    gslpp::vector<gslpp::complex>** ComputeCoeffsleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) const;
    
    
    
    gslpp::vector<gslpp::complex>** ComputeCoeffuleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) const;
    
    
    
    /**
     * @brief Computes the Wilson coefficient for the process \f$ B_d \to \tau \nu \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the Wilson coefficients for the process \f$ B_d \to \tau \nu \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffbtaunu(QCD::meson meson_i) const;

    gslpp::vector<gslpp::complex>** ComputeCoeffsnunu() const;

    gslpp::vector<gslpp::complex>** ComputeCoeffdnunu() const;

    /**
     * @brief Computes the Wilson coefficient for the process \f$ b \to s \gamma \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the Wilson coefficients for the process \f$ b \to s \gamma \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffsgamma(double mu, bool noSM = false, schemes scheme = NDR) const;

    /**
     * @brief Computes the chirality flipped Wilson coefficient for the process \f$ b \to s \gamma \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the chirality flipped Wilson coefficients for the process \f$ b \to s \gamma \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffprimesgamma(double mu, schemes scheme = NDR) const;

    /**
     * @brief Computes the Wilson coefficient for the process \f$ B \to V/P \ell^+ \ell^- \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the Wilson coefficients for the process \f$ B \to V/P \ell^+ \ell^- \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBMll(double mu, QCD::lepton lepton, bool noSM = false, schemes scheme = NDR) const;

    /**
     * @brief Computes the chirality flipped Wilson coefficient for the process \f$ B \to V/P \ell^+ \ell^- \f$.
     * @param[in] mu the lower matching scale for the process
     * @param[in] scheme the scheme in which the Wilson Coefficients need to be calculated
     * @return returns the chirality flipped Wilson coefficients for the process \f$ B \to V/P \ell^+ \ell^- \f$
     *
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffprimeBMll(double mu, QCD::lepton lepton, schemes scheme = NDR) const;

    /**
     * @brief Returns the initial and final state dependent object for \f$ B \to V \ell^+ \ell^- \f$.
     * @param[in] meson_i specifies the meson in the initial state
     * @param[in] vector_i specifies the vector in the final state
     * @param[in] lepton_i specifies the lepton in the final state
     * @return returns a pointer to the initial and final state dependent object for the process \f$ B \to V \ell^+ \ell^- \f$
     *
     */
    MVll& getMVll(QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) const;

    /**
     * @brief Returns the initial and final state dependent object for \f$ B \to P \ell^+ \ell^- \f$.
     * @param[in] meson_i specifies the meson in the initial state
     * @param[in] pseudoscalar_i specifies the vector in the final state
     * @param[in] lepton_i specifies the lepton in the final state
     * @return returns a pointer to the initial and final state dependent object for the process \f$ B \to P \ell^+ \ell^- \f$
     *
     */
    MPll& getMPll(QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) const;

    /**
     * @brief Returns the initial and final state dependent object for \f$ B \to V \gamma \f$.
     * @param[in] meson_i specifies the meson in the initial state
     * @param[in] vector_i specifies the vector in the final state
     * @return returns a pointer to the initial and final state dependent object for the process \f$ B \to V \gamma \f$
     *
     */
    MVgamma& getMVgamma(QCD::meson meson_i, QCD::meson vector_i) const;

    /**
     * @brief Returns the initial and final state dependent object for \f$ B \to V \ell \nu \f$.
     * @param[in] meson_i specifies the meson in the initial state
     * @param[in] vector_i specifies the vector in the final state
     * @param[in] lepton_i specifies the lepton in the final state
     * @return returns a pointer to the initial and final state dependent object for the process \f$ B \to V \ell^+ \ell^- \f$
     *
     */
    MVlnu& getMVlnu(QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) const;

    /**
     * @brief Returns the initial and final state dependent object for \f$ B \to P \ell \nu \f$.
     * @param[in] meson_i specifies the meson in the initial state
     * @param[in] pseudoscalar_i specifies the vector in the final state
     * @param[in] lepton_i specifies the lepton in the final state
     * @return returns a pointer to the initial and final state dependent object for the process \f$ B \to V \ell \nu \f$
     *
     */
    MPlnu& getMPlnu(QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) const;

    /**
     * @brief sets the update flag for the initial and final state dependent object for \f$ B \to V \ell^+ \ell^- \f$.
     * @param[in] meson_i specifies the meson in the initial state
     * @param[in] vector_i specifies the vector in the final state
     * @param[in] lepton_i specifies the lepton in the final state
     *
     */
    void setUpdateFlag(QCD::meson meson_i, QCD::meson meson_j, QCD::lepton lep_i, bool updated_i) const;

    /**
     * @brief gets the update flag for the initial and final state dependent object for \f$ B \to V \ell^+ \ell^- \f$.
     * @param[in] meson_i specifies the meson in the initial state
     * @param[in] vector_i specifies the vector in the final state
     * @param[in] lepton_i specifies the lepton in the final state
     *
     */
    bool getUpdateFlag(QCD::meson meson_i, QCD::meson meson_j, QCD::lepton lep_i) const;

    /**
     * @brief a member used for the caching for \f$ B \to V \ell^+ \ell^- \f$.
     *
     */
    void setSMupdated() const;
    
    bool setFlag(const std::string name, const bool value);

    bool setFlagUseDispersionRelation(bool dispersion)
    {
        return (this->dispersion = dispersion);
    }

    bool setFlagCLN(bool CLNflag)
    {
        return (this->CLNflag = CLNflag);
    }

    bool setFlagFixedWCbtos(bool FixedWCbtosflag)
    {
        return (this->FixedWCbtosflag = FixedWCbtosflag);
    }

    bool getFlagUseDispersionRelation() const
    {
        return dispersion;
    }

    bool getFlagCLN() const
    {
        return CLNflag;
    }

    bool getFlagFixedWCbtos() const
    {
        return FixedWCbtosflag;
    }

private:
    template<typename T, typename... Args> std::shared_ptr<T>& getPtr(std::shared_ptr<T>& x, Args... args) const;
    template <typename T, typename... Args> T& getM(std::map<std::vector<int>,std::shared_ptr<T> >& map, Args ... args) const;
    const StandardModel & mySM;
    mutable std::shared_ptr<HeffDF2> HDF2; ///< An Object for the Hamiltonian of the \f$ \Delta F = 2 \f$ processes.
    mutable std::shared_ptr<HeffDF1_Plepnu> HDF1_Plepnu; ///< An Object for the Hamiltonian of the \f$ \Delta F = 1 \f$ processes Pion leptonic decay.
    mutable std::shared_ptr<HeffDB1> HDB1; ///< An Object for the Hamiltonian of the \f$ \Delta B = 1 \f$ processes.
    mutable std::shared_ptr<HeffDC1> HDC1; ///< An Object for the Hamiltonian of the \f$ \Delta C = 1 \f$ processes.
    mutable std::shared_ptr<HeffDS1> HDS1; ///< An Object for the Hamiltonian of the \f$ \Delta S = 1 \f$ processes.
    mutable std::map<std::vector<int>, std::shared_ptr<MVll> > MVllMap;
    mutable std::map<std::vector<int>, std::shared_ptr<MVlnu> > MVlnuMap;
    mutable std::map<std::vector<int>, std::shared_ptr<MVgamma> > MVgammaMap;
    mutable std::map<std::vector<int>, std::shared_ptr<MPll> > MPllMap;
    mutable std::map<std::vector<int>, std::shared_ptr<MPlnu> > MPlnuMap;
    mutable std::map<std::vector<int>, bool> flagUpdateMap;
    
    mutable bool dispersion;
    mutable bool CLNflag;
    mutable bool FixedWCbtosflag;
};

#endif /* FLAVOUR_H */
