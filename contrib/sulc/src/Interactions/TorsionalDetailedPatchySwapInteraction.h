/*
 * TorsionalDetailedPatchySwapInteraction.h
 *
 *  Created on: 13/may/2021
 *      Author: lorenzo
 */

#ifndef TorsionalDetailedPatchySwapInteraction_H_
#define TorsionalDetailedPatchySwapInteraction_H_

#include "Interactions/BaseInteraction.h"
#include "../Particles/TorsionalPatchyParticle.h"

/**
 * @brief Manages the interaction between generic patchy particles.
 *
 * This interaction is selected with
 * interaction_type = TorsionalDetailedPatchySwapInteraction
 */

#define PLPATCH_VM1 0
#define PLPATCH_VM3 1

//#define PATCHY_POWER 200
//#define PATCHY_CUTOFF 0.18f

#define PLEXCL_S   1.0f
#define PLEXCL_R   0.9053f
#define PLEXCL_B   677.505671539f
#define PLEXCL_RC  0.99888f
#define PLEXCL_EPS 2.0f

//constants based on villars paper
#define PLEXCL_VS   1.0f
#define PLEXCL_VR   1.15f
#define PLEXCL_VB   (-0.3797564833966851)
#define PLEXCL_VRC  (2.757802724660435)
#define PLEXCL_VEPS 2.0f

//48-96 lj potnetial
#define PLEXCL_V48S   1.0
#define PLEXCL_V48R   1.15f
#define PLEXCL_V48B   (-2.11835)
#define PLEXCL_V48RC  1.19798
#define PLEXCL_V48EPS 2.0f

//there are two types of modulation
#define PLEXCL_NARROW_N 2




class TorsionalDetailedPatchySwapInteraction: public BaseInteraction {
protected:
	/// Number of particles of each species
	std::vector<int> _N_per_species;
	int _N = 0;
	int _N_patch_types = 0;

	/// Number of patches per particle
	std::vector<uint> _N_patches;
	/// Patch type for each particle species
	std::vector<std::vector<int> > _patch_types;
	/// Base position of the patches for each particle species
	std::vector<std::vector<LR_vector> > _base_patches;

    /// orientations of the patches for each particle species
	std::vector<std::vector<LR_vector> > _a1_patches;
    std::vector<std::vector<LR_vector> > _a2_patches;
    

	std::string _interaction_matrix_file;

	/// Repulsive interaction energy at the cut-off
	number _rep_rcut = 0.;
	number _sqr_rep_rcut = -1.;

	/// Patchy-related quantities
	std::vector<number> _patchy_eps;
	number _patch_rcut = -1.;
	number _sqr_patch_rcut = -1.;
	number _sigma_ss = 0.4;
	number _rcut_ss = -1.;
	number _lambda = 1.0;
	number _A_part = 0.;
	number _B_part = 0.;

	/// KF-related quantities
	/// Width of the patches
	bool _is_KF = false;


    bool _use_torsion = true;
    int _narrow_type = 0; 

	/// Exponent for the Gaussian-like potential well used for the patches
	int _patch_power = 30;
	number _patch_delta;
	/// Angular width of the patches
	number _patch_cosmax;
	/// _patch_alpha^10
	number _patch_pow_delta;
	/// _patch_cosmax^30
	number _patch_pow_cosmax;
	/// Angular cut-off for the patchy attraction
	number _patch_angular_cutoff;

	/// Optional spherical attraction
	number _spherical_attraction_strength = 0.;
	number _spherical_rcut = 2.5;
	number _sqr_spherical_rcut = 6.25;
	number _spherical_E_cut = 0.;

    //model constants
    //model constants
	number PLPATCHY_THETA_TC[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_T0[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_TS[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_A[PLEXCL_NARROW_N];
	number PLPATCHY_THETA_B[PLEXCL_NARROW_N];

    //custom functions
    number _V_mod(int type, number t);
    number _V_modD(int type, number t);
    number _V_modDsin(int type, number t);


	void _parse_interaction_matrix();
	//std::vector<LR_vector> _parse_base_patches(std::string filename, int N_patches);
    void _parse_base_patches(std::string filename, int N_patches,std::vector<LR_vector> &base_patches,  std::vector<LR_vector>  &a1_patches,  std::vector<LR_vector>  &a2_patches);
	number _patchy_two_body_point(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _patchy_two_body_KF(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	number _spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	number _three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces);

	inline std::vector<PatchyBond> &_particle_bonds(BaseParticle *p) {
		return static_cast<TorsionalPatchyParticle *>(p)->bonds;
	}

public:
	enum {
		PATCHY = 0,
		SPHERICAL = 1
	};

	bool no_three_body = false;

	TorsionalDetailedPatchySwapInteraction();
	virtual ~TorsionalDetailedPatchySwapInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	void begin_energy_computation() override;

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

extern "C" TorsionalDetailedPatchySwapInteraction *make_TorsionalDetailedPatchySwapInteraction();

#endif /* TorsionalDetailedPatchySwapInteraction_H_ */
