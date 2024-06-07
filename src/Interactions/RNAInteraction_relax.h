/*
 * RNAInteraction_relax.h
 *
 *  Created on: Nov 28, 2014
 *      Author: Ben Snodin
 */

#ifndef RNA_INTERACTION_RELAX_H
#define RNA_INTERACTION_RELAX_H

#include "BaseInteraction.h"
#include "RNAInteraction.h"

/**
 * @brief Modified version of RNAInteraction which modifies the bonded backbone-backbone potential so it does not diverge
 *
 * Replaces the bonded backbone-backbone FENE potential with a harmonic potential. This is to allow very stressed initial
 * configurations, which might otherwise cause the simulation to fail, to relax to a sensible structure
 *
 * This interaction takes 3 compulsory arguments:
 *
 * This interaction is selected with
 * interaction_type = RNA_relax
 *
 @verbatim
 relax_type = <string> (Possible values: constant_force, harmonic_force; Relaxation algorithm used)
 relax_strength = <float> (Force constant for the replacement of the FENE potential)
 soft_exc_vol = <bool> (True if we are using pass-through soft excluded volume)
 @endverbatim
 */


//constants for the soft repulsion nonbonded excluded volume (see excvol.nb)
 
 #define RNA_SOFT_EXCL_A1  2.1
 #define RNA_SOFT_EXCL_RS1 0.669
 #define RNA_SOFT_EXCL_B1   32.829
 #define RNA_SOFT_EXCL_RC1  0.711794

 #define RNA_SOFT_EXCL_A2 10.29
 #define RNA_SOFT_EXCL_RS2  0.29
 #define RNA_SOFT_EXCL_B2  66.1525
 #define RNA_SOFT_EXCL_RC2 0.335109


 #define RNA_SOFT_EXCL_A3 4.26
 #define RNA_SOFT_EXCL_RS3 0.45
 #define RNA_SOFT_EXCL_B3 26.7557
 #define RNA_SOFT_EXCL_RC3 0.521648

 #define RNA_SOFT_EXCL_A4 4.26
 #define RNA_SOFT_EXCL_RS4 0.45
 #define RNA_SOFT_EXCL_B4 26.7557
 #define RNA_SOFT_EXCL_RC4 0.521648





class RNAInteraction_relax: public RNAInteraction {
protected:
    inline number _soft_repulsion(const LR_vector &r, LR_vector &force, number a, number rstar, number b, number rc, number K, bool update_forces) ;

	inline virtual number _backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);
	inline virtual number _nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces);

	int _backbone_type;
	float _backbone_k;
	bool _soft_exc_vol;
	float _soft_exc_vol_K;

	int _constant_force;
	int _harmonic_force;

public:
	RNAInteraction_relax();
	virtual ~RNAInteraction_relax();



	void check_input_sanity(std::vector<BaseParticle *> &particles);
	void get_settings(input_file &inp);
};

#endif /*RNA_INTERACTION_RELAX_H*/
