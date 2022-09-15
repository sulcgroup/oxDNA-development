/*
 * TorsionalPatchyParticle.h
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#ifndef TTorsionalPatchyParticle_H_
#define TTorsionalPatchyParticle_H_

#include "../../../../src/Particles/BaseParticle.h"

struct PatchyBond {
	BaseParticle *other;
	number r_p;
	int p_patch, q_patch;
	number energy;
	LR_vector force;
	LR_vector p_torque, q_torque;

	PatchyBond(BaseParticle *o, number my_r_p, int pp, int qp, number e) :
					other(o),
					r_p(my_r_p),
					p_patch(pp),
					q_patch(qp),
					energy(e) {
	}
};

/**
 * @brief Incapsulates a patchy particle.
 */

class TorsionalPatchyParticle : public BaseParticle {
public:
	number _sigma;
	std::vector<LR_vector> _base_patches;
	std::vector<LR_vector> _a1_patches;  //a1 vector that is to be compared against the vector connecting the patches r_pp, needs to be parallel
    std::vector<LR_vector> _a2_patches; // vector which needs to be parallel with a2 on the complementary patch
	


public:
	//TorsionalPatchyParticle(std::vector<LR_vector> base_patches, int nt, number sigma);
	TorsionalPatchyParticle(std::vector<LR_vector> base_patches, std::vector<LR_vector> a1_patches, std::vector<LR_vector> a2_patches, int nt, number sigma);
	TorsionalPatchyParticle(int N_patches, int nt, number sigma);
	virtual ~TorsionalPatchyParticle();

	void set_positions();


	virtual bool is_rigid_body() {
		return true;
	}

	std::vector<LR_vector> base_patches() {
		return _base_patches;
	}

	std::vector<PatchyBond> bonds;
};

#endif /* TorsionalPatchyParticle_H_ */
