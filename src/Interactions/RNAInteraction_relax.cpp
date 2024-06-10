/*
 * RNAInteraction_relax.cpp
 *
 *  Created on: Nov 28, 2014
 *      Author: Ben Snodin
 */

#include <fstream>

#include "RNAInteraction_relax.h"
#include "../Particles/RNANucleotide.h"

RNAInteraction_relax::RNAInteraction_relax() :
				RNAInteraction() {
	OX_LOG(Logger::LOG_INFO, "Using unphysical backbone (RNA_relax interaction)");
	_constant_force = 0;
	_harmonic_force = 1;
}

RNAInteraction_relax::~RNAInteraction_relax() {
}

number RNAInteraction_relax::_soft_repulsion(const LR_vector &r, LR_vector &force, number a, number rstar, number b, number rc, number K, bool update_forces) {
	number rnorm = r.norm();
	
	number val = (number) 0.f; //energy
	
	if(rnorm < rc*rc) {
		number t = sqrt(rnorm); //rback.module();
		if(t > rstar) {
			// smoothing
			val = K * b * SQR(rc - t);
			if(update_forces)
				force = r * ((K * 2 * b) * (rc - t) / t);
		}
		else {
			val =  K*( (number) 1.f - a * SQR(t));
			if(update_forces)
				force = r * ( 2 * K * a );

		}
	}
	else if (update_forces) {
		force.x = force.y = force.z = (number) 0.f;
	}


return val;

}



number RNAInteraction_relax::_nonbonded_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_are_bonded(p, q)) {
		return (number) 0.f;
	}

	if(!_soft_exc_vol)
	{
		return RNAInteraction::_nonbonded_excluded_volume(p,q,compute_r,update_forces);
	}

	//here we calculate exc volume repulsion that is a quadratic

	LR_vector force(0, 0, 0);
	LR_vector torquep(0, 0, 0);
	LR_vector torqueq(0, 0, 0);

	// BASE-BASE
	LR_vector rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BASE];

	number energy = _soft_repulsion(rcenter, force, RNA_SOFT_EXCL_A2, RNA_SOFT_EXCL_RS2, RNA_SOFT_EXCL_B2, RNA_SOFT_EXCL_RC2, this->_soft_exc_vol_K,update_forces);
	
	if(update_forces) {
		torquep = -p->int_centers[RNANucleotide::BASE].cross(force);
		torqueq = q->int_centers[RNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BASE vs. Q-BACK
	rcenter = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BASE];
	energy += _soft_repulsion(rcenter, force, RNA_SOFT_EXCL_A3, RNA_SOFT_EXCL_RS3, RNA_SOFT_EXCL_B3, RNA_SOFT_EXCL_RC3, this->_soft_exc_vol_K,update_forces);// _soft_repulsion(rcenter, force, model->RNA_EXCL_S3, model->RNA_EXCL_R3, model->RNA_EXCL_B3, model->RNA_EXCL_RC3, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide::BASE].cross(force);
		torqueq += q->int_centers[RNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;
	}

	// P-BACK vs. Q-BASE
	rcenter = _computed_r + q->int_centers[RNANucleotide::BASE] - p->int_centers[RNANucleotide::BACK];
	energy += _soft_repulsion(rcenter, force, RNA_SOFT_EXCL_A4, RNA_SOFT_EXCL_RS4, RNA_SOFT_EXCL_B4, RNA_SOFT_EXCL_RC4, this->_soft_exc_vol_K,update_forces); //_soft_repulsion(rcenter, force, model->RNA_EXCL_S4, model->RNA_EXCL_R4, model->RNA_EXCL_B4, model->RNA_EXCL_RC4, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[RNANucleotide::BASE].cross(force);

		p->force -= force;
		q->force += force;
	}

	// BACK-BACK
	rcenter = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
	energy += _soft_repulsion(rcenter, force, RNA_SOFT_EXCL_A1, RNA_SOFT_EXCL_RS1, RNA_SOFT_EXCL_B1, RNA_SOFT_EXCL_RC1, this->_soft_exc_vol_K,update_forces); //_soft_repulsion(rcenter, force, model->RNA_EXCL_S1, model->RNA_EXCL_R1, model->RNA_EXCL_B1, model->RNA_EXCL_RC1, update_forces);

	if(update_forces) {
		torquep += -p->int_centers[RNANucleotide::BACK].cross(force);
		torqueq += q->int_centers[RNANucleotide::BACK].cross(force);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}





number RNAInteraction_relax::_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!this->_harmonic_spring)
	{
		return RNAInteraction::_backbone(p,q,compute_r,update_forces);
	}

	if(!_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	LR_vector rback = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
	number rbackmod = rback.module();
	number rbackr0 = rbackmod - model->RNA_FENE_R0;

	number energy = 0;
	if(_backbone_type == _constant_force) {
		energy = _backbone_k * fabs(rbackr0);
	}
	else if(_backbone_type == _harmonic_force) {
		energy = 0.5 * _backbone_k * rbackr0 * rbackr0;
	}

	if(update_forces) {
		LR_vector force;
		// this does not conserve energy very well -- but there is a discontinuity in the force, so...
		if(_backbone_type == _constant_force) {
			if(rbackr0 > 0) force = rback * (-_backbone_k / rbackmod);
			else force = rback * (_backbone_k / rbackmod);
		}
		// conserves energy about as well as the normal RNAInteraction
		else if(_backbone_type == _harmonic_force) force = rback * (-_backbone_k * rbackr0 / rbackmod);

		//printf("Force is %f %f %f, rback is %f %f %f, rbackr0 is %f, ene is %f , K is %f \n",force.x,force.y, force.z,rback.x,rback.y, rback.z, rbackr0,energy,_backbone_k);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		p->torque -= p->orientationT * p->int_centers[RNANucleotide::BACK].cross(force);
		q->torque += q->orientationT * q->int_centers[RNANucleotide::BACK].cross(force);
	}

	return energy;
}

void RNAInteraction_relax::check_input_sanity(std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		}
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
		}
	}
}

void RNAInteraction_relax::get_settings(input_file &inp) {
	RNAInteraction::get_settings(inp);

	char tmps[256];
	float ftmp;
	_harmonic_spring = true;

	if(getInputBool(&inp, "harmonic_backbone", &_harmonic_spring, 0) == KEY_FOUND)
	{
	
	 getInputString(&inp, "relax_type", tmps, 1);
	 if(strcmp(tmps, "constant_force") == 0) _backbone_type = _constant_force;
	 else if(strcmp(tmps, "harmonic_force") == 0) _backbone_type = _harmonic_force;
	 else throw oxDNAException("Error while parsing input file: relax_type '%s' not implemented; use constant_force or harmonic_force", tmps);
	
	if(getInputFloat(&inp, "relax_strength", &ftmp, 0) == KEY_FOUND) {
		_backbone_k = (number) ftmp;
		OX_LOG(Logger::LOG_INFO, "Using spring constant = %f for the RNA_relax interaction", _backbone_k);
	}
	else {
		if (_backbone_type == _harmonic_force) {
			_backbone_k = (number) 32.;
		}
		else {
			_backbone_k = (number) 1.;
		}
		OX_LOG(Logger::LOG_INFO, "Using default strength constant = %f for the RNA_relax interaction", _backbone_k);
	}

	}
	bool soft_exc_vol = false;
	if(getInputBool(&inp, "soft_exc_vol", &soft_exc_vol, 0) == KEY_FOUND)  
	{ 
		_soft_exc_vol = soft_exc_vol;
	} 
	else 
	{
		_soft_exc_vol = false;
	}
	if (soft_exc_vol)
	{
		if(getInputFloat(&inp, "exc_vol_strength", &ftmp, 0) == KEY_FOUND) {
			_soft_exc_vol_K = ftmp;
		}
		else 
		{
			_soft_exc_vol_K = 10.;
		}
		OX_LOG(Logger::LOG_INFO, "Using strength constant = %f for the RNA_relax soft non-bonded exclusion volume interaction",_soft_exc_vol_K);
	}


/*
 for(double x = 0.; fabs(x) < 1.0; x -= 0.005)
 {
  LR_vector force(0,0,0);
  LR_vector rcenter(0,0,x);
  number energy = _soft_repulsion(rcenter, force, RNA_SOFT_EXCL_A1, RNA_SOFT_EXCL_RS1, RNA_SOFT_EXCL_B1, RNA_SOFT_EXCL_RC1, 10.,true);
  printf("@@@ %f %f %f\n",x,energy,force.z);
 }
 */
}
