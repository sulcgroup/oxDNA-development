/*
 * TorsionalDetailedPatchySwapInteraction.cpp
 *
 *  Created on: 15/may/2021
 *      Author: lorenzo
 */

#include "TorsionalDetailedPatchySwapInteraction.h"
//#include "Particles/CustomParticle.h"
#include <Utilities/Utils.h>
#include "../Particles/TorsionalPatchyParticle.h"
//#include "Particles/CustomParticle.h"

#include <string>

using namespace std;


number TorsionalDetailedPatchySwapInteraction::_V_mod(int type, number t)
{
	number val = (number) 0;
	t -= PLPATCHY_THETA_T0[type];
	if (t < 0)
		t *= -1;

	if (t < PLPATCHY_THETA_TC[type]) {
		if (t > PLPATCHY_THETA_TS[type]) {
			// smoothing
			val = PLPATCHY_THETA_B[type] * SQR(PLPATCHY_THETA_TC[type] - t);
		} else
			val = (number) 1.f - PLPATCHY_THETA_A[type] * SQR(t);
	}

	return val;
}



number TorsionalDetailedPatchySwapInteraction::_V_modD(int type, number t)
{
	number val = (number) 0;
	number m = (number) 1;
	t -= PLPATCHY_THETA_T0[type];
	// this function is a parabola centered in t0. If t < 0 then the value of the function
	// is the same but the value of its derivative has the opposite sign, so m = -1
	if(t < 0) {
		t *= -1;
		m = (number) -1;
	}

	if(t < PLPATCHY_THETA_TC[type]) {
		if(t > PLPATCHY_THETA_TS[type]) {
			// smoothing
			val = m * 2 * PLPATCHY_THETA_B[type] * (t - PLPATCHY_THETA_TC[type]);
		}
		else val = -m * 2 * PLPATCHY_THETA_A[type] * t;
	}

	return val;
}



number TorsionalDetailedPatchySwapInteraction::_V_modDsin(int type, number t)
{
	    number val = (number) 0;
		number m = (number) 1;
		number tt0 = t - PLPATCHY_THETA_T0[type];
		// this function is a parabola centered in t0. If t < 0 then the value of the function
		// is the same but the value of its derivative has the opposite sign, so m = -1
		if(tt0 < 0) {
			tt0 *= -1;
			m = (number) -1;
		}

		if(tt0 < PLPATCHY_THETA_TC[type]) {
		    	number sint = sin(t);
			if(tt0 > PLPATCHY_THETA_TS[type]) {
				// smoothing
				val = m * 2 * PLPATCHY_THETA_B[type] * (tt0 - PLPATCHY_THETA_TC[type]) / sint;
			}
			else {
			    if(SQR(sint) > 1e-8) val = -m * 2 * PLPATCHY_THETA_A[type] * tt0 / sint;
			    else val = -m * 2 * PLPATCHY_THETA_A[type];
			}
		}

		return val;
}


TorsionalDetailedPatchySwapInteraction::TorsionalDetailedPatchySwapInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(SPHERICAL, _spherical_patchy_two_body);

}

TorsionalDetailedPatchySwapInteraction::~TorsionalDetailedPatchySwapInteraction() {

}

void TorsionalDetailedPatchySwapInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputNumber(&inp, "DPS_lambda", &_lambda, 0);
	getInputString(&inp, "DPS_interaction_matrix_file", _interaction_matrix_file, 1);

    getInputBool(&inp, "use_torsion", &_use_torsion, 0);
    getInputNumber(&inp, "narrow_type", &_narrow_type, 0);

	getInputBool(&inp, "DPS_is_KF", &_is_KF, 0);
	if(_is_KF) {
		getInputInt(&inp, "DPS_patch_power", &_patch_power, 0);
		getInputNumber(&inp, "DPS_KF_delta", &_patch_delta, 1);
		getInputNumber(&inp, "DPS_KF_cosmax", &_patch_cosmax, 1);

        if ( _use_torsion)
        {
            throw oxDNAException("use_torsion option is not supported with KF potential");
        }
	}
	else {
		getInputNumber(&inp, "DPS_sigma_ss", &_sigma_ss, 0);
	}

	getInputNumber(&inp, "DPS_spherical_attraction_strength", &_spherical_attraction_strength, 0.);
	if(_spherical_attraction_strength > 0.) {
		getInputNumber(&inp, "DPS_spherical_rcut", &_spherical_rcut, 1.);
	}
}

void TorsionalDetailedPatchySwapInteraction::init() {
	_rep_rcut = pow(2., 1. / 6.);
	_sqr_rep_rcut = SQR(_rep_rcut);

	if(_is_KF) {
		ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body_KF);

		_patch_rcut = 1.5 * _patch_delta;

		// the patch-patch radial attraction is centred around 0 (considering that we use the distance between the two surfaces as a variable)
		_sigma_ss = 0.;
		_patch_pow_delta = pow(_patch_delta, (number) 10.);
		_patch_pow_cosmax = pow(1. - _patch_cosmax, (number) _patch_power);
		// this makes sure that at the cutoff the angular modulation is 10^-2
		_patch_angular_cutoff = (1. - _patch_cosmax) * std::pow(4 * std::log(10), 1. / _patch_power);

		OX_LOG(Logger::LOG_INFO, "FS-KF parameters: lambda = %lf, patch_delta = %lf, patch_power = %d, patch_cosmax = %lf, patch_angular_cutoff = %lf, ", _lambda, _patch_delta, _patch_power, _patch_cosmax, _patch_angular_cutoff);
	}
	else {
		ADD_INTERACTION_TO_MAP(PATCHY, _patchy_two_body_point);

		_rcut_ss = 1.5 * _sigma_ss;
		_patch_rcut = _rcut_ss;

		number B_ss = 1. / (1. + 4. * SQR(1. - _rcut_ss / _sigma_ss));
		_A_part = -1. / (B_ss - 1.) / exp(1. / (1. - _rcut_ss / _sigma_ss));
		_B_part = B_ss * pow(_sigma_ss, 4.);

		OX_LOG(Logger::LOG_INFO, "FS parameters: lambda = %lf, A_part = %lf, B_part = %lf, narrow_type = %d, use_torsion=%d", _lambda, _A_part, _B_part,_narrow_type,_use_torsion);
	}
	if(_narrow_type == 0)
	{
 	  //default option, very wide!
	  PLPATCHY_THETA_T0[0] = PLPATCHY_THETA_T0[1] = 0.;
      PLPATCHY_THETA_TS[0] = PLPATCHY_THETA_TS[1] = 0.7;
      PLPATCHY_THETA_TC[0] = PLPATCHY_THETA_TC[1] = 3.10559;
      PLPATCHY_THETA_A[0] = PLPATCHY_THETA_A[1] = 0.46;
      PLPATCHY_THETA_B[0] = PLPATCHY_THETA_B[1] = 0.133855;
	}
	else if(_narrow_type == 1)
	{
		  //narrower!
		  PLPATCHY_THETA_T0[0] =  PLPATCHY_THETA_T0[1] = 0.;
		  PLPATCHY_THETA_TS[0] =  PLPATCHY_THETA_TS[1] = 0.2555;
		  PLPATCHY_THETA_TC[0] =  PLPATCHY_THETA_TC[1] = 1.304631441617743;
		  PLPATCHY_THETA_A[0]  =  PLPATCHY_THETA_A[1] = 3.;
		  PLPATCHY_THETA_B[0]  =  PLPATCHY_THETA_B[1] = 0.7306043547966398;
	}
	else if(_narrow_type == 2)
	{
	  //narrower!
	  PLPATCHY_THETA_T0[0] =  PLPATCHY_THETA_T0[1] = 0.;
	  PLPATCHY_THETA_TS[0] =  PLPATCHY_THETA_TS[1] = 0.2555;
	  PLPATCHY_THETA_TC[0] =  PLPATCHY_THETA_TC[1] = 0.782779;
	  PLPATCHY_THETA_A[0]  =  PLPATCHY_THETA_A[1] = 5.;
	  PLPATCHY_THETA_B[0]  =  PLPATCHY_THETA_B[1] = 2.42282;
	}
	else if (_narrow_type == 3)
	{
 	   //narrower:
       PLPATCHY_THETA_T0[0] = PLPATCHY_THETA_T0[1] = 0.;
       PLPATCHY_THETA_TS[0] = PLPATCHY_THETA_TS[1] = 0.17555;
       PLPATCHY_THETA_TC[0] = PLPATCHY_THETA_TC[1] = 0.4381832920710734;
       PLPATCHY_THETA_A[0] = PLPATCHY_THETA_A[1] = 13.;
       PLPATCHY_THETA_B[0] = PLPATCHY_THETA_B[1] = 8.68949241736805;
	}
	else if(_narrow_type == 4)
	{
      //narrowest:
      PLPATCHY_THETA_T0[0] = PLPATCHY_THETA_T0[1] = 0.f;
      PLPATCHY_THETA_TS[0] = PLPATCHY_THETA_TS[1] = 0.17555;
      PLPATCHY_THETA_TC[0] = PLPATCHY_THETA_TC[1] = 0.322741;
      PLPATCHY_THETA_A[0] = PLPATCHY_THETA_A[1] = 17.65;
      PLPATCHY_THETA_B[0] = PLPATCHY_THETA_B[1] = 21.0506;
	}
	else{
		 throw oxDNAException ("Invalid narrow_type option, has to be between 0 and 4");
	}


	_sqr_patch_rcut = SQR(_patch_rcut);
	_rcut = 1. + _patch_rcut;

	if(_spherical_attraction_strength > 0.) {
		_sqr_spherical_rcut = SQR(_spherical_rcut);
		_spherical_E_cut = 4. * _spherical_attraction_strength * (1. / pow(_sqr_spherical_rcut, 6) - 1. / pow(_sqr_spherical_rcut, 3));

		if(_spherical_rcut > _rcut) {
			_rcut = _spherical_rcut;
		}
	}

	_sqr_rcut = SQR(_rcut);
}

number TorsionalDetailedPatchySwapInteraction::_spherical_patchy_two_body(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number energy = (number) 0.f;

	// centre-centre
	if(sqr_r < _sqr_rep_rcut) {
		number ir2 = 1. / sqr_r;
		number lj_part = ir2 * ir2 * ir2;
		energy = 4 * (SQR(lj_part) - lj_part) + 1.0 - _spherical_attraction_strength - _spherical_E_cut;
		if(update_forces) {
			LR_vector force = _computed_r * (-24. * (lj_part - 2 * SQR(lj_part)) / sqr_r);
			p->force -= force;
			q->force += force;
		}
	}
	else {
		if(sqr_r < _sqr_spherical_rcut && _spherical_attraction_strength > 0.) {
			number ir2 = 1. / sqr_r;
			number lj_part = ir2 * ir2 * ir2;
			energy = 4 * _spherical_attraction_strength * (SQR(lj_part) - lj_part) - _spherical_E_cut;
			if(update_forces) {
				LR_vector force = _computed_r * (-24. * _spherical_attraction_strength * (lj_part - 2 * SQR(lj_part)) / sqr_r);
				p->force -= force;
				q->force += force;
			}
		}
	}

	return energy;
}

number TorsionalDetailedPatchySwapInteraction::_patchy_two_body_point(BaseParticle *pp, BaseParticle *qq, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	TorsionalPatchyParticle *p = dynamic_cast<TorsionalPatchyParticle *>(pp);
	TorsionalPatchyParticle *q = dynamic_cast<TorsionalPatchyParticle *>(qq);

    number  cosa1 = 1;
    number  cosb1 = 1;
    number  cosa2b2 = 1;
    number  ta1 = 0;
    number  tb1 = 0;
    number  ta2b2 = 0;
    number  fa1 =1;
    number  fb1 =1;
    number  fa2b2 = 1;
	number energy = (number) 0.f;

    number r_dist = sqrt(sqr_r);
    LR_vector r_dist_dir = _computed_r / r_dist;

	for(uint p_patch = 0; p_patch < p->N_int_centers(); p_patch++) {
		LR_vector p_patch_pos = p->int_centers[p_patch];
		for(uint q_patch = 0; q_patch < q->N_int_centers(); q_patch++) {
			LR_vector q_patch_pos = q->int_centers[q_patch];

			LR_vector patch_dist = _computed_r + q_patch_pos - p_patch_pos;
          
			number dist = patch_dist.norm();
			if(dist < _sqr_patch_rcut) {
				uint p_patch_type = _patch_types[p->type][p_patch];
				uint q_patch_type = _patch_types[q->type][q_patch];
				number epsilon = _patchy_eps[p_patch_type + _N_patch_types * q_patch_type];

				if(epsilon != 0.) {
					number r_p = sqrt(dist);
					number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
                    number f1 = epsilon * _A_part * exp_part * (_B_part / SQR(dist) - 1.);
					number tmp_energy = f1; //epsilon * _A_part * exp_part * (_B_part / SQR(dist) - 1.);

                    number angular_part = 1.;
                   
                    if(_use_torsion)
                    {
                        cosa1 = p->_a1_patches[p_patch] * r_dist_dir;
                        cosb1 = -q->_a1_patches[q_patch] * r_dist_dir;
                        cosa2b2 = p->_a2_patches[p_patch] * q->_a2_patches[q_patch];

                        ta1 = LRACOS(cosa1);
                        tb1 = LRACOS(cosb1);
                        ta2b2 = LRACOS(cosa2b2);
                        fa1 =  _V_mod(PLPATCH_VM1,ta1);
                        fb1 =  _V_mod(PLPATCH_VM1,tb1) ;

                        fa2b2 =   _V_mod(PLPATCH_VM3,ta2b2);
                        angular_part = fa1 * fb1 * fa2b2;
                        tmp_energy *= angular_part;
                    }
                

					energy += tmp_energy;

					number tb_energy = (r_p < _sigma_ss) ? epsilon : -tmp_energy;

					PatchyBond p_bond(q, r_p, p_patch, q_patch, tb_energy);
					PatchyBond q_bond(p, r_p, q_patch, p_patch, tb_energy);

					if(update_forces) {
						number force_mod = epsilon * _A_part * exp_part * (4. * _B_part / (SQR(dist) * r_p)) + _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
						LR_vector tmp_force = patch_dist * (force_mod / r_p);

                        LR_vector p_torque;
                        LR_vector q_torque;

                        if(_use_torsion)
                        {
                            //number f1D =  (5 * exp_part * r8b10);
                            //LR_vector<number> tmp_force = patch_dist * (f1D * angular_part);

                            tmp_force *= angular_part;
                            //printf("CRITICAL 1 Adding %f %f %f \n",tmp_force.x,tmp_force.y,tmp_force.z);

                            number fa1Dsin =  _V_modDsin(PLPATCH_VM1,ta1);
                            number fb1Dsin =  _V_modDsin(PLPATCH_VM1,tb1);

                        
                            //printf("CRITICAL 2 Adding %f %f %f \n",tmp_force.x,tmp_force.y,tmp_force.z);
                            //torque VM3
                            number fa2b2Dsin =  _V_modDsin(PLPATCH_VM3,ta2b2);
                            LR_vector dir = -p->_a2_patches[p_patch].cross(q->_a2_patches[q_patch]) *  (f1 * fa1 * fb1 * fa2b2Dsin );
                            LR_vector torqueq = dir;
                            LR_vector torquep = dir;


                            //torque VM1
                            dir = r_dist_dir.cross(p->_a1_patches[p_patch]);
                            torquep += dir * (f1 * fa1Dsin * fb1 );

                            //torque VM2
                            dir = r_dist_dir.cross(q->_a1_patches[q_patch]);
                            torqueq += dir * (f1 * fa1 * fb1Dsin );


                            torquep += p_patch_pos.cross(tmp_force);
                            torqueq += q_patch_pos.cross(tmp_force);

                            tmp_force += (p->_a1_patches[p_patch] -  r_dist_dir * cosa1) * (f1 * fa1Dsin *  fb1* fa2b2 / r_dist);
                            tmp_force += -(q->_a1_patches[q_patch] +  r_dist_dir * cosb1) * (f1 * fa1 *  fb1Dsin * fa2b2 / r_dist);

                            p_torque = p->orientationT * torquep;
                            q_torque = q->orientationT * torqueq;


                            p->torque -=  p_torque; // p->orientationT * torquep;
					        q->torque +=  q_torque; // q->orientationT * torqueq;

                            p->force -= tmp_force;
						    q->force += tmp_force;

                          
                        }

                        else 
						{ 
                          p_torque = p->orientationT * p_patch_pos.cross(tmp_force);
						  q_torque = q->orientationT * q_patch_pos.cross(tmp_force);
                          p->force -= tmp_force;
						  q->force += tmp_force;
						  p->torque -= p_torque;
						  q->torque += q_torque;
                        }

						
						if(r_p > _sigma_ss) {
							p_bond.force = tmp_force;
							p_bond.p_torque = p_torque;
							p_bond.q_torque = q_torque;

							q_bond.force = -tmp_force;
							q_bond.p_torque = -q_torque;
							q_bond.q_torque = -p_torque;
						}
					}

					_particle_bonds(p).emplace_back(p_bond);
					_particle_bonds(q).emplace_back(q_bond);

					if(!no_three_body) {
						energy += _three_body(p, p_bond, update_forces);
						energy += _three_body(q, q_bond, update_forces);

					}
				}
			}
		}
	}

	return energy;
}

// here we compute x^n as (x*x)^((n-1)/2) * x since we now that n is always an odd number
inline double _lr_pow(double x, size_t n){
    double res = x;
    x *= x;

    n = (n - 1) / 2;
    while(n-- > 0){
        res *= x;
    }

    return res;
}

number TorsionalDetailedPatchySwapInteraction::_patchy_two_body_KF(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) {
		return (number) 0.f;
	}

	number rmod = sqrt(sqr_r);
	LR_vector r_versor = _computed_r / rmod;

	number dist_surf = rmod - 1.;
	number dist_surf_sqr = SQR(dist_surf);
	number r8b10 = SQR(SQR(dist_surf_sqr)) / _patch_pow_delta;
	number exp_part = -1.001 * exp(-(number) 0.5 * r8b10 * dist_surf_sqr);

	number energy = (number) 0.f;
	for(uint p_patch = 0; p_patch < p->N_int_centers(); p_patch++) {
		LR_vector p_patch_pos = p->int_centers[p_patch] * 2;

		number cospr = p_patch_pos * r_versor;
		number cospr_minus_one = cospr - 1.;
		if(cospr_minus_one < _patch_angular_cutoff) {
			number cospr_base = _lr_pow(cospr_minus_one, _patch_power - 1);
			// we do this so that later we don't have to divide this number by (cospr - 1), which could be 0
			number cospr_part = cospr_base * cospr_minus_one;
			number p_mod = exp(-cospr_part / (2. * _patch_pow_cosmax));

			for(uint q_patch = 0; q_patch < q->N_int_centers(); q_patch++) {
				LR_vector q_patch_pos = q->int_centers[q_patch] * 2;

				number cosqr = -(q_patch_pos * r_versor);
				number cosqr_minus_one = cosqr - 1.;
				if(cosqr_minus_one < _patch_angular_cutoff) {
					uint p_patch_type = _patch_types[p->type][p_patch];
					uint q_patch_type = _patch_types[q->type][q_patch];
					number epsilon = _patchy_eps[p_patch_type + _N_patch_types * q_patch_type];

					if(epsilon != 0.) {
						number cosqr_base = _lr_pow(cosqr_minus_one, _patch_power - 1);
						number cosqr_part = cosqr_base * cosqr_minus_one;
						number q_mod = exp(-cosqr_part / (2. * _patch_pow_cosmax));

						number tmp_energy = exp_part * p_mod * q_mod;

						if(tmp_energy < 0.) {
							energy += tmp_energy;

							// when we do the swapping the radial part is the one that gets to one beyond the minimum, while the angular part doesn't change
							number tb_energy = (dist_surf < _sigma_ss) ? epsilon * p_mod * q_mod : -tmp_energy;
							PatchyBond p_bond(q, dist_surf, p_patch, q_patch, tb_energy);
							PatchyBond q_bond(p, dist_surf, q_patch, p_patch, tb_energy);

							if(update_forces) {
								// radial part
								LR_vector radial_force = r_versor * (p_mod * q_mod * 5. * (rmod - 1.) * exp_part * r8b10);

								// angular p part
								number der_p = exp_part * q_mod * (0.5 * _patch_power * p_mod * cospr_base / _patch_pow_cosmax);
								LR_vector p_ortho = p_patch_pos - cospr * r_versor;
								LR_vector angular_force = p_ortho * (der_p / rmod);

								// angular q part
								number der_q = exp_part * p_mod * (-0.5 * _patch_power * q_mod * cosqr_base / _patch_pow_cosmax);
								LR_vector q_ortho = q_patch_pos + cosqr * r_versor;
								angular_force += q_ortho * (der_q / rmod);

								LR_vector tot_force = radial_force + angular_force;

								LR_vector p_torque = p->orientationT * (r_versor.cross(p_patch_pos) * der_p);
								LR_vector q_torque = q->orientationT * (q_patch_pos.cross(r_versor) * der_q);

								p->force -= tot_force;
								q->force += tot_force;

								p->torque -= p_torque;
								q->torque += q_torque;

								p_bond.force = (dist_surf < _sigma_ss) ? angular_force : tot_force;
								p_bond.p_torque = p_torque;
								p_bond.q_torque = q_torque;

								q_bond.force = (dist_surf < _sigma_ss) ? -angular_force : -tot_force;
								q_bond.p_torque = -q_torque;
								q_bond.q_torque = -p_torque;
							}

							_particle_bonds(p).emplace_back(p_bond);
							_particle_bonds(q).emplace_back(q_bond);

							if(!no_three_body) {
								energy += _three_body(p, p_bond, update_forces);
								energy += _three_body(q, q_bond, update_forces);
							}
						}
					}
				}
			}
		}
	}

	return energy;
}

number TorsionalDetailedPatchySwapInteraction::_three_body(BaseParticle *p, PatchyBond &new_bond, bool update_forces) {
	number energy = 0.;

	number curr_energy = new_bond.energy;
	const auto &p_bonds = _particle_bonds(p);
	for(auto &other_bond : p_bonds) {
		// three-body interactions happen only when the same patch is involved in more than a bond
		if(other_bond.other != new_bond.other && other_bond.p_patch == new_bond.p_patch) {
			number other_energy = other_bond.energy;

			energy += _lambda * curr_energy * other_energy;

			if(update_forces) {
				{
					BaseParticle *other = new_bond.other;

					number factor = -_lambda * other_energy;
					LR_vector tmp_force = factor * new_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					p->torque -= factor * new_bond.p_torque;
					other->torque += factor * new_bond.q_torque;
				}

				{
					BaseParticle *other = other_bond.other;

					number factor = -_lambda * curr_energy;
					LR_vector tmp_force = factor * other_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					p->torque -= factor * other_bond.p_torque;
					other->torque += factor * other_bond.q_torque;
				}
			}
		}
	}

	return energy;
}

void TorsionalDetailedPatchySwapInteraction::begin_energy_computation() {
	BaseInteraction::begin_energy_computation();

	for(int i = 0; i < _N; i++) {
		_particle_bonds(CONFIG_INFO->particles()[i]).clear();
	}
}

number TorsionalDetailedPatchySwapInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		if(q != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = _box->min_image(p->pos, q->pos);
		}
	}

	number energy = pair_interaction_bonded(p, q, false, update_forces);
	energy += pair_interaction_nonbonded(p, q, false, update_forces);
	return energy;
}

number TorsionalDetailedPatchySwapInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.;

	return energy;
}

number TorsionalDetailedPatchySwapInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	number energy = _spherical_patchy_two_body(p, q, false, update_forces);

	if(_is_KF) {
		energy += _patchy_two_body_KF(p, q, false, update_forces);
	}
	else {
		energy += _patchy_two_body_point(p, q, false, update_forces);
	}

	return energy;
}

void TorsionalDetailedPatchySwapInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	int N = particles.size();
	int curr_limit = _N_per_species[0];
	int curr_species = 0;
	for(int i = 0; i < N; i++) {
		if(i == curr_limit) {
			curr_species++;
			curr_limit += _N_per_species[curr_species];
		}
		if(_N_patches[curr_species] > 0 && _base_patches[curr_species].size() > 0) {
			particles[i] = new TorsionalPatchyParticle(_base_patches[curr_species], _a1_patches[curr_species], _a2_patches[curr_species],curr_species, 1.);
		}
		else {
			auto new_particle = new TorsionalPatchyParticle(_N_patches[curr_species], curr_species, 1.);
			particles[i] = new_particle;
			// we need to save the base patches so that the CUDA backend has access to them
			_base_patches[curr_species] = new_particle->base_patches();
            _a1_patches[curr_species] = new_particle->_a1_patches; //base_patches();
            _a2_patches[curr_species] = new_particle->_a2_patches; //base_patches();
            
		}
		particles[i]->index = i;
		particles[i]->strand_id = i;
		particles[i]->type = particles[i]->btype = curr_species;
	}
}

void TorsionalDetailedPatchySwapInteraction::_parse_interaction_matrix() {
	// parse the interaction matrix file
	input_file inter_matrix_file;
	inter_matrix_file.init_from_filename(_interaction_matrix_file);
	if(inter_matrix_file.state == ERROR) {
		throw oxDNAException("Caught an error while opening the interaction matrix file '%s'", _interaction_matrix_file.c_str());
	}

	_patchy_eps.resize(_N_patch_types * _N_patch_types, 0.);

	for(int i = 0; i < _N_patch_types; i++) {
		for(int j = 0; j < _N_patch_types; j++) {
			number value;
			std::string key = Utils::sformat("patchy_eps[%d][%d]", i, j);
			if(getInputNumber(&inter_matrix_file, key.c_str(), &value, 0) == KEY_FOUND) {
				_patchy_eps[i + _N_patch_types * j] = _patchy_eps[j + _N_patch_types * i] = value;
			}
		}
	}
}

void TorsionalDetailedPatchySwapInteraction::_parse_base_patches(std::string filename, int N_patches,std::vector<LR_vector> &base_patches,  std::vector<LR_vector>  &a1_patches,  std::vector<LR_vector>  &a2_patches) {
	std::ifstream patch_file(filename);
	if(!patch_file.good()) {
		throw oxDNAException("Can't read patch file '%s'. Aborting", filename.c_str());
	}

    base_patches.resize(N_patches);
    a1_patches.resize(N_patches);
    a2_patches.resize(N_patches);
    
	
	string line;
	for(int i = 0; i < N_patches; i++) {
		if(!patch_file.good()) {
			throw oxDNAException("The patch file '%s' does not seem to contain enough lines (%d found, should be %d)", filename.c_str(), i, N_patches);
		}
		std::getline(patch_file, line);
		auto spl = Utils::split(line);
		if(spl.size() != 9) {
			throw oxDNAException("Patch file '%s': invalid line '%s'", filename.c_str(), line.c_str());
		}
		LR_vector v;
		v[0] = std::stof(spl[0]);
		v[1] = std::stof(spl[1]);
		v[2] = std::stof(spl[2]);

        LR_vector a1;
        a1[0] = std::stof(spl[3]);
		a1[1] = std::stof(spl[4]);
		a1[2] = std::stof(spl[5]);

        LR_vector a2;
        a2[0] = std::stof(spl[6]);
		a2[1] = std::stof(spl[7]);
		a2[2] = std::stof(spl[8]);





		base_patches[i] = v;
        a1_patches[i] = a1;
        a2_patches[i] = a2;
	}

	patch_file.close();

	
}

void TorsionalDetailedPatchySwapInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	_N = particles.size();
	*N_strands = _N;

	std::ifstream topology(_topology_filename, ios::in);
	if(!topology.good()) {
		throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	}
	string line;
	int N_species;
	std::getline(topology, line);
	sscanf(line.c_str(), "%*d %d\n", &N_species);

	if(N_species < 1) {
		throw oxDNAException("The number of species should be larger than 0");
	}

	_N_per_species.resize(N_species);
	_N_patches.resize(N_species);
	_base_patches.resize(N_species);
    _a1_patches.resize(N_species);
    _a2_patches.resize(N_species);

	_patch_types.resize(N_species);
	int N_tot = 0;
	for(int i = 0; i < N_species; i++) {
		std::getline(topology, line);
		auto spl = Utils::split(line);
		if(spl.size() < 3) {
			throw oxDNAException("The topology line '%s' is malformed, since it should contain at least two integer numbers (number of particles, number of patches and a comma-separated list of patch types)", line.c_str());
		}
		int N_s = std::stoi(spl[0]);
		_N_per_species[i] = N_s;
		N_tot += N_s;
		_N_patches[i] = std::stoi(spl[1]);

		auto patch_types_str = Utils::split(spl[2], ',');
		if(patch_types_str.size() != _N_patches[i]) {
			throw oxDNAException("Species n. %d should have %d patches, but %d patch types are specified in the topology", i, _N_patches[i], patch_types_str.size());
		}
		std::vector<int> patch_types;
		for(auto type_str : patch_types_str) {
			int type = std::stoi(type_str);
			if(type < 0) {
				throw oxDNAException("Invalid patch type %d", type);
			}
			if(type >= _N_patch_types) {
				_N_patch_types = type + 1;
			}
			patch_types.push_back(type);
		}
		_patch_types[i] = patch_types;

		if(spl.size() > 3) {
            std::vector<LR_vector> base_patches;
            std::vector<LR_vector> a1_patches;
            std::vector<LR_vector> a2_patches;
              
            _parse_base_patches(spl[3], _N_patches[i],base_patches,a1_patches,a2_patches);
			_base_patches[i] = base_patches;
            _a1_patches[i] = a1_patches;
            _a2_patches[i] = a2_patches;
            
		}

		OX_LOG(Logger::LOG_INFO, "TorsionalDetailedPatchySwapInteraction: species %d has %d particles and %d patches", i, N_s, _N_patches[i]);
	}

	topology.close();

	if(N_tot != _N) {
		throw oxDNAException("The sum of the particles belonging to each species (%d) is different from the number of particles found in the first line of the topology file (%d)", N_tot, _N);
	}

	allocate_particles(particles);
	_parse_interaction_matrix();
}

void TorsionalDetailedPatchySwapInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}

extern "C" TorsionalDetailedPatchySwapInteraction* make_TorsionalDetailedPatchySwapInteraction() {
	return new TorsionalDetailedPatchySwapInteraction();
}


//DELETE ME

TorsionalPatchyParticle::TorsionalPatchyParticle(std::vector<LR_vector> base_patches, std::vector<LR_vector> a1_patches, std::vector<LR_vector> a2_patches, int nt, number sigma) :
				BaseParticle(),
				_sigma(sigma),
				_base_patches(base_patches),
				_a1_patches(a1_patches),
				_a2_patches(a2_patches)  {
	type = btype = nt;
	int_centers.resize(base_patches.size());
	_a1_patches.resize(base_patches.size());
	_a2_patches.resize(base_patches.size());

	for(uint i = 0; i < N_int_centers(); i++) {
		_base_patches[i].normalize();
		_base_patches[i] *= 0.5;

		_a1_patches[i].normalize();
		_a2_patches[i].normalize();



	}
}

TorsionalPatchyParticle::TorsionalPatchyParticle(int N_patches, int nt, number sigma) :
				BaseParticle(),
				_sigma(sigma) {
	type = btype = nt;
	int_centers.resize(N_patches);
	_base_patches.resize(N_patches);
	_a1_patches.resize(N_patches);
	_a2_patches.resize(N_patches);


	switch(N_int_centers()) {
	case 0:
		break;
	case 1: {
		_base_patches[0] = LR_vector(1, 0, 0);
		break;
	}
	case 2: {
		_base_patches[0] = LR_vector(0, 1, 0);
		_base_patches[1] = LR_vector(0, -1, 0);

		break;
	}
	case 3: {
		number cos30 = cos(M_PI / 6.);
		number sin30 = sin(M_PI / 6.);

		_base_patches[0] = LR_vector(0, 1, 0);
		_base_patches[1] = LR_vector(cos30, -sin30, 0);
		_base_patches[2] = LR_vector(-cos30, -sin30, 0);

		break;
	}
	case 4: {
		_base_patches[0] = LR_vector(-HALF_ISQRT3, -HALF_ISQRT3, HALF_ISQRT3);
		_base_patches[1] = LR_vector( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3);
		_base_patches[2] = LR_vector( HALF_ISQRT3, HALF_ISQRT3, HALF_ISQRT3);
		_base_patches[3] = LR_vector(-HALF_ISQRT3, HALF_ISQRT3, -HALF_ISQRT3);
		break;
	}
	case 5: {
		_base_patches[0] = LR_vector(0.135000, -0.657372, -0.741375);
		_base_patches[1] = LR_vector(0.259200, 0.957408, -0.127224);
		_base_patches[2] = LR_vector(-0.394215, -0.300066, 0.868651);
		_base_patches[3] = LR_vector(-0.916202, 0.202077, -0.346033);
		_base_patches[4] = LR_vector(0.916225, -0.202059, 0.345982);
		break;
	}
	default:
		throw oxDNAException("Unsupported number of patches %d\n", N_int_centers());
	}

	for(uint i = 0; i < N_int_centers(); i++) {
		_base_patches[i].normalize();
		_base_patches[i] *= 0.5;
	}
}

TorsionalPatchyParticle::~TorsionalPatchyParticle() {

}

void TorsionalPatchyParticle::set_positions() {
	for(uint i = 0; i < N_int_centers(); i++) {
		int_centers[i] = (orientation * _base_patches[i]) * _sigma;
		_a1_patches[i] = (_a1_patches[i].x * this->orientationT.v1) + (_a1_patches[i].y * this->orientationT.v2) + (_a1_patches[i].z * this->orientationT.v3); //possibly can be accelerated
		_a2_patches[i] = (_a2_patches[i].x * this->orientationT.v1) + (_a2_patches[i].y * this->orientationT.v2) + (_a2_patches[i].z * this->orientationT.v3); //possibly can be accelerated

	}
}



