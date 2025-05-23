/*

   Copyright (c) 2006-2010, The Scripps Research Institute
   Copyright (c) 2015, The University of Georgia

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute
           
   Modifications for Vina-Carb 1.0 By: Anita K. Nivedha <nivedha@uga.edu>
                                       The Woods' Lab
                                       Complex Carbohydrate Research Center
                                       The University of Georgia

*/

#include "parallel.h"
#include "parallel_mc.h"
#include "coords.h"
#include "parallel_progress.h"

struct parallel_mc_task {
	model m;
	output_container out;
	rng generator;
	parallel_mc_task(const model& m_, int seed) : m(m_), generator(static_cast<rng::result_type>(seed)) {
		VINA_FOR_IN(i, m.atoms){
			explicit_atom_update(m_.atoms[i], m.atoms[i]);
		}
		VINA_FOR_IN(i, m.grid_atoms){
			explicit_atom_update(m_.grid_atoms[i], m.grid_atoms[i]);
		}
		m.weight_chpi = m_.weight_chpi;
		m.weight_repulsion = m_.weight_repulsion;
		m.weight_gauss1 = m_.weight_gauss1;
		m.weight_gauss2 = m_.weight_gauss2;
                m.chpi_explicit_hydrogen = m_.chpi_explicit_hydrogen;
                m.build_residue_info();
                m.build_ar_ring_info();
		m.prec_factor = m_.prec_factor;
		m.build_chpi_smooth();
	}
};

typedef boost::ptr_vector<parallel_mc_task> parallel_mc_task_container;

struct parallel_mc_aux {
	const monte_carlo* mc;
	const precalculate* p;
	const igrid* ig;
	const precalculate* p_widened;
	const igrid* ig_widened;
	const vec* corner1;
	const vec* corner2;
	const fl chi_coeff;
	const fl chi_cutoff;
	parallel_progress* pg;
	parallel_mc_aux(const monte_carlo* mc_, const precalculate* p_, const igrid* ig_, const precalculate* p_widened_, const igrid* ig_widened_, const vec* corner1_, const vec* corner2_, const fl chi_coeff_, const fl chi_cutoff_, parallel_progress* pg_)
		: mc(mc_), p(p_), ig(ig_), p_widened(p_widened_), ig_widened(ig_widened_), corner1(corner1_), corner2(corner2_), chi_coeff(chi_coeff_), chi_cutoff(chi_cutoff_), pg(pg_) {}
	void operator()(parallel_mc_task& t) const { //Here the Monte Carlo run starts
		(*mc)(t.m, t.out, *p, *ig, *p_widened, *ig_widened, *corner1, *corner2, pg, t.generator, chi_coeff, chi_cutoff);
	}
};

void merge_output_containers(const output_container& in, output_container& out, fl min_rmsd, sz max_size) {
	VINA_FOR_IN(i, in)
		add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_mc_task_container& many, output_container& out, fl min_rmsd, sz max_size) {
	min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
	VINA_FOR_IN(i, many)
		merge_output_containers(many[i].out, out, min_rmsd, max_size);
	out.sort();
}

void parallel_mc::operator()(const model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator,const fl chi_coeff, const fl chi_cutoff) const {
	parallel_progress pp;
	parallel_mc_aux parallel_mc_aux_instance(&mc, &p, &ig, &p_widened, &ig_widened, &corner1, &corner2, chi_coeff, chi_cutoff, (display_progress ? (&pp) : NULL));
	parallel_mc_task_container task_container;
	VINA_FOR(i, num_tasks)
		task_container.push_back(new parallel_mc_task(m, random_int(0, 1000000, generator)));
	if(display_progress) 
		pp.init(num_tasks * mc.num_steps);
	parallel_iter<parallel_mc_aux, parallel_mc_task_container, parallel_mc_task, true> parallel_iter_instance(&parallel_mc_aux_instance, num_threads);
	parallel_iter_instance.run(task_container);
	merge_output_containers(task_container, out, mc.min_rmsd, mc.num_saved_mins);
}
