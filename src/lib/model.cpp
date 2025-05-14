/*

   Copyright (c) 2006+6010, The Scripps Research Institute
   Copyright (c) 2015, The University of Georgia

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE+6.0

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

#include "model.h"
#include "file.h"
#include "curl.h"
#include "parse_pdbqt.h"
#include <boost/thread/mutex.hpp>
#include <algorithm>
#include <iterator> //std::distance

namespace
	{
	boost::mutex cout_mutex;
	}



template<typename T>
atom_range get_atom_range(const T& t) {
	atom_range tmp = t.node;
	VINA_FOR_IN(i, t.children) {
		atom_range r = get_atom_range(t.children[i]);
		if(tmp.begin > r.begin) tmp.begin = r.begin;
		if(tmp.end   < r.end  ) tmp.end   = r.end;
	}
	return tmp;
}

struct branch_metrics {
	sz length;
	sz corner2corner;
	branch_metrics() : length(0), corner2corner(0) {}
};

template<typename T>
branch_metrics get_branch_metrics(const T& t) {
	branch_metrics tmp;
	if(!t.children.empty()) {
		sz corner2corner_max = 0;
		szv lengths;
		VINA_FOR_IN(i, t.children) {
			branch_metrics res = get_branch_metrics(t.children[i]);
			if(corner2corner_max < res.corner2corner)
				corner2corner_max = res.corner2corner;
			lengths.push_back(res.length + 1); // FIXME? weird compiler warning (sz -> unsigned)
		}
		std::sort(lengths.begin(), lengths.end());

		tmp.length = lengths.back();

		tmp.corner2corner = tmp.length;
		if(lengths.size() >= 2)
			tmp.corner2corner += lengths[lengths.size() - 1];

		if(tmp.corner2corner < corner2corner_max)
			tmp.corner2corner = corner2corner_max;
	}
	return tmp;
}

sz model::ligand_longest_branch(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).length;
}

sz model::ligand_length(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).corner2corner;
}

void ligand::set_range() {
	atom_range tmp = get_atom_range(*this);
	begin = tmp.begin;
	end   = tmp.end;
}

/////////////////// begin MODEL::APPEND /////////////////////////

// FIXME hairy code - needs to be extensively commented, asserted, reviewed and tested

struct appender_info {
	sz grid_atoms_size;
	sz m_num_movable_atoms;
	sz atoms_size;

	appender_info(const model& m) : grid_atoms_size(m.grid_atoms.size()), m_num_movable_atoms(m.m_num_movable_atoms), atoms_size(m.atoms.size()) {}
};

class appender {
	appender_info a_info;
	appender_info b_info;
	sz new_grid_index(sz x) const {
		return (is_a ? x : (a_info.grid_atoms_size + x)); // a-grid_atoms spliced before b-grid_atoms
	}
public:
	bool is_a;

	appender(const model& a, const model& b) : a_info(a), b_info(b), is_a(true) {}

	sz operator()(sz x) const { // transform coord index
		if(is_a) {
			if(x < a_info.m_num_movable_atoms)  return x; // a-movable unchanged
			else                                return x + b_info.m_num_movable_atoms; // b-movable spliced before a-inflex
		}
		else {
			if(x < b_info.m_num_movable_atoms)  return x + a_info.m_num_movable_atoms; // a-movable spliced before b-movable
			else                                return x + a_info.atoms_size; // all a's spliced before b-inflex
		}
	}
	atom_index operator()(const atom_index& x) const { // transform atom_index
		atom_index tmp(x);
		if(tmp.in_grid) tmp.i = new_grid_index(tmp.i);
		else            tmp.i = operator()(tmp.i);
			return tmp;
	}

	// type-directed old -> new transformations
	void update(interacting_pair& ip) const {
		ip.a = operator()(ip.a);
		ip.b = operator()(ip.b);
	}
	void update(vec& v) const { // coordinates & forces - do nothing
	}
	void update(ligand& lig) const {
		lig.transform(*this); // ligand as an atom_range subclass
		transform_ranges(lig, *this);
		VINA_FOR_IN(i, lig.pairs)
			this->update(lig.pairs[i]);
		VINA_FOR_IN(i, lig.cont)
			this->update(lig.cont[i]); // parsed_line update, below
	}
	void update(residue_vc& r) const {
		transform_ranges(r, *this);
	}
	void update(parsed_line& p) const {
		if(p.second)
			p.second = operator()(p.second.get());
	}
	void update(atom_vc& a) const {
		VINA_FOR_IN(i, a.bonds) {
			bond_vc& b = a.bonds[i];
			b.connected_atom_index = operator()(b.connected_atom_index); // atom_index transformation, above
		}
	}

	// ligands, flex, flex_context, atoms; also used for other_pairs
	template<typename T>
	void append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbbbbbbb
		sz a_sz = a.size();
		vector_append(a, b);

		is_a = true;
		VINA_FOR(i, a_sz)
			update(a[i]);

		is_a = false;
		VINA_RANGE(i, a_sz, a.size())
			update(a[i]);
	}

	// internal_coords, coords, minus_forces, atoms
	template<typename T>
	void coords_append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbaab
		std::vector<T> b_copy(b); // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

		is_a = true;
		VINA_FOR_IN(i, a)
			update(a[i]);

		is_a = false;
		VINA_FOR_IN(i, b_copy)
			update(b_copy[i]);

		// interleave 
		typedef typename std::vector<T>::const_iterator const cci;
		cci b1 = b_copy.begin();
		cci b2 = b_copy.begin() + b_info.m_num_movable_atoms;
		cci b3 = b_copy.end();

		a.insert(a.begin() + a_info.m_num_movable_atoms , b1 , b2);
		a.insert(a.end()                                , b2 , b3);
	}
};

void model::append(const model& m) {
	VINA_CHECK(atom_typing_used() == m.atom_typing_used());

	appender t(*this, m);

	t.append(other_pairs, m.other_pairs);

	//Append ligand aromatic rings
	//this->lig_ar_rings.insert(this->lig_ar_rings.end(), m.lig_ar_rings.begin(), m.lig_ar_rings.end());

	VINA_FOR_IN(i, atoms)
		VINA_FOR_IN(j, m.atoms) {
			if(i >= m_num_movable_atoms && j >= m.m_num_movable_atoms) continue; // no need for inflex-inflex interactions

			const atom_vc& a =   atoms[i];
			const atom_vc& b = m.atoms[j];

			sz t1 = a.get(atom_typing_used());
			sz t2 = b.get(atom_typing_used());
			sz n = num_atom_types(atom_typing_used());

			if(t1 < n && t2 < n) {
				t.is_a =  true;
				sz new_i = t(i);
				t.is_a = false;
				sz new_j = t(j);
				sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
				other_pairs.push_back(interacting_pair(type_pair_index, new_i, new_j));
			}
		}

	VINA_CHECK(  minus_forces.size() ==   coords.size());
	VINA_CHECK(m.minus_forces.size() == m.coords.size());

	t.coords_append(internal_coords, m.internal_coords);
	t.coords_append(         coords, m         .coords);
	t.coords_append(   minus_forces, m   .minus_forces); // for now, minus_forces.size() == coords.size() (includes inflex)

	t.append(ligands,         m.ligands);
	t.append(flex,            m.flex);
	t.append(flex_context,    m.flex_context);

	t       .append(grid_atoms, m.grid_atoms);
	t.coords_append(     atoms, m     .atoms);

	m_num_movable_atoms += m.m_num_movable_atoms;

}

///////////////////  end  MODEL::APPEND /////////////////////////


/////////////////// begin MODEL::INITIALIZE /////////////////////////

atom_index model::sz_to_atom_index(sz i) const {
	if(i < grid_atoms.size()) return atom_index(i                    ,  true);
	else                      return atom_index(i - grid_atoms.size(), false);
}

distance_type model::distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const {
	if(i.in_grid && j.in_grid) return DISTANCE_FIXED;
	if(i.in_grid) return (j.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	if(j.in_grid) return (i.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	assert(!i.in_grid);
	assert(!j.in_grid);
	assert(i.i < atoms.size());
	assert(j.i < atoms.size());
	sz a = i.i;
	sz b = j.i;
	if(a == b) return DISTANCE_FIXED;
	return (a < b) ? mobility(a, b) : mobility(b, a);
}

const vec& model::atom_coords(const atom_index& i) const {
	return i.in_grid ? grid_atoms[i.i].coords : coords[i.i];
}

fl model::distance_sqr_between(const atom_index& a, const atom_index& b) const {
	return vec_distance_sqr(atom_coords(a), atom_coords(b));
}

struct bond_less { // FIXME rm!?
	bool operator()(const bond_vc& a, const bond_vc& b) const {
		return a.connected_atom_index.i < b.connected_atom_index.i;
	}
};


bool model::atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const { // there is an atom closer to both a and b then they are to each other and immobile relative to them
	fl r2 = distance_sqr_between(a, b);
	VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
		sz i = relevant_atoms[relevant_atoms_i];
		atom_index c = sz_to_atom_index(i);
		if(a == c || b == c) continue;
		distance_type ac = distance_type_between(mobility, a, c);
		distance_type bc = distance_type_between(mobility, b, c);
		if(ac != DISTANCE_VARIABLE &&
		   bc != DISTANCE_VARIABLE &&
		   distance_sqr_between(a, c) < r2 &&
		   distance_sqr_between(b, c) < r2){
			return true;
                   }
	}
	return false;
}

struct beads {
	fl radius_sqr;
	std::vector<std::pair<vec, szv> > data;
	beads(sz reserve_size, fl radius_sqr_) : radius_sqr(radius_sqr_) { data.reserve(reserve_size); }
	void add(sz index, const vec& coords) {
		VINA_FOR_IN(i, data) {
			if(vec_distance_sqr(coords, data[i].first) < radius_sqr) {
				data[i].second.push_back(index);
				return;
			}
		}
		// not found
		std::pair<vec, szv> tmp;
		tmp.first = coords;
		tmp.second.push_back(index);
		data.push_back(tmp);
	}
};

void model::assign_bonds(const distance_type_matrix& mobility) { // assign bonds based on relative mobility, distance and covalent length
	const fl bond_length_allowance_factor = 1.1;
	sz n = grid_atoms.size() + atoms.size();

	// construct beads
	const fl bead_radius = 15;
	beads beads_instance(n, sqr(bead_radius));
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		beads_instance.add(i, atom_coords(i_atom_index));
	}
	// assign bonds
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		const vec& i_atom_coords = atom_coords(i_atom_index);
		atom_vc& i_atom = get_atom(i_atom_index);
		const fl max_covalent_r = max_covalent_radius(); // FIXME mv to atom_constants
		fl i_atom_covalent_radius = max_covalent_r;
		if(i_atom.ad < AD_TYPE_SIZE)
			i_atom_covalent_radius = ad_type_property(i_atom.ad).covalent_radius;

		//find relevant atoms
		szv relevant_atoms;
		const fl bead_cutoff_sqr = sqr(bead_radius + bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r));
		VINA_FOR_IN(b, beads_instance.data) {
			if(vec_distance_sqr(beads_instance.data[b].first, i_atom_coords) > bead_cutoff_sqr) continue;
			const szv& bead_elements = beads_instance.data[b].second;
			VINA_FOR_IN(bead_elements_i, bead_elements) {
				sz j = bead_elements[bead_elements_i];
				atom_index j_atom_index = sz_to_atom_index(j);
				atom_vc& j_atom = get_atom(j_atom_index);
				const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
				distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
                                
       
				if(dt != DISTANCE_VARIABLE && i != j) {
					fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
					//if(r2 < sqr(bond_length_allowance_factor * bond_length))
					if(r2 < sqr(bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r)))
						relevant_atoms.push_back(j);
				}
			}
		}
		// find bonded atoms
		VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
			sz j = relevant_atoms[relevant_atoms_i];
			if(j <= i) continue; // already considered
			atom_index j_atom_index = sz_to_atom_index(j);
			atom_vc& j_atom = get_atom(j_atom_index);
			const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
			distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
			fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
			if(r2 < sqr(bond_length_allowance_factor * bond_length) && !atom_exists_between(mobility, i_atom_index, j_atom_index, relevant_atoms)) {
				bool rotatable = (dt == DISTANCE_ROTOR);
				fl length = std::sqrt(r2);
				i_atom.bonds.push_back(bond_vc(j_atom_index, length, rotatable));
				j_atom.bonds.push_back(bond_vc(i_atom_index, length, rotatable));
			}

		}
	}
}

bool model::bonded_to_HD(const atom_vc& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond_vc& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).ad == AD_TYPE_HD) 
			return true;
	}
	return false;
}

bool model::bonded_to_heteroatom(const atom_vc& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond_vc& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).is_heteroatom())
			return true;
	}
	return false;
}

bool model::bonded_to_hydrogen(const atom_vc& a) const {
        VINA_FOR_IN(i, a.bonds) {
                const bond_vc& b = a.bonds[i];
		const sz& ad = get_atom(b.connected_atom_index).ad;
                if(ad == AD_TYPE_HD || ad == AD_TYPE_H)
                        return true;
        }
        return false;
}

bool model::chpi_contact_possible(const atom_vc& central_carbon, vec& ring_centroid, vec& ring_normal, aptrv& ring, fl& ring_effective_radius, bool ring_is_receptor){
	if (!bonded_to_hydrogen(central_carbon)){
		return false;
	}

	VINA_FOR_IN(i, central_carbon.bonds) {
		const bond_vc& b = central_carbon.bonds[i];
		const atom_vc& neighbor = get_atom(b.connected_atom_index);
		const sz& ad = neighbor.ad;

		if (ad == AD_TYPE_HD || ad == AD_TYPE_H){
			//Check if this central carbon has at least one hydrogen within 3.5 A of the ring centroid
			const vec& h_coord = neighbor.coords;
			fl dx1 = h_coord.data[0] - ring_centroid.data[0];
			fl dy1 = h_coord.data[1] - ring_centroid.data[1];
			fl dz1 = h_coord.data[2] - ring_centroid.data[2];

			vec h_centroid(dx1, dy1, dz1);
			fl r1 = std::sqrt(sqr(dx1) + sqr(dy1) + sqr(dz1));

        		fl dot_product = h_centroid * ring_normal;
			fl hc_length = magnitude(h_centroid);
        		fl cosine = dot_product / hc_length;
			fl vo = std::abs(hc_length * cosine);
			fl ho = std::sqrt(hc_length*hc_length - vo*vo);

			if (vo < 4.5 && ho < ring_effective_radius){
				return true;
			}
			/*if (r1 < chpi_h_cutoff){
				return true;
			}

			//Or, at least one hydrogen within 3.5 A of any ring atom. 
			VINA_FOR_IN(j, ring){
				atom_vc* ra = ring[j];
				vec* ra_coord = NULL;
				if (ring_is_receptor){
					ra_coord = &(ra->coords);
				}
				else{
					sz nrp_index = get_nrp_atom_index(ra);
                        		ra_coord = &(this->coords[nrp_index]);
				}
				VINA_CHECK(ra_coord != NULL);

				fl dx2 = h_coord.data[0] - ra_coord->data[0];
				fl dy2 = h_coord.data[1] - ra_coord->data[1];
				fl dz2 = h_coord.data[2] - ra_coord->data[2];
				fl r1 = std::sqrt(sqr(dx1) + sqr(dy1) + sqr(dz1));
				fl r2 = std::sqrt(sqr(dx2) + sqr(dy2) + sqr(dz2));

				if (r2 < chpi_h_cutoff){
					//Comment this out and see what happens.
					//return true;
				}

			}*/

		}
	}

	return false;
}

void model::assign_types() {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		atom_vc& a = get_atom(ai);
		a.assign_el();
		sz& x = a.xs;

		bool acceptor   = (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA); // X-Score forumaltion apparently ignores SA
		bool donor_NorO = (a.el == EL_TYPE_Met || bonded_to_HD(a));

		switch(a.el) {
			case EL_TYPE_H    : break;
			//case EL_TYPE_C    : x = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H; break; //This is the original code
                        case EL_TYPE_C    : x = (a.ad == AD_TYPE_A ? (bonded_to_heteroatom(a) ? XS_TYPE_C_A_P : XS_TYPE_C_A_H) : (bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H));break;
			case EL_TYPE_N    : x = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break;
			case EL_TYPE_O    : x = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break;
			case EL_TYPE_S    : x = XS_TYPE_S_P; break;
			case EL_TYPE_P    : x = XS_TYPE_P_P; break;
			case EL_TYPE_F    : x = XS_TYPE_F_H; break;
			case EL_TYPE_Cl   : x = XS_TYPE_Cl_H; break;
			case EL_TYPE_Br   : x = XS_TYPE_Br_H; break;
			case EL_TYPE_I    : x = XS_TYPE_I_H; break;
			case EL_TYPE_Met  : x = XS_TYPE_Met_D; break;
			case EL_TYPE_NK   : x = XS_TYPE_N_K; break; //Amika
                        case EL_TYPE_NR   : x = XS_TYPE_N_R; break; //Amika
                        case EL_TYPE_ON   : x = XS_TYPE_O_N; break; //Amika
                        case EL_TYPE_OS   : x = XS_TYPE_O_S; break; //Amika
            		case EL_TYPE_SO   : x = XS_TYPE_S_O; break; //Amika
			case EL_TYPE_SIZE : break;
			default: VINA_CHECK(false);
		}
	}
}

sz model::find_ligand(sz a) const {
	VINA_FOR_IN(i, ligands)
		if(a >= ligands[i].begin && a < ligands[i].end)
			return i;
	return ligands.size();
}

void model::bonded_to(sz a, sz n, szv& out) const {
	if(!has(out, a)) { // not found
		out.push_back(a);
		if(n > 0) 
			VINA_FOR_IN(i, atoms[a].bonds) {
				const bond_vc& b = atoms[a].bonds[i];
				if(!b.connected_atom_index.in_grid)
					bonded_to(b.connected_atom_index.i, n-1, out);
			}
	}
}

szv model::bonded_to(sz a, sz n) const {
	szv tmp;
	bonded_to(a, n, tmp);
	return tmp;
}


void model::initialize_pairs(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, atoms) {
		sz i_lig = find_ligand(i);
		szv bonded_atoms = bonded_to(i, 3);
		VINA_RANGE(j, i+1, atoms.size()) {
			if(i >= m_num_movable_atoms && j >= m_num_movable_atoms) continue; // exclude inflex-inflex
			if(mobility(i, j) == DISTANCE_VARIABLE && !has(bonded_atoms, j)) {
				sz t1 = atoms[i].get  (atom_typing_used());
				sz t2 = atoms[j].get  (atom_typing_used());
				sz n  = num_atom_types(atom_typing_used());
				if(t1 < n && t2 < n) { // exclude, say, Hydrogens
					sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
					interacting_pair ip(type_pair_index, i, j);
					if(i_lig < ligands.size() && find_ligand(j) == i_lig)
						ligands[i_lig].pairs.push_back(ip);
					else
						other_pairs.push_back(ip);
				}
			}
		}
	}
}

void model::initialize(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, ligands)
		ligands[i].set_range();
	assign_bonds(mobility);
	assign_types();
	initialize_pairs(mobility);
}

///////////////////  end  MODEL::INITIALIZE /////////////////////////


sz model::num_internal_pairs() const {
	sz tmp = 0;
	VINA_FOR_IN(i, ligands)
		tmp += ligands[i].pairs.size();
	return tmp;
}

szv model::get_movable_atom_types(atom_type::t atom_typing_used_) const {
	szv tmp;
	sz n = num_atom_types(atom_typing_used_);
	VINA_FOR(i, m_num_movable_atoms) {
		const atom_vc& a = atoms[i];
		sz t = a.get(atom_typing_used_);
		if(t < n && !has(tmp, t))
			tmp.push_back(t);
	}
	return tmp;
}

conf_size model::get_size() const {
	conf_size tmp;
	tmp.ligands = ligands.count_torsions();
	tmp.flex    = flex   .count_torsions();
	return tmp;
}

conf model::get_initial_conf() const { // torsions = 0, orientations = identity, ligand positions = current
	conf_size cs = get_size();
	conf tmp(cs);
	tmp.set_to_null();
	VINA_FOR_IN(i, ligands)
		tmp.ligands[i].rigid.position = ligands[i].node.get_origin();
	return tmp;
}

grid_dims model::movable_atoms_box(fl add_to_each_dimension, fl granularity) const {
	vec corner1(0, 0, 0), corner2(0, 0, 0);
	VINA_FOR(i, num_movable_atoms()) {
		const vec& v = movable_coords(i);
		VINA_FOR_IN(j, v) {
			if(i == 0 || v[j] < corner1[j]) corner1[j] = v[j];
			if(i == 0 || v[j] > corner2[j]) corner2[j] = v[j];
		}
	}
	corner1 -= add_to_each_dimension / 2;
	corner2 += add_to_each_dimension / 2;

	grid_dims gd;
	{ // always doing this now FIXME ?
		vec center; center = 0.5 * (corner2 + corner1);
		VINA_FOR_IN(i, gd) {
			gd[i].n = sz(std::ceil((corner2[i] - corner1[i]) / granularity));
			fl real_span = granularity * gd[i].n;
			gd[i].begin = center[i] - real_span/2;
			gd[i].end = gd[i].begin + real_span;
		}
	}
	return gd;
}

void string_write_coord(sz i, fl x, std::string& str) {
	VINA_CHECK(i > 0);
	--i;
	std::ostringstream out;
	out.setf(std::ios::fixed, std::ios::floatfield);
	out.setf(std::ios::showpoint);
	out << std::setw(8) << std::setprecision(3) << x;
	VINA_CHECK(out.str().size() == 8); 
	VINA_CHECK(str.size() > i + 8);
	VINA_FOR(j, 8)
		str[i+j] = out.str()[j];
}
std::string coords_to_pdbqt_string(const vec& coords, const std::string& str) { 
	std::string tmp(str);
	string_write_coord(31, coords[0], tmp);
	string_write_coord(39, coords[1], tmp);
	string_write_coord(47, coords[2], tmp);
	return tmp;
}

void model::write_context(const context& c, ofile& out) const {
	verify_bond_lengths();
	VINA_FOR_IN(i, c) {
		const std::string& str = c[i].first;
		if(c[i].second) {
			out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
		}
		else
			out << str << '\n';
	}
}

void model::seti(const conf& c) {
	ligands.set_conf(atoms, internal_coords, c.ligands);
}

void model::sete(const conf& c) {
	VINA_FOR_IN(i, ligands)
		c.ligands[i].rigid.apply(internal_coords, coords, ligands[i].begin, ligands[i].end);
	flex.set_conf(atoms, coords, c.flex);
}

void model::set         (const conf& c) {
	ligands.set_conf(atoms, coords, c.ligands);
	//std::exit(1);
	flex   .set_conf(atoms, coords, c.flex);
}

fl model::gyration_radius(sz ligand_number) const {
	VINA_CHECK(ligand_number < ligands.size());
	const ligand& lig = ligands[ligand_number];
	fl acc = 0;
	unsigned counter = 0;
	VINA_RANGE(i, lig.begin, lig.end) {
		if(atoms[i].el != EL_TYPE_H) { // only heavy atoms are used
			acc += vec_distance_sqr(coords[i], lig.node.get_origin()); // FIXME? check!
			++counter;
		}
	}
	return (counter > 0) ? std::sqrt(acc/counter) : 0;
}


fl eval_interacting_pairs(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords) { // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		fl r2 = vec_distance_sqr(coords[ip.a], coords[ip.b]);
		if(r2 < cutoff_sqr) {
			fl tmp = p.eval_fast(ip.type_pair_index, r2);
			curl(tmp, v);
			e += tmp;
		}
	}
	return e;
}

fl eval_interacting_pairs_deriv(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces) { // adds to forces  // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	int count=0; 
	VINA_FOR_IN(i, pairs) {
		count++;
		const interacting_pair& ip = pairs[i];  
		vec r; r = coords[ip.b] - coords[ip.a]; // a -> b
		fl r2 = sqr(r);
		if(r2 < cutoff_sqr) {
			pr tmp = p.eval_deriv(ip.type_pair_index, r2); 
			vec force; force = tmp.second * r;
			curl(tmp.first, force, v);
			e += tmp.first;
			// FIXME inefficient, if using hard curl
			forces[ip.a] -= force; // we could omit forces on inflex here
			forces[ip.b] += force;
		}
	}
	return e;
}

fl model::evali(const precalculate& p,                                  const vec& v                          ) const { // clean up
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, internal_coords); // probably might was well use coords here
	return e;
}

fl model::evale(const precalculate& p, const igrid& ig, const vec& v                          ) const { // clean up
	fl e = ig.eval(*this, v[1]);
	e += eval_interacting_pairs(p, v[2], other_pairs, coords);
	return e;
}

double model::phi_alpha_energy(double phi_angle)
{
double LH = 2.97696467271672, Lc = -199.494365163839, LW = 677.808323900125, RH = 102.253303636096, Rc = 170.599580473404, RW = 1696.78443699429, aH = 10.7448005875571, ac = -105.313553566706, aW = 4724.58364072706, bH = 3.67344580413578, bc = 6.20116671874232, bW = 1347.72056251564, cH = 2.06094652655659, cc = 91.6553021324274, cW = 1500.02002601097, Off = 1.00501e-30, dH = 6.19388683252667, dc = -22.9786969888816, dW = 2122.27783139301, eH = -2.11153017593601, ec = 83.6019123356148, eW = 1254.13371108961, fH = -98.0013005657107, fc = 170.012289132741, fW = 1598.73272567307, Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=phi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
fx = fH * exp(-pow((x-(fc)),2.0)/fW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + fx;
return Totx;
}



double model::phi_beta_energy(double phi_angle)
{
double Lc = -330.769995527134, aH = 5.93533323829663, ac = -152.080139620062, aW = 6049.77220005964, bH = 22.467372096061, bc = -23.5159916173247, bW = 606.89715970453, cH = 10.0360057033439, cc = 120.962836525241, cW = 4037.89330459304, LH = 450.540038600828, LW = 4449.7622241787, RH = 23.7118506901333, Rc = 304.624980492529, RW = 8375.1929028027, /*Off = -2.27829251796721,*/ Off = -2.1283, dH = -18.1406478565247, dc = -24.2677756921736, dW = 543.050986049266, eH = 5.88226333077368, ec = 19.6321032903376, eW = 897.92664572344, Leftx, Rightx, ax, bx, cx, dx, ex, x, Totx;
x=phi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + Off;
return Totx;
}

double model::psi_2A3E_energy(double psi_angle)
{
double LH = 4.62366277694845, Lc = 5.045583934308, LW = 5005.75236060956, RH = 4.61387450239844, Rc = 362.487847702007, RW = 2090.63190217702, aH = 4.94191727813274, ac = 121.202321824468, aW = 2093.75214491931, bH = 0.402901504292045, bc = 241.428583877882, bW = 456.828754790442, cH = 0.798876573705798, cc = 68.425080241155, cW = 678.807178379645, Off = -0.125645118474882, dc = 192.925748017071, dW = 347.244734136509, dH = 0.222992242737354, Leftx, Rightx, ax, bx, cx, dx, x, Totx;
if(psi_angle<0)
{
psi_angle=360+psi_angle;}
x=psi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + Off;
return Totx;
}

double model::psi_2E3A_energy(double psi_angle)
{
double LH = 4.46811874171788, Lc = 1e-30, LW = 1279.58772056199, RH = 4.38204018882823, Rc = 357.770654336205, RW = 6050.14162479438, aH = 284.944948778136, ac = 146.644068129462, aW = 1551.75673776163, bH = 4.76134025362478, bc = 220.683118921686, bW = 5892.94143218231, cH = -169.197666368856, cc = 147.370828680332, cW = 1742.47541063603, Off = 1.0219924486158, dc = 146.05660843428, dW = 1359.82873591396, dH = -118.440552792375, Leftx, Rightx, ax, bx, cx, dx, x, Totx;
if(psi_angle<0){
psi_angle=360+psi_angle;}
x=psi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + Off;
return Totx;
}

double model::psi_6A_energy(double psi_angle)
{
double aH = 67.9431348410598, ac = -59.5393395706705, aW = 993.323581145538, bH = 6.13421142432396, bc = 10.4786088782815, bW = 945.770771330812, cH = 3.27628727235978, cc = 54.2960678151208, cW = 851.528141801851, dH = 0.727486729062442, dc = 131.067737803489, dW = 1037.41211378392, eH = 2.57362265878937, ec = 245.102891425541, eW = 2012.99451568206, fH = 5.75995973448166, fc = 359.999988549478, fW = 1153.3974275673, gH = 3.47492643928157, gc = 321.677942414686, gW = 2080.97053159226, hH = -0.741000462200939, hc = 199.106903524814, hW = 522.180434119001, ax, bx, cx, dx, ex, fx, gx, hx, x, Totx;

if(psi_angle<0)
{
psi_angle=360+psi_angle;
}
x=psi_angle;
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
fx = fH * exp(-pow((x-(fc)),2.0)/fW);
gx = gH * exp(-pow((x-(gc)),2.0)/gW);
hx = hH * exp(-pow((x-(hc)),2.0)/hW);
Totx = ax + bx + cx + dx + ex + fx + gx + hx;
return Totx;
}

double model::psi_6E_energy(double psi_angle)
{
double aH = 7.24858655753829, ac = 3.60600554520403, aW = 2459.23916629141, bH = 1.9, bc = 96.5930821702371, bW = 2683.88656516991, cH = 0.741022592342903, cc = 141.663521919709, cW = 1150.04756181103, dH = 0.2, dc = 162, dW = 400, eH = 0.287090039932611, ec = 228.171314273305, eW = 272.201363844744, fH = 1.22591971967808, fc = 292.206221787048, fW = 1134.52455512381, gH = 7.41063235334191, gc = 369.867701147817, gW = 3499.15994772992, hH = -0.61489499584011, hc = 271.319024293053, hW = 532.437194483944, iH = -0.35, ic = 183, iW = 100, ax, bx, cx, dx, ex, fx, gx, hx, ix, x, Totx;

if(psi_angle<0)
{
psi_angle=360+psi_angle;
}
x=psi_angle;
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
fx = fH * exp(-pow((x-(fc)),2.0)/fW);
gx = gH * exp(-pow((x-(gc)),2.0)/gW);
hx = hH * exp(-pow((x-(hc)),2.0)/hW);
ix = iH * exp(-pow((x-(ic)),2.0)/iW);
Totx = ax + bx + cx + dx + ex + fx + gx + hx + ix;
return Totx;
}

double model::omega_6A_energy(double omega_angle)
{
double x, energy, b, k=0.0025;
        if(omega_angle<0)
        {
        x=360+omega_angle;
        }
        else
        {
        x=omega_angle;
        }
        if((x>=0.0 && x<120.0)||(x>=360.0 && x<120.0))
        {
        b=0.0;
        energy=k*pow((x-60),2)+b;
        }
        else if(x>=120.0 && x<240.0)
        {
        b=0.3;
        energy=k*pow((x-180),2)+b;
        }
        else if(x>=240.0 && x<360.0)
        {
        b=1.0;
        energy=k*pow((x-300),2)+b;
        }
return energy;
}


double model::omega_6E_energy(double omega_angle)
{
double x, energy, b, k=0.0025;
        if(omega_angle<0)
        {
        x=360+omega_angle;
        }
        else
        {
        x=omega_angle;
        }
        if((x>=0.0 && x<120.0)||(x>=360.0 && x<120.0))
        {
        b=0.21;
        energy=k*pow((x-60),2)+b;
        }
        else if(x>=120.0 && x<240.0)
        {
        b=1.39;
        energy=k*pow((x-180),2)+b;
        }
        else if(x>=240.0 && x<360.0)
        {
        b=0.0;
        energy=k*pow((x-300),2)+b;
        }
return energy;
}



fl model::eval         (const precalculate& p, const igrid& ig, const vec& v, const conf& c           ) { // clean up
	set(c);
	fl e = evale(p, ig, v);
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords
	return e;
}

double model::get_torsion_coords_vecs_list(vec A, vec B, vec C, vec D)
{
double angle=0.0;
vec AB, BC, CD, ABCcross, BCDcross;
double ABC_BCD_dot, AB_BCD_dot, BCscalar_AB_BCD_dot;
AB=B-A;
BC=C-B;
CD=D-C;
ABCcross=cross_product(AB,BC);
BCDcross=cross_product(BC,CD);
ABC_BCD_dot=dot_product(ABCcross,BCDcross);
AB_BCD_dot=dot_product(AB,BCDcross);
BCscalar_AB_BCD_dot=magnitude(BC)*AB_BCD_dot;                        
angle=atan2(BCscalar_AB_BCD_dot,ABC_BCD_dot)*57.2957795; //angle in degrees
return angle;           
}                       

vecv model::get_flexible_coords(){
return coords;
}


fl model::eval_chi(const fl chi_coeff, const fl chi_cutoff)
{
fl e=0.0;
if(chi_coeff!=0){
	double phi_energy=0.0, psi_energy=0.0,omega_energy=0.0, phi=0.0, psi=0.0, omega=0.0, total=0.0,current_energy=0.0;
	std::vector< std::vector<size_t*> > glyco_info=glycan_info_func();
	VINA_FOR(i,glyco_info.size())
	{
	if(glyco_info[i][10][0]!=0){
	phi=get_torsion_coords_vecs_list(coords[glyco_info[i][0][0]],coords[glyco_info[i][1][0]],coords[glyco_info[i][2][0]],coords[glyco_info[i][3][0]]);
		if(glyco_info[i][5][0]==0)
		{//Alpha
			if(glyco_info[i][10][0]==1)
			{//D sugar (Alpha)
			current_energy=phi_alpha_energy(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}
			}
			else if(glyco_info[i][10][0]==2)
			{//L-Sugar (Alpha)
			current_energy=phi_alpha_energy(-phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
                                current_energy=0.0;
				}
			}
		}
		else
		{//Beta
			if(glyco_info[i][10][0]==1)
			{//D sugar (Beta)
			current_energy=phi_beta_energy(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}
			}
			else if(glyco_info[i][10][0]==2)
			{//L Sugar (Beta)
			current_energy=phi_beta_energy(-phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}
			}
		}
	}
	if(glyco_info[i][11][0]!=0){
	psi=get_torsion_coords_vecs_list(coords[glyco_info[i][1][0]],coords[glyco_info[i][2][0]],coords[glyco_info[i][3][0]],coords[glyco_info[i][4][0]]);
	omega=get_torsion_coords_vecs_list(coords[glyco_info[i][2][0]],coords[glyco_info[i][3][0]],coords[glyco_info[i][4][0]],coords[glyco_info[i][9][0]]);
		if(glyco_info[i][7][0]==2 || glyco_info[i][7][0]==4)
		{//2 or 4 linkage
			if(glyco_info[i][6][0]==0)
			{//Axial attachment
				if(glyco_info[i][11][0]==2)
				{//L sugar
				current_energy=psi_2A3E_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
					current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
				current_energy=psi_2A3E_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			}
			else
			{//Equatorial attachment
				if(glyco_info[i][11][0]==2)
				{//L sugar
				current_energy=psi_2E3A_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
					current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
				current_energy=psi_2E3A_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			
			}
		}
		else if(glyco_info[i][7][0]==3)
		{//3 linkage
			if(glyco_info[i][6][0]==0)
			{//Axial attachment
			
				if(glyco_info[i][11][0]==2)
				{//L sugar
				current_energy=psi_2E3A_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
                                current_energy=psi_2E3A_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			}
			else
			{//Equatorial attachment
				if(glyco_info[i][11][0]==2)
				{//L sugar
                                current_energy=psi_2A3E_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
				current_energy=psi_2A3E_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			
			}
		}
	else if(glyco_info[i][7][0]==6)
                {//6 linkage
                        if(glyco_info[i][8][0]==0)
                                {
                                current_energy=psi_6A_energy(psi);
                                        if(current_energy>chi_cutoff){
                                        psi_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                current_energy=omega_6A_energy(omega);
                                        if(current_energy>chi_cutoff){
                                        omega_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                }
                        else if(glyco_info[i][8][0]==1)
                                {
                                current_energy=psi_6E_energy(psi);
                                        if(current_energy>chi_cutoff){
                                        psi_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                current_energy=omega_6E_energy(omega);
                                        if(current_energy>chi_cutoff){
                                        omega_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                }
                }
	}
	
	}
	total=phi_energy+psi_energy+omega_energy;
	total=total*chi_coeff; 
	double combined_score=total+e;
	e=combined_score;
	}
	else
	e=0.0;
return e;
}


fl  model::eval_deriv  (const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g/*, std::vector< std::vector<int> > glyco_info*/,const fl chi_coeff, const fl chi_cutoff) { // clean up
	set(c);

	fl chi_energy = 0.0;
	chi_energy=eval_chi(chi_coeff,chi_cutoff);

	fl e = ig.eval_deriv(*this, v[1]); // sets minus_forces, except inflex

	e += this->eval_chpi(false);
	
        e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, minus_forces); // adds to minus_forces
        VINA_FOR_IN(i, ligands){
                e += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords, minus_forces); // adds to minus_forces
        }

        ligands.derivative(coords, minus_forces, g.ligands);
	flex.derivative(coords, minus_forces, g.flex); // inflex forces are ignored
	return (e+chi_energy);
}


fl model::eval_intramolecular(const precalculate& p, const vec& v, const conf& c) {
	set(c);
	fl e = 0;

	// internal for each ligand
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords

	sz nat = num_atom_types(atom_typing_used());
	const fl cutoff_sqr = p.cutoff_sqr();

	// flex-rigid
	VINA_FOR(i, num_movable_atoms()) {
		if(find_ligand(i) < ligands.size()) continue; // we only want flex-rigid interaction
		const atom_vc& a = atoms[i];
		sz t1 = a.get(atom_typing_used());
		if(t1 >= nat) continue;
		VINA_FOR_IN(j, grid_atoms) {
			const atom_vc& b = grid_atoms[j];
			sz t2 = b.get(atom_typing_used());
			if(t2 >= nat) continue;
			fl r2 = vec_distance_sqr(coords[i], b.coords);
			if(r2 < cutoff_sqr) {
				sz type_pair_index = triangular_matrix_index_permissive(nat, t1, t2);
				fl this_e = p.eval_fast(type_pair_index, r2);
				curl(this_e, v[1]);
				e += this_e;
			}
		}
	}

	// flex-flex
	VINA_FOR_IN(i, other_pairs) {
		const interacting_pair& pair = other_pairs[i];
		if(find_ligand(pair.a) < ligands.size() || find_ligand(pair.b) < ligands.size()) continue; // we only need flex-flex
		fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
		if(r2 < cutoff_sqr) {
			fl this_e = p.eval_fast(pair.type_pair_index, r2);
			curl(this_e, v[2]);
			e += this_e;
		}
	}
	return e;
}


fl model::eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy, bool score_in_place) {
	/*VINA_FOR_IN(i, this->coords){ //Yao 20231031: Is this necessary?
                this->adjusted_coords[i] = this->coords[i];
        }
	this->adjusted_atom_indices.clear();*/

	fl e = eval(p, ig, v, c); // sets c

	e += this->eval_chpi(score_in_place);
	return sf.conf_independent(*this, e - intramolecular_energy);
}

fl model::rmsd_lower_bound_asymmetric(const model& x, const model& y) const { // actually static
	sz n = x.m_num_movable_atoms; 
	VINA_CHECK(n == y.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, n) {
		const atom_vc& a =   x.atoms[i];
		if(a.el != EL_TYPE_H) {
			fl r2 = max_fl;
			VINA_FOR(j, n) {
				const atom_vc& b = y.atoms[j];
				if(a.same_element(b) && !b.is_hydrogen()) {
					fl this_r2 = vec_distance_sqr(x.coords[i], 
					                              y.coords[j]);
					if(this_r2 < r2)
						r2 = this_r2;
				}
			}
			assert(not_max(r2));
			sum += r2;
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_lower_bound(const model& m) const {
	return (std::max)(rmsd_lower_bound_asymmetric(*this,     m),
		            rmsd_lower_bound_asymmetric(    m, *this));
}

fl model::rmsd_upper_bound(const model& m) const {
	VINA_CHECK(m_num_movable_atoms == m.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, m_num_movable_atoms) {
		const atom_vc& a =   atoms[i];
		const atom_vc& b = m.atoms[i];
		assert(a.ad == b.ad);
		assert(a.xs == b.xs);
		if(a.el != EL_TYPE_H) {
			sum += vec_distance_sqr(coords[i], m.coords[i]);
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_ligands_upper_bound(const model& m) const {
	VINA_CHECK(ligands.size() == m.ligands.size());
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR_IN(ligand_i, ligands) {
		const ligand&   lig =   ligands[ligand_i];
		const ligand& m_lig = m.ligands[ligand_i];
		VINA_CHECK(lig.begin == m_lig.begin);
		VINA_CHECK(lig.end   == m_lig.end);
		VINA_RANGE(i, lig.begin, lig.end) {
			const atom_vc& a =   atoms[i];
			const atom_vc& b = m.atoms[i];
			assert(a.ad == b.ad);
			assert(a.xs == b.xs);
			if(a.el != EL_TYPE_H) {
				sum += vec_distance_sqr(coords[i], m.coords[i]);
				++counter;
			}
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}


void model::verify_bond_lengths() const {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		const atom_vc& a = get_atom(ai);
		VINA_FOR_IN(j, a.bonds) {
			const bond_vc& b = a.bonds[j];
			fl d = std::sqrt(distance_sqr_between(ai, b.connected_atom_index));
			bool ok = eq(d, b.length);
			if(!ok) {
				VINA_SHOW(d);
				VINA_SHOW(b.length);
			}
			VINA_CHECK(ok);
		}
	}
}

void model::check_internal_pairs() const {
	VINA_FOR_IN(i, ligands) {
		const ligand& lig = ligands[i];
		const interacting_pairs& pairs = lig.pairs;
		VINA_FOR_IN(j, pairs) {
			const interacting_pair& ip = pairs[j];
			VINA_CHECK(ip.a >= lig.begin);
			VINA_CHECK(ip.b  < lig.end);
		}
	}
}

void model::about() const {
	VINA_SHOW(atom_typing_used());
	VINA_SHOW(num_movable_atoms());
	VINA_SHOW(num_internal_pairs());
	VINA_SHOW(num_other_pairs());
	VINA_SHOW(num_ligands());
	VINA_SHOW(num_flex());
}

void model::print_stuff() const {
	std::cout << "coords:\n";
	VINA_FOR_IN(i, coords)
		printnl(coords[i]);

	std::cout << "internal_coords:\n";
	VINA_FOR_IN(i, internal_coords)
		printnl(internal_coords[i]);

	std::cout << "atoms:\n";
	VINA_FOR_IN(i, atoms) {
		const atom_vc& a = atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}

	std::cout << "grid_atoms:\n";
	VINA_FOR_IN(i, grid_atoms) {
		const atom_vc& a = grid_atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}
	about();
}

fl pairwise_clash_penalty(fl r, fl covalent_r) {
	// r = 0          -> max_penalty 
	// r = covalent_r -> 1
	// elsewhere      -> hyperbolic function
	assert(r >= 0);
	assert(covalent_r > epsilon_fl);
	const fl x = r / covalent_r;
	if(x > 2) return 0;
	return 1-x*x/4;
}

fl model::clash_penalty_aux(const interacting_pairs& pairs) const {
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));
		const fl covalent_r = atoms[ip.a].covalent_radius() + atoms[ip.b].covalent_radius();
		e += pairwise_clash_penalty(r, covalent_r);
	}
	return e;
}

fl model::clash_penalty() const {
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += clash_penalty_aux(ligands[i].pairs);
	e += clash_penalty_aux(other_pairs);
	return e;}

void model::build_ar_ring_info() {

	this->rec_ar_ring_info.clear();
	DetectAromaticCycles(this->receptor_residues, this->rec_ar_ring_info, false);

	/*VINA_FOR_IN(i,this->rec_ar_ring_info){
		aptrv& ring_atoms = this->rec_ar_ring_info[i].all_atom_ptrs;
		if (ring_atoms[0]->resname.find("TRP") != std::string::npos){
			VINA_FOR_IN(j, this->rec_ar_ring_info[i].fused_edges){
				sz& ei = this->rec_ar_ring_info[i].fused_edges[j].first;
				sz& ej = this->rec_ar_ring_info[i].fused_edges[j].second;
				std::cout << "Access fused edge: " << this->grid_atoms[ei].atomname << " and " << this->grid_atoms[ej].atomname << std::endl;
			}
		}
	}
	std::exit(1);*/

	this->lig_ar_ring_info.clear();
	DetectAromaticCycles(this->ligand_residues, this->lig_ar_ring_info, true);


	VINA_FOR_IN(i, this->lig_ar_ring_info){
		std::cout << "\nLigand aromatic ring " << i+1  << ":";
		aptrv& ring_atoms = this->lig_ar_ring_info[i].all_atom_ptrs;
		VINA_FOR_IN(j, ring_atoms){
			std::cout << ring_atoms[j]->atomname << ",";
		}
		std::cout << "\n";

		std::cout << "Aromatic carbons: ";
		aptrv& ar_carbons = this->lig_ar_ring_info[i].aromatic_carbon_ptrs;
		VINA_FOR_IN(j, ar_carbons){
			std::cout << ar_carbons[j]->atomname << ",";
		}
		std::cout << "\n";

		std::cout << "Heteroatoms: ";
		aptrv& heteros = this->lig_ar_ring_info[i].heteroatom_ptrs;
		VINA_FOR_IN(j, heteros){
			std::cout << heteros[j]->atomname << ",";
		}
		std::cout << "\n\n";
	}

	this->DetectAliphaticCarbons();
	
	this->num_r_ali_carbs = this->rec_ali_carb_info.size();
	this->num_l_ali_carbs = this->lig_ali_carb_info.size();
	this->num_r_aro_rings = this->rec_ar_ring_info.size();
	this->num_l_aro_rings = this->lig_ar_ring_info.size();

}

std::vector<aptrv> model::remove_redundant_cycles(std::vector<aptrv>& cycles, std::map<sz, szv>& subcycles, szv& returned_cycle_indices){
	//This function checks whether a cycle is a subset of another cycles. The TRP side chain for example, benzene is a subset of indole. 
	std::vector<aptrv> cycles_non_redundant;
	VINA_FOR_IN(i, cycles){
		aptrv& this_cycle = cycles[i];
		bool is_subset_of_another_cycle = false;

		VINA_FOR_IN(j, cycles){
			if (j != i){
				aptrv& another_cycle = cycles[j];
				bool is_subset_of_this_other_cycle = true;

				VINA_FOR_IN(k, this_cycle){
					atom_vc* this_cycle_atom = this_cycle[k];
					if (std::find(another_cycle.begin(), another_cycle.end(), this_cycle_atom) == another_cycle.end()){
						is_subset_of_this_other_cycle = false;
						break;
					}
				}

				if (is_subset_of_this_other_cycle){
					is_subset_of_another_cycle = true;
					subcycles[j].push_back(i);
					break;
				}
				
			}
		}

		if (!is_subset_of_another_cycle){
			cycles_non_redundant.push_back(this_cycle);
			returned_cycle_indices.push_back(i);
		}
	}

	//The above works if there is at least one ring that's different from any other ring. 
	//If all rings are the same, the non redundant vector is empty. If the input cycle vector isn't empty, then all rings are identical.
	//In this case, add the first input cycle to the non-redundant vector. 
	if (!cycles.empty() && cycles_non_redundant.empty()){
		cycles_non_redundant.push_back(cycles[0]);
	}

	return cycles_non_redundant;
}

std::vector<aptrv> model::remove_duplicate_cycles(std::vector<aptrv>& cycles){
	//The same ring may be detected many times, for example, once clockwise, and once again counter-clockwise.
	//With fused cycles, it's more complicated than that. 
	std::vector<aptrv> cycles_non_duplicate;
	VINA_FOR_IN(i, cycles){
		aptrv& this_cycle = cycles[i];
		bool is_identical_to_another_cycle = false;

		for (unsigned int j = i+1; j < cycles.size(); j++){
			aptrv& another_cycle = cycles[j];
			if (this_cycle.size() == another_cycle.size()){
				bool is_identical_to_this_other_cycle = true;
				//As long as there is at least one different atom, it is not a subset of this existing cycle.

				VINA_FOR_IN(k, this_cycle){
					atom_vc* this_cycle_atom = this_cycle[k];
					if (std::find(another_cycle.begin(), another_cycle.end(), this_cycle_atom) == another_cycle.end()){
						is_identical_to_this_other_cycle = false;
						break;
					}
				}

				if (is_identical_to_this_other_cycle){
					is_identical_to_another_cycle = true;
					break;
				}
			}
		}

		if (!is_identical_to_another_cycle){
			cycles_non_duplicate.push_back(this_cycle);
		}

	}
	return cycles_non_duplicate;
}

void model::DFSVisit(aptrv& atoms_to_search, std::map<sz, bool>& atom_visited_map, std::map<sz, int>& atom_parent_map, std::map<sz, atom_vc*>& index_atom_map, sz atom_index, int& counter, aptrv& current_path, std::vector<aptrv>& detected_cycles){
	atom_vc* atom = index_atom_map[atom_index];
    	atom_visited_map[atom_index] = true;
	current_path.push_back(atom);
	int parent_index = atom_parent_map[atom_index];

	//Determine if downstream branches exist.
	int num_neighbors_downstream = 0;
	VINA_FOR_IN(i, atom->bonds){
		const bond_vc& b = atom->bonds[i];
                atom_vc& neighbor = get_atom(b.connected_atom_index);
                atom_vc* n_ptr = &(neighbor);

		int n_index = -1; //-1 means that neighbor atom is in another residue. Don't go beyond current residue. 
                for (std::map<sz, atom_vc*>::iterator mapit = index_atom_map.begin(); mapit != index_atom_map.end(); mapit++){
                        if (mapit->second == n_ptr){
                                n_index = mapit->first;
                                break;
                        }
                }

		if (parent_index == -1){
			num_neighbors_downstream++;
		} 
		else if (n_index != parent_index){
			num_neighbors_downstream++;
		}
	}

	VINA_FOR_IN(i, atom->bonds){
		const bond_vc& b = atom->bonds[i];
		atom_vc& neighbor = get_atom(b.connected_atom_index);
		atom_vc* n_ptr = &(neighbor);

    		if (has(atoms_to_search, n_ptr)){
			int n_index = -1;
			for (std::map<sz, atom_vc*>::iterator mapit = index_atom_map.begin(); mapit != index_atom_map.end(); mapit++){
				if (mapit->second == n_ptr){
					n_index = mapit->first;
					break;
				}
			}
			VINA_CHECK(n_index != -1); //If has() returns true, this should never happen.

			//Check if neighbor is in the atoms to be searched, only continue searching if this is the case.
			//If ligand is covalently bonded to receptor, use this to block searching to the other molecule. 
			//If not visited
               		if(!atom_visited_map[n_index]){
               			atom_parent_map[n_index] = atom_index;
				//If there are branches downstream, duplicate tracker containers, and start recursion.
				if (num_neighbors_downstream > 1){
					std::map<sz, bool> child_atom_visited_map = atom_visited_map;
					std::map<sz, int> child_atom_parent_map = atom_parent_map;
					aptrv child_current_path = current_path;
              				DFSVisit(atoms_to_search, child_atom_visited_map, child_atom_parent_map, index_atom_map, n_index, counter, child_current_path, detected_cycles);
				}
				//If there is no branch, keep using containers parsed from upstream.
				else{
              				DFSVisit(atoms_to_search, atom_visited_map, atom_parent_map, index_atom_map, n_index, counter, current_path, detected_cycles);
				}
                	}
                	else{
				if (parent_index == -1) continue; //-1 means no parent.

                    		atom_vc* parent = index_atom_map[parent_index];
				//making sure we are not tracking back to the previous atom which is the parent of neigbor (current atom)
                    		if(n_ptr != parent){
					//VINA_CHECK(parent->atomname != "fake");
                        		counter++;

					//I'll find cycle a(neighbor)->b(atom) then b->a. Eliminate b->a by requiring that a is ahead of b on the path. 
					//In the case of b->a, b wil have been removed from current_path (visited)
					aptrv::iterator cycle_begin = std::find(current_path.begin(), current_path.end(), n_ptr);
					aptrv::iterator cycle_last = std::find(current_path.begin(), current_path.end(), atom);

					bool valid_cycle = true;
					//If a neighbor is visited, but already dequeued from current path, then I have b->a. False positive. 
					if (cycle_begin == current_path.end()){
						valid_cycle = false;
					}

					if (valid_cycle){
						bool reading_cycle = false;
						int num_ring_atoms = 0;
						aptrv new_cycle;

						VINA_FOR_IN(i, current_path){
							atom_vc* path_atom = current_path[i];
							if (path_atom == n_ptr){
								reading_cycle = true;
							}
	
							if(reading_cycle){
								new_cycle.push_back(path_atom);
								num_ring_atoms++;
							}
							if (path_atom == atom){
								reading_cycle = false;
								break;
							}

						}
						//std::cout << "At node " << atom->atomname << " found cycle from " << n_ptr->resname << " " << n_ptr->resnum<< "-" << n_ptr->atomname << " to " << atom->resname << "-" << atom->atomname << " Num ring atoms: " << num_ring_atoms << std::endl;

						detected_cycles.push_back(new_cycle);
					}
						
                    		}
                	}

        	}

	}

	current_path.pop_back();
}

std::vector<aptrv> model::DetectCyclesByDFS(std::map<sz, atom_vc*>& index_atom_map){
	std::vector<aptrv> cycles;
       	if(index_atom_map.empty()) return cycles;

	int counter = 0;
	//std::map<atom_vc*, bool> atom_visited_map;
	//std::map<atom_vc*, atom_vc*> atom_parent_map;
	std::map<sz, bool> atom_visited_map;
	std::map<sz, int> atom_parent_map;
	aptrv current_path; 
	aptrv atoms_to_search;

	atom_vc fake_atom;
	fake_atom.atomname = "fake";

	for(std::map<sz, atom_vc*>::iterator mapit = index_atom_map.begin(); mapit != index_atom_map.end(); mapit++){
		sz atom_index = mapit->first;
		atom_vc* aptr = mapit->second;
 		atom_visited_map[atom_index] = false;
       		atom_parent_map[atom_index] = -1; //Use -1 for fake atom.
		atoms_to_search.push_back(aptr);
	}       

	std::map<sz, atom_vc*>::iterator map_begin = index_atom_map.begin();
	sz first_atom_index = map_begin->first;
	
 	DFSVisit(atoms_to_search, atom_visited_map, atom_parent_map, index_atom_map, first_atom_index, counter, current_path, cycles);

    	return cycles;
		
}

void model::DetectAromaticCycles(resv& residues, ring_info& ring_attributes, bool is_ligand){
	aptrv temp_atom_ptrs = (is_ligand) ? this->atom_ptrs : this->grid_atom_ptrs;
	vptrv temp_coord_ptrs;
	if (is_ligand){
		VINA_FOR_IN(i, this->coords){
			temp_coord_ptrs.push_back(&(this->coords[i]));
		}
	}
	else{
		VINA_FOR_IN(i, this->grid_atoms){
			temp_coord_ptrs.push_back(&(this->grid_atoms[i].coords));
		}
	}

	VINA_FOR_IN(i, residues){
		std::map<sz, atom_vc*>& index_atom_map = residues[i].index_atom_map;
		//std::map<sz, vec*>& index_coord_map = residues[i].index_coord_map;
		std::vector<std::pair<aptrv, std::vector<aptrv> > > ar_cycles_this_residue = DetectAromaticCycles(index_atom_map);

		VINA_FOR_IN(j, ar_cycles_this_residue){
			aptrv& this_cycle = ar_cycles_this_residue[j].first;
			std::vector<aptrv>& subcycles = ar_cycles_this_residue[j].second;

			/*std::vector<szv> subcycle_atom_indices;
			VINA_FOR_IN(k, subcycles){
				aptrv& this_subcycle = subcycles[k];
				szv this_subcycle_atom_indices;

				VINA_FOR_IN(l, this_subcycle){
					atom_vc* cycle_atom = this_subcycle[l];
					int index = -1;
                                	for (std::map<sz, atom_vc*>::iterator mapit = index_atom_map.begin(); mapit != index_atom_map.end(); mapit++){
                                        	if (mapit->second == cycle_atom){
                                                	index = mapit->first;
                                                	break;
                                        	}
                                	}
                                	VINA_CHECK(index != -1);
                                	this_subcycle_atom_indices.push_back(sz(index));
				}
				subcycle_atom_indices.push_back(this_subcycle_atom_indices);
			}*/

			szv atom_indices;
			aptrv ring_hydrogens;
			szv ring_h_indices;
			vptrv ring_h_coords;
			szv ring_h_heavy_neighbor_indices;
			vptrv ring_h_heavy_neighbor_coords;

			VINA_FOR_IN(k, this_cycle){
				atom_vc* cycle_atom = this_cycle[k];
				int index = -1;
				for (std::map<sz, atom_vc*>::iterator mapit = index_atom_map.begin(); mapit != index_atom_map.end(); mapit++){
					if (mapit->second == cycle_atom){
						index = mapit->first;
						break;
					}
				}
				VINA_CHECK(index != -1);
				atom_indices.push_back(sz(index));

				VINA_FOR_IN(l, cycle_atom->bonds){
					atom_vc& neighbor = get_atom(cycle_atom->bonds[l].connected_atom_index);
					if (neighbor.ad == AD_TYPE_H){
						ring_hydrogens.push_back(&(neighbor));
						//int r_h_index = std::distance(temp_atom_ptrs.begin(), std::find(temp_atom_ptrs.begin(), temp_atom_ptrs.end(), &(neighbor)));
						int r_h_index = cycle_atom->bonds[l].connected_atom_index.i;
						ring_h_indices.push_back(r_h_index);
						ring_h_coords.push_back(temp_coord_ptrs[r_h_index]);
						ring_h_heavy_neighbor_indices.push_back(index);
						ring_h_heavy_neighbor_coords.push_back(temp_coord_ptrs[index]);
					}
				}
				
			}

			ring_attributes.push_back(ring_attribute(this, this_cycle, atom_indices, temp_coord_ptrs, is_ligand, ring_hydrogens, ring_h_indices, ring_h_coords, ring_h_heavy_neighbor_indices, ring_h_heavy_neighbor_coords, &(adjusted_coord_ptrs), false, temp_atom_ptrs, subcycles));
		}
	}
	return;
}

std::vector<std::pair<aptrv, std::vector<aptrv> > > model::DetectAromaticCycles(std::map<sz, atom_vc*>& index_atom_map){

	std::vector<aptrv> cycles = DetectCyclesByDFS(index_atom_map);
	std::vector<aptrv> reasonable_cycles;
        int max_num_atoms_allowed = 100;

        VINA_FOR_IN(i, cycles){
                aptrv& this_cycle = cycles[i];

                if (this_cycle.size() > max_num_atoms_allowed){
                        std::cout << "Warning, encountered a cycle > " << max_num_atoms_allowed << " atoms." << std::endl;
                        std::cout << "This cycle will be filtered out. Please check for bad input structures." << std::endl;
                        continue;
                }
                reasonable_cycles.push_back(this_cycle);
        }

        std::vector<aptrv> reasonable_cycles_non_duplicate =  remove_duplicate_cycles(reasonable_cycles);


	//Filter for rings that conform to aromaticity
	std::vector<aptrv> reasonable_aromatic_cycles;
	VINA_FOR_IN(i,  reasonable_cycles_non_duplicate){
		aptrv& this_cycle = reasonable_cycles_non_duplicate[i];
		if (is_aromatic_cycle(this_cycle)){
			reasonable_aromatic_cycles.push_back(this_cycle);
		}
	}

	std::map<sz, szv> subcycles;
	szv returned_cycle_indices;
        std::vector<aptrv> reasonable_cycles_non_redundant =  remove_redundant_cycles(reasonable_aromatic_cycles, subcycles, returned_cycle_indices);

	std::vector<std::pair<aptrv, std::vector<aptrv> > > cycle_subcycles;
	VINA_FOR_IN(i, returned_cycle_indices){
		sz& index = returned_cycle_indices[i];

		aptrv& this_cycle = reasonable_aromatic_cycles[index];
		szv& subcycle_indices = subcycles[index];
		std::vector<aptrv> subcycles;

		VINA_FOR_IN(j, subcycle_indices){
			aptrv& this_subcycle = reasonable_aromatic_cycles[subcycle_indices[j]];
			subcycles.push_back(this_subcycle);
		}
		
		/*std::string resname = this_cycle[0]->resname;
		if (resname.find("TRP") != std::string::npos){
		std::cout << "This " << resname << " cycle contains " << this_cycle.size() << " atoms " << " divided into " << subcycles.size() << " subcycles." << std::endl;
		VINA_FOR_IN(j, subcycles){
			std::cout << "This subcycle contains " << subcycles[j].size() <<  " atoms." << std::endl;
		}
		}*/
		
		cycle_subcycles.push_back(std::make_pair(this_cycle, subcycles));
	}

    	//return reasonable_cycles_non_redundant;
    	return cycle_subcycles;
}

fl model::eval_chpi_enthalpy_c(fl r){
	if (r >= chpi_dcut) return 0;
	return g_coeff * std::exp(-sqr(r-chpi_miu)/g_denominator);
}

fl model::eval_chpi_entropy(fl horizontal_offset){
	if (horizontal_offset >= chpi_ho_max) return 0;
	fl ho_effective = std::max(horizontal_offset, chpi_ho_epsilon);
	//return (-0.594 * std::log(ho_effective) + 1.0609);  //6.0A, using intra-protein interactions. 
	//return (-0.594 * std::log(ho_effective) + 1.0066);  //5.5A, using intra-protein interactions. 
	//return (-0.595 * std::log(ho_effective) + 0.9460);  //5.0A, using intra-protein interactions. 
	return (-0.595 * std::log(ho_effective) + 0.9087);  //4.7A, using intra-protein interactions.
	//return (-0.597 * std::log(ho_effective) + 0.8846);  //4.5A, using intra-protein interactions. 
	//return (-0.599 * std::log(ho_effective) + 0.8172);  //4.0A, using intra-protein interactions. 
	//return (-0.603 * std::log(ho_effective) + 0.7505);  3.5A
	//return (-0.605 * std::log(ho_effective) + 0.6564);  //3.0A
	//return (-0.607 * std::log(ho_effective) + 0.5472);  //2.5A
	//return (-0.611 * std::log(ho_effective) + 0.4192);  //2.0A
	//return (-0.613 * std::log(ho_effective) + 0.2108);  //1.4A, using intra-protein interactions. 
}

void model::eval_chpi_c_ring(aliphatic_carbon_attribute& c, ring_attribute& r, pr& dH_minusTdS, bool score_in_place){
	//std::cout << "Eval C " << std::endl;
	fl dH = 0.00;
	bool& ligand_aliphatic = c.is_ligand;

	sz& carbon_index = c.carbon_atom_index;
	vec* c_coord = c.c_coord;
	vec& centroid = r.centroid;

	vec centroid_c; centroid_c = *c_coord - centroid;
	fl ccl = magnitude(centroid_c); //C-centroid length
        if (ccl >= (chpi_dcut + r.max_atom_centroid_dist)) return;

	bool this_ring_interacting = false;
	vptrv& ar_coord_ptrs = r.aromatic_carbon_coord_ptrs;
	szv& ar_atom_indices = r.aromatic_carbon_indices;
	szv interacting_indices;
	//vec total_deriv(0,0,0);

	VINA_FOR(i, r.num_aromatic_carbons){
		//fl(&rd)[3] = ar_coord_ptrs[i]->data; //rd = ring atom coord data
		sz& ar_atom_index = ar_atom_indices[i];
		vec* rc = ar_coord_ptrs[i];
		//A vector of ligand coord - receptor coord;
		vec ra_to_la = (ligand_aliphatic) ? *c_coord - *rc : *rc - *c_coord;
		fl r = magnitude(ra_to_la);
		if (r >= chpi_dcut) continue;

                this_ring_interacting = true;
		interacting_indices.push_back(ar_atom_index);

		pr e_dor = (score_in_place) ? this->eval_chpi_enthalpy_c_fast(r): this->eval_chpi_enthalpy_c_deriv(r);

		dH += e_dor.first;
		//total_deriv += this_pair_deriv;
		vec this_pair_deriv = ra_to_la;  this_pair_deriv *= e_dor.second;

		if (ligand_aliphatic){
			this->minus_forces[carbon_index] += this_pair_deriv;
                }
		else{
			this->minus_forces[ar_atom_index] += this_pair_deriv; 
		}
	}

	if (!this_ring_interacting) return;
	dH_minusTdS.first += dH;

	//return; //no entropy

	//if (c.num_h_neighbors == 0) return; //Skipping entropy for C without hydrogen, mostly peptide bond carbonyl. 

	fl(&rnd)[3] = r.normal.data;
	fl ccndp = centroid_c[0]*rnd[0] + centroid_c[1]*rnd[1] + centroid_c[2]*rnd[2];
        fl vo = std::abs(ccndp);
        fl ho = std::sqrt(ccl*ccl - vo*vo);

	bool entropy_in_cutoff = (vo < chpi_vo_max_c && ho < chpi_ho_max);
	if (!entropy_in_cutoff) return;

	pr entropy_e_dor = (score_in_place) ? this->eval_chpi_entropy_fast(ho) : this->eval_chpi_entropy_deriv(ho);
	fl minusTdS = -this->weight_chpi * entropy_e_dor.first; fl tds_dor = -this->weight_chpi * entropy_e_dor.second; 
	dH_minusTdS.second += minusTdS;

	vec vertical(r.normal); vertical *= vo;
	if (ccndp < 0) vertical *= -1; //If carbon and ring normal on oppposite of ring, flip vertical to be on the same side. 
	vec p = centroid + vertical;
	vec entropy_deriv = (*c_coord - p);
	entropy_deriv *= tds_dor;

	if (ligand_aliphatic){
		//this->minus_forces[carbon_index] += total_deriv + entropy_deriv;
		this->minus_forces[carbon_index] += entropy_deriv;
		return;
	}

	entropy_deriv *= -1; //Reverse the direction of derivative so that it is still receptor->ligand
	vec entropy_deriv_per_ring_atom = entropy_deriv;
	//entropy_deriv_per_ring_atom /= (fl) interacting_indices.size(); 
	entropy_deriv_per_ring_atom /= (fl) ar_atom_indices.size(); 
	if (magnitude(entropy_deriv_per_ring_atom) < epsilon_fl) return;
	//VINA_FOR_IN(i, interacting_indices){
	VINA_FOR_IN(i, ar_atom_indices){
		sz& ar_atom_index = ar_atom_indices[i];
		//Apply total entropy evenly on each interacting ring atom
		this->minus_forces[ar_atom_index] += entropy_deriv_per_ring_atom; 
	}
	return;

	//return dH; //Testing without entropy
	//return (dH - eval_chpi_entropy(ho));
}

void model::eval_chpi_h_ring(aliphatic_carbon_attribute& c, ring_attribute& r, pr& dH_minusTdS, bool score_in_place){
	if (c.num_h_neighbors == 0) return;

	bool& ligand_aliphatic = c.is_ligand;
	sz& carbon_index = c.carbon_atom_index;

	szv& h_neighbor_indices = c.h_neighbor_indices;
	vptrv& h_coords = c.h_coords;
	fl dH = 0.00; 

	vec& centroid = r.centroid;
	vec& normal = r.normal;
       	vec* c_coord = c.c_coord;
	vec cc = *c_coord - centroid;
	fl ccl = magnitude(cc);
        fl ccdp = cc * normal;
        fl ccvo = std::abs(ccdp); //C-centroid vertical offset. All normal vectors are pre-normalized, length = 1.

	vptrv& ar_coord_ptrs = r.aromatic_carbon_coord_ptrs;
	szv& ar_atom_indices = r.aromatic_carbon_indices;
        bool this_ring_interacting = false;

	szv& ring_h_indices = r.ring_h_indices;
        vptrv& ring_h_coords = r.ring_h_coords;
        szv& ring_h_heavy_neighbor_indices = r.ring_h_heavy_neighbor_indices;
        vptrv& ring_h_heavy_neighbor_coords = r.ring_h_heavy_neighbor_coords;

	aptrv& all_atom_ptrs = r.all_atom_ptrs;
        vptrv& ring_atoms_coords = r.all_atom_coord_ptrs;
        szv& ring_atom_indices = r.all_atom_indices;

	std::map<sz, fl> h_centroid_dist;
	szv ip_hs;

	//std::vector<flv> h_ring_dists(c.num_h_neighbors, flv(r.num_ring_atoms, max_fl));

	VINA_FOR(i, c.num_h_neighbors){
                sz h_index = h_neighbor_indices[i];
		//if (ligand_aliphatic && has(this->adjusted_atom_indices, h_index)) continue;

		//vec* hv = (ligand_aliphatic) ? this->adjusted_coord_ptrs[h_index] : h_coords[i];
		vec* hv = h_coords[i];
                vec hc = (*hv) - centroid;
                fl hcl = magnitude(hc);

		//If H-centroid distance is greater than (hydrogen cutoff + ring effective radius), skip;
		if (hcl >= (chpi_hcut + r.max_atom_centroid_dist)) continue;

		fl hcdp = hc * normal;	
		//Step 1: check if the C and H atom is on the same face of ring. If not, their cosines have opposite sign.
		//if (ccdp * hcdp < 0) continue;

		fl hcvo = std::abs(hcdp); //H-centroid vertical offset
		//Step 2: Make sure CH bond doesn't point away from the ring, if so: H vertical offset >= C vertical offset
		//if (hcvo >= ccvo) continue;

		bool induction_polarization = (ccdp * hcdp > 0 && hcvo < ccvo);
                //if (!induction_polarization) continue;
		//continue;
		vec ch_bond = (*hv) - (*c_coord);

                fl chndp = ch_bond * normal; //CH-Normal dot product
		fl chn_cosine = chndp / magnitude(ch_bond);
		fl chn_angle = std::acos(chn_cosine);
		fl chn_cosine_abs = std::abs(chn_cosine);

		//VINA_FOR(j, r.num_aromatic_carbons)
                	//vec* rc = ar_coord_ptrs[j];
			//sz& ar_atom_index = ar_atom_indices[j];
		VINA_FOR(j, r.num_ring_atoms){
                	vec* rc = ring_atoms_coords[j];
			sz& ring_atom_index = ring_atom_indices[j];

			vec ra_to_la = (ligand_aliphatic) ? *hv - *rc : *rc - *hv; //This is the H-ring atom vector, but will be applied on C.
			//vec c_ra_to_la = (ligand_aliphatic) ? *c_coord - *rc : *rc - *c_coord; //Use C-ring instead?

			fl r = magnitude(ra_to_la);
			//h_ring_dists[i][j] = 1 / sqr(r);

			if (r >= chpi_hcut) continue;

			if (h_centroid_dist.find(i) == h_centroid_dist.end()){
				h_centroid_dist[i] = hcl;
				if (induction_polarization){
					ip_hs.push_back(i);
				}
			}
			
        	}

	}

	if (h_centroid_dist.empty()) return;

	int closest_h_i = this->choose_interacting_h(h_centroid_dist, ip_hs);
	if (closest_h_i == -1) return;
	vec* h_closest = h_coords[closest_h_i];
	vec* closest_centroid = this->choose_closest_centroid(h_closest, centroid, r.subcycle_centroids);
	sz h_closest_index = h_neighbor_indices[closest_h_i];

	vec hcc = *h_closest - *closest_centroid;
	fl hcc_dp = hcc * normal;
        fl hcc_vo = std::abs(hcc_dp); //H-centroid vertical offset
	fl hcc_ho = std::sqrt(sqr(magnitude(hcc))-sqr(hcc_vo));

	//VINA_FOR(j, r.num_aromatic_carbons)
		//vec* rc = ar_coord_ptrs[j];
                //sz& ar_atom_index = ar_atom_indices[j];
	VINA_FOR(j, r.num_ring_atoms){
                vec* rc = ring_atoms_coords[j];
                sz& ring_atom_index = ring_atom_indices[j];
		vec ra_to_la = (ligand_aliphatic) ? *h_closest - *rc : *rc - *h_closest;

		fl r = magnitude(ra_to_la);
                if (r >= chpi_hcut) continue;
		pr this_pair_e_dor(0,0);
			
		fl rep_distance = r - (xs_vdw_radii[all_atom_ptrs[j]->xs] + 1.1);
                if (rep_distance < 0){
                        fl repulsion_e = sqr(rep_distance); fl repulsion_deriv = 2.0 * rep_distance;
                        this_pair_e_dor.first -= repulsion_e; this_pair_e_dor.second -= repulsion_deriv;
                        dH -= repulsion_e;
                }

		fl e = eval_chpi_enthalpy_h(r); fl deriv = (chpi_miu_h - r) * inv_ssqr * e;
                fl e2 = eval_chpi_enthalpy_h2(r); fl deriv2 = (chpi_miu_h2 -r) * inv_ssqr2 * e2;

		//pr f_df = this->calc_horizontal_factor(hcc_ho);

		fl e_total = e + e2;
		fl e_deriv = deriv + deriv2;

		//std::cout << "Raw and factored deriv " << hcc_ho << "," << e_deriv << "," << deriv_total << std::endl;
                this_pair_e_dor.first += e_total; this_pair_e_dor.second += e_deriv;

		//fl dor = (r > 4) ? (this->weight_chpi * this_pair_e_dor.second) : (this->weight_chpi * this_pair_e_dor.second / r);
		fl dor = this->weight_chpi * this_pair_e_dor.second / r;
		//fl dor = this->weight_chpi * this_pair_e_dor.second;
		dH += e_total;
		vec this_pair_deriv = ra_to_la; this_pair_deriv *= dor;

		if (ligand_aliphatic){
			this->minus_forces[h_closest_index] += this_pair_deriv;
		}
		else{
			//this->minus_forces[ar_atom_index] += this_pair_deriv;
			this->minus_forces[ring_atom_index] += this_pair_deriv;
		}
	}

        dH_minusTdS.first += dH;
	return; //No entropy until I get at least 33% in performance

	//this->eval_chpi_entropy_set_force2(h_closest, closest_h_i, h_closest_index, h_ring_dists, r, dH_minusTdS, ligand_aliphatic);
	//this->eval_chpi_entropy_set_force(h_closest, closest_h_i, h_closest_index, h_ring_dists, r, dH_minusTdS, ligand_aliphatic);
	this->eval_chpi_entropy_set_force_old(h_closest, h_closest_index, r, dH_minusTdS, ligand_aliphatic);
	return;
}
void model::undo_attraction_scoring(fl& e_chpi, bool score_in_place){
	fl e_gau1 = 0, e_gau2 = 0, e_phobic = 0;
	VINA_FOR(i, this->num_l_aro_rings){
                ring_attribute& r_inf = this->lig_ar_ring_info[i];
                VINA_FOR(j, this->num_r_ali_carbs){
                        this->undo_attraction_scoring_each_pair(this->rec_ali_carb_info[j], r_inf, e_gau1, e_gau2, e_phobic, score_in_place);
                }
        }

        VINA_FOR(i, this->num_r_aro_rings){
                ring_attribute& r_inf = this->rec_ar_ring_info[i];
                VINA_FOR(j, this->num_l_ali_carbs){
                        this->undo_attraction_scoring_each_pair(this->lig_ali_carb_info[j], r_inf, e_gau1, e_gau2, e_phobic, score_in_place);
                }
        }

	if (score_in_place){
		e_gau1 *= this->weight_gauss1;
		e_gau2 *= this->weight_gauss2;
		e_phobic *= this->weight_hydrophobic;
	}
	//std::cout << "g1, g2, e phobic " << e_gau1 << "," << e_gau2 << "," << e_phobic << std::endl;
	e_chpi -= (e_gau1 + e_gau2 + e_phobic);

        return;	
}

void model::undo_attraction_scoring_each_pair(aliphatic_carbon_attribute& c, ring_attribute& r, fl& e_gau1, fl& e_gau2, fl& e_phobic, bool score_in_place){
	if (c.num_h_neighbors == 0) return;
        bool& ligand_aliphatic = c.is_ligand;
        sz& carbon_index = c.carbon_atom_index;

        vec* c_coord = c.c_coord;
	//std::cout << "C coord " << (*c_coord)[0] << "," << (*c_coord)[1] << "," << (*c_coord)[2] << "\n";
	//vec& centroid = r.centroid;
        //vec& normal = r.normal;
        //vec cc = *c_coord - centroid;
        //fl ccl = magnitude(cc);
        //fl ccdp = cc * normal;
        //fl ccvo = std::abs(ccdp); //C-centroid vertical offset. All normal vectors are pre-normalized, length = 1.

        vptrv& ar_coord_ptrs = r.aromatic_carbon_coord_ptrs;
        szv& ar_atom_indices = r.aromatic_carbon_indices;

        aptrv& all_atom_ptrs = r.all_atom_ptrs;
        vptrv& ring_atoms_coords = r.all_atom_coord_ptrs;
        szv& ring_atom_indices = r.all_atom_indices;

	VINA_FOR(j, r.num_ring_atoms){
                vec* rc = ring_atoms_coords[j];
                sz& ring_atom_index = ring_atom_indices[j];
		vec ra_to_la = (ligand_aliphatic) ? *c_coord - *rc : *rc - *c_coord;
		fl r = magnitude(ra_to_la);
                if (r >= chpi_dcut) continue;

		fl r_optimal = r - (xs_vdw_radii[all_atom_ptrs[j]->xs] + 1.9); 
		//std::cout << "R optimal " << r_optimal << "\n";
		//std::cout << "Other atom xs: " << all_atom_ptrs[j]->xs << "\n";
		//std::cout << "Other atom name: " << all_atom_ptrs[j]->atomname << "\n";
		//std::cout << "r is: " << r << std::endl;
		//std::cout << "Other atom xs vdw radii: " << xs_vdw_radii[all_atom_ptrs[j]->xs] << "\n";
		//fl r2 = sqr(r_optimal);
		fl r2 = sqr(r);

		if (!score_in_place){
			pr gau1_e_dor = this->eval_gauss1_deriv(r2);
			pr gau2_e_dor = this->eval_gauss2_deriv(r2);
			pr phobic_e_dor = this->eval_phobic_deriv(r2);
			e_gau1 += gau1_e_dor.first;
			e_gau2 += gau2_e_dor.first;
			e_phobic += phobic_e_dor.first;

			vec this_pair_deriv = ra_to_la;
			this_pair_deriv *= (gau1_e_dor.second + gau2_e_dor.second + phobic_e_dor.second);
                	if (ligand_aliphatic){
                        	this->minus_forces[carbon_index] -= this_pair_deriv;
                	}
                	else{
                        	this->minus_forces[ring_atom_index] -= this_pair_deriv;
                	}
		}
		else{
			/*pr gau1_e_fast = this->eval_gauss1_fast(r2);
			pr gau2_e_fast = this->eval_gauss2_fast(r2);
			pr phobic_e_fast = this->eval_phobic_fast(r2);
			e_gau1 += gau1_e_fast.first;
			e_gau2 += gau2_e_fast.first;
			e_phobic += phobic_e_fast.first;*/

			e_gau1 += gaussian(r_optimal, 0.5);
			e_gau2 += gaussian(r_optimal - 3, 2.0);
			e_phobic += slope_step(1.5, 0.5, r_optimal);
		}
		//e_gau1 += gaussian(r_optimal, 0.5);
		//e_gau2 += gaussian(r_optimal - 3, 2.0);
		//e_phobic += slope_step(1.5, 0.5, r_optimal);

		//How to get the derivative back?
		/*vec this_pair_deriv = ra_to_la;

		this_pair_deriv *= e_dor.second;
                if (ligand_aliphatic){
                        this->minus_forces[carbon_index] -= this_pair_deriv;
                }
                else{
                        this->minus_forces[ring_atom_index] -= this_pair_deriv;
                }*/
        }
	return;
}

fl model::get_chpi_scaling_factor(sz num_h_neighbors, sz num_ip_h){
	VINA_CHECK(num_ip_h <= num_h_neighbors);
	if (num_h_neighbors == 0) return 0; //Shouldn't happen
	if (num_h_neighbors == 1){
		if (num_ip_h == 0) return 0.008;
		if (num_ip_h == 1) return 1.000;
	}
	else if (num_h_neighbors == 2){
		if (num_ip_h == 0) return 0.009;
		if (num_ip_h == 1) return 1.000;
		if (num_ip_h == 2) return 0.120;
	}
	else if (num_h_neighbors == 3){
		if (num_ip_h == 0) return 0.006;
		if (num_ip_h == 1) return 1.000;
		else if (num_ip_h == 2) return 0.201;
		else if (num_ip_h == 3) return 0.007; //Necessary?
	}
	else if (num_h_neighbors == 4){
		return 1; //Shouldn't happen, unless it's methane, only good for score-in-place.
	}
	VINA_CHECK(false);
	return 0; //Shouldn't happen, just to make g++ happy. 
}

/*fl model::get_chpi_scaling_factor(sz num_h_neighbors, sz num_ip_h){
        VINA_CHECK(num_ip_h <= num_h_neighbors);
        if (num_h_neighbors == 0) return 0; //Shouldn't happen
        if (num_h_neighbors == 1){
                if (num_ip_h == 0) return 0.000;
                if (num_ip_h == 1) return 1.000;
        }
        else if (num_h_neighbors == 2){
                if (num_ip_h == 0) return 0.000;
                if (num_ip_h == 1) return 1.000;
                if (num_ip_h == 2) return 0.000;
        }
        else if (num_h_neighbors == 3){
                if (num_ip_h == 0) return 0.000;
                if (num_ip_h == 1) return 1.000;
                else if (num_ip_h == 2) return 0.000;
                else if (num_ip_h == 3) return 0.000; //Necessary?
        }
        else if (num_h_neighbors == 4){
                return 1; //Shouldn't happen, unless it's methane, only good for score-in-place.
        }
        VINA_CHECK(false);
        return 0; //Shouldn't happen, just to make g++ happy. 
}*/

int model::choose_interacting_h(std::map<sz, fl>& h_centroid_dist, szv& ip_hs){
	bool no_ip_h = ip_hs.empty();
	int closest_h_i = -1; fl min_hc_dist = max_fl;
        for (std::map<sz, fl>::iterator mapit = h_centroid_dist.begin(); mapit != h_centroid_dist.end(); mapit++){
                sz h_i = mapit->first;
                fl hc_dist = mapit->second;

		if (!no_ip_h){
			if (!has(ip_hs, h_i)) continue; //Temporary
		}
                if (hc_dist < min_hc_dist){
                        min_hc_dist = hc_dist;
                        closest_h_i = h_i;
                }
        }
	//VINA_CHECK(closest_h_i != -1);
	return closest_h_i;
}

fl model::eval_chpi(bool score_in_place){
	//return 0; //Skip temporarily

	pr dH_minusTdS(0.00, 0.00);
	VINA_FOR(i, this->num_l_aro_rings){
		this->lig_ar_ring_info[i].ComputeCentroidAndNormal(false);
		//this->lig_ar_ring_info[i].ComputeLargestAtomCentroidDistance(false);
	}

	if (this->chpi_explicit_hydrogen){
		this->eval_chpi_h(dH_minusTdS, score_in_place);
	}
	else{
		this->eval_chpi_c(dH_minusTdS, score_in_place);
	}

	fl e_chpi = this->weight_chpi * (dH_minusTdS.first - dH_minusTdS.second); 
	//this->undo_attraction_scoring(e_chpi, score_in_place);
	return e_chpi;

}

void model::eval_chpi_c(pr& dH_minusTdS, bool score_in_place){
	VINA_FOR(i, this->num_l_aro_rings){
		ring_attribute& r_inf = this->lig_ar_ring_info[i];
		VINA_FOR(j, this->num_r_ali_carbs){
			this->eval_chpi_c_ring(this->rec_ali_carb_info[j], r_inf, dH_minusTdS, score_in_place);
		}
	}

	VINA_FOR(i, this->num_r_aro_rings){
		ring_attribute& r_inf = this->rec_ar_ring_info[i];
		VINA_FOR(j, this->num_l_ali_carbs){
                        this->eval_chpi_c_ring(this->lig_ali_carb_info[j], r_inf, dH_minusTdS, score_in_place);
		}
	}

	return;
}

void model::eval_chpi_h(pr& dH_minusTdS, bool score_in_place){
	sz num_rings_to_save = 1;

        VINA_FOR(i, this->num_r_ali_carbs){
		aliphatic_carbon_attribute& c_attr = this->rec_ali_carb_info[i];
        	VINA_FOR(j, this->num_l_aro_rings){
                	ring_attribute& r_inf = this->lig_ar_ring_info[j];
                        this->eval_chpi_h_ring(c_attr, r_inf, dH_minusTdS, score_in_place);
                }
        }

        VINA_FOR(i, this->num_l_ali_carbs){
		aliphatic_carbon_attribute& c_attr = this->lig_ali_carb_info[i];
                VINA_FOR(j, this->num_r_aro_rings){
                	ring_attribute& r_inf = this->rec_ar_ring_info[j];
                        this->eval_chpi_h_ring(c_attr, r_inf, dH_minusTdS, score_in_place);
                }
        }
	//std::exit(1);

	return;
        //return (lig_rec_chpi + rec_lig_chpi);
}

fl model::eval_chpi_enthalpy_h(fl r){
	if (r >= chpi_hcut) return 0;
        return g_coeff_h * std::exp(-sqr(r-chpi_miu_h)/g_denominator_h);
}

fl model::eval_chpi_enthalpy_h2(fl r){
	//return 0;
        if (r >= chpi_hcut) return 0;
        return g_coeff_h2 * std::exp(-sqr(r-chpi_miu_h2)/g_denominator_h2);
}

void model::build_residues_from_atoms(resv& residues, aptrv& atoms, bool is_ligand){
	VINA_FOR_IN(i, atoms){
                atom_vc* atom = atoms[i];
                std::string& resname = atom->resname;
                std::string& resnum = atom->resnum;
                std::string& chainID = atom->chainID;
                //std::cout << "Target atom name, resnum, and chain ID " << resname << "," << resnum << "," << chainID << std::endl;

                if (residues.empty()){
                        residues.push_back(model_residue(resname, resnum, chainID));
                }

                model_residue& current_residue = residues.back();
                //If neither resnum and chain ID has changed, this atom still belongs to the current residue. 

		bool in_new_residue = true;
		for (int j = residues.size() -1; j >= 0; j--){
			model_residue& this_residue = residues[j];
			if (resnum == this_residue.number && chainID == this_residue.chainID){
				this_residue.index_atom_map[i] = atom;
				this_residue.index_coord_map[i] = (is_ligand) ? &(this->coords[i]) : &(atom->coords);
				in_new_residue = false;
				break;
			}
		}

		if (in_new_residue){
			residues.push_back(model_residue(resname, resnum, chainID));
			residues.back().index_atom_map[i] = atom;
			residues.back().index_coord_map[i] = (is_ligand) ? &(this->coords[i]) : &(atom->coords);
		}
        }
}

void model::build_residue_info(){
	//If there is one rigid receptor and one ligand. This will work. Haven't tested on multiple receptor conformations, or multiple ligands. 
	//Make sure to clear all the vectors. Because if this object is implicitly copied (new_m = old_m), these vectors won't be empty to begin with. 
	this->adjusted_atom_indices.clear();
	this->grid_atom_ptrs.clear();
	this->adjusted_coords = this->coords;
	this->adjusted_coord_ptrs.clear();
	VINA_FOR_IN(i, this->adjusted_coords){
		this->adjusted_coord_ptrs.push_back(&(this->adjusted_coords[i]));
	}

	this->all_atom_ptrs.clear();
	this->all_atom_coords.clear();
	VINA_FOR_IN(i, this->grid_atoms){
		atom_vc& atom = this->grid_atoms[i];
		this->grid_atom_ptrs.push_back(&atom);
		this->all_atom_ptrs.push_back(&atom);
	 	this->all_atom_coords.push_back(&(atom.coords));
	}

	this->receptor_residues.clear();
	this->build_residues_from_atoms(this->receptor_residues, this->grid_atom_ptrs, false);

	this->atom_ptrs.clear();
	for (unsigned int i = 0 ; i < this->atoms.size(); i++){
		atom_vc& atom = this->atoms[i];
		this->atom_ptrs.push_back(&atom);
		this->all_atom_ptrs.push_back(&atom);
		this->all_atom_coords.push_back(&(this->coords[i]));
	}

	this->ligand_residues.clear();
	this->build_residues_from_atoms(this->ligand_residues, this->atom_ptrs, true);
	std::cout << "Num ligand residue: " << this->ligand_residues.size() << std::endl;

}

bool model::is_trigonal_planar(atom_vc* atom){
	sz element = atom->el;
	if (atom->bonds.size() != 3) return false;

	atom_vc& neighbor1 = get_atom(atom->bonds[0].connected_atom_index);
	atom_vc& neighbor2 = get_atom(atom->bonds[1].connected_atom_index);
	atom_vc& neighbor3 = get_atom(atom->bonds[2].connected_atom_index);

	sz index0 = std::distance(this->all_atom_ptrs.begin(), std::find(this->all_atom_ptrs.begin(), this->all_atom_ptrs.end(), atom));
	sz index1 = std::distance(this->all_atom_ptrs.begin(), std::find(this->all_atom_ptrs.begin(), this->all_atom_ptrs.end(), &neighbor1));
	sz index2 = std::distance(this->all_atom_ptrs.begin(), std::find(this->all_atom_ptrs.begin(), this->all_atom_ptrs.end(), &neighbor2));
	sz index3 = std::distance(this->all_atom_ptrs.begin(), std::find(this->all_atom_ptrs.begin(), this->all_atom_ptrs.end(), &neighbor3));

	VINA_CHECK (index0 < this->all_atom_coords.size());
	VINA_CHECK (index1 < this->all_atom_coords.size());
	VINA_CHECK (index2 < this->all_atom_coords.size());
	VINA_CHECK (index3 < this->all_atom_coords.size());

	vec* coord0 = this->all_atom_coords[index0];
	vec* coord1 = this->all_atom_coords[index1];
	vec* coord2 = this->all_atom_coords[index2];
	vec* coord3 = this->all_atom_coords[index3];

	vec v1(coord1->data[0] - coord0->data[0], coord1->data[1] - coord0->data[1], coord1->data[2] - coord0->data[2]);
	vec v2(coord2->data[0] - coord1->data[0], coord2->data[1] - coord1->data[1], coord2->data[2] - coord1->data[2]);
	vec v3(coord3->data[0] - coord2->data[0], coord3->data[1] - coord2->data[1], coord3->data[2] - coord2->data[2]);

	vec normal1 = cross_product(v1, v2);
	vec normal2 = cross_product(v2, v3);
	normalize_vec_in_place(normal1);
	normalize_vec_in_place(normal2);

	fl improper_torsion_cosine_abs = std::abs(dot_product(normal1, normal2));
	if (improper_torsion_cosine_abs > trigonal_planar_improper_cosine_allowance){
		return true;
	}
	return false;
}

bool model::is_aromatic_heteroatom_element(sz& el){
	return (el == EL_TYPE_O || el == EL_TYPE_N || el == EL_TYPE_S);
}

bool model::is_aromatic_cycle(aptrv& ring){

	VINA_FOR_IN(i, ring){
		atom_vc* ring_atom = ring[i];
		sz element = ring_atom->el;
		sz num_neighbors = ring_atom->bonds.size();
		if (num_neighbors >= 4) return false;

		if (element == EL_TYPE_C){
			if (!is_trigonal_planar(ring_atom)) return false;
		}
		else if (!is_aromatic_heteroatom_element(element)) return false;
		else if (num_neighbors == 3){
			if (!is_trigonal_planar(ring_atom)) return false;
		}
	}

	return true;
}

void model::DetectAliphaticCarbons(){
	//std::string flexible_residues = "ARG_LYS_GLU_GLN_MET";
	this->rec_ali_carb_info.clear();
	aptrv all_receptor_aromatic_carbons;
	VINA_FOR_IN(i, this->rec_ar_ring_info){
		//all_receptor_aromatic_carbons.push_back(this->rec_ar_ring_info[i].aromatic_carbon_ptrs);
		aptrv& aromatic_carbons = this->rec_ar_ring_info[i].aromatic_carbon_ptrs;
		all_receptor_aromatic_carbons.insert(all_receptor_aromatic_carbons.end(), aromatic_carbons.begin(), aromatic_carbons.end());
	}
	VINA_FOR_IN(i, this->grid_atoms){
		atom_vc& atom = this->grid_atoms[i];
		sz element = atom.el;
		if (element != EL_TYPE_C) continue;
		if (has(all_receptor_aromatic_carbons, &atom)) continue;

		bool flexible = false;
		std::string resname = atom.resname;
		resname.erase(std::remove_if(resname.begin(), resname.end(), isspace), resname.end());
		std::string atomname = atom.atomname;
		atomname.erase(std::remove_if(atomname.begin(), atomname.end(), isspace), atomname.end());

		/*if (flexible_residues.find(resname) != std::string::npos){
			flexible = true;
		}
		//Skipping side chain carbon atoms of flexible residues (linear && >= 3 rotatable bonds)
		if (flexible && atomname != "CA"){
			continue; 
		}*/
		

		//vec* c_coord = this->all_atom_coord_ptr_map[&atom];
		vec* c_coord = &(atom.coords);

                aptrv h_neighbors;
                vptrv h_coords;
		szv h_neighbor_indices;
                VINA_FOR_IN(j, atom.bonds){
                        const bond_vc& b = atom.bonds[j];
                        atom_vc& h = get_atom(b.connected_atom_index);
                        if (h.ad != AD_TYPE_H) continue;

                        h_neighbors.push_back(&h);
			vec* h_coord = &(h.coords);
                        h_coords.push_back(h_coord);
			//sz h_index = std::distance(this->grid_atom_ptrs.begin(), std::find(this->grid_atom_ptrs.begin(), this->grid_atom_ptrs.end(), &h));
			sz h_index = b.connected_atom_index.i;
			h_neighbor_indices.push_back(h_index);
                }
                this->rec_ali_carb_info.push_back(aliphatic_carbon_attribute(&atom, i, c_coord, h_neighbors, h_neighbor_indices, h_coords, this->chpi_explicit_hydrogen, false));
	}

	this->lig_ali_carb_info.clear();
	aptrv all_ligand_aromatic_carbons;
	VINA_FOR_IN(i, this->lig_ar_ring_info){
		aptrv& aromatic_carbons = this->lig_ar_ring_info[i].aromatic_carbon_ptrs;
		all_ligand_aromatic_carbons.insert(all_ligand_aromatic_carbons.end(), aromatic_carbons.begin(), aromatic_carbons.end());
        }
	VINA_FOR_IN(i, this->atoms){
		atom_vc& atom = this->atoms[i];
		sz element = atom.el;
                if (element != EL_TYPE_C) continue;
		if (has(all_ligand_aromatic_carbons, &atom)) continue;

		//std::cout << "Found ligand aliphatic carbon: " << atom.atomname << std::endl;
                //vec* c_coord = this->all_atom_coord_ptr_map[&atom];
                vec* c_coord = &(this->coords[i]);
                aptrv h_neighbors;
                vptrv h_coords;
		szv h_neighbor_indices;
                VINA_FOR_IN(j, atom.bonds){
                        const bond_vc& b = atom.bonds[j];
                        atom_vc& h = get_atom(b.connected_atom_index);
                        if (h.ad != AD_TYPE_H) continue;

                        h_neighbors.push_back(&h);
			//sz h_index = std::distance(this->atom_ptrs.begin(), std::find(this->atom_ptrs.begin(), this->atom_ptrs.end(), &h));
			sz h_index = b.connected_atom_index.i;
			VINA_CHECK(h_index < this->atoms.size());
                        vec* h_coord = &(this->coords[h_index]);
                        h_coords.push_back(h_coord);
			h_neighbor_indices.push_back(h_index);
                }
                this->lig_ali_carb_info.push_back(aliphatic_carbon_attribute(&atom, i, c_coord, h_neighbors, h_neighbor_indices, h_coords, this->chpi_explicit_hydrogen, true));
	}
}

void model::build_chpi_smooth(){

	//flv rs = this->calculate_rs_dense();
	flv rs = this->calculate_rs(sqr(chpi_dcut));
	this->num_rs = rs.size();

	this->chpi_c_smooth.clear();
	this->chpi_h_smooth.clear();
	this->chpi_entropy_smooth.clear();
	this->gauss1_smooth.clear();
	this->gauss2_smooth.clear();
	this->phobic_smooth.clear();

	this->chpi_c_fast.clear();
	this->chpi_h_fast.clear();
	this->chpi_entropy_fast.clear();
	this->gauss1_fast.clear();
	this->gauss2_fast.clear();
	this->phobic_fast.clear();

	this->chpi_c_smooth.resize(this->num_rs, pr(0, 0));
	this->chpi_h_smooth.resize(this->num_rs, pr(0, 0));
	this->chpi_entropy_smooth.resize(this->num_rs, pr(0, 0));
	this->gauss1_smooth.resize(this->num_rs, pr(0, 0));
	this->gauss2_smooth.resize(this->num_rs, pr(0, 0));
	this->phobic_smooth.resize(this->num_rs, pr(0, 0));

	this->chpi_c_fast.resize(this->num_rs, 0);
	this->chpi_h_fast.resize(this->num_rs, 0);
	this->chpi_entropy_fast.resize(this->num_rs, 0);
	this->gauss1_fast.resize(this->num_rs, 0);
	this->gauss2_fast.resize(this->num_rs, 0);
	this->phobic_fast.resize(this->num_rs, 0);

	VINA_FOR(i, this->num_rs){
		fl r = rs[i];
		this->chpi_c_smooth[i].first = this->weight_chpi * this->eval_chpi_enthalpy_c(r);
		this->chpi_h_smooth[i].first = this->weight_chpi * (this->eval_chpi_enthalpy_h(r) + this->eval_chpi_enthalpy_h2(r));
		this->gauss1_smooth[i].first = this->weight_gauss1 * gaussian(r , 0.5);

		//std::cout << "Gauss1 smooth " << i << " is " << this->gauss1_smooth[i].first << "\n";
		this->gauss2_smooth[i].first = this->weight_gauss2 * gaussian(r-3 , 2);
		this->phobic_smooth[i].first = this->weight_hydrophobic * slope_step(1.5, 0.5, r);
		this->chpi_entropy_smooth[i].first = -this->weight_chpi * this->eval_chpi_entropy(r);
	}
	//std::exit(1);
	
	//Below is essentially a copy of precalculate_element.init_from_smooth_fst(const flv& rs). 
	VINA_FOR(i, this->num_rs){
		fl& dor_c = this->chpi_c_smooth[i].second;
		fl& dor_h = this->chpi_h_smooth[i].second;
                fl& dor_e = this->chpi_entropy_smooth[i].second;
		fl& dor_g1 = this->gauss1_smooth[i].second;
		fl& dor_g2 = this->gauss2_smooth[i].second;
		fl& dor_phobic = this->phobic_smooth[i].second;

		if(i == 0 || i == this->num_rs - 1){
			dor_c = 0; dor_h = 0; dor_e = 0, dor_g1 = 0, dor_g2 = 0, dor_phobic = 0;
		}
		else{
			fl delta = rs[i+1] - rs[i-1];
			fl r = rs[i];
			dor_c = (this->chpi_c_smooth[i+1].first - this->chpi_c_smooth[i-1].first) / (delta * r);
			dor_h = (this->chpi_h_smooth[i+1].first - this->chpi_h_smooth[i-1].first) / (delta * r);
			dor_e = (this->chpi_entropy_smooth[i+1].first - this->chpi_entropy_smooth[i-1].first) / (delta * r);
			dor_g1 = (this->gauss1_smooth[i+1].first - this->gauss1_smooth[i-1].first) / (delta * r);
			dor_g2 = (this->gauss2_smooth[i+1].first - this->gauss2_smooth[i-1].first) / (delta * r);
			dor_phobic = (this->phobic_smooth[i+1].first - this->phobic_smooth[i-1].first) / (delta * r);
		}

		fl f1c = this->chpi_c_smooth[i].first;
                fl f2c = (i+1 >= this->num_rs) ? 0 : this->chpi_c_smooth[i+1].first;
                this->chpi_c_fast[i] = (f2c + f1c) / 2;

		fl f1h = this->chpi_h_smooth[i].first;
                fl f2h = (i+1 >= this->num_rs) ? 0 : this->chpi_h_smooth[i+1].first;
                this->chpi_h_fast[i] = (f2h + f1h) / 2;

		fl f1e = this->chpi_entropy_smooth[i].first;
                fl f2e = (i+1 >= this->num_rs) ? 0 : this->chpi_entropy_smooth[i+1].first;
                this->chpi_entropy_fast[i] = (f2e + f1e) / 2;

		fl f1g1 = this->gauss1_smooth[i].first; 
                fl f2g1 = (i+1 >= this->num_rs) ? 0 : this->gauss1_smooth[i+1].first;
		this->gauss1_fast[i] = (f1g1 + f2g1) / 2;

		
		fl f1g2 = this->gauss2_smooth[i].first; 
                fl f2g2 = (i+1 >= this->num_rs) ? 0 : this->gauss2_smooth[i+1].first;
		this->gauss2_fast[i] = (f1g2 + f2g2) / 2;

		fl f1pho = this->phobic_smooth[i].first; 
                fl f2pho = (i+1 >= this->num_rs) ? 0 : this->phobic_smooth[i+1].first;
		this->phobic_fast[i] = (f1pho + f2pho) / 2;
		//std::cout << "Chpi c smooth dor " << rs_c[i] << " is " << dor << std::endl;
		//std::cout << "Chpi c fast " << rs_c[i] << " is " << this->chpi_c_fast[i] << std::endl;
	}

	return;
	//std::exit(1);
}

pr model::eval_chpi_enthalpy_c_deriv(fl r){
	sz i1, i2; fl rem;
        //this->calc_smooth_i1_i2_rem(i1, i2, rem, r2); 
        this->calc_smooth_i1_i2_rem_dense(i1, i2, rem, r); 
	if (i1 >= this->num_rs -1){
		std::cout << "Eval C deriv overflow " << "r: " << r << std::endl;
		std::exit(1);
		//return pr (0,0); //Overflow, ususally too far beyond cutoff, or, other bug. 
	}
        const pr& p1 = this->chpi_c_smooth[i1];
        const pr& p2 = this->chpi_c_smooth[i2];
        fl e   = p1.first  + rem * (p2.first  - p1.first);
        fl dor = p1.second + rem * (p2.second - p1.second);
        return pr(e, dor);
}

pr model::eval_chpi_enthalpy_h_deriv(fl r, fl chn_cos_abs){
	fl e1 = this->eval_chpi_enthalpy_h(r);
	fl deriv = (chpi_miu_h - r) * inv_ssqr * e1;
	return pr(e1, deriv);
	//return pr(e1*chn_cos_abs, deriv*chn_cos_abs);

	//fl e2 = this->eval_chpi_enthalpy_h2(r);
	//fl deriv = (chpi_miu_h - r) * inv_ssqr * e1 + (chpi_miu_h2 - r) * inv_ssqr2 * e2;
	//return pr((e1+e2)*chn_cos_abs, deriv*chn_cos_abs*r);
	//return pr((e1+e2), deriv*r);
	//return pr((e1+e2) * chn_cos_abs , deriv * chn_cos_abs * r);

	/*if (num_interacting_h > 1){
		return pr((e1+e2) * chn_cos_abs , deriv * chn_cos_abs);
	}
	return pr(e1+e2 , deriv);*/

	/*fl r_left = r - dr; fl r_right = r + dr;
	fl e1 = this->eval_chpi_enthalpy_h(r);
	fl e1_left = this->eval_chpi_enthalpy_h(r_left);
	fl e1_right = this->eval_chpi_enthalpy_h(r_right);

	fl e2 = this->eval_chpi_enthalpy_h2(r);
	fl e2_left = this->eval_chpi_enthalpy_h2(r_left);
	fl e2_right = this->eval_chpi_enthalpy_h2(r_right);

	fl e = e1 + e2;
	//fl dor = (e_right - e_left) / (r_right - r_left); dor /= r;
	fl dor1 = (e1_right - e1_left) / (r_right - r_left) ;
	fl dor2 = (e2_right - e2_left) / (r_right - r_left) ;
	//fl dor = (dor1 + dor2) / r;
	fl dor = (dor1 + dor2);

        return pr(e , dor);*/

	/*if (num_interacting_h > 1){
        	return pr(this->weight_chpi * e * chn_cos_abs, this->weight_chpi * dor * chn_cos_abs);
        }
        return pr(this->weight_chpi * e , this->weight_chpi * dor);*/
}

pr model::eval_chpi_entropy_deriv(fl rh){
        rh = std::max(rh, chpi_ho_epsilon);
        fl e = this->eval_chpi_entropy(rh);
	//fl dor = -0.616 / rh;
	fl dor = -0.616 / (rh * rh);
	
        /*fl rh_left = rh - dr; rh_left = std::max(rh_left, chpi_ho_epsilon);
        fl rh_right = rh + dr;

        fl e = this->eval_chpi_entropy(rh);
        fl e_left = this->eval_chpi_entropy(rh_left);
        fl e_right = this->eval_chpi_entropy(rh_right);
        fl dor = (e_right - e_left) / ((rh_right - rh_left)*rh);*/

        //fl e_final = -this->weight_chpi * e;
        //fl dor_final = -this->weight_chpi * dor;
        //curl(e_final, dor_final, chpi_entropy_curl_v);

        return pr(e , dor);
}

pr model::eval_repulsion_deriv(fl r, fl cutoff){
	fl e = this->eval_repulsion(r, cutoff);
	fl deriv = 2.0 * (r - cutoff); //f'(x-a)^2 = 2 * (x-a) 
	return pr(e, deriv);
        /*fl r_left = r - dr; fl r_right = r + dr;
        fl e = this->eval_repulsion(r, cutoff);
        fl e_left = this->eval_repulsion(r_left, cutoff);
        fl e_right = this->eval_repulsion(r_right, cutoff);
        //fl dor = (e_right - e_left) / ((r_right - r_left)*r);
        fl dor = (e_right - e_left) / (r_right - r_left);
        return pr(e, dor);*/

        //return pr(-this->weight_chpi * e, -this->weight_chpi * dor);
}

pr model::eval_gauss1_deriv(fl r2){
	sz i1, i2; fl rem;
        this->calc_smooth_i1_i2_rem(i1, i2, rem, r2); 
	if (i1 >= this->num_rs -1){
		std::cout << "Eval gauss1 overflow r: " << r2 << std::endl;
		std::exit(1);
	}
        const pr& p1 = this->gauss1_smooth[i1];
        const pr& p2 = this->gauss1_smooth[i2];
        fl e   = p1.first  + rem * (p2.first  - p1.first);
        fl dor = p1.second + rem * (p2.second - p1.second);
        return pr(e, dor);
}

pr model::eval_gauss2_deriv(fl r2){
        sz i1, i2; fl rem;
        this->calc_smooth_i1_i2_rem(i1, i2, rem, r2);
        if (i1 >= this->num_rs -1){
                std::cout << "Eval gauss2 overflow r: " << r2 << std::endl;
                std::exit(1);
        }
        const pr& p1 = this->gauss2_smooth[i1];
        const pr& p2 = this->gauss2_smooth[i2];
        fl e   = p1.first  + rem * (p2.first  - p1.first);
        fl dor = p1.second + rem * (p2.second - p1.second);
        return pr(e, dor);
}

pr model::eval_phobic_deriv(fl r2){
	sz i1, i2; fl rem;
        this->calc_smooth_i1_i2_rem(i1, i2, rem, r2);
        if (i1 >= this->num_rs -1){
                std::cout << "Eval phobic overflow r: " << r2 << std::endl;
                std::exit(1);
        }
        const pr& p1 = this->phobic_smooth[i1];
        const pr& p2 = this->phobic_smooth[i2];
        fl e   = p1.first  + rem * (p2.first  - p1.first);
        fl dor = p1.second + rem * (p2.second - p1.second);
        return pr(e, dor);
}

pr model::eval_gauss1_fast(fl r2){
        sz i = sz(this->prec_factor * r2);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
	return pr(this->gauss1_fast[i], 0.00);
}

pr model::eval_gauss2_fast(fl r2){
        sz i = sz(this->prec_factor * r2);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
	return pr(this->gauss2_fast[i], 0.00);
}

pr model::eval_phobic_fast(fl r2){
        sz i = sz(this->prec_factor * r2);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
	return pr(this->phobic_fast[i], 0.00);
}
void model::calc_smooth_i1_i2_rem(sz& i1, sz& i2, fl& rem, fl r2){
	fl r2_factored = this->prec_factor * r2;
        i1 = sz(r2_factored);
        i2 = i1 + 1;
        rem = r2_factored - i1;
	return;
}

void model::calc_smooth_i1_i2_rem_dense(sz& i1, sz& i2, fl& rem, fl r){
	i1 = sz(r/dr);
        i2 = i1 + 1;
        rem = fl(r - dr*i1);
        return;
}

pr model::eval_chpi_enthalpy_c_fast(fl r){
        //sz i = sz(this->prec_factor * r2);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
        sz i = sz(r/dr);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
        return pr(this->chpi_c_fast[i], 0.00);
}

pr model::eval_chpi_enthalpy_h_fast(fl r, fl chn_cosine_abs){

	return pr(this->eval_chpi_enthalpy_h(r) + this->eval_chpi_enthalpy_h2(r), 0.00);
	//return pr ((this->eval_chpi_enthalpy_h(r) + this->eval_chpi_enthalpy_h2(r)) * chn_cosine_abs, 0.00);


        //sz i = sz(this->prec_factor * r2 * sqr_chpi_c_h_cutoff_ratio);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
        sz i = sz(r/dr);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
       	return pr(this->chpi_h_fast[i] * chn_cosine_abs, 0.00);
}

pr model::eval_chpi_entropy_fast(fl rh){
	return pr(this->eval_chpi_entropy(rh), 0.00);
        //sz i = sz(this->prec_factor * r2h * sqr_chpi_c_e_cutoff_ratio);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
        sz i = sz(rh/dr);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
	fl e_final = this->chpi_entropy_fast[i];
        //curl(e_final, chpi_entropy_curl_v);
        return pr(e_final, 0.00);
}

pr model::eval_repulsion_fast(fl r, fl cutoff){
	return pr(this->eval_repulsion(r, cutoff), 0.00);
	//sz i = sz(this->prec_factor * r2 * sqr_chpi_c_r_cutoff_ratio);
        //sz i = sz(r/dr);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
	//return pr(this->chpi_hc_repul_fast[i], 0.00);
}

fl model::eval_repulsion(fl r, fl cutoff){
	fl d = r - cutoff;
	if (d > 0) return 0;
	return d*d;
}

flv model::calculate_rs(fl sqr_cutoff){
	sz n =  sz(this->prec_factor * sqr_cutoff) + 3;
	flv tmp(n, 0);
        VINA_FOR(i, n){
 	       tmp[i] = std::sqrt(i / this->prec_factor);
	}
        return tmp;
}

flv model::calculate_rs_dense(){
	sz n = sz(chpi_dcut / dr) + 3;
	flv tmp(n, 0);
	VINA_FOR(i, n){
		tmp[i] = (fl) i * dr;
	}
	return tmp;
}

void model::adjust_out_of_box_coords(){
	this->adjusted_atom_indices.clear();
	VINA_FOR_IN(i, this->coords){
		this->adjusted_coords[i] = this->coords[i];
	}

	VINA_FOR_IN(i, this->adjusted_coords){
		atom_vc& a = this->atoms[i];
		if (a.ad == AD_TYPE_H) continue;

		vec& old = this->coords[i];
		vec& adjusted = this->adjusted_coords[i];

		bool out_of_box = false;
		VINA_FOR(j, 3){
			if (adjusted[j] < gd_dims[j].first){
				out_of_box = true;
				adjusted[j] = gd_dims[j].first;
			}
			else if (adjusted[j] > gd_dims[j].second){
				out_of_box = true;
				adjusted[j] = gd_dims[j].second;	
			}
		}

		//If this heavy atom is adjusted, move its hydrogen neighbors along with it. 
		if (out_of_box){
			this->adjusted_atom_indices.push_back(i);
			vec translation_vector = adjusted - old;
			VINA_FOR_IN(j, a.bonds) {
                		const bond_vc& b = a.bonds[j];
                		const atom_vc& h = get_atom(b.connected_atom_index);
				if (h.ad != AD_TYPE_H) continue;

				sz h_index = b.connected_atom_index.i;
				this->adjusted_coords[h_index] += translation_vector;
				this->adjusted_atom_indices.push_back(h_index);				

			}
		}
	}
}

fl model::GetDihedral(vec* a1, vec* a2, vec* a3, vec* a4){
    /*GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);*/

	vec b1 = *a2 - *a1;
	vec b2 = *a3 - *a2;
	vec b3 = *a4 - *a3;
	vec b4 = b2;
	b4 *= -1;

    /*GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);*/
	vec b2xb3 = cross_product(b2, b3);
	vec b1_m_b2n(b1); b1_m_b2n *= magnitude(b2);
	vec b1xb2 = cross_product(b1, b2);

    //fl current_dihedral = std::atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    fl current_dihedral = std::atan2(dot_product(b1_m_b2n, b2xb3), dot_product(b1xb2, b2xb3));
    return current_dihedral;
}

void model::eval_chpi_entropy_set_force2(vec* h, sz closest_h_i, sz h_closest_index, std::vector<flv>& h_ring_dists, ring_attribute& r, pr& dH_minusTdS, bool ligand_aliphatic){
	szv& ring_atom_indices = r.all_atom_indices;
	vptrv& ring_atom_coords = r.all_atom_coord_ptrs;
	sz ring_size = r.num_ring_atoms;
	vec& centroid = r.centroid;
	vec& normal = r.normal;

	flv edge_dist_sum(ring_size, 0);		
	flv& closest_h_ring_dists = h_ring_dists[closest_h_i];
	//std::cout << "Ring size: " << ring_size << " and closest_h_ring_dists size: " << closest_h_ring_dists.size() << std::endl;
	//std::cout << "Closest h index is: " << closest_h_i << std::endl;

	VINA_FOR(i, ring_size){
		sz j = (i < ring_size -1) ? i+1 : 0;
		edge_dist_sum[i] = (closest_h_ring_dists[i] + closest_h_ring_dists[j]);
	}

	/*VINA_FOR_IN(i, ring_atom_coords){
		vec* ra = ring_atom_coords[i];
		vec* next = (i < ring_size -1) ? ring_atoms[i+1] : ring_atoms[0];
		fl rh = std::sqrt(sqr(h[0]-ra[0]) + sqr(h[1]-ra[1]) + sqr(h[2]-ra[2]));
		fl rh_next = std::sqrt(sqr(h[0]-next[0]) + sqr(h[1]-next[1]) + sqr(h[2]-next[2]));
		edge_dist_sum[i] = (rh + rh_next);
	}*/

	sz closest_edge_index = 0; fl lowest_sum_dist = max_fl;
	VINA_FOR_IN(i, edge_dist_sum){
		fl& dist_sum = edge_dist_sum[i];
		if (dist_sum < lowest_sum_dist){
			lowest_sum_dist = dist_sum;
			closest_edge_index = i;			
		}
	}

	//a1 and a2 are the two atoms forming the edge closest to the hydrogen atom.
	//std::cout << "Closest edge index is: " << closest_edge_index << std::endl;
	vec* a1 = ring_atom_coords[closest_edge_index];
	vec* a2 = (closest_edge_index < ring_size -1) ? ring_atom_coords[closest_edge_index+1] : ring_atom_coords[0];
	fl omega = GetDihedral(h, a1, a2, &centroid);
	//std::cout << "Omega degrees: " << omega / M_PI * 180.0000 << std::endl;
	fl omega_abs = std::abs(omega);

	//Project hydrogen onto the plane, calling it h_prime
	vec a1h = *h - *a1;
	fl ahdp = a1h * normal; //Dot product betwwen a1h and ring normal
	fl a1h_vo = std::abs(ahdp); //Vertical distance between h and ring plane.

	//Make a vector parallel to the ring normal, starting at h_prime, ending at hydrogen.
	vec vertical(normal);
        vertical *= a1h_vo;
        if (ahdp < 0) vertical *= -1; //If hydrogen and ring normal on oppposite of ring, flip vertical to be on the same side. 
	vec h_prime = *h - vertical;

	//In triangle h-a1-a2, draw the height on edge a1-a2, intercepting a1-a2 on point b
	vec a1a2 = *a2 - *a1;
	fl a1a2_length = magnitude(a1a2);
	fl a1b_length = std::abs(a1h * a1a2) / a1a2_length;
	vec a1b(a1a2); a1b *= (a1b_length / magnitude(a1a2));
	vec b = *a1 + a1b;
	vec b_h_prime = h_prime - b;

	//DEBUG: check b_hprime * a1a2 = 0
	fl hb_a1a2_dp = b_h_prime * a1a2;
	std::cout << "Hb a1a2 dp is : " << hb_a1a2_dp << std::endl;
	//std::cout << "Point b coord " << b[0] << "," << b[1] << "," << b[2] << std::endl;

	fl b_hp_length = magnitude(b_h_prime);
	//std::cout << "b hp length: " << b_hp_length << std::endl;

	vec midpoint_centroid_estimated(b_h_prime); normalize_vec_in_place(midpoint_centroid_estimated);
	midpoint_centroid_estimated *= chpi_mp_distance_estimated;
	//std::cout << "Midpoint centroid length: " << magnitude(midpoint_centroid_estimated) << std::endl;

	//std::cout << "mc estd b4 " << midpoint_centroid_estimated[0] << "," << midpoint_centroid_estimated[1] << "," << midpoint_centroid_estimated[2] << std::endl;
	if (omega_abs >= half_pi){
		//std::cout << "Omega beyond half pi." << std::endl;
		midpoint_centroid_estimated *= -1;
	}

	//Estimate centroid location based on edge a1a2, assuming it's benzene.
	vec midpoint = *a1 + *a2; midpoint *= 0.5;
	//std::cout << "Midpoint " << midpoint[0] << "," << midpoint[1] << "," << midpoint[2] << std::endl;
	//std::cout << "mc estd after " << midpoint_centroid_estimated[0] << "," << midpoint_centroid_estimated[1] << "," << midpoint_centroid_estimated[2] << std::endl;
	vec centroid_estimated = midpoint + midpoint_centroid_estimated;
	//std::cout << "Centroid estimated: " << centroid_estimated[0] << "," << centroid_estimated[1] << "," << centroid_estimated[2] << std::endl;
	vec h_prime_centroid = h_prime - centroid_estimated;
	//std::cout << "H prime xyz: " << h_prime[0] << "," << h_prime[1] << "," << h_prime[2] << std::endl;

	fl horizontal_offset = magnitude(h_prime_centroid);
	fl ho_effective = std::max(horizontal_offset, chpi_ho_epsilon);
	//std::cout << "Entropy ho effective: " << ho_effective << std::endl;
	if (ho_effective > chpi_ho_max) return;

	//TODO: Apply weight chpi on entropic force, vertical factor on both energy and force.
	fl vertical_factor = slope_step(chpi_vo_max_h, chpi_vertical_optimum_h, a1h_vo);
        fl minusTdS = vertical_factor * this->eval_chpi_entropy(ho_effective);
        fl deriv = vertical_factor * -this->weight_chpi * -0.616 / ho_effective;

	dH_minusTdS.second += minusTdS;

	vec entropic_force(h_prime_centroid); entropic_force *= (deriv / horizontal_offset);
	if (ligand_aliphatic){
		this->minus_forces[h_closest_index] += entropic_force;
	}
	else{
        	entropic_force *= -1; //Reverse the direction of derivative so that it is still receptor->ligand
		vec force_per_ring_atom = entropic_force;
        	force_per_ring_atom /= (fl) ring_size;

        	VINA_FOR(i, r.num_ring_atoms){
                	sz ra_index = ring_atom_indices[i];
                	//Apply total entropy evenly on each ring atom
                	this->minus_forces[ra_index] += force_per_ring_atom;
        	}
	}
	return;
	
}

void model::eval_chpi_entropy_set_force_old(vec* h_closest, sz h_closest_index, ring_attribute& r, pr& dH_minusTdS, bool ligand_aliphatic){
	//vec* closest_centroid = this->choose_closest_centroid(h_closest, r.centroid, r.subcycle_centroids);
	vec* closest_centroid = &(r.centroid);
	this->eval_chpi_entropy_set_force_old_each_centroid(h_closest, h_closest_index, r, *closest_centroid, dH_minusTdS, ligand_aliphatic);
	return;

	/*vecv& subcycle_centroids = r.subcycle_centroids;
	if (!subcycle_centroids.empty()){
		VINA_FOR_IN(i, subcycle_centroids){
			this->eval_chpi_entropy_set_force_old_each_centroid(h_closest, h_closest_index, r, subcycle_centroids[i], dH_minusTdS, ligand_aliphatic);
		}
	}
	else{
		this->eval_chpi_entropy_set_force_old_each_centroid(h_closest, h_closest_index, r, r.centroid, dH_minusTdS, ligand_aliphatic);
	}
	return;*/
}

void model::eval_chpi_entropy_set_force_old_each_centroid(vec* h_closest, sz h_closest_index, ring_attribute& r, vec& centroid, pr& dH_minusTdS, bool ligand_aliphatic){
	szv& ring_atom_indices = r.all_atom_indices;
        vptrv& ring_atom_coords = r.all_atom_coord_ptrs;
        sz ring_size = r.num_ring_atoms;
        vec& normal = r.normal;

	vec h_closest_c = (*h_closest) - centroid;
	fl h_closest_cl = magnitude(h_closest_c);
	fl h_closest_dp = h_closest_c * normal;
	fl hcvo_closest = std::abs(h_closest_dp);
	fl hcho_closest = std::sqrt(sqr(h_closest_cl) - sqr(hcvo_closest));
	hcho_closest = std::max(hcho_closest, chpi_ho_epsilon);
	//std::cout << "Hcho closest: " << hcho_closest << std::endl;

	bool entropy_in_cutoff = (hcvo_closest < chpi_vo_max_h && hcho_closest < chpi_ho_max);
	if (!entropy_in_cutoff) return;

	
	vec vertical(r.normal);
        vertical *= hcvo_closest;
        if (h_closest_dp < 0) vertical *= -1; //If hydrogen and ring normal on oppposite of ring, flip vertical to be on the same side. 
        vec p = centroid + vertical;
	vec horizontal = (*h_closest - p);

	//pr entropy_e_dor = (fast) ? this->eval_chpi_entropy_fast(ccho) : this->eval_chpi_entropy_deriv(ccho);
	//fl minusTdS = entropy_e_dor.first; fl entropy_dor = -this->weight_chpi * entropy_e_dor.second;
	
	//return (-0.594 * std::log(ho_effective) + 1.0609);  //6.0A, using intra-protein interactions. 
        //return (-0.594 * std::log(ho_effective) + 1.0066);  //5.5A, using intra-protein interactions. 
        //return (-0.595 * std::log(ho_effective) + 0.9460);  //5.0A, using intra-protein interactions. 
        //return (-0.595 * std::log(ho_effective) + 0.9087);  //4.7A, using intra-protein interactions. 
        //return (-0.597 * std::log(ho_effective) + 0.8846);  //4.5A, using intra-protein interactions. 
        //return (-0.599 * std::log(ho_effective) + 0.8172);  //4.0A, using intra-protein interactions. 
        //return (-0.603 * std::log(ho_effective) + 0.7505);  3.5A
        //return (-0.605 * std::log(ho_effective) + 0.6564);  //3.0A
        //return (-0.607 * std::log(ho_effective) + 0.5472);  //2.5A
        //return (-0.611 * std::log(ho_effective) + 0.4192);  //2.0A
        //return (-0.613 * std::log(ho_effective) + 0.2108);  //1.4A, using intra-protein interactions.

	fl vertical_factor = slope_step(chpi_vo_max_h, chpi_vertical_optimum_h, hcvo_closest);
	//fl vertical_factor = 1;
	if (hcho_closest > chpi_ho_max) return;
	fl minusTdS = vertical_factor * eval_chpi_entropy(hcho_closest);
	fl entropy_d = vertical_factor * -this->weight_chpi * -0.595 / hcho_closest;
	fl entropy_dor = entropy_d / hcho_closest; 

	//std::cout << "Entopy vertical offset " << hcvo_closest << " and factor: " << vertical_factor << std::endl;
	
        dH_minusTdS.second += minusTdS;
	//return;
        vec entropy_deriv = horizontal; entropy_deriv *= entropy_dor;

	if (ligand_aliphatic){
         	//this->minus_forces[carbon_index] += entropy_deriv;
         	this->minus_forces[h_closest_index] += entropy_deriv;
		return;
	}
	
	/*szv& ring_atom_indices = r.all_atom_indices;
        vptrv& ring_atom_coord_ptrs = r.all_atom_coord_ptrs;
        std::map<fl, sz> cc_dist_index_map;

        VINA_FOR(i, r.num_ring_atoms){
                vec* ra_coord = ring_atom_coord_ptrs[i];
                vec c_ra = (*ra_coord) - (*c_coord);
                fl r = magnitude(c_ra);
                sz ring_atom_index = ring_atom_indices[i];
                cc_dist_index_map[r] = ring_atom_index;
        }
	
	fl force_magnitude = magnitude(entropy_deriv) / (fl) r.num_ring_atoms;
	vec up_force(vertical); normalize_vec_in_place(up_force); up_force *= force_magnitude;
	vec down_force(up_force); down_force *= -1;
	
	sz num_ring_atoms_go_up = sz(r.num_ring_atoms / 2);
	int nth_ra = 0;
	for (std::map<fl, sz>::iterator mapit = cc_dist_index_map.begin(); mapit != cc_dist_index_map.end(); mapit++){
		nth_ra++;
		sz this_ra_index = mapit->second;
		if (nth_ra <= num_ring_atoms_go_up){
			this->minus_forces[this_ra_index] += up_force;		
		}
		else{
			this->minus_forces[this_ra_index] += down_force;
		}
	}*/
	

	//Method 2: Move the entire aromatic ring horizontally, just like what I do for the aliphatic carbon. Is it effective?
       	entropy_deriv *= -1; //Reverse the direction of derivative so that it is still receptor->ligand
	vec entropy_deriv_per_ring_atom = entropy_deriv;
       	entropy_deriv_per_ring_atom /= (fl) r.num_ring_atoms;

       	VINA_FOR(i, r.num_ring_atoms){
               	sz ra_index = ring_atom_indices[i];
               	//Apply total entropy evenly on each ring atom
               	this->minus_forces[ra_index] += entropy_deriv_per_ring_atom;
       	}
	return;

}

void model::eval_chpi_entropy_set_force(vec* h, sz closest_h_i, sz h_closest_index, std::vector<flv>& h_ring_dists, ring_attribute& r, pr& dH_minusTdS, bool ligand_aliphatic){
	szv& ring_atom_indices = r.all_atom_indices;
	vptrv& ring_atom_coords = r.all_atom_coord_ptrs;
	sz ring_size = r.num_ring_atoms;
	vec& centroid = r.centroid;
	vec& normal = r.normal;

	flv edge_dist_sum(ring_size, 0);		
	flv& closest_h_ring_dists = h_ring_dists[closest_h_i];
	//std::cout << "Ring size: " << ring_size << " and closest_h_ring_dists size: " << closest_h_ring_dists.size() << std::endl;
	//std::cout << "Closest h index is: " << closest_h_i << std::endl;

	VINA_FOR(i, ring_size){
		sz j = (i < ring_size -1) ? i+1 : 0;
		edge_dist_sum[i] = (closest_h_ring_dists[i] + closest_h_ring_dists[j]);
	}

	/*VINA_FOR_IN(i, ring_atom_coords){
		vec* ra = ring_atom_coords[i];
		vec* next = (i < ring_size -1) ? ring_atoms[i+1] : ring_atoms[0];
		fl rh = std::sqrt(sqr(h[0]-ra[0]) + sqr(h[1]-ra[1]) + sqr(h[2]-ra[2]));
		fl rh_next = std::sqrt(sqr(h[0]-next[0]) + sqr(h[1]-next[1]) + sqr(h[2]-next[2]));
		edge_dist_sum[i] = (rh + rh_next);
	}*/

	sz closest_edge_index = 0; fl lowest_sum_dist = max_fl;
	VINA_FOR_IN(i, edge_dist_sum){
		fl& dist_sum = edge_dist_sum[i];
		if (dist_sum < lowest_sum_dist){
			lowest_sum_dist = dist_sum;
			closest_edge_index = i;			
		}
	}

	//a1 and a2 are the two atoms forming the edge closest to the hydrogen atom.
	vec* a1 = ring_atom_coords[closest_edge_index];
	vec* a2 = (closest_edge_index < ring_size -1) ? ring_atom_coords[closest_edge_index+1] : ring_atom_coords[0];
	fl omega = GetDihedral(h, a1, a2, &centroid);
	//std::cout << "Omega degrees: " << omega / M_PI * 180.0000 << std::endl;
	fl omega_abs = std::abs(omega);

	//Project hydrogen onto the plane, calling it h_prime
	vec a1h = *h - *a1;
	fl ahdp = a1h * normal; //Dot product betwwen a1h and ring normal
	fl a1h_vo = std::abs(ahdp); //Vertical distance between h and ring plane.

	//Make a vector parallel to the ring normal, starting at h_prime, ending at hydrogen.
	vec vertical(normal);
        vertical *= a1h_vo;
        if (ahdp < 0) vertical *= -1; //If hydrogen and ring normal on oppposite of ring, flip vertical to be on the same side. 
	vec h_prime = *h - vertical;

	//In triangle h-a1-a2, draw the height on edge a1-a2, intercepting a1-a2 on point b
	vec a1a2 = *a2 - *a1;
	fl a1a2_length = magnitude(a1a2);
	fl a1b_length = std::abs(a1h * a1a2) / a1a2_length;
	vec a1b(a1a2); a1b *= (a1b_length / magnitude(a1a2));
	vec b = *a1 + a1b;
	vec b_h_prime = h_prime - b;

	//DEBUG:verify h prime is on the plane
	/*vec hp_centroid = h_prime - centroid;
	fl dp = hp_centroid * normal;
	std::cout << "Dp is: " << dp << std::endl;*/

	fl b_hp_length = magnitude(b_h_prime);
	fl ho_effective = 0;
	vec entropic_force(b_h_prime); normalize_vec_in_place(entropic_force);

	if (omega_abs >= half_pi){
		//std::cout << "Omega beyond half pi." << std::endl;
		ho_effective = chpi_ho_max - chpi_ho_out_of_edge + b_hp_length;
		if (ho_effective >= chpi_ho_max) return;
	}
	else{
		//std::cout << "Omega smaller than half pi." << std::endl;
		entropic_force *= -1;
		ho_effective = 	chpi_ho_max - chpi_ho_out_of_edge - b_hp_length;
	}

	ho_effective = std::max(ho_effective, chpi_ho_epsilon);
	//std::cout << "Entropy ho effective: " << ho_effective << std::endl;
	fl vertical_factor = slope_step(chpi_vo_max_h, chpi_vertical_optimum_h, a1h_vo);
	fl minusTdS = vertical_factor * this->eval_chpi_entropy(ho_effective);
	dH_minusTdS.second += minusTdS;
	//fl deriv = vertical_factor * -this->weight_chpi * -0.616 / ho_effective;
	fl deriv = vertical_factor * -0.616 / ho_effective;
	entropic_force *= deriv;

	if (ligand_aliphatic){
		this->minus_forces[h_closest_index] += entropic_force;
	}
	else{
        	entropic_force *= -1; //Reverse the direction of derivative so that it is still receptor->ligand
		vec force_per_ring_atom = entropic_force;
        	force_per_ring_atom /= (fl) ring_size;

        	VINA_FOR(i, r.num_ring_atoms){
                	sz ra_index = ring_atom_indices[i];
                	//Apply total entropy evenly on each ring atom
                	this->minus_forces[ra_index] += force_per_ring_atom;
        	}
	}
	return;
	
}

vec* model::choose_closest_centroid(vec* h, vec& centroid, vecv& centroids){
	if (centroids.empty()) return (&centroid);
	flv distances(centroids.size(), 0);

	VINA_FOR_IN(i, centroids){
		vec& c = centroids[i];
		vec hc = *h - c;
		fl r = magnitude(hc);
		distances[i] =  r;
	}
	
	flv::iterator it = std::min_element(distances.begin(), distances.end());
	int index = std::distance(distances.begin(), it);
	return &(centroids[index]);
}

pr model::calc_horizontal_factor(fl rh){
	if (rh <= chpi_ho_max) return pr(1,0);
	else if (rh < chpi_ho_end){
		fl a1 = rh - chpi_ho_end;
		fl a2 = chpi_ho_max - chpi_ho_end;
		fl z = a1/a2;

		//fl f = 1 - sqr(a1/a2);
		//fl deriv = -2.0 * a1 / sqr(a2);
		fl f = std::sqrt(z);
		fl deriv = 0.5 * std::pow(z, -1.5) / a2;
		return pr(f, deriv);
	}
	return pr(0,0);
}
