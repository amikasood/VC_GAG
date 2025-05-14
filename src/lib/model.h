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

#ifndef VINA_MODEL_H
#define VINA_MODEL_H

#include <boost/optional.hpp> // for context
#include <map>
#include <algorithm> //std::find
#include <cmath>

#include "file.h"
#include "tree.h"
#include "matrix.h"
#include "precalculate.h"
#include "igrid.h"
#include "grid_dim.h"
#include "atom_constants.h" //Yao 20231023

#include "atom.h"
struct interacting_pair {
	sz type_pair_index;
	sz a; 
	sz b;
	interacting_pair(sz type_pair_index_, sz a_, sz b_) : type_pair_index(type_pair_index_), a(a_), b(b_) {} 
};

typedef std::vector<interacting_pair> interacting_pairs; 

typedef std::pair<std::string, boost::optional<sz> > parsed_line; //Pair: This class couples together a pair of values, which may be of different types (T1 and T2). The individual values can be accessed through the public members first and second.
typedef std::vector<parsed_line> context; 

typedef std::vector<atom_vc*> aptrv; //Yao 20230604, vector of pointers to object atom_vc. Use pointers to avoid duplication of objects. 

struct ligand : public flexible_body, atom_range {
	unsigned degrees_of_freedom; // can be different from the apparent number of rotatable bonds, because of the disabled torsions
	interacting_pairs pairs;
	context cont;
	ligand(const flexible_body& f, unsigned degrees_of_freedom_) : flexible_body(f), atom_range(0, 0), degrees_of_freedom(degrees_of_freedom_) {}
	void set_range();
};

struct residue_vc : public main_branch {
	residue_vc(const main_branch& m) : main_branch(m) {}
};

//Yao added 20230610
struct model_residue {
	std::string name;
	std::string number;
	std::string chainID;
	std::map<sz, atom_vc*> index_atom_map;
	std::map<sz, vec*> index_coord_map;
	model_residue(std::string& resname, std::string& resnum, std::string& chainID_){
		name = resname;
		number = resnum;
		chainID = chainID_;
	}
};
typedef std::vector<model_residue> resv;

//Yao added 20230706
typedef std::vector<vec*> vptrv;

struct ring_attribute {
	ring_attribute(model* m_, aptrv& all_atoms_, szv& all_atom_indices_, vptrv& all_atom_coords, bool is_ligand_, aptrv& ring_hydrogens_, szv& ring_h_indices_, vptrv& ring_h_coords_, szv& ring_h_heavy_neighbor_indices_, vptrv& ring_h_heavy_neighbor_coords_, vptrv* adjusted_coords_, bool use_adjusted_coords, aptrv& model_atoms, std::vector<aptrv>& subcycles_){
		this->m = m_;
		this->adjusted_coords = adjusted_coords_;
		this->all_atom_ptrs = all_atoms_;
		this->all_atom_indices = all_atom_indices_;
		this->is_ligand = is_ligand_; 

		this->ring_hydrogens = ring_hydrogens_;
		this->ring_h_indices = ring_h_indices_;
		this->ring_h_coords  = ring_h_coords_;
		this->ring_h_heavy_neighbor_indices = ring_h_heavy_neighbor_indices_;
		this->ring_h_heavy_neighbor_coords = ring_h_heavy_neighbor_coords_;
		this->num_ring_hydrogens =  this->ring_hydrogens.size();
		this->subcycles = subcycles_;
		this->subcycle_centroids.resize(subcycles_.size(), vec(0,0,0));

		this->all_atom_coord_ptrs.clear(); //Unnecessary?
		this->aromatic_carbon_ptrs.clear();
                this->aromatic_carbon_indices.clear();
                this->aromatic_carbon_coord_ptrs.clear();
		this->heteroatom_ptrs.clear();
                this->heteroatom_indices.clear();
                this->heteroatom_coord_ptrs.clear();
		
		VINA_FOR_IN(i, this->all_atom_ptrs){
			atom_vc* ring_atom = this->all_atom_ptrs[i];

			sz atom_index = this->all_atom_indices[i];
			VINA_CHECK(atom_index < all_atom_coords.size());
			vec* coord_ptr = all_atom_coords[atom_index];
			this->all_atom_coord_ptrs.push_back(coord_ptr);

			sz element = ring_atom->el;
			if (element == EL_TYPE_C){
				this->aromatic_carbon_ptrs.push_back(ring_atom);
				this->aromatic_carbon_indices.push_back(atom_index);
				this->aromatic_carbon_coord_ptrs.push_back(coord_ptr);
			}
			else{
				this->heteroatom_ptrs.push_back(ring_atom);
				this->heteroatom_indices.push_back(atom_index);
				this->heteroatom_coord_ptrs.push_back(coord_ptr);
			}
			
		}

		this->subcycle_coords.clear(); //Index internal to this ring attribute object, not to the model. 
                VINA_FOR_IN(i, this->subcycles){
                        aptrv& this_cycle = this->subcycles[i];
			vptrv this_cycle_coords;

                        VINA_FOR_IN(j, this_cycle){
                                atom_vc* cycle_atom = this_cycle[j];
                                aptrv::iterator it = std::find(this->all_atom_ptrs.begin(), this->all_atom_ptrs.end(), cycle_atom);
				VINA_CHECK(it != this->all_atom_ptrs.end());

                                sz internal_index = std::distance(this->all_atom_ptrs.begin(), it);
				this_cycle_coords.push_back(this->all_atom_coord_ptrs[internal_index]);
				
                        }

			this->subcycle_coords.push_back(this_cycle_coords);
                }

		/*std::cout << "Num subcycle coords: " << this->subcycle_coords.size() << std::endl;

		//if (this->all_atom_ptrs[0]->resname.find("TRP") != std::string::npos){
			VINA_FOR_IN(i,this->subcycle_coords){
				std::cout << "Subcycle coords size: " << this->subcycle_coords[i].size() << std::endl;
			}
		//}*/

		this->num_ring_atoms = this->all_atom_ptrs.size();
		this->num_aromatic_carbons = this->aromatic_carbon_ptrs.size();
		this->num_heteroatoms = this->heteroatom_ptrs.size();

		this->ComputeRemainingAttributes(use_adjusted_coords);
		this->fused_edges.clear();
		this->DetectFusedEdges(model_atoms);
	}

	void DetectFusedEdges(aptrv& model_atoms){
		aptrv fused_atoms;
		szv fused_atom_indices;
		VINA_FOR_IN(i, this->all_atom_ptrs){
			atom_vc* a = this->all_atom_ptrs[i];
			sz num_ring_neighbors = 0;

			VINA_FOR_IN(j, a->bonds){
				const bond_vc& b = a->bonds[j];
				sz neighbor_index = b.connected_atom_index.i;
				atom_vc* a2 = model_atoms[neighbor_index];
				if (has(this->all_atom_ptrs, a2)){
					num_ring_neighbors++;
				}

			}

			if (num_ring_neighbors >= 3){
				fused_atoms.push_back(a);	
				sz atom_index = this->all_atom_indices[i];
				fused_atom_indices.push_back(atom_index);
			}
		}

		VINA_FOR_IN(i, fused_atoms){
			atom_vc* a = fused_atoms[i];
			sz a_index = fused_atom_indices[i];
			aptrv bonded_neighbors;
			
			VINA_FOR_IN(j, a->bonds){
                                const bond_vc& b = a->bonds[j];
				sz neighbor_index = b.connected_atom_index.i;
                                atom_vc* a2 = model_atoms[neighbor_index];
				bonded_neighbors.push_back(a2);
			}

			for (sz j = i+1; j < fused_atoms.size(); j++){
				atom_vc* a2 = fused_atoms[j];
				if (has(bonded_neighbors, a2)){
					sz a2_index = fused_atom_indices[j];
					//std::cout << "Detected fused edge " << a->atomname << " and " << a2->atomname << std::endl;
					this->fused_edges.push_back(std::pair<sz, sz>(a_index, a2_index));
				}
			}
		}
	}

	bool is_fused_edge(sz i, sz j){
		VINA_FOR_IN(a, this->fused_edges){
			std::pair<sz, sz>& edge = this->fused_edges[a];
			sz& ei = edge.first; sz& ej = edge.second;

			if ((i == ei && j == ej) || (i == ej && j == ei)){
				return true;
			}
		}
		return false;
	}

	void ComputeCentroid(bool use_adjusted_coords){
		fl avg_x = 0, avg_y = 0, avg_z = 0;
                VINA_FOR_IN(i, this->all_atom_coord_ptrs){
			sz this_atom_index = this->all_atom_indices[i];
                        vec* coord = (use_adjusted_coords) ? (*this->adjusted_coords)[this_atom_index] : this->all_atom_coord_ptrs[i];

                        fl& x = coord->data[0]; fl& y = coord->data[1]; fl& z = coord->data[2];
                        avg_x += x; avg_y += y; avg_z += z;
                }

		fl num_ring_atoms_fl = (fl) num_ring_atoms;
                avg_x /= num_ring_atoms_fl; avg_y /= num_ring_atoms_fl; avg_z /= num_ring_atoms_fl;
               	this->centroid.data[0] = avg_x; this->centroid.data[1] = avg_y; this->centroid.data[2] = avg_z;
	
		VINA_FOR_IN(i, this->subcycle_coords){
			vptrv& this_cycle_coords = this->subcycle_coords[i];
			fl avg_x = 0, avg_y = 0, avg_z = 0;

			VINA_FOR_IN(j, this_cycle_coords){
				vec* coord = this_cycle_coords[j];
				fl& x = coord->data[0]; fl& y = coord->data[1]; fl& z = coord->data[2];
				avg_x += x; avg_y += y; avg_z += z;
			}

			fl num_cycle_atoms = (fl) this_cycle_coords.size();
			avg_x /= num_cycle_atoms; avg_y /= num_cycle_atoms; avg_z /= num_cycle_atoms;
			vec& this_centroid = this->subcycle_centroids[i];

			this_centroid[0] = avg_x; this_centroid[1] = avg_y; this_centroid[2] = avg_z;
		}
		return;
	}

	void ComputeNormal(bool use_adjusted_coords){
		sz index0 = this->all_atom_indices[0];
		sz index1 = this->all_atom_indices[1];
		sz index2 = this->all_atom_indices[2];

		vec *a1c = (use_adjusted_coords) ? (*this->adjusted_coords)[index0] : this->all_atom_coord_ptrs[0];
		vec *a2c = (use_adjusted_coords) ? (*this->adjusted_coords)[index1] : this->all_atom_coord_ptrs[1];
		vec *a3c = (use_adjusted_coords) ? (*this->adjusted_coords)[index2] : this->all_atom_coord_ptrs[2];

                fl V1x = a2c->data[0] - a1c->data[0], V1y = a2c->data[1] - a1c->data[1], V1z = a2c->data[2] - a1c->data[2];
                fl V2x = a3c->data[0] - a2c->data[0], V2y = a3c->data[1] - a2c->data[1], V2z = a3c->data[2] - a2c->data[2];
                this->normal.data[0] = V1y*V2z - V1z*V2y; this->normal.data[1] = V1z*V2x - V1x*V2z; this->normal.data[2] = V1x*V2y - V1y*V2x;
		normalize_vec_in_place(this->normal);
		return;
	}

	void ComputeLargestAtomCentroidDistance(bool use_adjusted_coords){
		fl largest_distance = 0.00;
                VINA_FOR_IN(i, this->all_atom_coord_ptrs){
			sz atom_index = this->all_atom_indices[i];
                        vec* coord = (use_adjusted_coords) ? (*this->adjusted_coords)[atom_index] : this->all_atom_coord_ptrs[i];
                        fl dx = this->centroid.data[0] - coord->data[0];
                        fl dy = this->centroid.data[1] - coord->data[1];
                        fl dz = this->centroid.data[2] - coord->data[2];
                        fl r = std::sqrt(dx*dx + dy*dy + dz*dz);
			if (r > largest_distance) largest_distance = r; 
                }
                this->max_atom_centroid_dist =  largest_distance;
		return;
	}

	void ComputeRemainingAttributes(bool use_adjusted_coords){
		
		this->ComputeCentroid(use_adjusted_coords);	
		this->ComputeNormal(use_adjusted_coords);
		this->ComputeLargestAtomCentroidDistance(use_adjusted_coords);
		return;
	}
	
	void ComputeCentroidAndNormal(bool use_adjusted_coords){
		this->ComputeCentroid(use_adjusted_coords);
		this->ComputeNormal(use_adjusted_coords);
		return;
	}

	aptrv all_atom_ptrs;
	szv all_atom_indices;
	vptrv all_atom_coord_ptrs; 

	aptrv aromatic_carbon_ptrs;
	szv aromatic_carbon_indices;
	vptrv aromatic_carbon_coord_ptrs;

	aptrv heteroatom_ptrs;
	szv heteroatom_indices;
	vptrv heteroatom_coord_ptrs;

	sz num_ring_atoms;
	sz num_aromatic_carbons;
	sz num_heteroatoms;

	vec centroid;
	vec normal;
	fl max_atom_centroid_dist; //Largest atom-centroid distance. 
	bool is_ligand;

 	aptrv ring_hydrogens;
        szv ring_h_indices;
        vptrv ring_h_coords;
        szv ring_h_heavy_neighbor_indices;
	vptrv ring_h_heavy_neighbor_coords;
	sz num_ring_hydrogens;

	vptrv* adjusted_coords;
	std::vector<std::pair<sz,sz> > fused_edges;
	model* m;

	std::vector<aptrv> subcycles;
	std::vector<vptrv> subcycle_coords;
	vecv subcycle_centroids;
};
typedef std::vector<ring_attribute> ring_info;

struct aliphatic_carbon_attribute {
	aliphatic_carbon_attribute(atom_vc* carbon_, sz& carbon_index_, vec* c_coord_, aptrv& h_neighbors_, szv& h_neighbor_indices_, vptrv& h_coords_, bool chpi_explicit_hydrogen, bool is_ligand_){
		VINA_CHECK(h_neighbors_.size() == h_coords_.size());

		this->carbon = carbon_;
		this->carbon_atom_index = carbon_index_;
		this->c_coord = c_coord_;
		this->h_neighbors = h_neighbors_;
		this->h_neighbor_indices = h_neighbor_indices_;
		this->h_coords = h_coords_;
		this->num_h_neighbors = h_coords_.size();
		this->explicit_hydrogen = chpi_explicit_hydrogen;
		this->is_ligand = is_ligand_;
	}

	void ComputeCHBondVectors(sz index, vec& bond){
		VINA_CHECK(this->explicit_hydrogen);
		VINA_CHECK(index < this->h_coords.size());

		vec* h_coord = this->h_coords[index];
		bond.data[0] = h_coord->data[0] - this->c_coord->data[0];
		bond.data[1] = h_coord->data[1] - this->c_coord->data[1];
		bond.data[2] = h_coord->data[2] - this->c_coord->data[2];
		return;
	}

	sz carbon_atom_index;
	atom_vc* carbon;
	vec* c_coord;
	aptrv h_neighbors;
	szv h_neighbor_indices;
	vptrv h_coords;
	sz num_h_neighbors;
	bool explicit_hydrogen;
	bool is_ligand;
};
typedef std::vector<aliphatic_carbon_attribute> aliphatic_carbon_info;

enum distance_type {DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_VARIABLE};
typedef strictly_triangular_matrix<distance_type> distance_type_matrix;

struct non_cache; // forward declaration
struct naive_non_cache; // forward declaration
struct cache; // forward declaration
struct szv_grid; // forward declaration
struct terms; // forward declaration
struct conf_independent_inputs; // forward declaration
struct pdbqt_initializer; // forward declaration - only declared in parse_pdbqt.cpp
struct model_test;

struct model {
	void append(const model& m);
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }

	sz num_movable_atoms() const { return m_num_movable_atoms; }
	sz num_internal_pairs() const;
	sz num_other_pairs() const { return other_pairs.size(); }
	sz num_ligands() const { return ligands.size(); }
	sz num_flex() const { return flex.size(); }
	sz ligand_degrees_of_freedom(sz ligand_number) const { return ligands[ligand_number].degrees_of_freedom; }
	sz ligand_longest_branch(sz ligand_number) const;
	sz ligand_length(sz ligand_number) const;

	szv get_movable_atom_types(atom_type::t atom_typing_used_) const;

	conf_size get_size() const;
	conf get_initial_conf() const; // torsions = 0, orientations = identity, ligand positions = current

	grid_dims movable_atoms_box(fl add_to_each_dimension, fl granularity = 0.375) const;

	void write_flex  (                  const path& name, const std::string& remark) const { write_context(flex_context, name, remark); }
	void write_ligand(sz ligand_number, const path& name, const std::string& remark) const { VINA_CHECK(ligand_number < ligands.size()); write_context(ligands[ligand_number].cont, name, remark); }
	void write_structure(ofile& out) const {
		VINA_FOR_IN(i, ligands)
			write_context(ligands[i].cont, out);
		if(num_flex() > 0) // otherwise remark is written in vain
			write_context(flex_context, out);
	}
	void write_structure(ofile& out, const std::string& remark) const {
		out << remark;
		write_structure(out);
	}
	void write_structure(const path& name) const { ofile out(name); write_structure(out); }
	void write_model(ofile& out, sz model_number, const std::string& remark) const {
		out << "MODEL " << model_number << '\n';
		write_structure(out, remark);
		out << "ENDMDL\n";
	}
	void seti(const conf& c);
	void sete(const conf& c);
	void set (const conf& c);

	fl gyration_radius(sz ligand_number) const; // uses coords

	const atom_base& movable_atom  (sz i) const { assert(i < m_num_movable_atoms); return  atoms[i]; }
	const vec&       movable_coords(sz i) const { assert(i < m_num_movable_atoms); return coords[i]; }

	const vec& atom_coords(const atom_index& i) const;
	fl distance_sqr_between(const atom_index& a, const atom_index& b) const;
	bool atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const; // there is an atom closer to both a and b then they are to each other and immobile relative to them
	bool is_aromatic_heteroatom_element(sz& el);
	bool is_aromatic_cycle(aptrv& ring);
	bool is_trigonal_planar(atom_vc* atom);

	distance_type distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const;

	// clean up
	fl evali     (const precalculate& p,                  const vec& v                          ) const;
	fl evale     (const precalculate& p, const igrid& ig, const vec& v                          ) const;
	fl eval      (const precalculate& p, const igrid& ig, const vec& v, const conf& c           );
	double get_torsion_coords_vecs_list(vec A, vec B, vec C, vec D);
	double phi_alpha_energy(double phi_angle);
	double phi_beta_energy(double phi_angle);
	double psi_2A3E_energy(double psi_angle);
	double psi_2E3A_energy(double psi_angle);
	double psi_6A_energy(double psi_angle);
        double psi_6E_energy(double psi_angle);
        double omega_6A_energy(double omega_angle);
        double omega_6E_energy(double omega_angle);

	fl eval_deriv(const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g, const fl chi_coeff, const fl chi_cutoff);

	fl eval_intramolecular(                            const precalculate& p,                  const vec& v, const conf& c);
	fl eval_internal_torsional(const precalculate& p, const vec& v, const conf& c);
	fl eval_chi(const fl chi_coeff, const fl chi_cutoff);
	fl eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy, bool score_in_place);


	fl rmsd_lower_bound(const model& m) const; // uses coords
	fl rmsd_upper_bound(const model& m) const; // uses coords
	fl rmsd_ligands_upper_bound(const model& m) const; // uses coords

	void verify_bond_lengths() const;
	void about() const;

	vecv get_ligand_internal_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			tmp.push_back(internal_coords[i]);
		return tmp;
	}
	vecv get_flexible_coords();
	vecv get_ligand_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			{
			tmp.push_back(coords[i]);
			}
		return tmp;
	}
	vecv get_heavy_atom_movable_coords() const { // FIXME mv
		vecv tmp;
		VINA_FOR(i, num_movable_atoms())
			if(atoms[i].el != EL_TYPE_H)
				tmp.push_back(coords[i]);
		return tmp;
	}
	void check_internal_pairs() const;
	void print_stuff() const; // FIXME rm

	fl clash_penalty() const;
 
        bool bonded_to_HD(const atom_vc& a) const; //Yao: make this public, used to be private.
        bool bonded_to_heteroatom(const atom_vc& a) const; //Yao: make this public, used to be private. 
	bool bonded_to_hydrogen(const atom_vc& a) const;
	bool chpi_contact_possible(const atom_vc& central_carbon, vec& ring_centroid, vec& ring_normal, aptrv& ring, fl& ring_effective_radius, bool ring_is_receptor);
	bool is_atom_aromatic(atom_vc* atom);

	void build_ar_ring_info(); //Yao added 2023602, find receptor aromatic rings, pre-compute ring centroids and normals. 
        void DFSVisit(aptrv& atoms, std::map<sz, bool>& atom_visited_map, std::map<sz, int>& atom_parent_map, std::map<sz, atom_vc*>& index_atom_map, sz atom_index, int& counter, aptrv& current_path, std::vector<aptrv>& detected_cycles);

	sz num_r_ali_carbs;
	sz num_l_ali_carbs;
	sz num_r_aro_rings;
	sz num_l_aro_rings;

	fl eval_chpi(bool score_in_place);
	void eval_chpi_c(pr& dH_minusTdS, bool score_in_place);
	void eval_chpi_h(pr& dH_minusTdS, bool score_in_place);
	void eval_chpi_c_ring(aliphatic_carbon_attribute& c_attr, ring_attribute& r_attr, pr& dH_minusTdS, bool score_in_place);
	void eval_chpi_h_ring(aliphatic_carbon_attribute& c_attr, ring_attribute& r_attr, pr& dH_minusTdS, bool score_in_place);
	int choose_interacting_h(std::map<sz, fl>& h_centroid_dist, szv& ip_hs);
	fl get_chpi_scaling_factor(sz num_h_neghbors, sz num_ip_h);

	fl eval_chpi_enthalpy_c(fl r);
	fl eval_chpi_enthalpy_h(fl r);
	fl eval_chpi_enthalpy_h2(fl r);
	pr eval_chpi_enthalpy_c_deriv(fl r2);
	pr eval_chpi_enthalpy_h_deriv(fl r2, fl chn_cos_abs);
	//pr eval_chpi_enthalpy_h_deriv_r(fl r, fl chn_cos_abs, sz num_interacting_h);
	fl eval_chpi_entropy(fl horizontal_offset);
	pr eval_chpi_entropy_deriv(fl r2h);
	pr eval_chpi_entropy_deriv_r(fl rh);
	fl eval_chpi_entropy_log_deriv(fl rh);
	pr eval_repulsion_deriv(fl r2, fl cutoff);
	fl eval_chpi_enthalpy_gaussian_deriv(fl r);
	fl eval_chpi_enthalpy_gaussian_deriv2(fl r);
	fl eval_gaussian_deriv(fl r, fl miu, fl inv_sigma_sqr, fl gaussian_value);
	
	fl eval_repulsion(fl r, fl cutoff);

	pr eval_chpi_enthalpy_c_fast(fl r2); //For all fast's, pr.second = 0 (no derivative).
	pr eval_chpi_enthalpy_h_fast(fl r2, fl chn_cosine_abs);
	//pr eval_chpi_enthalpy_h_fast_r(fl r, fl chn_cosine_abs, sz num_interacting_h);
	pr eval_chpi_entropy_fast(fl r2h); //Notice that I'm using vertical distance itself, but the sqr of horizontal distance. 
	pr eval_chpi_entropy_fast_r(fl r2h); //Notice that I'm using vertical distance itself, but the sqr of horizontal distance. 
	pr eval_repulsion_fast(fl r2, fl cutoff);

	pr eval_gauss1_deriv(fl r2);
	pr eval_gauss2_deriv(fl r2);
	pr eval_phobic_deriv(fl r2);
	pr eval_gauss1_fast(fl r2);
	pr eval_gauss2_fast(fl r2);
	pr eval_phobic_fast(fl r2);

	std::vector<aptrv> remove_duplicate_cycles(std::vector<aptrv>& cycles);
	std::vector<aptrv> remove_redundant_cycles(std::vector<aptrv>& cycles, std::map<sz, szv>& subcycles, szv& returned_cycle_indices);

	void DetectAliphaticCarbons();
        std::vector<aptrv> DetectCyclesByDFS(std::map<sz, atom_vc*>& index_atom_map);
        std::vector<std::pair<aptrv, std::vector<aptrv> > > DetectAromaticCycles(std::map<sz, atom_vc*>& index_atom_map);
	void DetectAromaticCycles(resv& residues, ring_info& ring_attributes, bool is_ligand);
	void build_residue_info();
	void build_residues_from_atoms(resv& residues, aptrv& atoms, bool is_ligand);
	void calc_smooth_i1_i2_rem(sz& i1, sz& i2, fl& rem, fl r2);
	void calc_smooth_i1_i2_rem_dense(sz& i1, sz& i2, fl& rem, fl r2);
	void adjust_out_of_box_coords();
	fl GetDihedral(vec* p1, vec* p2, vec* p3, vec* p4);//Get dihedral 1-2-3-4. 
	void eval_chpi_entropy_set_force(vec* h_closest, sz closest_h_i, sz h_closest_index, std::vector<flv>& h_ring_dists, ring_attribute& r, pr& dH_minusTdS, bool ligand_aliphatic);
	void eval_chpi_entropy_set_force2(vec* h, sz closest_h_i, sz h_closest_index, std::vector<flv>& h_ring_dists, ring_attribute& r, pr& dH_minusTdS, bool ligand_aliphatic);

	void eval_chpi_entropy_set_force_old(vec* h_closest, sz h_closest_index, ring_attribute& r, pr& dH_minusTdS, bool ligand_aliphatic);
	void eval_chpi_entropy_set_force_old_each_centroid(vec* h_closest, sz h_closest_index, ring_attribute& r, vec& centroid, pr& dH_minusTdS, bool ligand_aliphatic);
	vec* choose_closest_centroid(vec* h, vec& centroid, vecv& centroids);
	pr calc_horizontal_factor(fl rh);

	//Yao added 20230602
        //Could contain multiple ligands. Each vector of aptrv is all the aromatic rings of that ligand.
	aliphatic_carbon_info lig_ali_carb_info;//cleared
	aliphatic_carbon_info rec_ali_carb_info;//cleared
	ring_info lig_ar_ring_info;//cleared
	ring_info rec_ar_ring_info;//cleared

	resv receptor_residues;//cleared
	resv ligand_residues;//cleared

	aptrv grid_atom_ptrs;//cleared
	aptrv atom_ptrs;//cleared
	aptrv all_atom_ptrs;
	vptrv all_atom_coords;
	fl weight_chpi;
	fl weight_repulsion;
	fl weight_gauss1;
	fl weight_gauss2;
	fl weight_hydrophobic;
	bool chpi_explicit_hydrogen;
	void undo_attraction_scoring(fl& score, bool score_in_place);
	void undo_attraction_scoring_each_pair(aliphatic_carbon_attribute& c, ring_attribute& r, fl& e_gau1, fl& e_gau2, fl& e_phobic, bool score_in_place);

	void build_chpi_smooth();
	prv chpi_c_smooth;
	prv chpi_h_smooth;
	prv chpi_entropy_smooth;
	prv gauss1_smooth;
	prv gauss2_smooth;
	prv phobic_smooth;

	flv chpi_c_fast;
	flv chpi_h_fast;
	flv chpi_entropy_fast;
	flv gauss1_fast;
	flv gauss2_fast;
	flv phobic_fast;

	fl prec_factor;
	prv gd_dims;

private:
	friend struct non_cache;
	friend struct naive_non_cache;
	friend struct cache;
	friend struct szv_grid;
	friend struct terms;
	friend struct conf_independent_inputs;
	friend struct appender_info;
	friend struct pdbqt_initializer;
	friend struct model_test;
	friend struct parallel_mc_task; 

	model() : m_num_movable_atoms(0), m_atom_typing_used(atom_type::XS) {};

	const atom_vc& get_atom(const atom_index& i) const { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }
	      atom_vc& get_atom(const atom_index& i)       { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }

	void write_context(const context& c, ofile& out) const;
	void write_context(const context& c, ofile& out, const std::string& remark) const {
		out << remark;
	}
	void write_context(const context& c, const path& name) const {
		ofile out(name);
		write_context(c, out);
	}
	void write_context(const context& c, const path& name, const std::string& remark) const {
		ofile out(name);
		write_context(c, out, remark);
	}
	fl rmsd_lower_bound_asymmetric(const model& x, const model& y) const; // actually static
	
	atom_index sz_to_atom_index(sz i) const; // grid_atoms, atoms
	//bool bonded_to_HD(const atom_vc& a) const; //Yao: make it public for now
	//bool bonded_to_heteroatom(const atom_vc& a) const; //Yao: make it public for now
	sz find_ligand(sz a) const;
	void bonded_to(sz a, sz n, szv& out) const;
	szv bonded_to(sz a, sz n) const;

	void assign_bonds(const distance_type_matrix& mobility); // assign bonds based on relative mobility, distance and covalent length
	void assign_types();
	void initialize_pairs(const distance_type_matrix& mobility);
	void initialize(const distance_type_matrix& mobility);
	fl clash_penalty_aux(const interacting_pairs& pairs) const;

	flv calculate_rs(fl sqr_cutoff);
	flv calculate_rs_dense();

	vecv internal_coords;
	vecv coords;
	vecv adjusted_coords; //Yao added 20231029
	vptrv adjusted_coord_ptrs;
	vecv minus_forces;

	atomv grid_atoms;
	atomv atoms; // movable, inflex
	vector_mutable<ligand> ligands;  //type declarations (?) 
	vector_mutable<residue_vc> flex;
	context flex_context;
	interacting_pairs other_pairs; // all except internal to one ligand: ligand-other ligands; ligand-flex/inflex; flex-flex/inflex

	sz num_rs;
	sz num_rs_c, num_rs_h, num_rs_e, num_rs_r, num_rs_r_h;
	sz imax;
	szv adjusted_atom_indices;

	sz m_num_movable_atoms;
	atom_type::t m_atom_typing_used;
};

#endif
