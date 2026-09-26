// bdry_gen — btree's membrane shell for one boundary directive, from btree's own boundary_surface.cpp / vec_quat_manipulator.cpp
// (compiled with btree's flags so the std::unordered_set vertex order and the float ops are the same as in btree_chromo).
//   usage: bdry_gen <directive> <r_bdry> <out.bin>
//     directive: spherical_bdry:R,x0,y0,z0  |  overlapping_spheres_bdry:h,r,x0,y0,z0,u,v,w
//   writes n×3 float64 (btree order, translated as LAMMPS_sys::generate_*_bdry) to out.bin and prints on stdout the three
//   bbox lines exactly as LAMMPS_sys::write_data prints them for a system made of these atoms (calc_bbox, s = 1.2):
//     "<xlo>\t<xhi>" etc. with ostream default formatting.
#include <boundary_surface.hpp>
#include <vec_quat_manipulator.hpp>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
  if (argc != 4) { std::cerr << "usage: bdry_gen <directive> <r_bdry> <out.bin>" << std::endl; return 2; }
  std::string d = argv[1]; double r_bdry = std::stod(argv[2]);
  size_t c = d.find(':'); std::string kind = d.substr(0, c); std::vector<double> p;
  { std::stringstream ss(d.substr(c + 1)); std::string tok; while (std::getline(ss, tok, ',')) p.push_back(std::stod(tok)); }
  vec_quat_manipulator vqm; boundary_surface b; vec r0;
  if (kind == "spherical_bdry" && p.size() == 4) { b.generate_sphere(p[0], r_bdry); r0 = vqm.v_new(p[1], p[2], p[3]); }
  else if (kind == "overlapping_spheres_bdry" && p.size() == 8) { b.generate_overlapping_spheres(p[0], p[1], r_bdry, p[5], p[6], p[7]); r0 = vqm.v_new(p[2], p[3], p[4]); }
  else { std::cerr << "unsupported directive: " << d << std::endl; return 3; }
  std::vector<vec> co = b.get_coords();
  for (size_t i = 0; i < co.size(); i++) co[i] = vqm.v_xpy(co[i], r0);
  { std::ofstream o(argv[3], std::ios::binary); for (const vec &v : co) { double t[3] = {v.x, v.y, v.z}; o.write((const char *)t, sizeof t); } }
  // LAMMPS_sys::calc_bbox over these atoms
  vec mn = co[0], mx = co[0];
  for (size_t i = 1; i < co.size(); i++) {
    const vec &a = co[i];
    if (a.x > mx.x) mx.x = a.x; if (a.y > mx.y) mx.y = a.y; if (a.z > mx.z) mx.z = a.z;
    if (a.x < mn.x) mn.x = a.x; if (a.y < mn.y) mn.y = a.y; if (a.z < mn.z) mn.z = a.z;
  }
  vec r_mid = vqm.v_linterp(0.5, mn, mx); vec dr = vqm.v_xpy(mx, vqm.v_inv(mn)); dr = vqm.v_ax(0.5, dr);
  double s = 1.2; vec lo = vqm.v_axpy(-s, dr, r_mid), hi = vqm.v_axpy(s, dr, r_mid);
  std::cout << co.size() << "\n" << lo.x << "\t" << hi.x << "\n" << lo.y << "\t" << hi.y << "\n" << lo.z << "\t" << hi.z << std::endl;
  return 0;
}
