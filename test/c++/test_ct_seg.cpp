/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by M. Ferrero, O. Parcollet, H.Hafermann, T.Ayral
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#include <triqs/gf/imtime.hpp>
#include <triqs/gf/block.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/parameters/parameters.hpp>
//#include <triqs/python_tools/array_interface.hpp>
#include <boost/mpi.hpp>
//#include <triqs/arrays.hpp>
namespace triqs { namespace applications { namespace impurity_solvers { namespace ctqmc_seg {

 //using namespace gf; //?
 // THIS IS UGLY ! 
 using triqs::gfs::gf; using triqs::gfs::gf_view; using triqs::gfs::imtime; using triqs::gfs::block;
 using triqs::gfs::make_gf; using triqs::gfs::Fermion; using triqs::gfs::Boson; using triqs::gfs::half_bins;

}}}}

#include "applications/impurity_solvers/ctqmc_seg/operator_maps.hpp"
#include "applications/impurity_solvers/ctqmc_seg/configuration.hpp"
#include "applications/impurity_solvers/ctqmc_seg/qmc_data.hpp"
#include "applications/impurity_solvers/ctqmc_seg/ctqmc_seg.hpp"
#include "applications/impurity_solvers/ctqmc_seg/move_insert_segment.hpp"
#include "applications/impurity_solvers/ctqmc_seg/move_insert_spin_segment.hpp"
#include "applications/impurity_solvers/ctqmc_seg/move_remove_segment.hpp"
#include "applications/impurity_solvers/ctqmc_seg/move_remove_spin_segment.hpp"
using namespace triqs::applications::impurity_solvers::ctqmc_seg;
using namespace triqs;
//using namespace triqs::gf;



void test_operator_maps(){


 //create operator_maps from scratch

 operator_maps c(2);
 operator_maps::flavored_const_iterator it3;// bool ok;
 try{ it3 = c.insert_operator(utility::time_pt(3.0,10.0),0,true, c.end(0));}
 catch(insertion_error const & e){ std::cerr << "Insertion error " << std::endl;}
 /*std::cout << "ok = " << ok << std::endl;
   std::tie(it3, ok) =c.insert_operator(3.1,0,true, c.end(0));
   std::cout << "ok = " << ok << std::endl;
   std::tie(it3, ok) =c.insert_operator(2.0,0,true, c.end(0));
   std::cout << "ok = " << ok << std::endl;
   std::tie(it3, ok) =c.insert_operator(3.0,1,false, c.end(1));
   std::cout << "ok = " << ok << std::endl;//should fail*/
 /*std::tie(it3, ok) =c.insert_operator(1.0,1,false);
   std::cout << "ok = " << ok << std::endl;
   std::tie(it3, ok) = c.insert_operator(2.5,1,false);
   std::cout << "ok = " << ok << std::endl;*/
 //std::cout << "printing" <<std::endl;
 //std::cout << c << std::endl;
 //std::cout << "erasing" <<std::endl;
 //c.erase_operator(it3);
 //std::cout << "printing" <<std::endl;
 std::cout << c << std::endl;

}
void test_configuration(){

 //test configuration
 double beta=10;
 configuration config(2, beta) ;

 //test insert segment

 qmc_time_t tau_1(3.0,beta),tau_2(2.0,beta);
 try{  auto seg = config.insert_segment(tau_1,tau_2,0,config.begin(0)); }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 


  qmc_time_t tau_3(4.5,beta),tau_4(3.5,beta);
 try{  auto seg = config.insert_segment(tau_3,tau_4,0,config.begin(0)); }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 
/*
 try{
  seg =  config.insert_segment(3.1,2.0,1,config.operator_maps().begin(1));//fails
  std::cout << "2. try insert 3.1 2.0 line 1 ok ="  << " [ko]" <<std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 

 try{
  seg =  config.insert_segment(3.0,2.0,0,config.operator_maps().begin(0));//fails
  std::cout << "3. try insert 3.0 2.0 line 0 ok ="  << "[ko]"<<std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 

 try{
  seg =  config.insert_segment(5.0,2.5,1,config.operator_maps().begin(0));
  std::cout << "4. try insert 5.0 2.5 line 1 ok ="  << std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 


 try{
  seg =  config.insert_segment(8.0,5.5,0,config.operator_maps().begin(0));
  std::cout << "5. try insert 8.0 5.5 line 0 ok ="  << std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 

 std::cout << "Removing 8.0 5.5"<<std::endl;
 config.remove_segment(seg);

 try{
  seg =  config.insert_segment(1.0,6.0,1,config.operator_maps().begin(0));
  std::cout << "6. Try insert 1.0 6.0 line 1  ok ="  << std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 

 try{
  seg2 =  config.insert_segment(.5,7.0,1,config.operator_maps().begin(0));
  std::cout << "7. try insert .5 7.0 line 1 ok ="  << std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 



 //test remove segment
 std::cout << "8. test remove" << std::endl;
 config.remove_segment(seg2);
 std::cout << config << std::endl; 

 //test overlap
 std::cout << "8. test overlap" << std::endl;
 double length;triqs::arrays::vector<double> overlap;
 length=config.length(seg);
 overlap=config.overlaps(seg);
 std::cout << "\t length ["<< seg.fit_l->tau <<","<< config.operator_maps().cyclic_right(seg.fit_l)->tau<<"] = " << length << std::endl;
 std::cout << "\t overlap = " << overlap << std::endl;


 try{
  seg =  config.insert_segment(8.0,5.5,0,config.operator_maps().begin(0));
  std::cout << "ok ="  << std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 


 length=config.length(seg);
 overlap=config.overlaps(seg);
 std::cout << "length ["<< seg.fit_l->tau <<","<< config.operator_maps().cyclic_right(seg.fit_l)->tau<<"] = " << length << std::endl;
 std::cout << "overlap = " << overlap << std::endl;


 try{
  seg =  config.insert_segment(.5,7.0,1,config.operator_maps().begin(1));//hint is annihilation->antiseg
  std::cout << "ok ="  << std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 
 length=config.length(seg);
 overlap=config.overlaps(seg);
 std::cout << "length ["<< seg.fit_l->tau <<","<< config.operator_maps().cyclic_right(seg.fit_l)->tau<<"] = " << length << std::endl;
 std::cout << "overlap = " << overlap << std::endl;
 //test maximal_length

 double lmax;operator_maps::flavored_const_iterator it4;
 std::tie(lmax,it4) =  config.maximal_length(1.5, 0);
 std::cout << "lmax = " << lmax <<", right_op = " << it4->tau << " on line " << it4->flavor << std::endl;


 try{
  seg =  config.insert_segment(1.5,9.5,0,it4);
  std::cout << "ok ="  << std::endl;
 }
 catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
 std::cout << config << std::endl; 
 length=config.length(seg);
 overlap=config.overlaps(seg);
 std::cout << "length ["<< seg.fit_l->tau <<","<< config.operator_maps().cyclic_right(seg.fit_l)->tau<<"] = " << length << std::endl;
 std::cout << "overlap = " << overlap << std::endl;
 */
}
/*
void test_ordered_links(){


  ordered_links links;
  links.insert_raw(std::make_pair(0,0));
  links.insert_raw(std::make_pair(1,1));
  links.insert_raw(std::make_pair(2,2));
  links.insert_raw(std::make_pair(3,3));

  std::cout <<links << std::endl;

  std::cout << "now, inserting (1,1)" << std::endl;
  links.insert(std::make_pair(1,1));
  std::cout <<links << std::endl;
  std::cout << "[(0, 0) - (1, 1) - (2, 2) - (3, 3) - (4, 4) -]"<<std::endl;

  std::cout << "now, removing (1,1)" << std::endl;
  links.remove(std::make_pair(1,1));
  std::cout <<links << std::endl;
  std::cout << "[(0, 0) - (1, 1) - (2, 2) - (3, 3) -]"<<std::endl;



  ordered_links links2;
  links2.insert_raw(std::make_pair(1,0));
  links2.insert_raw(std::make_pair(0,1));
  links2.insert_raw(std::make_pair(2,2));
  links2.insert_raw(std::make_pair(3,3));

  std::cout <<links2 << std::endl;

  std::cout << "now, inserting (1,1)" << std::endl;

  links2.insert(std::make_pair(1,1));

  std::cout <<links2 << std::endl;
  std::cout << "[(0, 2) - (1, 1) - (2, 0) - (3, 3) - (4, 4) -]"<<std::endl;

  
}*/

void test_qmc_data(){

 //test qmc_data
 int n_flavor_=2;
 double beta_=10;
 arrays::matrix<double> U(2,2);U(0,0)=0.;U(1,1)=0.0;U(0,1)=1.0;U(1,0)=1.0;
 arrays::vector<double> mu(2);mu()=.5;

 auto delta_block = triqs::gfs::make_gf<triqs::gfs::imtime> (beta_, triqs::gfs::Fermion, arrays::make_shape(1,1));
 triqs::clef::placeholder<0> tau_;
 delta_block(tau_) << (-.5);//works
 //delta_block(tau_) << (-.5)+(tau_*(beta_-tau_));//works but wrong tail...
 //delta_block(tau_) << (-.5)*(exp(-1.0*tau_)+exp(-1.0*beta_+tau_));//does not work
 std::vector<triqs::gfs::gf_view<triqs::gfs::imtime> >  D; 
 D.push_back(delta_block); D.push_back(delta_block); 

 auto Delta = triqs::gfs::make_gf_view<triqs::gfs::block<triqs::gfs::imtime>> (D);
 bool dynamical_U_=false;
 triqs::gfs::gf<triqs::gfs::imtime> K_ ;//empty
 auto J_perp = triqs::gfs::make_gf<triqs::gfs::imtime> (beta_, triqs::gfs::Boson, arrays::make_shape(1,1));
 J_perp(tau_) << (-.5);
 qmc_data data(n_flavor_, beta_, U, mu, Delta, dynamical_U_, K_, J_perp); 

 //test move_insert
 triqs::mc_tools::random_generator RNG("mt19937", 2345);
 move_insert_segment move_insert(&data,RNG);
 move_remove_segment move_remove(&data,RNG);

 double prob,p;
 for(int n=0;n<200;n++){
  std::cout << "-------------------------------" << std::endl;
  prob = move_insert.attempt();
  std::cout << "prob = " << prob << std::endl;
  p = RNG(1.0);
  if(abs(prob)>p) move_insert.accept();
  else move_insert.reject();


  std::cout << "Config = " <<data.config << std::endl;


  std::cout << "-------------------------------" << std::endl;
  prob = move_remove.attempt();
  std::cout << "prob = " << prob << std::endl;
  p = RNG(1.0);
  if(abs(prob)>p) move_remove.accept();
  else move_remove.reject();
  std::cout << "Config = " <<data.config << std::endl;

 }


 std::cout << "######### Test move insert spin ########## " << std::endl;

 move_insert_spin_segment move_insert_spin(&data,RNG);
 move_remove_spin_segment move_remove_spin(&data,RNG);
 for(int n=0;n<1000;n++){
   std::cout << "--------insert spin seg-----------------------" << std::endl;
   prob = move_insert_spin.attempt();
   std::cout << "prob = " << prob << std::endl;
   p = RNG(1.0);
   if(abs(prob)>p) {std::cout << "ACCEPT " <<std::endl;move_insert_spin.accept();}
   else{std::cout << "REJECT " << std::endl; move_insert_spin.reject();}


   std::cout << "Config = " <<data.config << std::endl;

   std::cout << "----------remove spin seg---------------------" << std::endl;
   prob = move_remove_spin.attempt();
   std::cout << "prob = " << prob << std::endl;
   p = RNG(1.0);
   if(abs(prob)>p) {std::cout << "ACCEPT " <<std::endl;move_remove_spin.accept();}
   else{std::cout << "REJECT " << std::endl; move_remove_spin.reject();}
 }


}


void test_solver(int argc, char ** argv){


  boost::mpi::environment env(argc, argv);
  double beta_=10.;
  utility::parameters param;
  param["beta"]=beta_;


  ctqmc_seg solver(param);

  //solver.solve(param);//throws exception: missing required keys.

  arrays::array<double,1> mu(2);mu()=.5;
  //arrays::vector<double> mu(2);mu()=.5;
  arrays::array<double,2> U(2,2);U(0,1)=1.0;U(1,0)=1.0;
  //arrays::matrix<double> U(2,2);U(0,1)=1.0;U(1,0)=1.0;
  param["mu"] = mu;
  param["N_Cycles"]=100;
  param["Length_Cycle"]=1;
  param["U"] = U;

  //preparing hybridization function
  auto delta_block = triqs::gfs::make_gf<triqs::gfs::imtime> (beta_, triqs::gfs::Fermion, arrays::make_shape(1,1));
  triqs::clef::placeholder<0> tau_;
  delta_block(tau_) << (-.5);//works
  //delta_block(tau_) << (-.5)+(tau_*(beta_-tau_));//works but wrong tail...
  //delta_block(tau_) << (-.5)*(exp(-1.0*tau_)+exp(-1.0*beta_+tau_));//does not work
  std::vector<triqs::gfs::gf_view<triqs::gfs::imtime> >  D; 
  D.push_back(delta_block); D.push_back(delta_block); 

  auto sha1 = arrays::make_shape(1,1);      //fermionic quantities
  std::vector<std::string> block_names;
  block_names.push_back("up");
  block_names.push_back("dn");
  auto Delta = make_gf<triqs::gfs::block<triqs::gfs::imtime>>(block_names, make_gf<triqs::gfs::imtime>(beta_, triqs::gfs::Fermion, sha1, 10000, triqs::gfs::half_bins) );
  Delta = triqs::gfs::make_gf_view<triqs::gfs::block<triqs::gfs::imtime>> (D);
  //solver.deltat_view() = Delta;//crashes
  //param["deltat"]=Delta;//causes segfault


  solver.solve(param);
}


void test_spin_configuration(){

  //test configuration
  double beta=10;
  configuration config(2, beta) ;

  //test insert segment

  std::cout << "##1. Trying to insert" << std::endl;
  qmc_time_t tau_1(3.5,beta),tau_2(2.0,beta);
  configuration::spin_segment_const_iterator seg;
  try{  seg = config.insert_spin_segment(tau_1,tau_2,config.begin(0),config.begin(1)); }
  catch(insertion_error const & e){std::cerr << e.what() << std::endl;}
  std::cout << "c = " <<config << std::endl; 


  std::cout << "##2. Trying to insert" << std::endl;

  config.full_lines[0]=true;
  try{  seg = config.insert_spin_segment(tau_1,tau_2,config.begin(0),config.begin(1)); }
  catch(insertion_error const & e){std::cerr << e.what() << std::endl;}
  std::cout << "c = "<<config << std::endl; 
  
  //test remove segment
  std::cout << "##3. Trying to remove" << std::endl;
  config.remove_spin_segment(seg);
  std::cout << "c = " <<config << std::endl; 

  std::cout << "Full_lines = " << config.full_lines[0] << ", " << config.full_lines[1] << std::endl;


  std::cout <<"## insert segment on down" << std::endl;
  qmc_time_t tau_3(9.0,beta),tau_4(1.0,beta);
  configuration::hyb_segment rseg;
  try{  rseg = config.insert_segment(tau_3,tau_4,1,config.begin(1)); }
  catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
  std::cout << "rseg.{l,r} = " << double(rseg.l->tau) << "," << double(rseg.r->tau) << std::endl;
  std::cout << "c = " << config << std::endl<<std::endl; 


  std::cout << "## trying to insert spin" << std::endl;
  try{  seg = config.insert_spin_segment(tau_1,tau_2,config.begin(0),config.unhyb(rseg.l)); }
  catch(insertion_error const & e){std::cerr << e.what() << std::endl;}
  std::cout << config << std::endl;
  std::cout << "[(9,1,1:0,0,) - (3.5,0,1:1,0,) - (3.5,1,0:1,1,) - (2,1,1:1,0,) - (2,0,0:0,0,) - (1,1,0:0,1,) -]" << std::endl; 

  std::cout << "*** removing spin segment (3.5, 2)" <<std::endl;

  config.remove_spin_segment(seg);
  std::cout << "c = " <<config << std::endl;
  std::cout << "*** removing segment (9, 1)" <<std::endl;
  config.remove_segment(rseg);
  std::cout << "c = " <<config << std::endl;

  std::cout << "\n-- same thing, reversed -- " << std::endl;

  config.full_lines[0]=false;
  std::cout <<"## insert segment on down" << std::endl;
  //configuration::segment rseg;
  try{  rseg = config.insert_segment(tau_4,tau_3,1,config.begin(1)); }
  catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
  std::cout << "rseg.{l,r} = " << double(rseg.l->tau) << "," << double(rseg.r->tau) << std::endl;
  std::cout << "c=" <<config << std::endl<<std::endl; 


  std::cout << "## trying to insert spin" << std::endl;
  try{  seg = config.insert_spin_segment(tau_1,tau_2,config.begin(0),config.unhyb(rseg.r)); }
  catch(insertion_error const & e){std::cerr << e.what() << std::endl;}
  std::cout << config << std::endl;
  std::cout << "[(9,1,0:0,1,) - (3.5,1,1:0,0,) - (3.5,0,0:1,0,) - (2,0,1:0,0,) - (2,0,0:0,1,) - (1,1,1:0,0,) -]" << std::endl; 

  std::cout << "Removing spin " << std::endl;
  config.remove_spin_segment(seg);
  std::cout << "c =" <<config << std::endl;
  std::cout << "Removing segment " << std::endl;
  config.remove_segment(rseg);
  std::cout << "c ="<< config << std::endl;


  std::cout << "--- test when already operators on both lines ---" << std::endl;

  std::cout <<"## insert segment on down ({9,1})" << std::endl;
  configuration::hyb_segment dseg,useg;
  try{  dseg = config.insert_segment(tau_4,tau_3,1,config.begin(1)); }
  catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
  std::cout << "dseg.{l,r} = " << double(dseg.l->tau) << "," << double(dseg.r->tau) << std::endl;
  std::cout << "c ="<< config << std::endl;


  std::cout <<"## insert segment on up ({9.5,6})" << std::endl;
  tau_4=6.0;
  tau_3=9.5;
  try{  useg = config.insert_segment(tau_4,tau_3,0,config.begin(0)); }
  catch(insertion_error const & e){std::cerr << "Insertion error" << std::endl;}
  std::cout << "useg.{l,r} = " << double(useg.l->tau) << "," << double(useg.r->tau) << std::endl;
  std::cout << config << std::endl<<std::endl;

  std::cout << "## trying to insert spin" << std::endl;
  tau_1 = 3.;tau_2=2.;
  try{  seg = config.insert_spin_segment(tau_1,tau_2,config.unhyb(useg.l),config.unhyb(dseg.r)); }
  catch(insertion_error const & e){std::cerr << e.what() << std::endl;}
  std::cout << config << std::endl;
  std::cout << "[(9.5,0,0:1,0,) - (9,1,0:1,1,) - (6,0,1:0,1,) - (3,1,1:0,0,) - (3,0,0:1,0,) - (2,0,1:0,0,) - (2,1,0:0,1,) - (1,1,1:0,0,) -]"<<std::endl;


  config.remove_segment(dseg);
  config.remove_segment(useg);//leads to spurious right_occupations because no longer makes sense as a "segment" (->operators in between!)

  std::cout << config << std::endl;

  config.remove_spin_segment(seg);
  std::cout << config << std::endl;
}

int main(int argc, char ** argv){

  std::cout << "*********test operator_maps****************"<<std::endl;
  test_operator_maps();

  std::cout << "**********test configuration***************"<<std::endl;
  test_configuration();
  
  //std::cout << "**********test ordered links***************"<<std::endl;
  //test_ordered_links();

  std::cout << "************test spin insertion****************" << std::endl; 
  test_spin_configuration();

  std::cout << "***************test qmc_data*************" << std::endl; 
  //test_qmc_data();

  std::cout << "***********************test solver*****" << std::endl; 
  //test solver;
  test_solver(argc,argv);

}
