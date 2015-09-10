#include "./analysis.hpp"

// -------------------------------------------------------

int main() {

 calc C;

 auto c_dag_0_up = ( c_dag("00-up",0) + c_dag("10-up",0) + c_dag("01-up",0) + c_dag("11-up",0))/2;
 auto c_dag_1_up = ( c_dag("00-up",0) - c_dag("10-up",0) + c_dag("01-up",0) - c_dag("11-up",0))/2;
 auto c_dag_2_up = ( c_dag("00-up",0) - c_dag("10-up",0) - c_dag("01-up",0) + c_dag("11-up",0))/2;
 auto c_dag_3_up = ( c_dag("00-up",0) + c_dag("10-up",0) - c_dag("01-up",0) - c_dag("11-up",0))/2;

 auto c_dag_0_down = ( c_dag("00-down",0) + c_dag("10-down",0) + c_dag("01-down",0) + c_dag("11-down",0))/2;
 auto c_dag_1_down = ( c_dag("00-down",0) - c_dag("10-down",0) + c_dag("01-down",0) - c_dag("11-down",0))/2;
 auto c_dag_2_down = ( c_dag("00-down",0) - c_dag("10-down",0) - c_dag("01-down",0) + c_dag("11-down",0))/2;
 auto c_dag_3_down = ( c_dag("00-down",0) + c_dag("10-down",0) - c_dag("01-down",0) - c_dag("11-down",0))/2;

 auto ops_0hole = std::vector<many_body_operator<double>>{ 
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up)/2,
(c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up)/2
};
 
auto ops_1hole = std::vector<many_body_operator<double>>{ 
c_dag_0_up*(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up),
c_dag_0_down*(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up),
(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up)*c_dag_3_up,
(c_dag_1_up*c_dag_2_down-c_dag_1_down*c_dag_2_up)*c_dag_3_down,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_2_up,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_2_down,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_1_up,
(c_dag_0_up*c_dag_3_down-c_dag_0_down*c_dag_3_up)*c_dag_1_down,
 c_dag_0_up*(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up),
 c_dag_0_down*(c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up),
 (c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up)*c_dag_1_up,
 (c_dag_2_up*c_dag_3_down-c_dag_2_down*c_dag_3_up)*c_dag_1_down,
 (c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_2_up,
 (c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_2_down,
 (c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_3_up,
 (c_dag_0_up*c_dag_1_down-c_dag_0_down*c_dag_1_up)*c_dag_3_down
};

// n in space
 auto n0up = c_dag_0_up * dagger(c_dag_0_up);
 auto n1up = c_dag_1_up * dagger(c_dag_1_up);
 auto n2up = c_dag_2_up * dagger(c_dag_2_up);
 auto n3up = c_dag_3_up * dagger(c_dag_3_up);
 std::cout << "n0up = "<< C.average(n0up) << std::endl;
 std::cout << "n1up = "<< C.average(n1up) << std::endl;
 std::cout << "n2up = "<< C.average(n2up) << std::endl;
 std::cout << "n3up = "<< C.average(n3up) << std::endl;

 auto n0down = c_dag_0_down * dagger(c_dag_0_down);
 auto n1down = c_dag_1_down * dagger(c_dag_1_down);
 auto n2down = c_dag_2_down * dagger(c_dag_2_down);
 auto n3down = c_dag_3_down * dagger(c_dag_3_down);
 std::cout << "n0down = "<< C.average(n0down) << std::endl;
 std::cout << "n1down = "<< C.average(n1down) << std::endl;
 std::cout << "n2down = "<< C.average(n2down) << std::endl;
 std::cout << "n3down = "<< C.average(n3down) << std::endl;

 // n in k space (to compare with gf)
 auto n_00 =  triqs::utility::n("00-up", 0);
 auto n_10 =  triqs::utility::n("10-up", 0);
 auto n_01 =  triqs::utility::n("01-up", 0);
 auto n_11 =  triqs::utility::n("11-up", 0);
 std::cout << "<n00> = "<< C.average(n_00) << std::endl;
 std::cout << "<n10> = "<< C.average(n_10) << std::endl;
 std::cout << "<n01> = "<< C.average(n_01) << std::endl;
 std::cout << "<n11> = "<< C.average(n_11) << std::endl;

 // dimer states are NOT orthogonals.
 std::cout << "dimer prod" << dot_product(C.act_on_vacuum(ops_0hole[0]), C.act_on_vacuum(ops_0hole[0])) << std::endl;
 std::cout << "dimer prod" << dot_product(C.act_on_vacuum(ops_0hole[0]), C.act_on_vacuum(ops_0hole[1])) << std::endl;
 std::cout << "dimer prod" << dot_product(C.act_on_vacuum(ops_0hole[1]), C.act_on_vacuum(ops_0hole[1])) << std::endl;
 std::cout << "dimer prod" << dot_product(C.act_on_vacuum(ops_1hole[1]), C.act_on_vacuum(ops_0hole[1])) << std::endl;
 for (auto const & op : ops_1hole) 
 std::cout << "<op0|op>" << dot_product(C.act_on_vacuum(op), C.act_on_vacuum(ops_1hole[2])) << std::endl;

std::vector<double> av_0hole;
for (auto const & op : ops_0hole) av_0hole.push_back( C.average_projector(C.act_on_vacuum(op)));

std::vector<double> av_1hole;
for (auto const & op : ops_1hole) av_1hole.push_back( C.average_projector(C.act_on_vacuum(op/std::sqrt(2))));

std::cout<< std::endl;
for (int i = 0; i <2 ; ++i) std::cout<< i << "  "<< av_0hole[i] <<std::endl;
std::cout<< std::endl;
for (int i = 0; i <16 ; ++i) std::cout<< i << "  "<< av_1hole[i] <<std::endl;
}
