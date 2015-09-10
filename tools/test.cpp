
#include "./analysis.hpp"

// -------------------------------------------------------

int main() {

 calc C;

 // here comes the enumeration of the state we are interested in ...
 auto rvb_op = c_dag("00-up", 0) * c_dag("01-up", 0) * c_dag("11-up", 0) + c_dag("01-down", 0);

 auto rvb_state = C.act_on_vacuum(rvb_op);
 std::cout << "RVB state: " << rvb_state << std::endl;

 auto av = C.average_projector(rvb_state);

}
