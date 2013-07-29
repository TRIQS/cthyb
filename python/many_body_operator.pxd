from libcpp.string cimport string as string

cdef extern from "<sstream>" namespace "std":
    cdef cppclass stringstream:
        stringstream()
        string str()

cdef extern from "c++/operator.hpp" namespace "triqs::utility":
        
    cdef cppclass many_body_operator[double,string,string]:
        
        many_body_operator() except +

        many_body_operator operator+(const many_body_operator&)
        many_body_operator operator+(double)
        many_body_operator operator-(const many_body_operator&)
        many_body_operator operator-(double)
        many_body_operator operator-()
        many_body_operator operator*(const many_body_operator&)
        many_body_operator operator*(double)

    cdef many_body_operator c(string,string)
    cdef many_body_operator c_dag(string,string)
    cdef many_body_operator n(string,string)
    
    cdef void operator<<(stringstream&, many_body_operator&)
    
cdef class Operator:
    """
    Many-body operator
    """
    cdef many_body_operator[double,string,string] * _c