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
    cdef many_body_operator[double,string,string] _c
        
    def __cinit__(self):
        self._c = many_body_operator[double,string,string]()
    
    def __repr__(self):
        cdef stringstream ss
        ss << self._c
        return ss.str()
    def __str__(self):
        return self.__repr__()
        
    # Addition            
    def __iadd__(self, other):
        if isinstance(other,Operator):
            self._c += (<Operator?> other)._c
        else:
            self._c += <double?> other
        return self
        
    def __add__(A, B):
        O = Operator()
        for x in (A,B):
            if isinstance(x,Operator):
                (<Operator?> O)._c += (<Operator?> x)._c
            else:
                (<Operator?> O)._c += <double?> x
        return O
    
    # Subtraction
    def __isub__(self, other):
        if isinstance(other,Operator):
            self._c -= (<Operator?> other)._c
        else:
            self._c -= <double?> other
        return self
        
    def __sub__(A, B):
        O = Operator()
        if isinstance(A,Operator):
            (<Operator?> O)._c += (<Operator?> A)._c
        else:
            (<Operator?> O)._c += <double?> A
        if isinstance(B,Operator):
            (<Operator?> O)._c -= (<Operator?> B)._c
        else:
            (<Operator?> O)._c -= <double?> B
        return O
    
    def __neg__(self):
        O = Operator()
        O._c = -self._c
        return O
    
    # Multiplication
    def __imul__(self, other):
        if isinstance(other,Operator):
            self._c *= (<Operator?> other)._c
        else:
            self._c *= <double?> other
        return self
        
    def __mul__(A, B):
        O = Operator()
        if isinstance(A,Operator):
            if isinstance(B,Operator):
                (<Operator?> O)._c = (<Operator?> A)._c * (<Operator?> B)._c
            else:
                (<Operator?> O)._c = (<Operator?> A)._c * (<double?> B)
        else:
            (<Operator?> O)._c = (<Operator?> B)._c * (<double?> A)
        return O
        
def C(block,index):
    O = Operator()
    O._c = c(str(block),str(index))
    return O

def C_dag(block,index):
    O = Operator()
    O._c = c_dag(str(block),str(index))
    return O

def N(block,index):
    O = Operator()
    O._c = n(str(block),str(index))
    return O
