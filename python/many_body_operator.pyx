from many_body_operator cimport many_body_operator

cdef class Operator:
    """
    Many-body operator
    """
    #cdef many_body_operator[double,string,string] * _c
        
    def __cinit__(self):
        self._c = new many_body_operator[double,string,string]()
    def __dealloc__(self):
        del self._c
    
    def __repr__(self):
        cdef stringstream ss
        ss << self._c[0]
        return ss.str()
    def __str__(self):
        return self.__repr__()
        
    # Addition            
    def __iadd__(self, other):
        if isinstance(other,Operator):
            self._c[0] += (<Operator?> other)._c[0]
        else:
            self._c[0] += <double?> other
        return self
        
    def __add__(A, B):
        O = Operator()
        for x in (A,B):
            if isinstance(x,Operator):
                (<Operator?> O)._c[0] += (<Operator?> x)._c[0]
            else:
                (<Operator?> O)._c[0] += <double?> x
        return O
    
    # Subtraction
    def __isub__(self, other):
        if isinstance(other,Operator):
            self._c[0] -= (<Operator?> other)._c[0]
        else:
            self._c[0] -= <double?> other
        return self
        
    def __sub__(A, B):
        O = Operator()
        if isinstance(A,Operator):
            (<Operator?> O)._c[0] += (<Operator?> A)._c[0]
        else:
            (<Operator?> O)._c[0] += <double?> A
        if isinstance(B,Operator):
            (<Operator?> O)._c[0] -= (<Operator?> B)._c[0]
        else:
            (<Operator?> O)._c[0] -= <double?> B
        return O
    
    def __neg__(self):
        O = Operator()
        O._c[0] = -self._c[0]
        return O
    
    # Multiplication
    def __imul__(self, other):
        if isinstance(other,Operator):
            self._c[0] *= (<Operator?> other)._c[0]
        else:
            self._c[0] *= <double?> other
        return self
        
    def __mul__(A, B):
        O = Operator()
        if isinstance(A,Operator):
            if isinstance(B,Operator):
                (<Operator?> O)._c[0] = (<Operator?> A)._c[0] * (<Operator?> B)._c[0]
            else:
                (<Operator?> O)._c[0] = (<Operator?> A)._c[0] * (<double?> B)
        else:
            (<Operator?> O)._c[0] = (<Operator?> B)._c[0] * (<double?> A)
        return O
        
def C(block,index):
    O = Operator()
    O._c[0] = c(str(block),str(index))
    return O

def C_dag(block,index):
    O = Operator()
    O._c[0] = c_dag(str(block),str(index))
    return O

def N(block,index):
    O = Operator()
    O._c[0] = n(str(block),str(index))
    return O