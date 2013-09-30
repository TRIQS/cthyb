Generic many-body operators
===========================

.. Warning::

    Add information about C++/Cython header files and namespaces (once the header files have landed to the lib).

``many_body_operator`` is a template class, which implements the algebra of fermion operators.
An object of this class represents a general fermion operator and supports all standard algebraic operations (sums, products, multiplication by a scalar).
It allows to write readable and clean C++ code involving various operators, such as Hamiltonians and observables of many-body systems.

.. note::
   
    The internal storage of a ``many_body_operator`` object is not based on a matrix representation.
    Instead of that the object stores a list of normally-ordered monomials in basis elements
    (creation and annihilation operators), accompanied by scalar coefficients. This approach allows
    to minimize the required storage space.

``many_body_operator`` is defined as ::

    template<typename scalar_t, LessThanComparable... IndexTypes> class many_body_operator;

Template parameters
-------------------

    * **scalar_t** determines the type of a scalar to construct the algebra
    * **IndexTypes** is a parameter pack, which determines types of indices of a single creation/annihilation operator.
      The number of types may vary from 0 to the maximum supported by the compiler. Each type must model LessThanComparable.
      
.. note::

    **scalar_t** is often chosen once for a whole project. It is convenient to shorten declarations by using the template aliasing feature of C++11: ::
        
        template<typename ...IndexTypes> using op_with_double_coeff = many_body_operator<double, IndexTypes...>;
      
Construction/factories
----------------------

``many_body_operator`` provides a minimal set of constructors: default, copy and move constructor. The default constructor
creates a zero operator.

Three factory functions can be used to construct nontrivial operators: ::

     many_body_operator<scalar_t,IndexTypes...> Op1 = c(IndexTypes... indices); // annihilation
     many_body_operator<scalar_t,IndexTypes...> Op2 = c_dag(IndexTypes... indices); // creation
     many_body_operator<scalar_t,IndexTypes...> Op3 = n(IndexTypes... indices); // number of particles

Creation and annihilation operators obey the canonical anticommutation relation 

.. math::
    \hat C^\dagger_i \hat C_j + \hat C_j \hat C^\dagger_i = \delta_{ij},

and the number of particle is defined as

.. math::
    \hat n_i = \hat C^\dagger_i \hat C_i.
     
There is no need to preregister valid values of ``indices`` before they are used to create an elementary operator.
That means, that an algebra can be extended with new basis elements on-the-fly, after you have created some expressions.
     
Overloaded operators
--------------------

``many_body_operator`` class defines a number of arithmetic operations with objects of the class and constants of type ``scalar_t``.
If ``A`` and ``B`` are objects of class ``many_body_operator`` (instantiated with *the same template parameters*) and ``x`` is an
instance of ``scalar_t``, then the following expressions are valid: ::
    
    // Addition
    A + B
    A + x
    x + A
    A += B
    A += x
    
    // Subtraction
    A - B
    A - x
    x - A
    A -= B
    -A
    
    // Multiplication
    x*A
    A*x
    A*B

The result of any of the defined operations is guaranteed to preserve its normally ordered form.

Within the current implementation, ``many_body_operator`` provides no type conversions between objects with
different scalar types or index types.
For example, one cannot mix operators with ``scalar_t=double`` and ``scalar_t=std::complex<double>`` in a single expression.
Nevertheless, one operator can be assigned from another if their scalar types are compatible.

An instance of ``many_body_operator`` can be inserted into an output stream, provided that ``scalar_t`` and all index types support
the stream insertion. ::
    
    many_body_operator<double,int> x = c(0);
    many_body_operator<double,int> y = c_dag(1);
    
    std::cout << (x + y)*(x - y) << std::endl; // prints "2*C^+(1)C(0)"
    
Serialization & HDF5
--------------------

Objects of ``many_body_operator`` are ready to be serialized/deserialized with Boost.Serialization.
It also allows to transparently send/receive them through Boost.MPI calls.

There is no special code to support HDF5-storage of operators. However, the core HDF5 library of TRIQS
automatically performs text-serialization of an operator and puts the resulting string into the HDF5-archive.

Iteration over monomials
------------------------

The aim of ``many_body_operator`` is to have a class, which allows to encode different operator expressions in C++ in the form closest to the mathematical notation.
But at the same time, one would like to explicitly extract the structure of a defined operator (to calculate its matrix elements, for example).
For this purpose ``many_body_operator`` exposes the following part of its interface:

- ``struct canonical_ops_t``
    This structure represents an elementary operator (basis element of the algebra).
    ::

        struct canonical_ops_t { 
            bool dagger;    // true = creation, false = annihilation
            std::tuple<IndexTypes...> indices; // values of indices
        };

- ``typedef ... monomial_t;``
    An ordered sequence of elementary operators (monomial).

- ``typedef ... const_iterator;``
    A bidirectional constant iterator to the list of monomials.
    It can be dereferenced into a special proxy object, which carries two data members: ``coef`` and ``monomial``.

- ``begin()``/``cbegin()``
    Returns ``const_iterator`` pointing at the first monomial.

- ``end()``/``cend()``
    Returns ``const_iterator`` pointing past the end.

Here is an example of use: ::
    
    typedef many_body_operator<double,int> Op;
    Op H = -0.5*(n(0) + n(1)) + n(0)*n(1);
    
    for(Op::const_iterator it = H.begin(); it != H.end(); ++it){
        double coef = it->coef;
        Op::monomial_t monomial = it->monomial;
        
        std::cout << "Coefficient: " << coef << std::endl;
        std::cout << "Monomial: " << std::endl;
        for(auto const& o : monomial){
            std::cout << "dagger: " << o.dagger << " index: " << std::get<0>(o.indices) << " "; // only 1 index per elementary operator 
        }
        std::cout << std::endl;
    }
    
The output should be ::

    Coefficient: -0.5
    Monomial: 
    dagger: 1 index: 0 dagger: 0 index: 0 
    Coefficient: -0.5
    Monomial: 
    dagger: 1 index: 1 dagger: 0 index: 1 
    Coefficient: 1
    Monomial: 
    dagger: 1 index: 0 dagger: 1 index: 1 dagger: 0 index: 1 dagger: 0 index: 0 

Cythonization
-------------

This class can be cimported from Cython: ::
    
    from many_body_operator cimport many_body_operator
    
It corresponds to a specialized version of ``many_body_operator``: ``double`` as the scalar type and two indices of type ``std::string``.
There are also an extension type ``Operator`` and three factory function to be imported from Python: ::
    
    from many_body_operator import Operator, C, C_dag, N

All arithmetic operations implemented in C++ are also available in Python as well as special methods ``__repr__()`` and ``__str__()``.
The factory functions accept two arguments of any types and convert them into strings using Python operator ``str()``.

.. Warning::

    ``many_body_operator.pyx`` must be split into *.pxd* and *.pyx* parts when it is moved into the library.
    That was the initial intention, but there is an issue with how Cython compiler looks for header files.
    In short words, *.pxd*-files must reside at the same place of build and install directory hierarchies.
    But this is impossible to achieve with the standard directory layout of a TRIQS application.