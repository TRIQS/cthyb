#include <pomerol/Misc.h>
#include <pomerol/Logger.h>
#include <pomerol/Lattice.h>
#include <pomerol/LatticePresets.h>
#include <pomerol/Index.h>
#include <pomerol/IndexClassification.h>
#include <pomerol/Operator.h>
#include <pomerol/OperatorPresets.h>
#include <pomerol/IndexHamiltonian.h>
#include <pomerol/Symmetrizer.h>
#include <pomerol/StatesClassification.h>
#include <pomerol/HamiltonianPart.h>
#include <pomerol/Hamiltonian.h>
#include <pomerol/FieldOperatorContainer.h>
#include <pomerol/GFContainer.h>

#include <triqs/gf/local/fourier_matsubara.hpp>
#include <triqs/gf/block.hpp>
#include <triqs/gf/imtime.hpp>
#include <triqs/gf/imfreq.hpp>

#include <cstdlib>

using triqs::gf::gf;
using triqs::gf::block_index;
using triqs::gf::imfreq;
using triqs::gf::imtime;
using triqs::gf::make_gf;
using triqs::gf::Fermion;

using namespace Pomerol;

RealType beta = 10;
RealType U = 2.0;
RealType mu = 1.0;
RealType h = 0.0;
RealType V = 1.0;
RealType epsilon = 2.3;

int main(int argc, char* argv[])
{
    Lattice L;
    // Correlated site
    L.addSite(new Lattice::Site("C",1,2));
    // Bath sites
    L.addSite(new Lattice::Site("0",1,2));
    L.addSite(new Lattice::Site("1",1,2));
    
    LatticePresets::addCoulombS(&L, "C", U, -mu);
    LatticePresets::addMagnetization(&L, "C", 2*h);
    LatticePresets::addLevel(&L, "0", -epsilon);
    LatticePresets::addLevel(&L, "1", epsilon);
    LatticePresets::addHopping(&L, "C", "0", V);
    LatticePresets::addHopping(&L, "C", "1", V);

    IndexClassification IndexInfo(L.getSiteMap());
    IndexInfo.prepare();
    
    IndexHamiltonian HStorage(&L,IndexInfo);
    HStorage.prepare();

    Symmetrizer Symm(IndexInfo, HStorage);
    Symm.compute();

    StatesClassification S(IndexInfo,Symm);
    S.compute();

    Hamiltonian H(IndexInfo, HStorage, S);
    H.prepare();
    H.diagonalize();

    DensityMatrix rho(S,H,beta);
    rho.prepare();
    rho.compute();
    
    FieldOperatorContainer Operators(IndexInfo, S, H);
    Operators.prepare();

    ParticleIndex down_index = IndexInfo.getIndex("C",0,down);
    ParticleIndex up_index = IndexInfo.getIndex("C",0,up);

    GreensFunction GF_down(S,H,
    	Operators.getAnnihilationOperator(down_index),
    	Operators.getCreationOperator(down_index),
	rho);
    
    GreensFunction GF_up(S,H,
    	Operators.getAnnihilationOperator(up_index),
    	Operators.getCreationOperator(up_index),
	rho);

    GF_down.prepare(); GF_up.prepare();
    GF_down.compute(1000); GF_up.compute(1000);

    // TRIQS part
    std::vector<std::string> block_names;
    block_names.push_back("up");
    block_names.push_back("down");

    auto sha = triqs::arrays::make_shape(1,1);
    auto TRIQS_G_w = make_gf<block_index, gf<imfreq>>(block_names, make_gf<imfreq>(beta, Fermion, sha, 1000));
    auto TRIQS_G_tau = make_gf<block_index, gf<imtime>>(block_names, make_gf<imtime>(beta, Fermion, sha, 1000));

    for(int n=0; n<1000; ++n){
    	TRIQS_G_w[0][n](0,0) = GF_up(n);
	TRIQS_G_w[1][n](0,0) = GF_down(n);
    }

    TRIQS_G_w[0].singularity()(1)(0,0) = 1;
    TRIQS_G_w[1].singularity()(1)(0,0) = 1;

    for(int block=0; block<2; block++) TRIQS_G_tau()[block] = triqs::gf::lazy_inverse_fourier(TRIQS_G_w[block]);

    H5::H5File G_file("anderson.ed.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_up",TRIQS_G_tau[0]);
    h5_write(G_file,"G_down",TRIQS_G_tau[1]);

    return EXIT_SUCCESS;
}
