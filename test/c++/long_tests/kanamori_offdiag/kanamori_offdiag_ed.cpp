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

RealType beta = 20;
int num_o = 2;
RealType U = 2.0;
RealType J = 0.2;
RealType mu = 1.0;
RealType h = 0.0;
RealType V = 1.0;
RealType Vp = 0.1;
RealType epsilon = 2.3;

int main(int argc, char* argv[])
{
    Lattice L;
    // Correlated site
    L.addSite(new Lattice::Site("C",2,2));
    
    // Bath sites
    L.addSite(new Lattice::Site("B0p", 1, 2));
    L.addSite(new Lattice::Site("B1p", 1, 2));
    L.addSite(new Lattice::Site("B0m", 1, 2));
    L.addSite(new Lattice::Site("B1m", 1, 2));
    
    LatticePresets::addCoulombP(&L, "C", U, J, -mu);
    LatticePresets::addMagnetization(&L, "C", 2*h);

    LatticePresets::addLevel(&L, "B0p", -epsilon);
    LatticePresets::addLevel(&L, "B1p", -epsilon);
    LatticePresets::addLevel(&L, "B0m", epsilon);
    LatticePresets::addLevel(&L, "B1m", epsilon);
    
    LatticePresets::addHopping(&L, "C", "B0p", V, 0, 0);
    LatticePresets::addHopping(&L, "C", "B0p", Vp, 1, 0);
    LatticePresets::addHopping(&L, "C", "B1p", Vp, 0, 0);
    LatticePresets::addHopping(&L, "C", "B1p", V, 1, 0);
    LatticePresets::addHopping(&L, "C", "B0m", V, 0, 0);
    LatticePresets::addHopping(&L, "C", "B0m", Vp, 1, 0);
    LatticePresets::addHopping(&L, "C", "B1m", Vp, 0, 0);
    LatticePresets::addHopping(&L, "C", "B1m", V, 1, 0);
    
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

    // TRIQS part
    std::vector<std::string> block_names;
    block_names.push_back("up");
    block_names.push_back("down");

    auto sha = triqs::arrays::make_shape(num_o,num_o);
    auto TRIQS_G_w = make_gf<block_index, gf<imfreq>>(block_names, make_gf<imfreq>(beta, Fermion, sha, 1000));
    auto TRIQS_G_tau = make_gf<block_index, gf<imtime>>(block_names, make_gf<imtime>(beta, Fermion, sha, 1000));

    for(int o1=0; o1<num_o; ++o1)
    for(int o2=0; o2<num_o; ++o2){
        
        GreensFunction GF_up(S,H,
	    Operators.getAnnihilationOperator(IndexInfo.getIndex("C",o1,up)),
	    Operators.getCreationOperator(IndexInfo.getIndex("C",o2,up)),
	    rho);
	
        GreensFunction GF_down(S,H,
	    Operators.getAnnihilationOperator(IndexInfo.getIndex("C",o1,down)),
	    Operators.getCreationOperator(IndexInfo.getIndex("C",o2,down)),
	    rho);

        GF_down.prepare(); GF_up.prepare();
        GF_down.compute(1000); GF_up.compute(1000);

        for(int n=0; n<1000; ++n){
	    TRIQS_G_w[0][n](o1,o2) = GF_up(n);
	    TRIQS_G_w[1][n](o1,o2) = GF_down(n);
        }
	TRIQS_G_w[0].singularity()(1)(o1,o2) = (o1==o2) ? 1 : 0;
	TRIQS_G_w[1].singularity()(1)(o1,o2) = (o1==o2) ? 1 : 0;
    }
    
    for(int block=0; block<2; block++) TRIQS_G_tau()[block] = triqs::gf::lazy_inverse_fourier(TRIQS_G_w[block]);

    H5::H5File G_file("kanamori_offdiag.ed.h5",H5F_ACC_TRUNC);
    for(int o=0; o<num_o; ++o){
        h5_write(G_file,"G_up",TRIQS_G_tau[0]);
        h5_write(G_file,"G_down",TRIQS_G_tau[1]);
    }

    return EXIT_SUCCESS;
}
