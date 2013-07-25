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
RealType t = 1.0;
RealType epsilon = 2.3;

int main(int argc, char* argv[])
{
    Lattice L;
    // Correlated site
    L.addSite(new Lattice::Site("1",1,1));
    L.addSite(new Lattice::Site("2",1,1));
    
    // Bath sites
    L.addSite(new Lattice::Site("1+",1,1));
    L.addSite(new Lattice::Site("2+",1,1));
    L.addSite(new Lattice::Site("1-",1,1));
    L.addSite(new Lattice::Site("2-",1,1));

    LatticePresets::addLevel(&L,"1",-mu);
    LatticePresets::addLevel(&L,"2",-mu);
    LatticePresets::addHopping(&L,"1","2",-t,0,0,0,0);
    L.addTerm(Lattice::Term::Presets::NupNdown("1","2",U,0,0,0,0));

    LatticePresets::addHopping(&L,"1","1+",1.0,0,0,0,0);
    LatticePresets::addHopping(&L,"1","1-",1.0,0,0,0,0);
    LatticePresets::addHopping(&L,"2","2+",1.0,0,0,0,0);
    LatticePresets::addHopping(&L,"2","2-",1.0,0,0,0,0);

    LatticePresets::addLevel(&L,"1+",-epsilon);
    LatticePresets::addLevel(&L,"2+",-epsilon);
    LatticePresets::addHopping(&L,"1+","2+",-t,0,0,0,0);
    LatticePresets::addLevel(&L,"1-",epsilon);
    LatticePresets::addLevel(&L,"2-",epsilon);
    LatticePresets::addHopping(&L,"1-","2-",-t,0,0,0,0);

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

    ParticleIndex indices[2];
    indices[0] = IndexInfo.getIndex("1",0,0);
    indices[1] = IndexInfo.getIndex("2",0,0);

    // TRIQS part
    std::vector<std::string> block_names;
    block_names.push_back("tot");

    auto sha = triqs::arrays::make_shape(2,2);
    auto TRIQS_G_w = make_gf<block_index, gf<imfreq>>(block_names, make_gf<imfreq>(beta, Fermion, sha, 1000));
    auto TRIQS_G_tau = make_gf<block_index, gf<imtime>>(block_names, make_gf<imtime>(beta, Fermion, sha, 1000));

    for(int i1=0; i1<2; ++i1)
    for(int i2=0; i2<2; ++i2){
    	GreensFunction GF(S,H,
		Operators.getAnnihilationOperator(indices[i1]),
		Operators.getCreationOperator(indices[i2]),
		rho);
	GF.prepare();
	GF.compute(1000);

	for(int n=0; n<1000; ++n) TRIQS_G_w[0][n](i1,i2) = GF(n);
    }

    TRIQS_G_w[0].singularity()(1)(0,0) = 1;
    TRIQS_G_w[0].singularity()(1)(1,0) = 0;
    TRIQS_G_w[0].singularity()(1)(0,1) = 0;
    TRIQS_G_w[0].singularity()(1)(1,1) = 1;

    TRIQS_G_tau()[0] = triqs::gf::lazy_inverse_fourier(TRIQS_G_w[0]);

    H5::H5File G_file("spinless.ed.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_tau",TRIQS_G_tau[0]);
    
    return EXIT_SUCCESS;
}
