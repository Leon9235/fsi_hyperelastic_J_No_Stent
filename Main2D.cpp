/**
 * @file main.cpp
 * @brief 
 * @author AAA group
 * @version 
 * @date 2015-01-16
 *
 * Modify by Shihua Gong 30-04-2016
 * 
 */

#include <dolfin.h>
#include "solverbase.h"

using namespace dolfin;


class TestFSI: public SolverBase
{
    public:

        TestFSI(std::shared_ptr<Mesh> mesh, std::shared_ptr<MeshFunction<std::size_t>> cell_mark,
                std::shared_ptr<MeshFunction<std::size_t>> bd_mark, std::string para_file_fsi) 
            : SolverBase(mesh, cell_mark, bd_mark, para_file_fsi)
    {


        bcs_fix_Domain = new std::vector<std::size_t>;
        bcs_fix_Domain->push_back(FacetType::GAMMA_FI); 
        bcs_fix_Domain->push_back(FacetType::GAMMA_FO); 

        class Inflow : public Expression 
        {
            private:
                double* __t;
            public:
                Inflow() : Expression(2) {}
                Inflow(double* t) : Expression(2), __t(t) {}
            public:
                void eval(Array<double>& values, const Array<double>& x, const
                        ufc::cell& ufc_cell) const {
                    // aaa
                    //values[0] = 0.0;
                    //values[1] = -1000.0 *(-x[0]*x[0]+ 0.005*0.005)*std::abs(sin(*__t));

                    // pi
                    //values[1] = 0.0;
                    //values[0] =  (-x[1]*x[1]+ 0.5*0.5)*std::abs(sin(*__t));

                    // aaalong
                    values[0] = 0.0;
                    //double t = ( *__t)*10e5;
                    //long n = t/0.8;
                    //double dt = t-n*0.8;
                    int n = *__t/0.8;
                    double dt = *__t-n*0.8;
                    double v= 0;

                    if (dt<0.4)
                        v = 0.5 * (sin( 2*DOLFIN_PI/0.4 *dt - DOLFIN_PI/2) +1 )/2 ;
                    else 
                        v = 0.01 * (sin( 2*DOLFIN_PI/0.4 *(dt-0.4) - DOLFIN_PI/2) +1 )/2 ;
                    values[0] = 0.0;
                    values[1] = v* (x[0] - 0) *(x[0] - 0.02)* 10000 ;


                }

        };
        class OutflowStress : public Expression
        {
            public:

                OutflowStress(const Mesh& mesh) : Expression(2), mesh(mesh) {}

                void eval(Array<double>& values, const Array<double>& x,
                        const ufc::cell& ufc_cell) const
                {
                    dolfin_assert(ufc_cell.local_facet >= 0);

                    Cell cell(mesh, ufc_cell.index);
                    Point n = cell.normal(ufc_cell.local_facet);

                    //const double g = -13330.0;
                    const double g = 0.0;
                    values[0] = g*n[0];
                    values[1] = g*n[1];
                    //values[2] = g*n[1];
                }

            private:

                const Mesh& mesh;

        };

        sigmaN.reset(new OutflowStress(*mesh));
        std::shared_ptr<Inflow> inflow = std::make_shared<Inflow>(&t);
        std::shared_ptr<Constant> zero = std::make_shared<Constant>(0.0,0.0);

        bdfuncs_FSI = new std::vector<std::shared_ptr<GenericFunction>>;
        bcsubDomain_FSI = new std::vector<std::size_t>;

        bdfuncs_FSI->push_back(zero);
        bcsubDomain_FSI->push_back(FacetType::GAMMA_0);
        //bdfuncs_FSI->push_back(zero);
        //bcsubDomain_FSI->push_back(FacetType::GAMMA_S);
        //bdfuncs_FSI->push_back(zero);
        //bcsubDomain_FSI->push_back(FacetType::GAMMA_I);

        bdfuncs_FSI->push_back(inflow);
        bcsubDomain_FSI->push_back(FacetType::GAMMA_FI);

    }

        ~TestFSI()
        {
        };

};

/**
 * @brief 
 *
 * @param argc
 * @param argv[]
 *
 * @return 
 */
int main(int argc, char * argv[])
{
    double t, end_T, dt;

    HDF5File filer(MPI_COMM_WORLD,"data/mesh/aaalong.h5","r");
    /// Load test mesh 
    std::shared_ptr<Mesh> mesh=std::make_shared<Mesh>();
    filer.read(*mesh,std::string("mesh"),false);

    /// Load mesh-function to specify the subdomains
    std::shared_ptr<MeshFunction<std::size_t> > cell_marker= std::make_shared<MeshFunction<std::size_t>>
        (mesh, mesh->topology().dim());
    std::shared_ptr<MeshFunction<std::size_t> > bd_marker = std::make_shared<
            MeshFunction<std::size_t>>(mesh, mesh->topology().dim()-1);

    filer.read(*cell_marker,std::string("subdomains_mark"));
    filer.read(*bd_marker,std::string("facet_mark"));
    filer.close();


    std::vector<double>& coord = mesh->coordinates();
    for(int i = 0; i < coord.size(); i++)
    {
        coord[i]*= 0.001;

    }

    dolfin::info("Number of Global Vertices: %d", mesh->size_global(0));
    dolfin::info("Number of Mesh Vertices: %d", mesh->num_vertices());
    dolfin::info("Number of Mesh Elements: %d", mesh->num_cells());

    /// Load parameter file 
    std::string para_file_fsi("./data/parameters_fsi.xml");


    TestFSI fsisolver(mesh, cell_marker, bd_marker, para_file_fsi);

    fsisolver.solve();




    return 0;
}

/**
 * end of file 
 *
 */
