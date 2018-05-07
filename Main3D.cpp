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

//LZ std::size_t means unsigned integer
//LZ std::shared_ptr a smart pointer
//LZ std::shared_ptr<GenericFunction> A common base class for functions. 
//LZ This functionality is implemented by sub-classes that implement the eval(what value) and restrict(which boundary) functions.
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
        bcs_fix_Domain->push_back(FacetType::GAMMA_F2I);
        bcs_fix_Domain->push_back(FacetType::GAMMA_VI);
        bcs_fix_Domain->push_back(FacetType::GAMMA_VO);
        bcs_fix_Domain->push_back(FacetType::GAMMA_SI);
        //LZ bcs_fix_Domain->push_back(FacetType::GAMMA_SO);
        //LZ bcs_fix_Domain->push_back(FacetType::GAMMA_IFS);
        class Inflow_center_norm : public Expression
        {
            private:
                const Mesh& mesh;

            public:

                Inflow_center_norm(const Mesh& mesh) : Expression(3), mesh(mesh) {}

                void eval(Array<double>& values, const Array<double>& x,
                        const ufc::cell& ufc_cell) const
                {
                    dolfin_assert(ufc_cell.local_facet >= 0);

                    Cell cell(mesh, ufc_cell.index);
                    Point n = cell.normal(ufc_cell.local_facet);
                    std::cout << "Point"<<"   "<<x[0]<<"   "<<x[1]<<"   "<<x[2]<<"   "<<"Norm"<<"   "<<n[0]<<"   "<<n[1]<<"   "<<n[2]<<std::endl;
                }
        };
        class Inflow : public Expression 
        {
            private:
                double* __t;
                const Mesh& __mesh; //LZ
            public:
                Inflow(double* t, const Mesh& mesh) : Expression(3), __t(t), __mesh(mesh){}//LZ
            public:
                void eval(Array<double>& values, const Array<double>& x, const
                        ufc::cell& ufc_cell) const {

                    int n = *__t/0.8;
                    double dt = *__t-n*0.8;
                    double v= 0;
                    // if (dt<0.4)
                    //     v = 0.5 * (sin( 2*DOLFIN_PI/0.4 *dt - DOLFIN_PI/2.0) +1.0 )/2.0 ;
                    // else 
                    //     v = 0.01 * (sin( 2*DOLFIN_PI/0.4 *(dt-0.4) - DOLFIN_PI/2.0) +1.0 )/2.0 ;
                    
                    if (dt<0.2)
                        v = 0.5 * (sin( 2*DOLFIN_PI/0.4 *dt - DOLFIN_PI/2.0) +1.0 )/2.0 ;
                    else 
                        v = 0.5 * (sin( 2*DOLFIN_PI/0.4 *0.2 - DOLFIN_PI/2.0) +1.0 )/2.0 ;
                    
                    dolfin_assert(ufc_cell.local_facet >= 0);
                    Cell cell(__mesh, ufc_cell.index);
                    Point norm = cell.normal(ufc_cell.local_facet);
                    //LZ the inner normal of the inflow
                    // double n0 = -0.000268409;
                    // double n1 = 0.0481119;
                    // double n2 = -0.998842;
                    double n0 = -norm[0];
                    double n1 = -norm[1];
                    double n2 = -norm[2];

                    //LZ the Center point of inflow 
                    double x0 = 0.0189097394;//LZ
                    double x1 = 0.1421571318;//LZ
                    double x2 = -0.2198169226;//LZ

                    double d = 0.00004470652;//LZ radius^2                 
                    double r = (x[0]-x0)*(x[0]-x0)+(x[1]-x1)*(x[1]-x1)+(x[2]-x2)*(x[2]-x2);//LZ r^2
                    double dis = std::abs(d-r)/d;
                    //double dis = 1.0;

                    values[0] = n0*v*2.0*dis; //LZ
                    values[1] = n1*v*2.0*dis; //LZ
                    values[2] = n2*v*2.0*dis; //LZ

                }

        };

        class beta1Function : public Expression 
        {
            private:
                double* __t;
                const Mesh& __mesh; //LZ
            public:
                beta1Function(double* t, const Mesh& mesh) : Expression(), __t(t), __mesh(mesh){}//LZ
            public:
                void eval(Array<double>& values, const Array<double>& x, const
                        ufc::cell& ufc_cell) const {
                    values[0] = 17.4e+4;
                    //values[0] = 1.0e+07 * (6.0183 * x[2] * x[2] + 3.1837  * x[2] + 0.4212); //LZ beta1*0.01
                    //LZ values[0] = 1.0e+07 * (6.0731 * x[2] * x[2] + 3.2126 * x[2] + 0.4249); //LZ beta1*0.001

                }

        };
        class beta2Function : public Expression 
        {
            private:
                double* __t;
                const Mesh& __mesh; //LZ
            public:
                beta2Function(double* t, const Mesh& mesh) : Expression(), __t(t), __mesh(mesh){}//LZ
            public:
                void eval(Array<double>& values, const Array<double>& x, const
                        ufc::cell& ufc_cell) const {
                    values[0] = 188.1e+4; //LZ Const
                    //values[0] = 1.0e+08 * (6.5060 * x[2] * x[2] +3.4417 * x[2] + 0.4554); //LZ beta2*0.01
                    //LZ values[0] = 1.0e+08 * (6.5652 * x[2] * x[2] +3.4730 * x[2] + 0.4593); //LZ beta2*0.001

                }

        };
        class OutflowStress : public Expression
        {
            public:

                OutflowStress(const Mesh& mesh) : Expression(3), mesh(mesh) {}

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
                    values[2] = g*n[2];

                }

            private:

                const Mesh& mesh;

        };
        beta1Func.reset(new beta1Function(&t,*mesh));
        beta2Func.reset(new beta2Function(&t,*mesh));
        sigmaN.reset(new OutflowStress(*mesh));
        std::shared_ptr<Inflow> inflow = std::make_shared<Inflow>(&t,*mesh);
        std::shared_ptr<Constant> zero = std::make_shared<Constant>(0.0,0.0,0.0);


        bdfuncs_FSI = new std::vector<std::shared_ptr<GenericFunction>>;
        bcsubDomain_FSI = new std::vector<std::size_t>;

        //LZ Find the norm and center point of the inflow
        /*{
            std::shared_ptr<Inflow_center_norm> inflow_center_norm = std::make_shared<Inflow_center_norm>(*mesh);
            bdfuncs_FSI->push_back(inflow_center_norm);
            bcsubDomain_FSI->push_back(FacetType::GAMMA_FI);
            bdfuncs_FSI->push_back(inflow_center_norm);
            bcsubDomain_FSI->push_back(FacetType::GAMMA_F2I);
            bdfuncs_FSI->push_back(inflow_center_norm);
            bcsubDomain_FSI->push_back(FacetType::GAMMA_SI);
        }*/

        //LZ the inflow boundary condition must be imposed first, then the other conditions.
        bdfuncs_FSI->push_back(inflow);
        bcsubDomain_FSI->push_back(FacetType::GAMMA_FI);
        bdfuncs_FSI->push_back(inflow);
        bcsubDomain_FSI->push_back(FacetType::GAMMA_F2I);
        bdfuncs_FSI->push_back(inflow);
        bcsubDomain_FSI->push_back(FacetType::GAMMA_SI);

        bdfuncs_FSI->push_back(zero);
        bcsubDomain_FSI->push_back(FacetType::GAMMA_VI);
        bdfuncs_FSI->push_back(zero);
        bcsubDomain_FSI->push_back(FacetType::GAMMA_VO);
        //LZ bdfuncs_FSI->push_back(zero);
        //LZ bcsubDomain_FSI->push_back(FacetType::GAMMA_SI);

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

    // #ifdef HAS_PETSC
    // parameters["linear_algebra_backend"] = "PETSc";
    // PetscInitialize(&argc,&argv,NULL,NULL);
    // PetscOptionsInsertFile(PETSC_COMM_WORLD,NULL,"../data/petsc_options",PETSC_TRUE);
    // if(dolfin::MPI::rank(MPI_COMM_WORLD) == 0){
    //     info("has PETSc");
    //     //list_petsc_snes_methods();
    //     //list_petsc_ksp_methods();
    //     //list_petsc_pre_methods();
    // }
    // #endif
    
    double t, end_T, dt;

    // HDF5File filer(MPI_COMM_WORLD,"data/mesh/pipe.h5","r");
    HDF5File filer(MPI_COMM_WORLD,"../data/mesh/merge.h5","r");
    /// Load test mesh 
    std::shared_ptr<Mesh> mesh=std::make_shared<Mesh>();
    filer.read(*mesh,std::string("mesh"),false);

    /// Load mesh-function to specify the subdomains
    std::shared_ptr<MeshFunction<std::size_t> > cell_marker= std::make_shared<MeshFunction<std::size_t>>
        (mesh, mesh->topology().dim());
    std::shared_ptr<MeshFunction<std::size_t> > bd_marker = std::make_shared<MeshFunction<std::size_t>>
        (mesh, mesh->topology().dim()-1);

    filer.read(*cell_marker,std::string("subdomains_mark"));//LZ Read all the marks in domain
    filer.read(*bd_marker,std::string("facet_mark"));//LZ Read all the marks on facet: Boundary, interface?
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
    std::string para_file_fsi("../data/parameters_fsi.xml");

    ///  Set up Boundary-Initial condition
    TestFSI fsisolver(mesh, cell_marker, bd_marker, para_file_fsi);

    /// Solving procedure
    fsisolver.solve();




    return 0;
}

/**
 * end of file 
 *
 */
