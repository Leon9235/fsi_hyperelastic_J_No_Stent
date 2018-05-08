/** \file    solverbase.cpp
 *  \brief   solver base for FSI problems
 *  \author Lu WANG
 *  \date   03/12/2014
 *  
 * Modify by Shihua 05/13/2016
 *
 */


#include "solverbase.h"
#include "FSI-weakform/Fluid2DP1P1.h"
#include "FSI-weakform/Elastic2DP1P1.h"
#include "FSI-weakform/LinearElastic2DP1P1.h"

#include "FSI-weakform/Fluid3DP1P1.h"
#include "FSI-weakform/Elastic3DP1P1.h"
#include "FSI-weakform/LinearElastic3DP1P1.h"

using namespace std;
using namespace dolfin;

#define _TIMER
#define _DEBUG

/**
 *  \fn SolverBase(ProblemBase& prob)
 *
 *	\brief	Class Constructor
 *
 *  \param  prob:  problem setup info for fsi system
 *  \author Lu WANG
 *  \date   03/12/2014
 */
SolverBase::SolverBase(std::shared_ptr<Mesh> _mesh0, 
        std::shared_ptr<MeshFunction<std::size_t>> _cell_mark0,
        std::shared_ptr<MeshFunction<std::size_t>> _bd_mark0, 
        std::string pfile_fsi) 
: mesh0(_mesh0), cell_mark0(_cell_mark0), bd_mark0(_bd_mark0)
{

    std::cout << "begin Constructor" << std::endl;
    mesh.reset(new Mesh(*mesh0));
    cell_mark.reset( new MeshFunction<std::size_t>(mesh, mesh->topology().dim()));
    bd_mark.reset( new MeshFunction<std::size_t>(mesh, mesh->topology().dim()-1));
    std::copy(cell_mark0->values(), cell_mark0->values() + cell_mark0->size(), cell_mark->values());
    std::copy(bd_mark0->values(), bd_mark0->values() + bd_mark0->size(), bd_mark->values());
    std::cout<< "mesh load"<<std::endl;



    read_parameters(pfile_fsi);

    t = 0;
    dt = (double)(_parameters["timestep"]);
    T  = _parameters["endtime"];
    int N = (_parameters["sample_number"]);
    sample_period  = T/N;
    last_sample = 0;

    tol   = _parameters["iteration_tolerance"];
    iter_max = _parameters["maximum_iterations"];
    iter_FS_max = _parameters["FSnewtonitrmax"];

    out_velocity = new File("./result/velocity.pvd");
    out_pressure = new File("./result/pressure.pvd");
    out_S = new File("./result/S.pvd");
    out_displacement = new File("./result/displacement.pvd");


    // LU Solver
    //lusolver = new LUSolver;
    //PETSc solver
    //FASPsolver = new PETScFASP("gmres","amg");
    /*
       FASPsolver->parameters["relative_tolerance"] = 1.0e-6;
       FASPsolver->parameters["absolute_tolerance"] = 1.0e-10;
       FASPsolver->parameters["divergence_limit"] = 1.0e1;
       FASPsolver->parameters["maximum_iterations"] = 10;
       FASPsolver->parameters["report"] = true;
       FASPsolver->parameters["error_on_nonconvergence"] = true;
       FASPsolver->parameters["nonzero_initial_guess"] = false;
       FASPsolver->parameters["monitor_convergence"] = true;

    //Krylov Solver
    const std::string prec(has_krylov_solver_preconditioner("amg")
    ? "amg" : "default");
    ksolver  = new KrylovSolver("gmres", prec);
    ksolver->parameters["relative_tolerance"] = 1.0e-12;
    ksolver->parameters["absolute_tolerance"] = 1.0e-15;
    ksolver->parameters["divergence_limit"] = 1.0e1;
    ksolver->parameters["maximum_iterations"] = 500;
    ksolver->parameters["report"] = true;
    ksolver->parameters["error_on_nonconvergence"] = true;
    ksolver->parameters["nonzero_initial_guess"] = true;
    ksolver->parameters["monitor_convergence"] = true;

    //FASP solver
    fsolver  = new FaspSolver();

    M.reset(new Matrix);
    preM.reset(new Matrix);
    b.reset(new Vector);
    */


    //assembler = new FSIAssembler;
    //assembler->init(_mesh0, elemark);
}

/**
 * @brief Read parameters 
 *
 * @param para_file
 */
void SolverBase::read_parameters(std::string para_file) 
{
    std::cout << "Reading parameters ... ";
    File file(para_file);
    if (!file.exists(para_file)) {
        dolfin_error("fsisolver.cpp",
                "read parameters",
                "Can not find the xml file");
    }
    //Parameters new_para("problem_parameters");
    file >> _parameters;
    //_parameters.add(new_para);

    std::cout << "Done!" << std::endl;
    dolfin::info(_parameters, true);
}

/**
 *  \fn init()
 *
 *	\brief	Initialization of the spaces and coefficients
 *  \author Lu WANG
 *  \date   03/12/2014
 *
 */
void SolverBase::init()
{

    // FunctionSpace
    //LZ In ufl, we only define a function space without mesh. 
    //LZ Here, we take the space out and then discrete it to form a discrete function space.
    //LZ Then let it be the space of v , p...
    //LZ And the ufl even do not know the discrete.
    //LZ We then tell the ufl what is dx and ds which transfer mesh information to ufl
    size_t d = mesh->geometry().dim();
    if (_parameters["FEspace"].value_str() == "P1P1")
    {
        // if (d == 2) {
        //     VQQ.reset(new Elastic2DP1P1::FunctionSpace(mesh));
        //     DG1.reset(new Elastic2DP1P1::CoefficientSpace_beta1(mesh));
        //     DG.reset(new Fluid2DP1P1::CoefficientSpace_delta_SUPG(mesh));
        // }
        //else if (d == 3) {
            VQQ.reset(new Elastic3DP1P1::FunctionSpace(mesh));
            DG1.reset(new Elastic3DP1P1::CoefficientSpace_rho_s(mesh));
            DG.reset(new Fluid3DP1P1::CoefficientSpace_delta_SUPG(mesh));
        //}
    }
    else
    {
        dolfin_error("fsisolver.cpp",
                "identify FEM",
                "support P1 at present");
        //break;
    }

    _sub_V  = VQQ->sub(0);
    _sub_Q  = VQQ->sub(1);
    _sub_Qs = VQQ->sub(2);
    V  = _sub_V->collapse(); //LZ Why we need to collapse?
    Q  = _sub_Q->collapse(); //LZ after collapse, which should we use? V or _sub_V
    Qs = _sub_Qs->collapse(); //LZ after collapse, which should we use? V or _sub_V



    //Generate bilinear form and linear form
    if (_parameters["FEspace"].value_str() == ("P1P1"))
    {
        // if (d == 2) {
        //     J_f.reset(new Fluid2DP1P1::JacobianForm(VQQ,VQQ));

        //     F_f.reset(new Fluid2DP1P1::ResidualForm(VQQ));

        //     if(_parameters["ElasticModel"].value_str() =="Linear")
        //     {
        //             info("Use linear elastic model\n");
        //             J_s.reset(new LinearElastic2DP1P1::JacobianForm(VQ,VQ));
        //             F_s.reset(new LinearElastic2DP1P1::ResidualForm(VQ));
        //     }
        //     else{
        //             info("Use hyper elastic model\n");
        //             J_s.reset(new Elastic2DP1P1::JacobianForm(VQQ,VQQ));
        //             F_s.reset(new Elastic2DP1P1::ResidualForm(VQQ));
                
        //     }
        // }
        // else if (d == 3) {

            J_f.reset(new Fluid3DP1P1::JacobianForm(VQQ,VQQ));

            F_f.reset(new Fluid3DP1P1::ResidualForm(VQQ));

            // if(_parameters["ElasticModel"].value_str() =="Linear")
            // {
            //         info("Use linear elastic model\n");
            //         J_s.reset(new LinearElastic3DP1P1::JacobianForm(VQ,VQ));
            //         F_s.reset(new LinearElastic3DP1P1::ResidualForm(VQ));
            // }
            //else{
                    info("Use hyper elastic model\n");
                    J_s.reset(new Elastic3DP1P1::JacobianForm(VQQ,VQQ));
                    F_s.reset(new Elastic3DP1P1::ResidualForm(VQQ));
                
            //}
        //}
    }
    else
    {
        dolfin_error("fsisolver.cpp",
                "identify FEM",
                "support P1 at present");
        //break;
    }

    // variable
    v_p_p.reset(new Function(VQQ));  //LZ velocity and pressure in mixed form
    //Function ini_p(Q);
    //ini_p = Constant(13330.0);
    //assign(reference_to_no_delete_pointer( (*v_p)[1]), reference_to_no_delete_pointer(ini_p));

    u.reset(new Function(V));  //LZ displacement
    u0.reset(new Function(V)); //LZ displacement at last time-step 

    uA.reset(new Function(V));  //LZ displacement for smoothing at current mesh iteration
    uA0.reset(new Function(V)); //LZ displacement for smoothing at last time-step???
    uA1.reset(new Function(V)); //LZ displacement for smoothing at last mesh iteration
    w.reset(new Function(V));   //LZ (uA-uA0)/dt velocity for smoothing
    dv.reset(new Function(V));  //LZ dv =uA-uA1 The movement of mesh from last mesh iteration to current mesh iteration
    a0.reset(new Function(V));  //LZ For Newmark scheme
                        

    v.reset(new Function(V));  //LZ velocity
    v0.reset(new Function(V)); //LZ velocity at last time-step
    v1.reset(new Function(V)); //LZ ???
 
    p.reset(new Function(Q));  //LZ pressure
    ps.reset(new Function(Qs));  //LZ pressure

    // Common Coefficient
    _dt.reset(new Constant( double(_parameters["timestep"]) ));

    // Fluid Coefficient
    _rho_f.reset(new Constant(double(_parameters["rho_f"])));
    _mu_f.reset( new Constant(double(_parameters["mu_f"])));
    _alpha.reset(new Constant(double(_parameters["alpha"])));
    _alpha2.reset(new Constant(double(_parameters["alpha2"]))); //LZ add parameter for structure stable term
    delta_SUPG.reset(new Function(DG));


    // Elastic Coeffcient
    _rho_s.reset(new Constant(double(_parameters["rho_s"])));
    double nu_s = double(_parameters["nu_s"]);
    double mu_s = double(_parameters["mu_s"]);
    //double lamb_s = 2*mu_s*nu_s/(1-2*nu_s);
    _mu_s.reset(new Function(DG1));
    *_mu_s = Constant(mu_s);
    _nu_s.reset(new Constant(nu_s));
    _p_s.reset(new Constant(double(_parameters["p_s"])));
    _beta1.reset(new Constant(double(_parameters["beta1"])));
    _beta2.reset(new Constant(double(_parameters["beta2"])));


    std::map<std::string, std::shared_ptr<const GenericFunction>> coef_list
        = { {"v_p_p",         v_p_p},
            {"v0",          v0},
            {"u0",          u0},
            {"v1",          v1},
            {"a0",          a0},
            {"rho_s",       _rho_s},
            {"p_s",         _p_s},
            {"beta1",       _beta1},
            {"beta2",       _beta2},
            {"beta1Func",   beta1Func},
            {"beta2Func",   beta2Func},
            {"dt",          _dt},
            {"sigmaN",      sigmaN},
            {"mu_s",        _mu_s},
            {"nu_s",        _nu_s},
            {"vm",          w},
            {"delta_SUPG",  delta_SUPG},
            {"rho_f",       _rho_f},
            {"mu_f",        _mu_f},
            {"alpha1",      _alpha},
            {"alpha2",      _alpha2}};
    J_f->set_some_coefficients(coef_list);
    F_f->set_some_coefficients(coef_list);
    J_s->set_some_coefficients(coef_list);
    F_s->set_some_coefficients(coef_list);

    J_s->dx = cell_mark; J_s->ds = bd_mark; J_s->dS = bd_mark;
    F_s->dx = cell_mark; F_s->ds = bd_mark; F_s->dS = bd_mark;
    J_f->dx = cell_mark; J_f->ds = bd_mark; J_f->dS = bd_mark;
    F_f->dx = cell_mark; F_f->ds = bd_mark; F_f->dS = bd_mark;

    //std::string smooth_method("harmonicsmoothing"); //elasticsmoothing
    std::string smooth_method("elasticsmoothing"); //elasticsmoothing
    meshsmoother.reset(new MeshSmoothing(mesh0, cell_mark0, bd_mark0, *bcs_fix_Domain, smooth_method));

    //LZ Impose all Dirichlet Boundary Condition
    for(int i = 0; i < bcsubDomain_FSI->size(); ++i)
    {   
        //LZ Why use _sub_V? not V? cllopse?
        _bcsFSI.push_back( std::shared_ptr<const DirichletBC> (
                    new DirichletBC(_sub_V, ((*bdfuncs_FSI)[i]), bd_mark,(*bcsubDomain_FSI)[i] )));
    }
    //Constant * pressure_outflow =new Constant(133322);
    //_bcsFSI.push_back( std::shared_ptr<const DirichletBC> (
    //            new DirichletBC(*_sub_Q, *pressure_outflow, *bd_mark, FacetType::GAMMA_FO)));


    pressure_dirichlet();
    
    nonlinear_problem = std::shared_ptr<NonlinearDiscreteProblem>(new NonlinearDiscreteProblem(
                reference_to_no_delete_pointer(*this)));

    // The Newton solver
    newton_solver = std::shared_ptr<NewtonSolver>(new NewtonSolver()); 
    newton_solver->parameters["maximum_iterations"] = int(_parameters["FSnewtonitrmax"]);
    newton_solver->parameters["relative_tolerance"] =1e-6* double(_parameters["FSnewtonsoltol"]);
    newton_solver->parameters["absolute_tolerance"] = double(_parameters["FSnewtonsoltol"]);
    newton_solver->parameters["relaxation_parameter"] = double(_parameters["FSrelaxation_parameter"]);//1.0;
    newton_solver->parameters["convergence_criterion"] = "residual"; //"incremental"
    newton_solver->parameters["report"] = true;
    newton_solver->parameters["error_on_nonconvergence"] = false;


}

/**
 *  \fn ~SolverBase()
 *
 *	\brief	Class destructor
 *
 *  \author Lu WANG
 *  \date   03/12/2014
 */
SolverBase::~SolverBase()
{
    /*    delete lusolver;
          delete ksolver;

    //delete M;
    //delete b;
    delete smooth;

    delete out_velocity;
    delete out_pressure;
    delete out_theta;
    delete out_displacement;

    delete dx;
    delete dA;
    delete assembler;
    */
}



/**
 *  \fn step()
 *
 *	\brief	Step from t^n to t^{n+1}
 *
 *  \author Lu WANG
 *  \date   03/12/2014
 *
 */
void SolverBase::step()
{
    double iter_res = 0.;
    std::size_t itercounter = 0;

    std::cout << "begin stepping" <<std::endl;
    /// loop for fixed-point iteration
    for(int i=0; i<iter_max; i++)
    {
        // prepareOuterIteration();
        *(uA1->vector()) = *(uA->vector());


        // Newton solver
        update_delta_SUPG();
        std::cout << "begin Newtonsolver" <<std::endl;
        newton_solver->solve(*nonlinear_problem, *(v_p_p->vector()));
        std::cout << "Newtonsolver done" <<std::endl;


        assign(v, reference_to_no_delete_pointer((*v_p_p)[0]));
        assign(p, reference_to_no_delete_pointer((*v_p_p)[1]));
        assign(ps, reference_to_no_delete_pointer((*v_p_p)[2]));//LZ added



        //LZ calculate displacement u 
        *(u->vector()) = *(v0->vector());
        *(u->vector()) += *(v->vector());
        *(u->vector()) *= dt/2;
        *(u->vector()) += *(u0->vector());
        *mesh = *mesh0;
        
        // update mesh
        std::cout << "begin smoothing" <<std::endl;
        meshsmoother->smoothing(mesh, *u, *uA); //LZ !!!!!!!!!!! u can not change!!!!!!!!!!
        *(w->vector()) = *(uA->vector());
        *(w->vector()) -= *(uA0->vector());
        *(w->vector()) /= dt;

        #ifdef _DEBUG
        std::cout << "mesh smoothing done" <<std::endl;
        std::cout << "norm 4 u  is " << u->vector()->norm("l2") << "\n";
        std::cout << "norm 4 v  is " << v->vector()->norm("l2") << "\n";
        std::cout << "norm 4 p  is " << p->vector()->norm("l2") << "\n";
        std::cout << "norm 4 uA  is " << uA->vector()->norm("l2") << "\n";
        std::cout << "norm 4 uA0  is " << uA0->vector()->norm("l2") << "\n";
        std::cout << "norm 4 u0  is " << u0->vector()->norm("l2") << "\n";
        std::cout << "norm 4 v0  is " << v0->vector()->norm("l2") << "\n";
        #endif
        

        *(dv->vector()) = *(uA->vector());
        *(dv->vector()) -= *(uA1->vector());
        iter_res = dv->vector()->norm("l2");
        #ifdef _DEBUG
        std::cout << "iter_res  is " << iter_res << "\n";
        #endif

        itercounter++;
        //if (dolfin::MPI::process_number() == 0) {
        std::cout << "\x1b[35;1m error of outer iteration "<< i << ": " << iter_res << "\x1b[0m\n";
        //}
        if(iter_res < tol) break;
    }

    std::cout << "\x1b[32;1m=======================================================\x1b[0m\n";
    std::cout << "\x1b[32;1m # of iteration : " << itercounter << "\x1b[0m\n";
    //}
}


/**
 *  \fn update_delta_SUPG()
 *
 */
void SolverBase::update_delta_SUPG()
{
    //*delta_SUPG = Constant(0.0);

    int D = mesh->topology().dim();
    int d = mesh->geometry().dim();

    const std::size_t num_vertices = mesh->num_vertices();

    std::size_t num_cell_dofs = 0;
    std::vector<dolfin::la_index> cell_dofs;
    std::vector<double> cell_values;

    ufc::cell ufc_cell;

    std::vector<double> v_values, w_values;
    v0->compute_vertex_values(v_values, *mesh);
    w->compute_vertex_values(w_values, *mesh);

    for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
        cell->get_cell_data(ufc_cell);
        double h = cell->circumradius();

        // Tabulate dofs on cell
        std::vector <double > v_mean(d);
        for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
        {
            for (std::size_t i = 0; i < d; ++i)
            {
                const std::size_t local_index  = vertex.pos()*d + i;
                const std::size_t global_index = i*num_vertices+vertex->index();

                v_mean[i] +=  (v_values[global_index] - w_values[global_index])/(d+1);
            }
        }

        double vbar;
        if(d==3)
            vbar = sqrt(v_mean[0]*v_mean[0] + v_mean[1]*v_mean[1] + v_mean[2]*v_mean[2]);
        else
            vbar = sqrt(v_mean[0]*v_mean[0] + v_mean[1]*v_mean[1] );
        double reynold = (double) (*_rho_f) * vbar * h / (2 * (double) (*_mu_f));
        double delta;

        double delta_SD = _parameters["delta_SD"];
        if (reynold > 1)
            //delta = h / (2 * vbar)*(1-1./reynold) ;
            delta = reynold*h*vbar* delta_SD;
        else
            delta = 0. ;

        /*if (reynold > 3)
          delta = h / (2 * vbar) ;
          else
          delta = h*h/12*(double)(*_rho_f) /(double)(*_mu_f) ;
          */
        //if(delta>0)
        //std::cout << "h :" <<h<<"  vbar:  "<<vbar<<"   delta: "<<delta<<std::endl;
        //LZ std::cout << "reynold :" << reynold<<std::endl;

        ArrayView<const dolfin::la_index>  delta_dof = delta_SUPG->function_space()->dofmap()->cell_dofs(cell->index());
        cell_dofs.push_back(delta_dof[0]);
        cell_values.push_back(delta);
        num_cell_dofs++;
    }

    delta_SUPG->vector()->zero();
    delta_SUPG->vector()->set(cell_values.data(),
            num_cell_dofs, cell_dofs.data());
    delta_SUPG->vector()->apply("insert");
    delta_SUPG->update();


}

/**
 *  \fn save(Function& velocity,
 *           Function& pressure,
 *           Function& s_displacement,
 *           Function& somethingelse
 *           double t)
 *
 *	\brief	Output solution
 *	\param  U: solution
 *	        t: time
 *
 *  \author Lu WANG
 *  \date   03/12/2014
 *
 */
void SolverBase::save(const Function& velocity,
        const Function& pressure,
        const Function& displacement,
        const Function& S,
        double t)
{
    if(t < DOLFIN_EPS)
    {
        *out_velocity << velocity;
        *out_pressure << pressure;
        *out_displacement << displacement;
        *out_S << S;
    }

    while(last_sample + sample_period <= t)
    {
        last_sample = std::min(t, last_sample + sample_period);

        *out_velocity << velocity;
        *out_pressure << pressure;
        *out_displacement << displacement;
        *out_S << S;
    }
}


/**
 *  \fn solve()
 *
 *  \brief	Solve the whole FSI problem (in general a mixed system)
 *
 *  \author Lu WANG
 *  \date   03/12/2014
 */
int SolverBase::solve()
{
    init();
    //if (dolfin::MPI::process_number() == 0) {
    std::cout << "System solving begin ... \n";
    //}

    #ifdef _TIMER
    double pde_timer, pde_timer2;
    double step_timer, step_timer2;
    pde_timer = time();
    #endif
    save(*v, *p, *u, *ps,t);
    while(t < T)
    {
        //if (dolfin::MPI::process_number() == 0) {
        std::cout << "\n\x1b[33;1m=======================\x1b[0m\n";
        std::cout << "\x1b[33;1m solve for t = " << t << "\x1b[0m\n";
        std::cout << "\x1b[33;1m=======================\x1b[0m\n";
        //}
        #ifdef _TIMER
        step_timer = time();
        #endif

        prepareStep();
        step();
        save(*v, *p, *u, *ps,t);

        t += dt;
        #ifdef _TIMER
        step_timer2 = time();
        //if (dolfin::MPI::process_number() == 0) {
        std::cout << "Step timer for t="<< t << ": " << step_timer2 - step_timer << "\n";
        // }
        #endif
    }

    #ifdef _TIMER
    pde_timer2 = time();
    std::cout << "Total PDE timer: " << pde_timer2 - pde_timer << "\n";
    #endif
    std::cout << "\nSystem solved!\n";


    return 0;
}

/**
 *  \fn prepareStep()
 *
 *	\brief	Prepare for step()
 *
 *  \author Lu WANG
 *  \date   03/12/2014
 *
 */
void SolverBase::prepareStep()
{
    //LZ std::cout << "prepares step" <<std::endl;
    *(a0->vector()) *= -1.0;           //LZ For Newmark
    *(u0->vector()) = *(v->vector());  //LZ For temporary use 2/dt*(v-v0)
    *(u0->vector()) *= 2.0/dt;         //LZ For temporary use 2/dt*(v-v0)
    *(a0->vector()) += *(u0->vector());//LZ For temporary use 2/dt*(v-v0)
    *(u0->vector()) = *(v0->vector()); //LZ For temporary use 2/dt*(v-v0)
    *(u0->vector()) *= -2.0/dt;        //LZ For temporary use 2/dt*(v-v0)
    *(a0->vector()) += *(u0->vector());//LZ For temporary use 2/dt*(v-v0)

    *(v1->vector()) = *(v0->vector()); //LZ add for second order time scheme
    *(v0->vector()) = *(v->vector());
    *(u0->vector())  = *(u->vector());
    *(uA0->vector()) = *(uA->vector());

    std::cout << "prepares step" <<std::endl;

}

void SolverBase::pressure_dirichlet()
{
    //std::shared_ptr<const FunctionSpace> _sub_Q = (*VQ)[1];
    const GenericDofMap& dofmap = *(VQQ->dofmap());
    const GenericDofMap& dofmap0 = *(_sub_V->dofmap());
    const GenericDofMap& dofmap1 = *(_sub_Q->dofmap());
    const GenericDofMap& dofmap2 = *(_sub_Qs->dofmap()); //LZ added for J=1

    // std::size_t d = _euler_mesh->geometry().dim();
    std::size_t num_dof_v=0,num_dof_p_blood=0;

    std::vector<bool> is_dirichlet_bc(dofmap.global_dimension(), true);//LZ for pressure

    // Iterate over cells
    // Create map from cells attached to boundary to local dofs.
    for (CellIterator cell(*mesh); !cell.end(); ++cell)
    {
        std::size_t dim_marker = (*cell_mark)[cell->index()];
        if ((dim_marker == DomainType::BLOOD) || (dim_marker == DomainType::BLOOD2) || (dim_marker == DomainType::STENT)) { 
            const ArrayView<const dolfin::la_index> cell_dofs1 =
                dofmap1.cell_dofs(cell->index());
            for (std::size_t i=0; i<cell_dofs1.size(); ++i) {
                is_dirichlet_bc[cell_dofs1[i]] = false; 
                ++num_dof_p_blood;
            }
        }
        //LZ added pressure dirichlet for hyperelastic
        if (dim_marker == DomainType::VESSEL) { 
            const ArrayView<const dolfin::la_index> cell_dofs2 =
                dofmap2.cell_dofs(cell->index());
            for (std::size_t i=0; i<cell_dofs2.size(); ++i) {
                is_dirichlet_bc[cell_dofs2[i]] = false; 
                ++num_dof_p_blood;
            }
        }

        const ArrayView<const dolfin::la_index> cell_dofs0 =
            dofmap0.cell_dofs(cell->index());
        for (std::size_t i=0; i<cell_dofs0.size(); ++i) {
            is_dirichlet_bc[cell_dofs0[i]] = false; 
            ++num_dof_v;
        }
    }

    for (dolfin::la_index i=0; i<is_dirichlet_bc.size(); ++i) {
        if (is_dirichlet_bc[i]) {
            pressure_dirichlet_dofs.push_back(i);
        }
    }
    pressure_dirichlet_values.resize(pressure_dirichlet_dofs.size(), 0.0);
    /*std::cout   << "num_dof_v: " << num_dof_v<<std::endl
      << "num_dof_p_blood: " << num_dof_p_blood<<std::endl
      << "num_dof_total: " << dofmap.global_dimension()<<std::endl
      << "num_dof_p_vessel: " << pressure_dirichlet_dofs.size()<<std::endl;
      std::getchar();
      */
}

//-----------------------------------------------------------------------------
// Implementation of NonlinearDiscreteProblem
//-----------------------------------------------------------------------------
SolverBase::NonlinearDiscreteProblem::
    NonlinearDiscreteProblem(std::shared_ptr<SolverBase> problem)
: _problem(problem)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
SolverBase::NonlinearDiscreteProblem::~NonlinearDiscreteProblem()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void SolverBase::
NonlinearDiscreteProblem::F(GenericVector& b, const GenericVector& x)
{
    // Get problem data
    //dolfin_assert(_problem);
    std::cout << "assemble residual vector" <<std::endl;
    std::shared_ptr<const Form> F_s(_problem->F_s);
    std::shared_ptr<const Form> F_f(_problem->F_f);
    std::vector<std::shared_ptr<const DirichletBC>>& bcs(_problem->_bcsFSI);
    std::vector<dolfin::la_index>& dirichlet_dofs = _problem->pressure_dirichlet_dofs;
    std::vector<double>& dirichlet_values= _problem->pressure_dirichlet_values;
    Vector tempb;
    std::vector<double>& coord0 = (_problem->mesh0->coordinates());
    std::vector<double>& coord = (_problem->mesh->coordinates());
    std::vector<double> backup_coord = (_problem->mesh->coordinates());

    // Assemble right-hand side
    //dolfin_assert(F);
    Assembler assembler;
    assembler.finalize_tensor = false;
    assembler.assemble(b, *F_f);

    assembler.finalize_tensor = true;
    assembler.add_values = true;
    coord = coord0;
    assembler.assemble(b, *F_s);
    coord = backup_coord;
    //b += tempb;




    // Apply boundary conditions
    for (std::size_t i = 0; i < bcs.size(); i++)
    {
        dolfin_assert(bcs[i]);
        bcs[i]->apply(b, x);
    }

    /// Modify RHS vector (b[i] = value) and apply changes
    b.set(dirichlet_values.data(), dirichlet_dofs.size(), dirichlet_dofs.data());
    b.apply("insert");

}
//-----------------------------------------------------------------------------
void SolverBase::NonlinearDiscreteProblem::J(GenericMatrix& A,
        const GenericVector& x)
{

    std::cout << "assemble momentum matrix" <<std::endl;
    // Get problem data
    //dolfin_assert(_problem);
    std::shared_ptr<const Form> J_s(_problem->J_s);
    std::shared_ptr<const Form> J_f(_problem->J_f);
    std::vector<std::shared_ptr<const DirichletBC>>& bcs(_problem->_bcsFSI);
    std::vector<dolfin::la_index>& dirichlet_dofs = _problem->pressure_dirichlet_dofs;
    std::vector<double>& dirichlet_values= _problem->pressure_dirichlet_values;
    //Matrix tempA, tempB;
    std::vector<double>& coord0 = (_problem->mesh0->coordinates());
    std::vector<double>& coord = (_problem->mesh->coordinates());
    std::vector<double> backup_coord = (_problem->mesh->coordinates());

    // Assemble right-hand side
    //dolfin_assert(F);
    Assembler assembler;
    assembler.finalize_tensor = false;
    assembler.assemble(A, *J_f);

    assembler.finalize_tensor = true;
    assembler.add_values = true;
    coord = coord0;
    assembler.assemble(A, *J_s);
    coord = backup_coord;


    // Apply boundary conditions
    for (std::size_t i = 0; i < bcs.size(); i++)
    {
        dolfin_assert(bcs[i]);
        bcs[i]->apply(A);
    }

    /// Modify linear system (A_ii = 1) and apply changes, diag identity
    A.ident(dirichlet_dofs.size(), dirichlet_dofs.data());
    A.apply("insert");

    //   info(A, true);
}
//-----------------------------------------------------------------------------

/**
 * end of file 
 *
 */
