/** \file solverbase.h
 *  \brief the class for FSI problem
 *  \author Feiteng Huang
 *  \date 09/17/2013
 *  \copyright GNU Public License.
 *
 *  Modify by Lu Wang
 *  Modify by Shihua Gong
 *
 */

#ifndef __SOLVER_BASE_H__
#define __SOLVER_BASE_H__

#include <dolfin.h>
#include "meshsmoothing.h"

using namespace dolfin;

/**
 * \class SolverBase
 *	@brief	A class for FSI problem
 */
class SolverBase// : public ProblemBase
{
    public:

        /**
         *  \fn SolverBase::SolverBase(Mesh& mesh, MeshFunction<std::size_t>& mark, std::string pfile)
         *
         *	\brief	Class Constructor
         *
         */
        SolverBase(std::shared_ptr<Mesh> mesh0, std::shared_ptr<MeshFunction<std::size_t>> cell_mark0,
                std::shared_ptr<MeshFunction<std::size_t>> bd_mark0, std::string pfile_fsi);

        /**
         *  \fn ~SolverBase()
         *
         *	\brief	Destructor
         *
         */
        virtual ~SolverBase();

        /**
         *  \fn init()
         *
         *	\brief	Initialization of the spaces and coefficients
         *
         */
        void init();

        /**
         *  \fn prepareStep()
         *
         *	\brief	Prepare for step()
         *
         */
        virtual void prepareStep();

        /**
         *  \fn step()
         *
         *	\brief	Step from t^n to t^{n+1}
         *
         */
        virtual void step();

        /**
         *  \fn solve()
         *
         *	\brief	Solve the whole FSI problem (in general a mixed system)
         */
        virtual int solve();

        /**
         *   \fn save(Function& f_velocity,
         *            Function& f_pressure,
         *            Function& s_velocity,
         *            Function& s_displacement,
         *            double t)
         *
         *	\brief	Output solution
         *
         */
        virtual void save(const Function& velocity,
                const Function& pressure,
                const Function& displacement,
                const Function& S,
                double t);

        /**
         * @brief Read parameters
         *
         * @param para_file
         */
        void read_parameters(std::string para_file);


    protected:


        // mesh
        std::shared_ptr<Mesh> mesh0;
        std::shared_ptr<Mesh> mesh;
        std::shared_ptr<MeshFunction<std::size_t>> bd_mark0, cell_mark0, bd_mark, cell_mark; 

        // FunctionSpace
        std::shared_ptr<FunctionSpace> VQQ, V, Q, Qs, DG, DG1;
        std::shared_ptr<FunctionSpace> _sub_V, _sub_Q,_sub_Qs;

    public:
        // Weakform
        std::shared_ptr<Form> J_f;
        std::shared_ptr<Form> J_s;
        std::shared_ptr<Form> F_f;
        std::shared_ptr<Form> F_s;

    protected:

        // variable
        std::shared_ptr<Function> v_p_p;  //!< velocity for fluid mesh

        std::shared_ptr<Function> u;    //!< displacement
        std::shared_ptr<Function> u0;   //!< displacement at last time-step

        std::shared_ptr<Function> uA;   //!< displacement for smoothing
        std::shared_ptr<Function> uA0;  //!< displacement for smoothing at last time-step
        std::shared_ptr<Function> uA1;  //!< displacement for smoothing at last time-step
        std::shared_ptr<Function> w;    //!< displacement for smoothing at last time-step

        std::shared_ptr<Function> v;    //!< velocity
        std::shared_ptr<Function> v0;   //!< velocity at last time-step
        std::shared_ptr<Function> v1;   //!< velocity at last inner-iteration-step
        std::shared_ptr<Function> dv;

        std::shared_ptr<Function> p;    //!< pressure
        std::shared_ptr<Function> ps;    //!< pressure
        std::shared_ptr<Function> a0;


        // Common Coefficient
        std::shared_ptr<Constant> _dt;

        // Fluid Coefficient
        std::shared_ptr<Constant> _rho_f;
        std::shared_ptr<Constant> _mu_f;
        std::shared_ptr<Constant> _alpha;
        std::shared_ptr<Constant> _alpha2;
        std::shared_ptr<Function> delta_SUPG;


        // Elastic Coeffcient
        std::shared_ptr<Constant> _rho_s;
        std::shared_ptr<Constant> _lam_s;  //lamb_s = 2*mu_s*nu_s/(1-2*nu_s);
        std::shared_ptr<Constant> _nu_s;
        std::shared_ptr<Function> _mu_s;
        std::shared_ptr<Constant> _p_s;
        std::shared_ptr<Constant> _beta1;
        std::shared_ptr<Constant> _beta2;


        // 
        std::shared_ptr<MeshSmoothing> meshsmoother;


        // boundary conditions

    public:
        std::vector<std::shared_ptr<const DirichletBC>> _bcsFSI;

    protected:

        // Need to initialize outside
        std::vector<std::shared_ptr<GenericFunction>>* bdfuncs_FSI;
        std::vector<std::size_t>* bcsubDomain_FSI;
        std::vector<std::size_t>* bcs_fix_Domain;
        std::shared_ptr<GenericFunction> sigmaN;
        std::shared_ptr<GenericFunction> beta1Func;
        std::shared_ptr<GenericFunction> beta2Func;




        // Solver
        //KrylovSolver* ksolver;
        //LUSolver* lusolver;
        //PETScFASP* FASPsolver;
        //FaspSolver* fsolver;

        //std::string st;
        //std::string name;


        // auxiliary data
        Parameters _parameters;
        double dt, t, T;
        int iter_max, iter_FS_max;
        double tol;



        double sample_period;
        double last_sample;
        File* out_velocity;
        File* out_pressure;
        File* out_S;
        File* out_displacement;

        /// Record solid pressure dirichlet bcs 
    public:
        std::vector<dolfin::la_index> pressure_dirichlet_dofs;
        std::vector<double> pressure_dirichlet_values;

    protected:
        void pressure_dirichlet();
        void update_delta_SUPG();

        // Nonlinear (algebraic) problem
        class NonlinearDiscreteProblem : public NonlinearProblem
    {
        public:

            // Constructor
            NonlinearDiscreteProblem(
                    std::shared_ptr<SolverBase> problem);

            // Destructor
            ~NonlinearDiscreteProblem();

            // Compute F at current point x
            virtual void F(GenericVector& b, const GenericVector& x);

            // Compute J = F' at current point x
            virtual void J(GenericMatrix& A, const GenericVector& x);

        private:

            // Problem and solver objects
            std::shared_ptr<SolverBase> _problem;

    };

        std::shared_ptr<NonlinearDiscreteProblem> nonlinear_problem;

        // The Newton solver
        std::shared_ptr<NewtonSolver> newton_solver;

};

#endif
