/** \file meshsmoothing.h
 * \brief the class for mesh smoothing
 * \author Lu WANG
 * \date 09/13/2013
 * \copyright GNU Public License.
 *
 * Modify by Shihua Gong on 05/05/2016
 *
 */

#ifndef __MESHSMOOTHING_H__
#define __MESHSMOOTHING_H__

#include <dolfin.h>
#include "common.h"
//#include <boost/timer.hpp>

using namespace dolfin;
//namespace FASPFSI {

/**
 * \class MeshSmoothing
 *	@brief	A class for Mesh smoothing
 */
class MeshSmoothing

{
    public:
        /**
         *  \fn MeshSmoothing()
         *
         *	\brief	Class Constructor
         *
         */
        MeshSmoothing(std::shared_ptr<Mesh> mesh0, std::shared_ptr<MeshFunction<std::size_t>> cell_marker0,
                std::shared_ptr<MeshFunction<std::size_t>> bd_marker0, std::vector<size_t> fixed_bd,
                std::string smooth);

        /**
         * \fn ~ProblemBase()
         *
         *	\brief	Destructor
         *
         */
        ~MeshSmoothing() {};

        /**
         * \fn void meshmoothing(Mesh& Fmesh, Mesh& Smesh,const Function& u_s,std::vector<std::vector<std::size_t> > FSmapping)
         *
         * \brief Move the coordinates of mesh according to displacement function and smoothing the mesh position by 2D harmonic equation.
         *
         * \param mesh          A Mesh object need to be modified.
         * \param Mmark         mesh mark for fluid and solid
         * \param u_s           A Function discribe the boundary moves.
         * \param uA            GenericVector of the smoothed result
         *
         */
       /* void smoothing(std::shared_ptr<Mesh> mesh,
                const Function& u_s,
                Function& uA);*/ //LZ
        void smoothing(std::shared_ptr<Mesh> mesh,
                Function& u_s,
                Function& uA);
        /**
         * \fn void faspsolve(GenericMatrix &A,GenericVector &b,GenericVector &u)
         *
         * \brief Solve a mesh problem by fasp
         *
         * \param A          matrix
         * \param b          right hand side
         * \param u          the solution
         *
         */
        //void faspsolve(GenericMatrix &A,GenericVector &b,GenericVector &u);
    protected:

        std::shared_ptr<Mesh> _mesh0;
        std::shared_ptr<MeshFunction<std::size_t>> _cell_marker0;
        std::shared_ptr<MeshFunction<std::size_t>> _bd_marker0;
        std::shared_ptr<FunctionSpace> V;
        std::shared_ptr<Form> form;
        Assembler assembler;
        std::shared_ptr<Matrix> A;
        Function *rhs;
        std::vector<size_t> _fixed_bd;

        std::vector<const DirichletBC*>* _bcsMesh;

        std::string smooth_method;
        std::size_t D;
        std::size_t d;
        std::shared_ptr<Constant> fixbd;
        std::shared_ptr<Constant> Cmaxt;
        std::shared_ptr<Constant> Cmint;
        double maxt;
        double mint;
        double mu_m;
        double lambda_m;

        //input_param     inparam;  // parameters from input files
        //itsolver_param  itparam;  // parameters for itsolver
        //AMG_param       amgparam; // parameters for AMG
        //ILU_param       iluparam; // parameters for ILU
        //Schwarz_param   schparam; // parameters for Schwarz
};
//}
#endif
