/** \file meshsmoothing.cpp
 * \brief impliment of the meshing smoothing class.
 *  \author Lu WANG
 *  \date   03/12/2014
 *
 * Modify by Shihua Gong 30-04-2016
 *  
 * \copyright GNU Public License.
 */

#include <dolfin.h>
//#include <iostream>
//#include <fstream>
//#include <string>
#include "dolfin/ale/ALE.h"

#include "meshsmoothing.h"
#include "ALE/VectorPoisson.h"
#include "ALE/Elastic.h"
#include "ALE/VectorPoisson3D.h"
#include "ALE/Elastic3D.h"

using namespace std;
using namespace dolfin;
//using namespace FASPFSI;

//#define _DEBUG

/**
 *  \fn MeshSmoothing()
 *
 *	\brief	Class Constructor
 *
 *  \author Lu WANG
 *  \date   03/22/2014
 *
 */
MeshSmoothing::MeshSmoothing(std::shared_ptr<Mesh> mesh0, std::shared_ptr<MeshFunction<std::size_t>> cell_marker0,
        std::shared_ptr<MeshFunction<std::size_t>> bd_marker0, std::vector<size_t> fixed_bd,
        std::string smooth): _mesh0(mesh0), _cell_marker0(cell_marker0), _bd_marker0(bd_marker0),
    _fixed_bd(fixed_bd), smooth_method(smooth) 
{
    A.reset(new Matrix);
    D = _mesh0->topology().dim();
    d = _mesh0->geometry().dim();

    if (smooth_method == "elasticsmoothing")
    {
        maxt = -1.0;
        mint = 1e+6;

        for (CellIterator c(*_mesh0); !c.end(); ++c)
        {
            if (c->volume() > maxt) maxt = c->volume();
            if (c->volume() < mint) mint = c->volume();
        }

        Cmaxt = std::make_shared<Constant>(maxt);
        Cmint = std::make_shared<Constant>(mint);

        if (d == 2) {
            V    = std::make_shared<Elastic::FunctionSpace>(_mesh0);
            form = std::make_shared<Elastic::BilinearForm>(V, V);
            form->set_coefficient("max_t",Cmaxt);
            form->set_coefficient("min_t",Cmint);
            fixbd = std::make_shared<Constant>(0.0,0.0);
        }
        else if (d == 3) {
            V    = std::make_shared<Elastic3D::FunctionSpace>(_mesh0);
            form = std::make_shared<Elastic3D::BilinearForm>(V, V);
            form->set_coefficient("max_t",Cmaxt);
            form->set_coefficient("min_t",Cmint);
            fixbd = std::make_shared<Constant>(0.0,0.0,0.0);
        }
    }
    else if (smooth_method == "harmonicsmoothing")
    {
        if (d == 2) {
            V     =  std::make_shared<VectorPoisson::FunctionSpace>(_mesh0);
            form  = std::make_shared<VectorPoisson::BilinearForm>(V, V);
            fixbd = std::make_shared<Constant>(0.0,0.0);
        }
        else if (d == 3) {
            V     = std::make_shared<VectorPoisson3D::FunctionSpace>(_mesh0);
            form  = std::make_shared<VectorPoisson3D::BilinearForm>(V, V);
            fixbd = std::make_shared<Constant>(0.0,0.0,0.0);
        }
    }
    rhs = new Function(V);
    rhs->vector()->zero();
    assembler.assemble(*A, *form);

    _bcsMesh = new std::vector<const DirichletBC*>(_fixed_bd.size());
    for (int i=0; i<_fixed_bd.size(); i++)
    {
        (*_bcsMesh)[i] = new DirichletBC(V, fixbd, _bd_marker0, _fixed_bd[i]);
        (*_bcsMesh)[i]->apply(*A,*(rhs->vector()));
    }

    //char input[] = "ini/input.dat";
    //fasp_param_input(input,&inparam);
    //fasp_param_init(&inparam, &itparam, &amgparam, &iluparam, &schparam);

}

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
 *  \author Lu WANG
 *  \date   03/12/2014
 *
 */


/*void MeshSmoothing::smoothing(std::shared_ptr<Mesh> mesh,
        const Function& u_s,
        Function& uA)*///LZ
void MeshSmoothing::smoothing(std::shared_ptr<Mesh> mesh,
        Function& u_s,
        Function& uA)
{

    // RHS vector
    const std::size_t num_vertices = _mesh0->num_vertices();

    // Dof range
    const dolfin::la_index n0 = V->dofmap()->ownership_range().first;
    const dolfin::la_index n1 = V->dofmap()->ownership_range().second;
    const dolfin::la_index num_dofs = n1 - n0;

    // Mapping of mesh vertex numbers to dofs (including ghost dofs)
    const std::vector<dolfin::la_index> dof_to_vertex_map =
        vertex_to_dof_map(*V);//V->dofmap()->dof_to_vertex_map(*_mesh);

    // Array of all dofs (including ghosts) with global numbering
    std::vector<dolfin::la_index> all_global_dofs(d*num_vertices);
    for (std::size_t i = 0; i < num_vertices; i++)
        for (std::size_t j = 0;j < d;j++)
            all_global_dofs[d*i+j] = dof_to_vertex_map[d*i+j] + n0;

    // Create arrays for setting bcs.
    // Their indexing does not matter - same ordering does.
    std::size_t num_boundary_dofs = 0;
    std::vector<dolfin::la_index> boundary_dofs;
    std::vector<double> boundary_values;

    // Update ghosts dofs
    u_s.update();

    const FiniteElement element = *(u_s.function_space()->element());
    // Get restriction if any
    //std::shared_ptr<const Restriction> restriction = u_s.function_space()->dofmap()->restriction();

    // Local data for interpolation on each cell
    const std::size_t num_cell_vertices = _mesh0->type().num_vertices(D);

    // Compute in tensor (one for scalar function, . . .)
    const std::size_t value_size_loc = u_s.value_size();

    // Create vector to hold cell vertex values
    std::vector<double> cell_vertex_values(value_size_loc*num_cell_vertices);

    // Create vector for expansion coefficients
    std::vector<double> coefficients(element.space_dimension());

    ufc::cell ufc_cell;
    std::vector<double> vertex_coordinates;

    std::vector<double> u_s_values;
    u_s.compute_vertex_values(u_s_values, *mesh);
    u_s.set_allow_extrapolation(true);//LZ In order to fix an unknown error

    for (CellIterator cell(*_mesh0); !cell.end(); ++cell)
    {
        // Get cell vertex coordinates
        cell->get_vertex_coordinates(vertex_coordinates);

        //if (restriction && !restriction->contains(*cell))
        //    continue;

        std::size_t meshval = (*_cell_marker0)[*cell];
        //LZ added BLOOD2 for Stent
        if ((meshval != DomainType::BLOOD) && (meshval != DomainType::BLOOD2) && (meshval != DomainType::STENT)) {
            cell->get_cell_data(ufc_cell);

            // Pick values from global vector
            u_s.restrict(&coefficients[0], element, *cell,vertex_coordinates.data(),ufc_cell);
            const int cell_orientation = 0;
            element.interpolate_vertex_values(&cell_vertex_values[0],
                    &coefficients[0],
                    vertex_coordinates.data(),
                    cell_orientation,
                    ufc_cell);

            // Tabulate dofs on cell
            for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
            {
                for (std::size_t i = 0; i < d; ++i)
                {
                    const std::size_t local_index  = vertex.pos()*d + i;
                    const std::size_t global_index = i*num_vertices+vertex->index();

                    const dolfin::la_index dof = dof_to_vertex_map[d*vertex->index()+i];
                    // Skip ghosts
                    if (dof >= 0 && dof < num_dofs)
                    {
                        // Global dof numbers
                        boundary_dofs.push_back(dof + n0);
                        boundary_values.push_back(u_s_values[global_index]);
                        num_boundary_dofs++;
                    }
                }
            }

        }
    }

    // Modify matrix (insert 1 on diagonal)
    A->ident(num_boundary_dofs, boundary_dofs.data());
    A->apply("insert");

    // Store bc into RHS and solution so that CG solver can be used
    rhs->vector()->set(boundary_values.data(), num_boundary_dofs, boundary_dofs.data());
    //rhs.vector()->set(boundary_values0.data(), num_boundary_dofs0, boundary_dofs0.data());
    rhs->vector()->apply("insert");
    rhs->update();

    //File testrhs("rhs.pvd");
    //testrhs << rhs;

    uA.vector()->zero();

    // Solve system
    // Pick amg as preconditioner if available
    const std::string prec(has_krylov_solver_preconditioner("amg")
            ? "amg" : "default");
    KrylovSolver solver(A,"gmres",prec);

    solver.parameters["relative_tolerance"] = 1.0e-4;
    solver.parameters["absolute_tolerance"] = 1.0e-15;
    solver.parameters["divergence_limit"] = 1.0e1;
    solver.parameters["maximum_iterations"] = 500;
    solver.parameters["error_on_nonconvergence"] = false;
    solver.parameters["nonzero_initial_guess"] = false;

    solver.solve(*(uA.vector()), *(rhs->vector()));

    //if (MPI::num_processes() > 1)
    uA.update();

    ALE::move(*mesh, uA);

    /*
    // Modify mesh coordinates
    std::vector<double> vertex_values;
    uA.compute_vertex_values(vertex_values, *_mesh0);

    //std::vector<double>& geometry = mesh.geometry().x();
    MeshGeometry& geometry  = mesh->geometry();
    MeshGeometry& geometry0 = _mesh0->geometry();

    std::vector<double> coord(d);
    for(VertexIterator v(*mesh); !v.end(); ++v)
    {
        for(unsigned int i = 0; i < d; i++)
        {
            Vertex& vertex = *v;
            geometry.x(vertex.index(), i) += vertex_values[i*num_vertices + vertex.index()];
        }
    }
    */

    //exit(1);
    //dolfin::MPI::barrier();
}

