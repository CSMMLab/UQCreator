#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <chrono>
#include <fstream>
#include <omp.h>
#include <spdlog/spdlog.h>

#include "closure.h"
#include "mesh.h"
#include "polynomial.h"
#include "problem.h"
#include "settings.h"
#include "timesolver.h"
#include "typedefs.h"

class MomentSolver
{
  private:
    Settings* _settings;                       // settings class
    Closure* _closure;                         // closure class for computing duals, defines u(Lambda), du(Lambda)
    Mesh* _mesh;                               // specified mesh
    TimeSolver* _time;                         // time solver for evolving moment vector in time
    MatTens _lambda;                           // dual matrix of dimension nCells x states x total number of moments
    Problem* _problem;                         // specified problem defines right hand side and initial condition
    double _dt, _tStart, _tEnd;                // timestep, start and end time
    unsigned _nCells;                          // number of spatial cells
    unsigned _nStates;                         // number of states of the original system
    unsigned _nMultiElements;                  // total number of multi elements
    unsigned _nQuadPoints;                     // number of moments in one uncertain dimension
    unsigned _nQTotal;                         // total number of quad points
    VectorU _nQTotalForRef;                    // total number of quad points for each refinement level
    unsigned _nTotal;                          // total number of moments
    std::shared_ptr<spdlog::logger> _log;      // log file writer
    std::vector<Vector> _referenceSolution;    // reference solution stores expected value and variance for all states
    VectorU _nTotalForRef;                     // nTotal for different refinement levels
    std::vector<unsigned> _cellIndexPE;        // gives cell index on each PE

    /**
     * numerical flux for moment system
     * @param output moment matrix
     * @param values deterministic fluxes
     * @param refinement level
     */
    void numFlux( Matrix& out, const Matrix& g, unsigned level );
    /**
     * source term for moment system
     * @param Moments for in- and output
     * @param solution at quadrature points
     * @param time step size
     * @param current time
     * @param refinement level
     */
    void Source( MatTens& uQNew, const MatTens& uQ, double dt, double t, const VectorU& refLevel ) const;
    /**
     * sets up moments for specified initial condition
     * @return output moment matrix
     */
    MatTens SetupIC() const;
    /**
     * exports moments and duals
     * @param moment matrix for export
     * @param dual matrix for export
     */
    void Export( const MatTens& u, const MatTens& lambda ) const;
    /**
     * imports previous settings from restart configfile if specified
     * @return pointer to previous settings
     */
    Settings* ImportPrevSettings() const;
    /**
     * imports previous moments from restart momentfile if specified
     * @param number of moments for previous settings
     * @return previous moment matrix
     */
    MatTens ImportPrevMoments( unsigned nPrevTotal ) const;
    /**
     * imports previous duals from restart momentfile if specified
     * @param number of moments/duals for previous settings
     * @return previous dual matrix
     */
    MatTens ImportPrevDuals( unsigned nPrevTotal );
    /**
     * compute distance of computed mean and variance to reference solution
     * @param export matrix after computation
     * @param norm used for error (can be 1 or 2)
     * @param defines left bottom point of rectangle in which error is computed
     * @param defines right top point of rectangle in which error is computed
     * @return error variance for different states
     */
    Vector CalculateError( const Matrix& solution, unsigned LNorm, const Vector& a, const Vector& b ) const;
    /**
     * loads moments from restart file or computes moments fron initial condtion
     * @param number of moments
     * @return moment matrix
     */
    MatTens DetermineMoments( unsigned nTotal ) const;
    /**
     * writes duals on _lambda and recomputes moment matrix for specified duals
     * @param previous settings
     * @param previous closure
     * @param previous moment matrix
     */
    void SetDuals( Settings* prevSettings, Closure* prevClosure, MatTens& u );
    /**
     * loads previous settings or current settings if no restart file specified
     * @return previous settings
     */
    Settings* DeterminePreviousSettings() const;
    /**
     * loads previous closure or current closure if no restart file specified
     * @return previous closure
     */
    Closure* DeterminePreviousClosure( Settings* prevSettings ) const;
    /**
     * computes distance of computed mean and variance to reference solution on every cell
     * @return error mean and variance on every spatial cell
     */
    Matrix CalculateErrorField( const Matrix& solution, unsigned LNorm ) const;

    void WriteErrors( const VectorU& refinementLevel );

    Matrix WriteMeanAndVar( const VectorU& refinementLevel, double t, bool writeExact ) const;
    void ExportRefinementIndicator( const VectorU& refinementLevel, const MatTens& u, unsigned index ) const;
    double ComputeRefIndicator( const VectorU& refinementLevel, const Tensor& u, unsigned refLevel ) const;
    void PerformInitialStep( const VectorU& refinementLevel, MatTens& u );
    void DetermineGradients( MatTens& duQx, MatTens& duQy, const MatTens& uQ, const VectorU& refLevel ) const;
    void DetermineGradientsScalarField( Matrix& dux, Matrix& duy, const Matrix& u ) const;
    void WriteGradientsScalarField( const Matrix& u ) const;
    void Write2ndDerMeanAndVar( const Matrix& meanAndVar ) const;

  public:
    /**
     * constructor
     * @param settings class
     * @param specified mesh
     * @param specified problem
     */
    MomentSolver( Settings* settings, Mesh* mesh, Problem* problem );
    /**
     * destructor
     */
    ~MomentSolver();
    /**
     * computed time evolution of moment system for chosen problem until specified residual or end time is reached
     */
    void Solve();
};

#endif    // MOMENTSOLVER_H
