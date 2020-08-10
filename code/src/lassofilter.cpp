#include "lassofilter.h"

#include "polynomial.h"
#include "quadraturegrid.h"

LassoFilter::LassoFilter( Settings* settings ) : Filter( settings ) {
    unsigned maxDegree = _settings->GetMaxDegree();
    _filterParam       = Vector( _settings->GetNTotal(), 1.0 );
    _l1Norms           = Vector( _settings->GetNTotal(), 0.0 );

    std::vector<Polynomial*> quad;
    QuadratureGrid* quadGrid;
    unsigned nQuadPoints = _settings->GetNQuadPoints();

    // initialize classes: Basis Functions and Quadrature rules are defined for Legendre and Hermite
    quad.resize( 2 );
    quad[0] = Polynomial::Create( _settings, nQuadPoints, DistributionType::D_LEGENDRE );
    quad[1] = Polynomial::Create( _settings, nQuadPoints, DistributionType::D_HERMITE );

    // get number of uncertain dimensions
    unsigned numDimXi = _settings->GetNDimXi();

    // get index vector for basis function calculation
    std::vector<std::vector<unsigned>> indices = _settings->GetPolyIndices();
    unsigned nTotal                            = _settings->GetNTotal();

    std::vector<Vector> xiGrid;
    std::vector<Vector> wGrid;
    std::vector<Vector> quadNodes;
    VectorU nQTotalForRef;
    VectorU nTotalForRef;
    Matrix phiTildeTrans;    // stores scaled basis functions evaluated at quadrature points
    Matrix phiTildeF;        // stores scaled basis functions evaluated at quadrature points times pdf
    Matrix phiTildeWf;       // stores scaled basis functions evaluated at quadrature points times weight and pdf
    std::vector<Vector> phiTildeVec;

    // define quadrature
    wGrid.resize( _settings->GetNRefinementLevels() );
    nQTotalForRef.resize( _settings->GetNRefinementLevels() );
    VectorU quadLevel  = _settings->GetQuadLevel();
    unsigned oldQLevel = 0;
    quadGrid           = nullptr;
    for( unsigned rlevel = 0; rlevel < _settings->GetNRefinementLevels(); ++rlevel ) {
        if( oldQLevel != quadLevel[rlevel] ) {    // new quadrature weights needed
            if( quadGrid ) delete quadGrid;
            quadGrid = QuadratureGrid::Create( _settings, quadLevel[rlevel] );
        }
        wGrid[rlevel]         = quadGrid->GetWeights();
        nQTotalForRef[rlevel] = wGrid[rlevel].size();
        oldQLevel             = quadLevel[rlevel];
    }
    xiGrid = quadGrid->GetNodes();
    // give number of quad points per level to settings. Quadrature variables are computed here.
    _settings->SetNQTotalForRef( nQTotalForRef );

    // set total number of quadrature points
    unsigned nQTotal = quadGrid->GetNodeCount();

    // compute basis functions evaluated at the quadrature points
    Matrix phiTilde = Matrix( nQTotal, nTotal, 1.0 );

    unsigned n = 0;
    for( unsigned k = 0; k < nQTotal; ++k ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned l = 0; l < numDimXi; ++l ) {
                if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                phiTilde( k, i ) *= quad[n]->Evaluate( indices[i][l], xiGrid[k][l] ) / quad[n]->L2Norm( indices[i][l] );    // sqrt( 2.0 * i + 1.0 );
                phiTildeF( k, i ) *=
                    quad[n]->Evaluate( indices[i][l], xiGrid[k][l] ) / quad[n]->L2Norm( indices[i][l] ) * quad[n]->fXi( xiGrid[k][l] );
            }
            // P(xi \in I_l) probability that xi lies in multi element I_l
            // double P            = 1.0 / _nMultiElements;
            //_phiTildeF( k, i )  = _phiTildeF( k, i ) / P;    // modify pdf to multielement ansatz, only valid for uniform distributions
            phiTildeWf( k, i ) = phiTildeF( k, i ) * wGrid[_settings->GetNRefinementLevels() - 1][k];
            phiTildeVec[k][i]  = phiTilde( k, i );
        }
    }

    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned k = 0; k < _settings->GetNQTotal(); ++k ) {
            _l1Norms[i] += std::fabs( phiTildeWf( k, i ) );
        }
    }

    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
            _filterParam[i] *= index * ( index + 1 );
        }
    }
}

void LassoFilter::SetupFilter() {}

LassoFilter::~LassoFilter() {}

void LassoFilter::FilterMoments( Tensor& u ) const {
    double scL1, filterStrength;
    unsigned nMax = _settings->GetNTotal() - 1;
    for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
        for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
            double uLastMoment = u( s, l, nMax );
            filterStrength     = std::fabs( uLastMoment ) / ( _filterParam[nMax] * _l1Norms[nMax] );
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                scL1 = 1.0 - filterStrength * _filterParam[i] * _l1Norms[i] / std::fabs( u( s, l, i ) );
                if( scL1 < 0 || std::fabs( u( s, l, i ) ) < 1e-7 ) scL1 = 0.0;
                u( s, l, i ) = scL1 * u( s, l, i );
            }
        }
    }
}

void LassoFilter::FilterMoments( Tensor& v, const Tensor& u ) const {
    double scL1, filterStrength;
    unsigned nMax = _settings->GetNTotal() - 1;
    for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
        for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
            double uLastMoment = u( s, l, nMax );
            filterStrength     = std::fabs( uLastMoment ) / ( _filterParam[nMax] * _l1Norms[nMax] );
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                scL1 = 1.0 - filterStrength * _filterParam[i] * _l1Norms[i] / std::fabs( u( s, l, i ) );
                if( scL1 < 0 || std::fabs( u( s, l, i ) ) < 1e-7 ) scL1 = 0.0;
                v( s, l, i ) = scL1 * u( s, l, i );
            }
        }
    }
}

void LassoFilter::FilterMoments( Matrix& v, const Tensor& u, unsigned l ) const {
    double scL1, filterStrength;
    unsigned nMax = _settings->GetNTotal() - 1;
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        double uLastMoment = u( s, l, nMax );
        filterStrength     = std::fabs( uLastMoment ) / ( _filterParam[nMax] * _l1Norms[nMax] );
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            scL1 = 1.0 - filterStrength * _filterParam[i] * _l1Norms[i] / std::fabs( u( s, l, i ) );
            if( scL1 < 0 || std::fabs( u( s, l, i ) ) < 1e-7 ) scL1 = 0.0;
            v( s, i ) = scL1 * u( s, l, i );
        }
    }
}
