/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_gsg_driver.cpp
    \brief A driver program for PECOS */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "CombinedSparseGridDriver.hpp"
#include "SharedProjectOrthogPolyApproxData.hpp"
#include "ProjectOrthogPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "CubatureDriver.hpp"
#include "pecos_data_types.hpp"
//#include "LocalRefinableDriver.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

#include "TestFunctions.hpp"

using namespace std;

#define VERB

#define MAX_CHARS_PER_LINE 1000
#define BTYPE              GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL
#define NUMVARS            3
#define STARTLEV           1
#define NITER              10
#define MAXORD             5
#define NQOI               1
#define GRIDFILE           "savedgrid.dat"
#define FCNFILE            "savedfcn.dat"


void restartGSGdriver(const char *grid, const char *fcnvals, 
                      Pecos::RealMatrix &storedSets, Pecos::RealVector
                      &storedVals) ;
int checkSetsInStoredSet(const Pecos::RealMatrix &storedSets, const Pecos::RealMatrix &variable_sets, 
			 std::vector<bool> &computedGridIDs);
void addNewSets(Pecos::RealMatrix &storedSets,RealVector &storedVals,
		Pecos::RealMatrix &newSets, RealVector &fev, std::vector<bool> &computedGridIDs);

void write_USAS(std::ostream& s, const Pecos::UShortArraySet &a);
void write_US2A(std::ostream& s, const Pecos::UShort2DArray  &a);
RealMatrix feval(const RealMatrix &dataMat, const int nQoI, std::vector<bool> &computedGridIDs, void *funInfo)  ;

/// A driver program for PECOS.
int main(int argc, char* argv[])
{

  using namespace Pecos;

  // Define defaults
  int nQoI = NQOI;
  short btype = (short) BTYPE;

  // Command-line arguments
  // empty for now

#ifdef VERB
  std::cout << "Instantiating CombinedSparseGridDriver:\n";
#endif
  unsigned short level = STARTLEV;  // reference grid level
  RealVector dimension_pref;        // empty -> isotropic
  short growth_rate = UNRESTRICTED_GROWTH;
  short refine_cntl = DIMENSION_ADAPTIVE_CONTROL_GENERALIZED;

  // Store grid and model evaluations 
  RealMatrix storedSets;
  RealVector storedVals;   

  // Restart is data available 
  // restartGSGdriver((char *)GRIDFILE, (char *)FCNFILE, storedSets, storedVals);

  // Start
  // Can either use IntegrationDriver(driver_type) and then assign data or
  // use IntegrationDriver() and then assign letter.
  IntegrationDriver int_driver; // empty envelope
  // assign letter using assign_rep()
  CombinedSparseGridDriver* csg_driver
    = new CombinedSparseGridDriver(level, dimension_pref, growth_rate,
				   refine_cntl);
  int_driver.assign_rep(csg_driver, false); // don't increment ref count

#ifdef VERB
  std::cout << "Instantiating basis...\n";
#endif

  size_t num_vars = NUMVARS;
  std::vector<BasisPolynomial> poly_basis(num_vars); // array of envelopes
  for (int i=0; i<num_vars; ++i)
    poly_basis[i] = BasisPolynomial(LEGENDRE_ORTHOG); // envelope operator=
  csg_driver->initialize_grid(poly_basis);

#ifdef VERB
  std::cout << "  - done\n";
#endif

  // Instantiate Pecos Objects
#ifdef VERB
  std::cout << "Instantiating pecos objects...\n";
#endif

  ExpansionConfigOptions expcfgopt(COMBINED_SPARSE_GRID, // expsolnapp
                                   DEFAULT_BASIS,        // expbassus
                                   SILENT_OUTPUT,        // output level
                                   true,                 // vbd flag
                                   2,                    // vbd order
                                   refine_cntl,          // refinement control
                                   100,                  // max iter
                                   1.e-5,                // conv tol
                                   2);                   // soft conv limit
  BasisConfigOptions bcopt;
  UShortArray aord(MAXORD,NUMVARS);
  SharedBasisApproxData shared_data;                          // Envelope
  SharedProjectOrthogPolyApproxData* shared_poly_data = new   // Letter
    SharedProjectOrthogPolyApproxData(BTYPE,aord,NUMVARS,expcfgopt,bcopt);
  shared_data.assign_rep(shared_poly_data, false); // don't increment ref count
  shared_poly_data->integration_driver_rep(csg_driver);
  shared_poly_data->polynomial_basis(poly_basis);

  // Instantiate Project poly approx
  std::vector<BasisApproximation> poly_approx(nQoI); // array of envelopes
  for ( int iQoI=0; iQoI<nQoI; iQoI++)
    poly_approx[iQoI].assign_rep(new 
      ProjectOrthogPolyApproximation(shared_data), false); // assign letter

#ifdef VERB
  std::cout << "  - done\n";
#endif
 
  // initial grid and compute reference approximation
  RealMatrix var_sets;
  csg_driver->compute_grid(var_sets);
  int numPts = var_sets.numCols();
  assert(num_vars==var_sets.numRows());
  
#ifdef VERB
  std::cout << "Evaluate function on reference grid, instantiate SurrogateData and compute coefficients ...\n"; 
#endif
  // Create SurrogateData instances and assign to ProjectOrthogPolyApproximation instances
  std::vector<bool> computedGridIDs(numPts,true) ; 
  RealMatrix fev0 = feval(var_sets,nQoI,computedGridIDs,NULL);
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    SurrogateDataVars sdv(num_vars,0,0);
    SurrogateDataResp sdr(1,num_vars); // no gradient or hessian
    SurrogateData     sdi;
    for( int jCol = 0; jCol < numPts; jCol++) {
      sdv.continuous_variables(Teuchos::getCol<int,double>(Teuchos::Copy,var_sets,jCol));
      sdr.response_function(fev0(jCol,iQoI));
      sdi.push_back(sdv,sdr);
    }
    poly_approx[iQoI].surrogate_data(sdi);
  } 

  shared_poly_data->allocate_data();    
  std::cout << "  - done\n";
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    std::cout<<"QoI="<<iQoI<<std::endl;
    poly_approx[iQoI].compute_coefficients();
  }
  
#ifdef VERB
  std::cout << "  - done\n";
#endif

  // ----------------- Comment from now restarts/etc------------------------
  // // if restart, check if grid is in the restart
  // std::vector<bool> computedGridIDs ;
  // if (storedVals.length() > 0) {
  //   int numHits = checkSetsInStoredSet(storedSets,var_sets,computedGridIDs);
  //   if ( numHits != var_sets.numCols() )
  //   {
  //     std::cout<<"main() error: initial sets not found in restart"
  //              <<std::endl;
  //     std::terminate();
  //   }
  // } 
  // else
  // {
  //   computedGridIDs.resize(var_sets.numCols(),true);
  //   RealVector fev = feval(var_sets,computedGridIDs,NULL);
  //   addNewSets(storedSets,storedVals,var_sets,fev,computedGridIDs);
  // }

#ifdef DEBUG
  write_data(std::cout, var_sets, true, true, true);
  write_US2A(std::cout, csg_driver->smolyak_multi_index());
#endif

  // start refinement
  csg_driver->initialize_sets();
  UShortArraySet a;
  for ( int iter = 0; iter<NITER; iter++) {

    /* Compute base variance */
    RealVector respVariance(nQoI,0.0);  
    for ( int iQoI=0; iQoI<nQoI; iQoI++) {
      PolynomialApproximation* poly_approx_rep =
	(PolynomialApproximation *) poly_approx[iQoI].approx_rep();
      respVariance[iQoI] = poly_approx_rep->variance() ;
    }
    Real deltaVar = 0.0;
    std::vector<short unsigned int> asave;

    a = csg_driver->active_multi_index();
    std::cout<<"Refine, iteration: "<<iter+1<<'\n';

    int choose = 0;
    for (UShortArraySet::iterator it=a.begin(); it!=a.end(); ++it) {

      // int pick = std::rand();
      // if ( pick > choose) {
      //   asave  = *it;
      //   choose = pick;
      // }

      csg_driver->push_trial_set(*it);

      // Update surrogate data
      numPts = 0;
      if (shared_poly_data->restore_available()) {

        // Set available -> restore in csg and the rest
	csg_driver->restore_set();

        size_t idxRestore = shared_poly_data->restoration_index();
	shared_poly_data->pre_restore_data();
	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  poly_approx[iQoI].restore_coefficients();
          // Also restore the corresponding surrogate data
	  SurrogateData sdi = poly_approx[iQoI].surrogate_data();
	  numPts = sdi.restore(idxRestore,true);
	}
	shared_poly_data->post_restore_data();

      }
      else {

	// New set -> compute
	// Create SurrogateData instances and assign to ProjectOrthogPolyApproximation instances
	csg_driver->compute_trial_grid(var_sets);
        numPts = var_sets.numCols();
	computedGridIDs.resize(numPts,true) ; 
        RealMatrix fev = feval(var_sets,nQoI,computedGridIDs,NULL);

	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  SurrogateDataVars sdv(num_vars,0,0);
	  SurrogateDataResp sdr(1,num_vars); // no gradient or hessian
	  SurrogateData     sdi = poly_approx[iQoI].surrogate_data();
	  for( int jCol = 0; jCol < numPts; jCol++) {
	    sdv.continuous_variables(Teuchos::getCol<int,double>(Teuchos::Copy,var_sets,jCol));
	    sdr.response_function(fev(jCol,iQoI));
	    sdi.push_back(sdv,sdr);
	  } // done loop over number of points
	  //poly_approx[iQoI].surrogate_data(sdi);
	} // done loop over QoIs
        
	shared_poly_data->increment_data();
	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  poly_approx[iQoI].increment_coefficients();
	}
      }


      /* Compute (normalized) change in variance */
      RealVector respVarianceNew(nQoI,0.0);  
      for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	PolynomialApproximation* poly_approx_rep =
	  (PolynomialApproximation *) poly_approx[iQoI].approx_rep();
	respVarianceNew[iQoI] = poly_approx_rep->variance() ;
      }
      respVarianceNew -= respVariance;
      Real normChange = respVarianceNew.normFrobenius()/csg_driver->unique_trial_points();
      if (normChange > deltaVar) {
        deltaVar = normChange;
        asave = *it;
      }

      // ----------------- Comment from now restarts/etc------------------------
      // int numHits = checkSetsInStoredSet(storedSets,var_sets,computedGridIDs);
      // if (numHits<var_sets.numCols())
      // {
      //   RealVector fev = feval(var_sets,computedGridIDs,NULL);
      //   addNewSets(storedSets,storedVals,var_sets,fev,computedGridIDs);
      // }

      csg_driver->pop_trial_set();

      shared_poly_data->decrement_data();
      for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	poly_approx[iQoI].decrement_coefficients();
	// Also restore the corresponding surrogate data
	SurrogateData sdi = poly_approx[iQoI].surrogate_data();
	sdi.pop(numPts,true);
      }

    }

    std::cout<<"Choosing "<<asave<<std::endl ;

    csg_driver->update_sets(asave);
    csg_driver->update_reference();

     //need to restore the data
    size_t idxRestore = shared_poly_data->restoration_index();
    shared_poly_data->pre_restore_data();
    for ( int iQoI=0; iQoI<nQoI; iQoI++) {
      poly_approx[iQoI].restore_coefficients();
      SurrogateData sdi = poly_approx[iQoI].surrogate_data();
      int numPts = sdi.restore(idxRestore,true);
    }
    shared_poly_data->post_restore_data();

  }

  csg_driver->finalize_sets(true, false); // use embedded output option

  // sequence from ApproximationInterface::finalize_approximation():

  // shared pre-finalize:
  shared_poly_data->pre_finalize_data();

  // per-approximation finalize:
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    // from Approximation::finalize() called from PecosApproximation::finalize()
    SurrogateData sdi = poly_approx[iQoI].surrogate_data();
    size_t i, num_restore = sdi.saved_trials(); // # of saved trial sets
    for (i=0; i<num_restore; ++i)
      sdi.restore(shared_poly_data->finalization_index(i),false);
    sdi.clear_saved();
    // from PecosApproximation::finalize()
    poly_approx[iQoI].finalize_coefficients();
  }

  // shared post-finalize:
  shared_poly_data->post_finalize_data();

  for ( int iQoI=0; iQoI<nQoI; iQoI++)
    poly_approx[iQoI].print_coefficients(std::cout,false);


  // Print final sets
  //std::cout<<"Final set:\n";
  //write_US2A(std::cout, csg_driver->smolyak_multi_index());
  
  return (0);

}

void write_US2A(std::ostream& s, const Pecos::UShort2DArray &a)
{
  s << "-----------------------------------------\n";
  size_t i, j, num_a = a.size();
  for (i=0; i<num_a; ++i) {
    const Pecos::UShortArray& aa = a[i];
    for (j=0; j < aa.size(); ++j)
      s<<std::setw(5)<<aa[j];
    s<<"\n";
  }
  s << "-----------------------------------------\n";
  return ;
}

void write_USAS(std::ostream& s, const Pecos::UShortArraySet &a)
{
  s << "-----------------------------------------\n";
  Pecos::UShortArraySet::const_iterator cit;
  for (cit=a.begin(); cit!=a.end(); ++cit) {
    const Pecos::UShortArray& aa = *cit;
    for (size_t j=0; j < aa.size(); ++j)
      s<<std::setw(5)<<aa[j];
    s<<'\n';
  }
  s << "-----------------------------------------\n";
  return ;
}

RealMatrix feval(const RealMatrix &dataMat, const int nQoI, std::vector<bool> &computedGridIDs, void *funInfo) 
{

  assert(nQoI==1);

  int i, j ;

  int numDim = dataMat.numRows(); // Dimensionality
  int numPts = dataMat.numCols(); // Number of function evaluations

  /* Count the number of function evaluations; */
  RealMatrix fev;

  int nEval=0;
  for (i=0; i<numPts; ++i) 
    if (computedGridIDs[i]) nEval++;
  if (nEval==0) return fev;

  fev.shape(nEval,nQoI);
  int ieval=0;
  for (i=0; i<numPts; ++i) {
    if (!computedGridIDs[i]) continue;
    RealVector xIn(numDim);
    for (j=0; j<numDim; ++j)
      xIn[j] = dataMat(j,i);
    fev(ieval,0)=genz(String("cp1"), xIn);
    ieval++;   
  }
  
  return fev;

}

// Restart if data available 
void restartGSGdriver(const char *grid, const char *fcnvals, 
                      Pecos::RealMatrix &storedSets, Pecos::RealVector &storedVals) 
{

  std::ifstream fin ;
  fin.open(grid);
  if (!fin.good()) return; // No restart file available

  std::vector<String> ftext;
  while (!fin.eof())
  {
    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
    String bufstr(buf) ;
    ftext.push_back(bufstr);
  }
  fin.close();

  // get the number of sets and dimensionality 
  int numSets = ftext.size();
  int numVars = 0 ;
  double tmp; 
  std::stringstream parseDbl(ftext[0]);
  while( parseDbl >> tmp ) numVars++ ;

  // get sets from stored text
  storedSets.reshape(numSets,numVars);
  for (int i=0; i<numSets; i++) {
    parseDbl.str(ftext[i]);
    int j = 0;
    while( parseDbl >> storedSets(i,j) ) j++ ;
    assert (j==numVars);
  }
  
  fin.open(fcnvals);
  if (!fin.good()) {
    std::cout<<"restartGSGdriver(): Error: found sets but no function values !" <<std::endl;
    std::terminate();
  }

  int i=0;
  for(std::string line; std::getline(fin, line); )   //read stream line by line
  {
    std::istringstream in(line);      //make a stream for the line itself
    in >> storedVals[i];
    i++;
  }
  assert(i==numSets);

  return ;

}

int checkSetsInStoredSet(const Pecos::RealMatrix &storedSets, const Pecos::RealMatrix &newSets, 
			 std::vector<bool> &computedGridIDs)
{
  // warning: in storedSets each line is a grid point, while in variable_sets each column is a point
  // need to reconcile this at some point
  assert (storedSets.numCols() == newSets.numRows());

  int nHits;
  computedGridIDs.resize(newSets.numCols(),false);
  for (int i=0; i<newSets.numRows(); i++ )
  {
    int numMatch;
    for (int j=0;j<storedSets.numCols(); j++ )
    {
      numMatch = 0;
      for (int k=0;k<storedSets.numRows(); k++ )
      {
	if (storedSets(k,j) != newSets(i,k) )
          break;
        else
	  numMatch++ ;
      }
      if (numMatch == storedSets.numRows()) {
        computedGridIDs[j] = true;
        nHits++;
        break;
      } 
    }
  }
  return (nHits);
}

void addNewSets(Pecos::RealMatrix &storedSets,RealVector &storedVals,
		Pecos::RealMatrix &newSets, RealVector &fev,
		std::vector<bool> &computedGridIDs)
{
    using namespace Pecos;

    
    RealMatrix storedSetsBkp = storedSets;
   
    // resize
    storedSets.shape(storedSetsBkp.numRows()+fev.length(),storedSetsBkp.numCols());

    // copy old data
    for (int i=0; i<storedSetsBkp.numRows(); i++ )
      for (int j=0;j<storedSetsBkp.numCols(); j++ )
	storedSets(i,j) = storedSetsBkp(i,j);

    // add new data
    int irow = storedSetsBkp.numRows();
    for  (int i=0; i<computedGridIDs.size(); i++)
    {
      if (computedGridIDs[i]) {
        for (int j=0;j<storedSets.numCols(); j++ )
	  storedSets(irow,j) = newSets(j,i);
	irow++;
      }
    }

    storedVals.resize(storedVals.length()+fev.length());
    for  (int i=0; i<fev.length(); i++) 
      storedVals[storedSetsBkp.numRows()+i]=fev[i];

    return ;

}
