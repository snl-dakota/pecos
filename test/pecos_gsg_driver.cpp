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

//#define CHGPROTFUNCS

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

  // Command-line arguments
  
  std::cout << "Instantiating CombinedSparseGridDriver:\n";
  unsigned short level = STARTLEV;  // reference grid level
  RealVector dimension_pref;        // empty -> isotropic
  short growth_rate = UNRESTRICTED_GROWTH;
  short refine_cntl = DIMENSION_ADAPTIVE_CONTROL_GENERALIZED;

  // Store grid and model evaluations 
  RealMatrix storedSets;
  RealVector storedVals;   

  // Restart is data available 
  restartGSGdriver((char *)GRIDFILE, (char *)FCNFILE, storedSets, storedVals);

  // Start 
  CombinedSparseGridDriver
    csg_driver(level, dimension_pref, growth_rate, refine_cntl);

  std::cout << "Instantiating basis:\n";
  size_t num_vars = NUMVARS;
  std::vector<BasisPolynomial> poly_basis(num_vars);
  for (int i=0; i<num_vars; ++i)
    poly_basis[i] = BasisPolynomial(LEGENDRE_ORTHOG);
  csg_driver.initialize_grid(poly_basis);

  // Instantiate Pecos Objects
  UShortArray aord(MAXORD,NUMVARS);
  SharedProjectOrthogPolyApproxData srdPolyApprox(BTYPE,aord,NUMVARS);
  srdPolyApprox.integration_driver_rep(&csg_driver);

  // Instantiate Project poly approx
  vector<ProjectOrthogPolyApproximation> polyProjApproxVec;
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    ProjectOrthogPolyApproximation polyProjApprox(srdPolyApprox);
    polyProjApproxVec.push_back(polyProjApprox);
  }

  //initial grid and compute reference approximation
  RealMatrix variable_sets;
  csg_driver.compute_grid(variable_sets);

  // Create SurrogateData instances and assign to ProjectOrthogPolyApproximation instances
  std::vector<bool> computedGridIDs(variable_sets.numCols(),true) ; 
  RealMatrix fev = feval(variable_sets,nQoI,computedGridIDs,NULL);
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    SurrogateDataVars sdv(variable_sets.numRows(),0,0);
    SurrogateDataResp sdr(1,variable_sets.numRows()); // no gradient or hessian
    SurrogateData     sdi;
    for( int jCol=0; jCol<variable_sets.numCols(); jCol++) {
      sdv.continuous_variables(Teuchos::getCol<int,double>(Teuchos::View,variable_sets,jCol));
      sdr.response_function(fev(jCol,iQoI));
      sdi.push_back(sdv,sdr);
    }
    polyProjApproxVec[iQoI].surrogate_data(sdi);
  } 

  //return(0);
#ifdef CHGPROTFUNCS
  srdPolyApprox.allocate_data(); // Cosmin both here and in the loop
                                 // below functions are protected   
  for ( int iQoI=0; iQoI<nQoI; iQoI++) 
     polyProjApproxVec[iQoI].compute_coefficients();
#endif
  
  // ----------------- Comment from now restarts/etc------------------------
  // // if restart, check if grid is in the restart
  // std::vector<bool> computedGridIDs ;
  // if (storedVals.length() > 0) {
  //   int numHits = checkSetsInStoredSet(storedSets,variable_sets,computedGridIDs);
  //   if ( numHits != variable_sets.numCols() )
  //   {
  //     std::cout<<"main() error: initial sets not found in restart"
  //              <<std::endl;
  //     std::terminate();
  //   }
  // } 
  // else
  // {
  //   computedGridIDs.resize(variable_sets.numCols(),true);
  //   RealVector fev = feval(variable_sets,computedGridIDs,NULL);
  //   addNewSets(storedSets,storedVals,variable_sets,fev,computedGridIDs);
  // }

#define DEBUG
#ifdef DEBUG
  write_data(std::cout, variable_sets, true, true, true);
  write_US2A(std::cout, csg_driver.smolyak_multi_index());
#endif

  // start refinement
  csg_driver.initialize_sets();

  UShortArraySet a;
  for ( int iter = 0; iter<NITER; iter++) {

    a = csg_driver.active_multi_index();
    std::cout<<"Refine, iteration: "<<iter+1<<'\n';
    write_USAS(std::cout, a) ;

    std::vector<short unsigned int> asave;
    RealMatrix vsets1;
    int choose = 0;
    for (UShortArraySet::iterator it=a.begin(); it!=a.end(); ++it) {

      int pick = std::rand();
      if ( pick > choose) {
        asave  = *it;
        choose = pick;
      }
      csg_driver.push_trial_set(*it);

      // Surrogate data needs to be updated 

#ifdef CHGPROTFUNCS //need to bring surrogate data up2date: restoration index
      int numPts = 0;
      if (srdPolyApprox.restore_available()) {
        // Set available -> restore
        size_t idxRestore = srdPolyApprox.restoration_index();
	srdPolyApprox.pre_restore_data();
	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  polyProjApproxVec[iQoI].restore_coefficients();
          // Also restore the corresponding surrogate data
	  SurrogateData sdi = polyProjApproxVec[iQoI].surrogate_data();
	  numPts = sdi.restore(idxRestore,true);
	}
	srdPolyApprox.post_restore_data();

      }
      else {
	// New set -> compute
	// Create SurrogateData instances and assign to ProjectOrthogPolyApproximation instances
	csg_driver.compute_trial_grid(vsets1);
        numPts = vsets1.numCols();
	computedGridIDs.resize(vsets1.numCols(),false) ; 
        fev = feval(vsets1,nQoI,computedGridIDs,NULL);
	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  SurrogateDataVars sdv(variable_sets.numRows(),0,0);
	  SurrogateDataResp sdr(1,variable_sets.numRows()); // no gradient or hessian
	  SurrogateData     sdi = polyProjApproxVec[iQoI].surrogate_data();
	  for( int jCol=0; jCol<numPts; jCol++) {
	    sdv.continuous_variables(Teuchos::getCol<int,double>(Teuchos::Copy,variable_sets,jCol));
	    sdr.response_function(fev(jCol,iQoI));
	    std::cout<<fev(jCol,iQoI)<<std::endl;
	    sdi.push_back(sdv,sdr);
	  }
	  //polyProjApproxVec[iQoI].surrogate_data(sdi);
	}
        
	srdPolyApprox.increment_data();
	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  polyProjApproxVec[iQoI].increment_coefficients();
	}
      }

#endif

      // ----------------- Comment from now restarts/etc------------------------
      // int numHits = checkSetsInStoredSet(storedSets,vsets1,computedGridIDs);
      // if (numHits<vsets1.numCols())
      // {
      //   RealVector fev = feval(vsets1,computedGridIDs,NULL);
      //   addNewSets(storedSets,storedVals,vsets1,fev,computedGridIDs);
      // }

      //RealVector fev = feval(vsets1,NULL);
      //write_data(std::cout, vsets1, false, true, true);
      //write_data(std::cout, fev, false, true, true);

      csg_driver.pop_trial_set();
#ifdef NOTPFUNC
      srdPolyApprox.decrement_data();
      for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	polyProjApproxVec[iQoI].decrement_coefficients();
	// Also restore the corresponding surrogate data
	SurrogateData sdi = polyProjApproxVec[iQoI].surrogate_data();
	sdi.pop(numPts,True);
      }
#endif
    
    }

    std::cout<<asave<<std::endl ;

    csg_driver.update_sets(asave);
    csg_driver.update_reference();

#ifdef NOTPFUNC //need to restore the data
    size_t idxRestore = srdPolyApprox.restoration_index();
    srdPolyApprox.pre_restore_data();
    for ( int iQoI=0; iQoI<nQoI; iQoI++) {
      polyProjApproxVec[iQoI].restore_coefficients();
      SurrogateData sdi = polyProjApproxVec[iQoI].surrogate_data();
      numPts = sdi.restore(idxRestore,True);
    }
    srdPolyApprox.post_restore_data();
#endif

  }

  csg_driver.finalize_sets(true, false); // use embedded output option

#ifdef NOTPFUNC //DakotaApprox look at finalize
  srdPolyApprox.pre_finalize_data();
  for ( int iQoI=0; iQoI<nQoI; iQoI++)
    polyProjApproxVec[iQoI].finalize_coefficients();
  srdPolyApprox.post_finalize_data();

  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    SurrogateData sdi = polyProjApproxVec[iQoI].surrogate_data();
    size_t i, num_restore = sdi.saved_trials(); // # of saved trial sets
    for (i=0; i<num_restore; ++i)
      sdi.restore(polyProjApproxVec[iQoI]->finalization_index(i),false);
    sdi.clear_saved();
  }
  
#endif

  // Print final sets
  //std::cout<<"Final set:\n";
  //write_US2A(std::cout, csg_driver.smolyak_multi_index());
  
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

  int i, j, numPts, numDim ;

  numDim = dataMat.numRows(); // Dimensionality
  numPts = dataMat.numCols(); // Number of function evaluations

  // Count the number of function evaluations;
  RealMatrix fev;

  int nEval=0;
  for (i=0; i<numPts; ++i) 
    if (computedGridIDs[i]) nEval++;
  if (nEval==0)
    return fev;

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
