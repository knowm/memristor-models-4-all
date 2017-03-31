//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Additions/Modifications Copyright 2017 Tim Molter
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2016 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Purpose        : Implementation of the MSS memristor model.
//                  
// Creator        : Tim Molter, Knowm Inc.
//
// Creation Date  : 12/26/2016
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MemristorKnowm_h
#define Xyce_N_DEV_MemristorKnowm_h

#include <N_DEV_fwd.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>
#include <N_UTL_RandomNumbers.h>
#include <N_DEV_MemristorTEAM.h>

namespace Xyce {
namespace Device {
namespace MemristorKnowm {

class Model;
class Instance;

// sensitivity functor
// not yet implemented.
class MemristorKnowmSensitivity : public baseSensitivity
{
  public:
  MemristorKnowmSensitivity() :
    baseSensitivity() {};

  virtual ~MemristorKnowmSensitivity() {};

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &name,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const ;
};

static MemristorKnowmSensitivity memrSens;

struct Traits : public DeviceTraits<Model, Instance, MemristorTEAM::Traits>
{
  static const char *name() {return "MemristorKnowm";}
  static const char *deviceTypeName() {return "YMEMRISTOR level 5";}
  static int numNodes() {return 2;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &p);
  static void loadInstanceParameters(ParametricData<Instance> &p);
};

//-----------------------------------------------------------------------------
// Class         : Xyce::Device::MemristorKnowm::Instance
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// MemristorKnowm device instance class.
//
// An instance is created for each occurance of the device in the netlist.
//
// It contains "unique" device information - ie stuff that will be
// true of only one memristor in the circuit, such as the nodes to
// which it is connected.  A memristor is connected to only two
// circuit nodes.
//
// This class does not directly contain information about its node
// indices. It contains indices into the 5 parts (dFdx, dQdx, dx, F,
// and Q) of the matrix problem A*dx = b, and also the solution
// vector x.  A is the Jacobian matrix that will be formed from dFdx
// and d(dQ/dt)dx, dx is the update to x, and b is right hand side
// function vector that will be formed from F and dQ/dt.  These
// indices are global, and determined by topology during the
// initialization stage of execution.
//
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;              
  friend class Model;
  friend class Traits;
  friend class Master;
  friend class MemristorKnowmSensitivity;

public:
  Instance(
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::MemristorKnowm::Instance::~Instance
  // Purpose       : destructor
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Destroys this instance
  //
  // @author Eric Keiter, SNL
  // @date   3/16/00
  ~Instance() {}

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::MemristorKnowm::Instance::getModel
  // Purpose       : destructor
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Gets the resistor model that owns this instance.
  //
  // @return reference to the owning MemristorKnowm::Model
  //
  // @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  // @date   Mon Aug 12 08:36:37 2013
  Model &getModel()   { return model_;  }

  virtual void registerLIDs(const std::vector<int> & intLIDVecRef, const std::vector<int> & extLIDVecRef) /* override */;
  virtual void registerStateLIDs(const std::vector<int> & staLIDVecRef) /* override */;
  virtual void registerStoreLIDs(const std::vector<int> & stoLIDVecRef) /* override */;
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef) /* override */;
  virtual void registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec) /* override */;

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  virtual bool processParams() /* override */;
  virtual bool updateTemperature(const double & temp_tmp) /* override */;
  virtual bool updateIntermediateVars() /* override */;
  virtual bool updatePrimaryState() /* override */;

//---------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//---------------------------------------------------------------------------
//
// Return Jacobian stamp that informs topology of the layout of the
// resistor jacobian.
//
// The Jacobian stamp describes the shape of the Jacobian to the
// Topology subsystem.  The Topology subsystem, in turn, returns
// the offsets into the matrix and solution vectors where this
// instance data is located.
//
// @return const reference to a std::vector of std::vector of
// integers describing Jacobian stamp shape
//
// @author Robert Hoekstra
// @date 8/20/2001
  virtual const std::vector< std::vector<int> > &jacobianStamp() const  /* override */ {
    return jacStamp;
  }

  virtual bool loadDAEQVector() /* override */;
  virtual bool loadDAEFVector() /* override */;
  virtual bool loadDAEdQdx() /* override */;
  virtual bool loadDAEdFdx() /* override */;


  virtual void setupPointers() /* override */;

private:
  static std::vector< std::vector<int> >  jacStamp; //< All MemristorKnowm instances have a common Jacobian Stamp
  static void initializeJacobianStamp();

  Model &     model_;                 //< Owning model

  // User-specified parameters:
  double      R_init_;
  double      Temp_; // instance temperature (TEMP)

  // Derived parameters:
  double      G;                      //< Conductance(1.0/ohms)
  double      dReffdvpos;             //< derivative of Reff with respect to Vpos
  double      dReffdvneg;             //< derivative of Reff with respect to Vneg
  double      dReffdx;                //< derivative of Reff with respect to x  
  double      dIdx;                   //< derivative of I with respect to x
  double      i0;                     //< Current(ohms)
  double      xVarFContribution;      //< x, internal variable for thickness of conductive layer, F-vector contribution
  double      dxFEqdVpos;             //< derivative of X F equation with respect to Vpos
  double      dxFEqdVneg;             //< derivative of X F equation with respect to Vneg
  double      dxFEqdx;                //< derivative of X F equation with respect to X

  int         li_Pos;                 //< Index for Positive Node
  int         li_Neg;                 //< Index for Negative Node
  int         li_x;                   //< Index for internal x, thickness, variable
  int         li_store_R;             //< Index to store resistance value
  int         li_store_tdt;           //< Index to store for next RTN time time delta t 
  int         li_branch_data;         //< Index for Lead Current and junction voltage (for power calculations)

  // Offset variables corresponding to the above declared indices.
  int         APosEquPosNodeOffset;   //< Column index into matrix of Pos/Pos conductance
  int         APosEquNegNodeOffset;   //< Column index into matrix of Pos/Neg conductance
  int         APosEquXNodeOffset;     //< Column index into matrix for internal varaible x, layer thickness
  int         ANegEquPosNodeOffset;   //< Column index into matrix of Neg/Pos conductance
  int         ANegEquNegNodeOffset;   //< Column index into matrix of Neg/Neg conductance
  int         ANegEquXNodeOffset;     //< Column index into matrix for internal varaible x, layer thickness
  int         XEquVPosOffset;         //< Thickness governing equation, VPos dependence
  int         XEquVNegOffset;         //< Thickness governing equation, VNeg dependence
  int         XEquXOffset;            //< Thickness variable, in thickness governing equation equation

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for Jacobian
  double *    f_PosEquPosNodePtr;
  double *    f_PosEquNegNodePtr;
  double *    f_PosEquXNodePtr;
  double *    f_NegEquPosNodePtr;
  double *    f_NegEquNegNodePtr;
  double *    f_NegEquXNodePtr;
  double *    f_XEquPosNodePtr;
  double *    f_XEquNegNodePtr;
  double *    f_XEquXNodePtr;
  double *    q_XEquXNodePtr;
#
#endif

};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::MemristorKnowm::Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// MemristorKnowm model class
//
class Model : public DeviceModel
{
  friend class ParametricData<Model>;               //< Allow ParametricData to changes member values
  friend class Instance;                            //< Don't force a lot of pointless getters
  friend class Traits;
  friend class Master;                              //< Don't force a lot of pointless getters

public:
  typedef std::vector<Instance *> InstanceVector;

  Model(
     const Configuration &     configuration,
     const ModelBlock &        model_block,
     const FactoryBlock &      factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::MemristorKnowm::Model::addInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014 
  //---------------------------------------------------------------------------
  //
  // Add an instance to the list of instances associated with this model
  //
  // @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  // @date   8/12/2013
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

  virtual bool processParams() /* override */;
  virtual bool processInstanceParams() /* override */;

private:
  InstanceVector      instanceContainer;            //< List of owned resistor instances

  // model parameters for Knowm model

  double      Roff_;
  double      Ron_;
  double      Voff_;
  double      Von_;
  double      Tau_;
  
  // model parameters for generalized Knowm model

  double      Phi_;
  double      SchottkyForwardAlpha_;
  double      SchottkyForwardBeta_;
  double      SchottkyReverseAlpha_;
  double      SchottkyReverseBeta_;
};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::MemristorKnowm::Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014 
//-----------------------------------------------------------------------------
//
// MemristorKnowm master
//
// The "master" class is the one that contains the updateState, loadDAEVectors
// and loadDAEMatrices methods that are actually called when it is time to
// compute and load device contributions.
//
// The default implementations of these methods in the DeviceMaster
// template class simply loops over all instances and calls their
// updatePrimaryState, loadDAEFVector/loadDAEQVector, and
// loadDAEdFdx/loadDAEdQdx methods, respectively.
//
// For efficiency, the MemristorKnowm class reimplements these methods to do the
// work directly, instead of calling instance-level functions.
//
class Master : public DeviceMaster<Traits>
{
  friend class Instance;                           
  friend class Model;                             

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::MemristorKnowm::Master::Master
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Construct a MemristorKnowm Device.
  //
  // @param configuration
  // @param factory_block
  // @param solver_state
  // @param device_options
  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options)
  {}

  virtual bool updateState(double * solVec, double * staVec, double * stoVec);
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx);

};

void registerDevice();

} // namespace MemristorKnowm
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_MemristorKnowm_h
