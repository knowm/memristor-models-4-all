// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//
// Purpose        : Implementation of the Knowm MKnowm memristor model.
//
// Creator        : Tim Molter, Knowm Inc.
//
// Creation Date  : 01/04/17
//
//
//----------------------------------------------------------------------------
#include <Xyce_config.h>
#include <N_DEV_Const.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <Sacado.hpp>
#include "../include/N_DEV_MemristorKnowm.h"

namespace Xyce {
namespace Device {
namespace MemristorKnowm {







template <typename ScalarT>
ScalarT p0ff2on( const ScalarT & V1, const ScalarT & V2, double VON, double TC, double VT )
{
	ScalarT exponent = -1 * ((V1-V2) - VON) / VT;
  double alpha = 1 / TC;
  ScalarT fval = alpha / (1.0 + exp(exponent));
  return fval;
}

template <typename ScalarT>
ScalarT pOn2Off( const ScalarT & V1, const ScalarT & V2, double VOFF, double TC, double VT )
{

	ScalarT exponent = -1 * ((V1-V2) + VOFF) / VT;
  double alpha = 1 / TC;
  ScalarT fval = alpha * (1.0 - 1.0 / (1.0 + exp(exponent)));
  return fval;
}

template <typename ScalarT>
ScalarT dXdt( const ScalarT & V1, const ScalarT & V2, const ScalarT & X, double RON, double ROFF, double VON, double VOFF, double TC, double VT)
{

//	if (DEBUG_DEVICE){
//		Xyce::dout()  << "  V1 = " <<  V1 << std::endl;
//		Xyce::dout()  << "  V2 = " <<  V2 << std::endl;
//		Xyce::dout()  << "  X = " <<  X << std::endl;
//	}

  // Probabilities
	ScalarT p0ff2onVal = p0ff2on(V1, V2, VON, TC, VT);
	ScalarT pOn2OffVal = pOn2Off(V1, V2, VOFF, TC, VT);

  // Number of switches making a transition
	ScalarT n0ff2on = (1 - X) * p0ff2onVal;
	ScalarT nOn2Off = X  * pOn2OffVal;

	ScalarT fval = n0ff2on - nOn2Off;

	return fval;
}






template <typename ScalarT>
ScalarT SchottkyCurrent( const ScalarT & V1, const ScalarT & V2, double SchottkyForwardAlpha, double SchottkyForwardBeta, double SchottkyReverseAlpha, double SchottkyReverseBeta )
{
	return SchottkyReverseAlpha * (-1 * exp(-1 * SchottkyReverseBeta * (V1-V2))) + SchottkyForwardAlpha * (exp(SchottkyForwardBeta * (V1-V2)));
}

template <typename ScalarT>
ScalarT Geff( const ScalarT & X, double RON, double ROFF )
{
	return  X / RON + (1 - X) / ROFF;
}

template <typename ScalarT>
ScalarT I_V( const ScalarT & V1, const ScalarT & V2, const ScalarT & X, double RON, double ROFF, double PHI, double SchottkyForwardAlpha, double SchottkyForwardBeta, double SchottkyReverseAlpha, double SchottkyReverseBeta ){

	ScalarT Gval=	Geff( X, RON, ROFF );
	ScalarT KnowmCurrentval = (V1-V2)*Gval;

	ScalarT SchottkyCurrentval = SchottkyCurrent(V1, V2, SchottkyForwardAlpha, SchottkyForwardBeta, SchottkyReverseAlpha, SchottkyReverseBeta);

	ScalarT fval = PHI * KnowmCurrentval + (1 - PHI) * SchottkyCurrentval;

//	if (DEBUG_DEVICE){
//		Xyce::dout()  << "  Geff = " <<  Gval << std::endl;
//		Xyce::dout()  << "  RON = " <<  RON << std::endl;
//		Xyce::dout()  << "  ROFF = " <<  ROFF << std::endl;
//		Xyce::dout()  << "  fval = " <<  fval << std::endl;
//	}

	return fval;
}






//
// Common Jacobian Stamp for all MemristorKnowm devices.
// Because all memristors have identical Jacobian stamps, this data is
// declared static and is shared by all memristor instances.
// 
std::vector<std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::initializeJacobianStamp
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 2/11/2014
//-----------------------------------------------------------------------------
//
// @brief Common Jacobian stamp initializer for all MemristorKnowm devices.
//
// The Jacobian stamp is a sparse-matrix representation of the pattern
// of non-zero elements that a device will put into the Jacobian matrix.
//
// The Jacobian stamp is used by the Topology package to determine indices
// into the full Jacobian matrix for elements that correspond to this 
// device.
//
// There is one row of the Jacobian stamp for each equation associated with
// a device.  The number of elements in a row are the number of non-zero
// elements in that row of the device's contribution to the Jacobian.
// The values of the elements are numbers local to the device that
// represent the column in which the non-zero element appears.
//
// For this memristor, there are two external nodes (the positive and negative
// terminals of the device).  The positive node is referred to as the 0th
// node of the device, and the negative node the 1st node.
//
void Instance::initializeJacobianStamp()
{
  if (jacStamp.empty())
  {
    jacStamp.resize(3);
    jacStamp[0].resize(3);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[0][2] = 2;
    jacStamp[1].resize(3);
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
    jacStamp[1][2] = 2;
    jacStamp[2].resize(3);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
    jacStamp[2][2] = 2;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Traits::loadInstanceParameters
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Loads the parameter definition into the instance parameter map.
//
// @param p     instance parameter map
//
// Each parameter supported by the memristor device instance is
// defined via the addPar call.  The minimum information required is
// the name of the parameter, the default value, and a member pointer
// to a variable in the instance class to which the value will be
// assigned.
//
// Additional information such as documentation for the parameter, units
// (used only in output of tables for documentation), whether a special
// boolean should be set if the parameter is given, etc. may be specified
// using the various set* calls of the Xyce::Device::Descriptor class.
//
// It is important to note that since parameters defined by addPar are
// defined as metadata separate from any instance, it is not possible to
// establish interrelationships between parameter defaults here.  Parameter
// defaults specified in addPar are always constants.  If a device requires
// that parameter defaults depend on values of other parameters, this has to
// be done in the instance constructor.  Examples of such parameters in this 
// device are the "DTEMP" and "W" parameters, which are set to special defaults
// at device instantiation time.  Defaults for those parameters in the addPar
// calls (and hence in the LaTeX tables in the reference guide produced from
// this data) are misleading.
//
void Traits::loadInstanceParameters(ParametricData<MemristorKnowm::Instance> &p)
{
  p.addPar("Rinit", 0.0, &MemristorKnowm::Instance::rInit_)
    .setUnit(U_OHM)
    .setDescription("Initial value for resistance");
  p.addPar("TEMP", 300.0, &MemristorKnowm::Instance::Temp_)
    .setUnit(U_DEGK)
    .setDescription("Device Temperature");

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Traits::loadModelParameters
// Purpose       : 
// Special Notes : The addPar calls here were refactored and moved here
//                 from the model constructor.  Those addPars had been
//                 in place from 2/4/2005.
// Scope         : private
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Loads the parameter definition into the model parameter map.
//
// @param p     model parameter map
//
// @see Xyce::Device::MemristorKnowm::Traits::loadInstanceParameters
//
//
void Traits::loadModelParameters(ParametricData<MemristorKnowm::Model> &p)
{
  // Create parameter definitions for parameter member variables

  // NOTE: The first string arg need to be ALL CAPS!

  p.addPar("ROFF", 16000.0, &MemristorKnowm::Model::Roff_)
    .setUnit(U_OHM)
    .setDescription("Off resistance.");
  p.addPar("RON", 100.0, &MemristorKnowm::Model::Ron_)
    .setUnit(U_OHM)
    .setDescription("On resistance");
  p.addPar("VOFF", 1.0, &MemristorKnowm::Model::Voff_)
    .setUnit(U_VOLT)
    .setDescription("Threshold voltage to turn device off");
  p.addPar("VON", 1.0, &MemristorKnowm::Model::Von_)
    .setUnit(U_VOLT)
    .setDescription("Threshold voltage to turn device on");
  p.addPar("TAU", 10.0e-9, &MemristorKnowm::Model::Tau_)
    .setUnit(U_SECOND)
    .setDescription("Time constant");

  p.addPar("PHI", 1.0, &MemristorKnowm::Model::Phi_)
    .setUnit(U_NONE)
    .setDescription("Schottky Forward Alpha_");
  p.addPar("SFA", 0.0, &MemristorKnowm::Model::SchottkyForwardAlpha_)
    .setUnit(U_NONE)
    .setDescription("Schottky Forward Alpha_");
  p.addPar("SFB", 0.0, &MemristorKnowm::Model::SchottkyForwardBeta_)
    .setUnit(U_NONE)
    .setDescription("Schottky Forward Beta");
  p.addPar("SRA", 0.0, &MemristorKnowm::Model::SchottkyReverseAlpha_)
    .setUnit(U_NONE)
    .setDescription("Schottky Reverse Alpha");
  p.addPar("SRB", 0.0, &MemristorKnowm::Model::SchottkyReverseBeta_)
    .setUnit(U_NONE)
    .setDescription("Schottky Reverse Beta");


}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::Instance
// Purpose       : Instance constructor
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Construct a memristor instance.
//
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
	  rInit_(0.0),
		Temp_(300.0),
    G(0.0),
    i0(0.0),
    li_Pos(-1),
    li_Neg(-1),
    li_x(-1),
    li_store_tdt(-1), 
    li_branch_data(0),


    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    APosEquXNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ANegEquXNodeOffset(-1),
    XEquVPosOffset(-1),
    XEquVNegOffset(-1),
    XEquXOffset(-1)

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    ,
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_PosEquXNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0),
    f_NegEquXNodePtr(0),
    f_XEquPosNodePtr(0),
    f_XEquNegNodePtr(0),
    f_XEquXNodePtr(0),
    q_XEquXNodePtr(0)
#endif

{
  // Initialize DeviceInstance values
  numIntVars   = 1;    // Initialize number if internal nodes in DeviceInstance
  numExtVars   = 2;    // Initialize number if external nodes in DeviceInstance
  numStateVars = 0;    // Initialize number if state variables in DeviceInstance
  setNumStoreVars(2);  // Initialize number if store variables in DeviceInstance 

  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

  initializeJacobianStamp();

  // Set params to constant default values from parameter definition
  setDefaultParams();

  // Set params according to instance line and constant defaults from metadata
  setParams(instance_block.params);

  // Calculate any parameters specified as expressions
  updateDependentParameters();

  // Process the parameters to complete initialization
  processParams();
  
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::processParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
// Process parameters.
//
// @return true on success
//
// In general, the processParams method is intended as a place for complex
// manipulation of parameters that must happen if temperature or other
// external changes happen.
//
bool Instance::processParams()
{

  // initialize X state (width between 0 and 1 and scaled to D)
  if (!given("Rinit")){
  	rInit_ = model_.Roff_;
  }

//	if (DEBUG_DEVICE){
//		Xyce::dout()  << "----------Instance::processParams"  << std::endl;
//        Xyce::dout()  << " rInit_  = " << rInit_ << std::endl;
//	}

  if (!given("TEMP")){
  	Temp_ = getDeviceOptions().temp.getImmutableValue<double>();
  }

  updateTemperature(Temp_);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::registerLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Register local IDs
//
// Register the local internal and external node IDs.
//
// @param intLIDVecRef internal local IDs from topology package
// @param extLIDVecRef external local IDs from topology package
//
void Instance::registerLIDs(
  const std::vector<int> & intLIDVecRef,
  const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // Copy the local ID lists.
  intLIDVec = intLIDVecRef;                           // Set the internal local IDs in DeviceInstance
  extLIDVec = extLIDVecRef;                           // Set the external local IDs in DeviceInstance

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  // This lid is for the internal variable for layer thickness that determines resistance
  li_x   = intLIDVec[0];
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::registerStateLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Register the local state IDs
//
// @note The memristor does not have any state vars, so this function
// does nothing.
//
// @param staLIDVecRef State variable local IDs
//
// In general, devices may declare at construction that they require storage
// locations in the "state vector."  Topology assigns locations in the 
// state vector and returns "Local IDs" (LIDs) for devices to use for their
// state vector entries.  If a device uses state vector entries, it
// uses the registerStateLIDs method to save the local IDs for later use.
// 
// @author Robert Hoekstra, SNL, Parallel Computational Sciences
// @date   06/12/02
//
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::registerStoreLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
// Register the local store IDs
//
// In addition to state vector, Xyce maintains a separate datastructure
// called a "store" vector.  As with other such vectors, the device
// declares at construction time how many store vector entries it needs,
// and later Topology assigns locations for devices, returning LIDs.
//
//
//
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef)
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());
  li_store_R = stoLIDVecRef[0];
  li_store_tdt = stoLIDVecRef[1];
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
// Register the local store IDs
//
// In addition to state vector, Xyce maintains a separate datastructure
// called a "branch data" vector.  As with other such vectors, the device
// declares at construction time how many branch vector entries it needs,
// and later Topology assigns locations for devices, returning LIDs.
//
// These LIDs are stored in this method for later use.
//
// The memristor device uses exactly one "branch data vector" element, where
// it keeps the "lead current" that may be used on .PRINT lines as
// "I(ymemristor)" for the current through ymemristor. and a junction voltage.
//
// @param stoLIDVecRef Store variable local IDs
//
// @author Richard Schiek, Electrical Systems Modeling
// @date   12/18/2012
void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Sysetms Modeling 
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  addInternalNode(symbol_table, li_x, getName(), "x");

  addStoreNode(symbol_table, li_store_R, getName(), "R");
  addStoreNode(symbol_table, li_store_tdt, getName(), "TDT");

  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::registerJacLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Register the Jacobian local IDs
//
// @param jacLIDVec Jacobian local Ids
//
// @see Xyce::Device::MemristorKnowm::Instance::initializeJacobianStamp
//
// Having established local IDs for the solution variables, Topology must
// also assign local IDs for the elements of the Jacobian matrix.
//
// For each non-zero element that was identified in the jacobianStamp,
// Topology will assign a Jacobian LID.  The jacLIDVec will have the 
// same structure as the JacStamp, but the values in the jacLIDVec will
// be offsets into the row of the sparse Jacobian matrix corresponding
// to the non-zero identified in the stamp.
// 
// These offsets are stored in this method for use later when we load
// the Jacobian.
//
// @author Robert Hoekstra, SNL, Parallel Computational Sciences
// @date   08/27/01
void Instance::registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec)
{
  // Let DeviceInstance do its work.
  DeviceInstance::registerJacLIDs(jacLIDVec);

  // Store offsets of the components of the Jacobian of this instance
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  APosEquXNodeOffset   = jacLIDVec[0][2];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
  ANegEquXNodeOffset   = jacLIDVec[1][2];
  XEquVPosOffset       = jacLIDVec[2][0]; 
  XEquVNegOffset       = jacLIDVec[2][1];
  XEquXOffset          = jacLIDVec[2][2];
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::setupPointers
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Setup direct access pointer to solution matrix and vectors.
//
// @see Xyce::Device::MemristorKnowm::Instance::registerJacLIDs
//
// As an alternative to the row offsets defined in registerJacLIDs, it 
// is also possible to obtain direct pointers of the Jacobian elements.
//
// This method uses the offsets obtained in registerJacLIDs to retrieve
// the pointers.
//
// In this device the pointers to the matrix are only saved
// (and are only used for matrix access) if
// Xyce_NONPOINTER_MATRIX_LOAD is NOT defined at compile time.  For
// some devices the use of pointers instead of array indexing can be
// a performance enhancement.
//
// Use of pointers in this device is disabled by defining
// Xyce_NONPOINTER_MATRIX_LOAD at compile time.  When disabled, array
// indexing with the offsets from registerJacLIDs is used in
// the matrix load methods.
//
// @author Eric Keiter, SNL
// @date   11/30/08
void Instance::setupPointers()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  Linear::Matrix & dQdx = *(extData.dQdxMatrixPtr);
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_PosEquXNodePtr   = &(dFdx[li_Pos][APosEquXNodeOffset]);
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
  f_NegEquXNodePtr   = &(dFdx[li_Neg][ANegEquXNodeOffset]);
  f_XEquPosNodePtr   = &(dFdx[li_x][XEquVPosOffset]);
  f_XEquNegNodePtr   = &(dFdx[li_x][XEquVNegOffset]);
  f_XEquXNodePtr     = &(dFdx[li_x][XEquXOffset]);
  q_XEquXNodePtr     = &(dQdx[li_x][XEquXOffset]);
#endif
}

// The following 6 methods can just simply return true becuase the methods in the Master implementation cover this functionality.

bool Instance::loadDAEQVector()
{
  return true;
}

bool Instance::loadDAEFVector()
{
  return true;
}

bool Instance::loadDAEdQdx()
{
  return true;
}

bool Instance::loadDAEdFdx()
{
  return true;
}

bool Instance::updatePrimaryState()
{
  return true;
}

bool Instance::updateIntermediateVars()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Update the parameters that depend on the temperature of the device
//
// @param temp_tmp temperature
//
// Xyce has a number of mechanisms that allow temperature to be changed
// after a device has been instantiated.  These include .STEP loops over
// temperature.  When temperature is changed, any device that has parameters
// that depend on temperature must be updated.  That updating happens here.
//
// The MemristorKnowm device supports temperature-dependent resistance through its
// TC1 (linear dependence) and TC2 (quadratic dependence) parameters.
// If these parameters are specified, the resistance must be updated.
//
// @return true on success
//
// @author Tom Russo, Component Information and Models
// @date   02/27/01
bool Instance::updateTemperature(const double & temp_tmp)
{
//	TEMP;
	if (DEBUG_DEVICE){
		Xyce::dout()  << "  temp_tmp = " <<  temp_tmp << std::endl;
	}
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Process model parameters
//
// @return true on success
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   6/03/02
bool Model::processParams()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Model::processInstanceParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//----------------------------------------------------------------------------
//
// Process the instance parameters of instance owned by this model
//
// This method simply loops over all instances associated with this
// model and calls their processParams method.
//
// @return true
//
// @author Dave Shirely, PSSI
// @date   03/23/06
bool Model::processInstanceParams()
{
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    (*it)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Construct a memristor model from a "model block" that was created
// by the netlist parser.
//
// @param configuration
// @param model_block
// @param factory_block
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   5/16/00
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block),
	Roff_(0.0),
	Ron_(0.0),
	Von_(0.0),
	Voff_(0.0),
	Tau_(0.0)
{
  // Set params to constant default values.
  setDefaultParams();

  // Set params according to .model line and constant defaults from metadata.
  setModParams(model_block.params);

  // Set any non-constant parameter defaults.
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions.
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors.
  processParams();


}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Destroy this model.
//
// Also destroys all instances that use this model.
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   3/16/00
Model::~Model()
{
  // Destory all owned instances
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    delete (*it);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_MemristorKnowmModel::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Print instances associated with this model.
//
// printOutInstances was intended as a debugging tool, but a refactor some time back took out all the places where this function would be called.
//
// @param os output stream
//
// @return reference to output stream
//
// @author Eric Keiter, SNL, Parallel Computational Sciences
// @date   4/03/00
//
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  os << std::endl;
  os << "Number of MemristorKnowm Instances: " << instanceContainer.size() << std::endl;
  os << "    name     model name  Parameters" << std::endl;

  int i = 0;
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    os << "  " << i << ": " << (*it)->getName() << "\t";
    os << getName();
    os << "\tG(T) = " << (*it)->G;
    os << std::endl;
    ++i;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
// Apply a device instance "op" to all instances associated with this
// model
// 
// @param[in] op Operator to apply to all instances.
// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}

double initialX( double rOn, double rOff, double rInit){

	double x = (rOn * (rInit - rOff )) / (rInit * (rOn - rOff));
	if(x== -0){
		x = 0.0;
	}

  return x;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Master::updateState
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Update state for all memristor instances, regardless of model.
//
// @param solVec solution vector
// @param staVec state vector
// @param stoVec store vector
//
// @return true on success
//
// @see Xyce::Device::MemristorKnowm::Instance::updatePrimaryState
// @author Eric Keiter, SNL
// @date   11/26/08
bool Master::updateState(double * solVec, double * staVec, double * stoVec)
{

//	if (DEBUG_DEVICE){
//		Xyce::dout()  << " updateState "  << std::endl;
//	}

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);
    
    double v_pos    = solVec[ri.li_Pos];
    double v_neg    = solVec[ri.li_Neg];
    double x        = solVec[ri.li_x];
//		if (DEBUG_DEVICE){
//			Xyce::dout()  << "  x_before = " <<  x << std::endl;
//		}


		{
				Sacado::Fad::SFad<double,3> varV1( 3, 0, v_pos );
				Sacado::Fad::SFad<double,3> varV2( 3, 1, v_neg );
				Sacado::Fad::SFad<double,3> varX( 3, 2, x );
				Sacado::Fad::SFad<double,3> resultFad;
				resultFad = I_V( varV1, varV2, varX, ri.model_.Ron_, ri.model_.Roff_, ri.model_.Phi_, ri.model_.SchottkyForwardAlpha_, ri.model_.SchottkyForwardBeta_, ri.model_.SchottkyReverseAlpha_, ri.model_.SchottkyReverseBeta_ );

				ri.i0 = resultFad.val(); // current
				ri.G  = resultFad.dx(0); // di/dv = conductance
				ri.dIdx = resultFad.dx(2); // di/dx

//				if (DEBUG_DEVICE){
//					Xyce::dout()  << "  ri.i0 = " <<  ri.i0 << std::endl;
//					Xyce::dout()  << "  ri.G = " <<  ri.G << std::endl;
//					Xyce::dout()  << "  ri.dIdx = " <<  ri.dIdx << std::endl;
//				}
			}

			{
				// evaluate the state variable equation
				Sacado::Fad::SFad<double,3> varV1( 3, 0, v_pos );
				Sacado::Fad::SFad<double,3> varV2( 3, 1, v_neg );
				Sacado::Fad::SFad<double,3> varX( 3, 2, x );
				Sacado::Fad::SFad<double,3> resultFad;

				double vT = 0.026;
				vT = ri.Temp_ * CONSTKoverQ;
				if (DEBUG_DEVICE){
					Xyce::dout()  << "  vT = " <<  vT << std::endl;
				}

				resultFad = dXdt( varV1, varV2, varX, ri.model_.Ron_, ri.model_.Roff_, ri.model_.Von_, ri.model_.Voff_, ri.model_.Tau_, vT );

				ri.xVarFContribution = resultFad.val();
				if( getSolverState().dcopFlag )
				{
					ri.xVarFContribution = 0;
				}
				ri.dxFEqdVpos = resultFad.dx(0);
				ri.dxFEqdVneg = resultFad.dx(1);
				ri.dxFEqdx = resultFad.dx(2);

//				if (DEBUG_DEVICE){
//					Xyce::dout()  << "  ri.xVarFContribution = " <<  ri.xVarFContribution << std::endl;
//					Xyce::dout()  << "  ri.dxFEqdVpos = " <<  ri.dxFEqdVpos << std::endl;
//					Xyce::dout()  << "  ri.dxFEqdVneg = " <<  ri.dxFEqdVneg << std::endl;
//					Xyce::dout()  << "  ri.dxFEqdx = " <<  ri.dxFEqdx << std::endl;
//				}

			}

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Master::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load DAE vectors of all memristor instances, regardless of model
//
// @param solVec solution vector
// @param fVec f vector
// @param qVec q vector
// @param leadF store lead current f vector
// @param leadQ store lead current q vector
//
// @return true on success
//
// @note Because the memristor device re-implements the base-class
// Master::loadDAEVectors, the Instance::loadDAEFVector method is
// never called.  This method replaces those, and does the same work
// but inside a loop over all memristor instances.
//
// @see Xyce::Device::MemristorYakopcic::Instance::loadDAEFVector
//
// @author Eric Keiter, SNL
// @date   11/26/08

bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
//  	if (DEBUG_DEVICE){
//  		Xyce::dout()  << "-loadDAEVectors"  << std::endl;
//  	}

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);
    fVec[ri.li_Pos] += -ri.i0;
    fVec[ri.li_Neg] += ri.i0;
    fVec[ri.li_x]   += ri.xVarFContribution;
    qVec[ri.li_x]   -= solVec[ri.li_x];
    if( getSolverState().dcopFlag )
		{
    	double x = initialX(ri.model_.Ron_, ri.model_.Roff_, ri.rInit_);
      qVec[ri.li_x] -= x;
//			if (DEBUG_DEVICE){
//				Xyce::dout()  << "DCOP x = "  << x  << std::endl;
//			}
    }
    if( ri.G != 0 )
    {
      double * storeVec = ri.extData.nextStoVectorRawPtr;
      storeVec[ri.li_store_R] = 1.0/ri.G;
    }

    if( ri.loadLeadCurrent )
    {
      leadF[ri.li_branch_data] = ri.i0;
      junctionV[ri.li_branch_data] = solVec[ri.li_Pos] - solVec[ri.li_Neg];
    }

  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorYakopcic::Master::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Load DAE matrices for all memristor instances, regardless of model
//
// @param dFdx matrix of derivatives of F vector with respect to solution
// @param dQdx matrix of derivatives of Q vector with respect to solution
//
// @return true on success
//
// @note Because the memristor device re-implements the base-class
// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
// never called.  This method replaces those, and does the same work
// but inside a loop over all memristor instances.
//
// @see Xyce::Device::MemristorYakopcic::Instance::loadDAEdFdx
//
// @author Eric Keiter, SNL
// @date   11/26/08

bool Master::loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{

//  	if (DEBUG_DEVICE){
//  		Xyce::dout()  << "-loadDAEMatrices"  << std::endl;
//  	}

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *(ri.f_PosEquPosNodePtr) += ri.G;
    *(ri.f_PosEquNegNodePtr) -= ri.G;
    *(ri.f_NegEquPosNodePtr) -= ri.G;
    *(ri.f_NegEquNegNodePtr) += ri.G;

    *(ri.f_PosEquXNodePtr)   += ri.dIdx;
    *(ri.f_NegEquXNodePtr )  += ri.dIdx;

    *(ri.f_XEquPosNodePtr )  += ri.dxFEqdVpos;
    *(ri.f_XEquNegNodePtr )  += ri.dxFEqdVneg;
    *(ri.f_XEquXNodePtr )    += ri.dxFEqdx;

    *(ri.q_XEquXNodePtr )    = -1.0;

#else
    dFdx[ri.li_Pos][ri.APosEquPosNodeOffset] += ri.G;
    dFdx[ri.li_Pos][ri.APosEquNegNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquPosNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquNegNodeOffset] += ri.G;
    dFdx[ri.li_Pos][ri.APosEquXNodeOffset]   += ri.dIdx;
    dFdx[ri.li_Neg][ri.ANegEquXNodeOffset]   += ri.dIdx;
    dFdx[ri.li_x][ri.XEquVPosOffset]         += ri.dxFEqdVpos;
    dFdx[ri.li_x][ri.XEquVNegOffset]         += ri.dxFEqdVneg;
    dFdx[ri.li_x][ri.XEquXOffset]            += ri.dxFEqdx;
    dQdx[ri.li_x][ri.XEquXOffset] = -1.0;
#endif
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::Traits::factory
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Create a new instance of the MemristorKnowm device.
//
// @param configuration
// @param factory_block
//
Device * Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::MemristorKnowm::registerDevice
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Define how to use the device in a netlist.
//
// This method is called from the Xyce::Device::registerOpenDevices
// function, which in turn is called by the device manager.
//
// The device is declared here to be an "memristor" device, which must 
// take a model card of type "memristor".  This device will correspond to model
// level 1 of memristor models.
void
registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("Knowm", 1)
    .registerModelType("Knowm", 1);
}

//-----------------------------------------------------------------------------
// Function      : MemristorKnowmSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p=R.  
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
void MemristorKnowmSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  double * solVec = in->extData.nextSolVectorRawPtr;
  double v_pos = solVec[in->li_Pos];
  double v_neg = solVec[in->li_Neg];

  double dfdpLoc = -(v_pos-v_neg)*in->G*in->G;

  dfdp.resize(2);
  dfdp[0] = +dfdpLoc;
  dfdp[1] = -dfdpLoc;

  Findices.resize(2);
  Findices[0] = in->li_Pos;
  Findices[1] = in->li_Neg;
}

} // namespace MemristorKnowm
} // namespace Device
} // namespace Xyce
