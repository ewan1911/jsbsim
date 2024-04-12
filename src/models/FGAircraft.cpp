/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Module:       FGAircraft.cpp
 Author:       Jon S. Berndt
 Date started: 12/12/98
 Purpose:      Encapsulates an aircraft
 Called by:    FGFDMExec

 ------------- Copyright (C) 1999  Jon S. Berndt (jon@jsbsim.org) -------------

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; either version 2 of the License, or (at your option) any
 later version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License along
 with this program; if not, write to the Free Software Foundation, Inc., 59
 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

 Further information about the GNU Lesser General Public License can also be
 found on the world wide web at http://www.gnu.org.

FUNCTIONAL DESCRIPTION
--------------------------------------------------------------------------------
Models the aircraft reactions and forces. This class is instantiated by the
FGFDMExec class and scheduled as an FDM entry.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMMENTS, REFERENCES,  and NOTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INCLUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "FGAircraft.h"
#include "FGFCS.h"
#include "input_output/FGXMLElement.h"

#include "FGFDMExec.h"

//INCLUDE DE AUXILIARY
#include "models/FGAuxiliary.h"

#define PI 3.14159

using namespace std;

namespace JSBSim {

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

FGAircraft::FGAircraft(FGFDMExec* fdmex) : FGModel(fdmex)
{
  /* std::shared_ptr<JSBSim::FGAuxiliary>auxiliaryPtr = fdmex->GetAuxiliary(); // VOIR CHATGPT JSBSIM CODE EXPL
  Auxiliary = auxiliaryPtr.get(); */

  Name = "FGAircraft";
  WingSpan = 0.0;
  WingArea = 0.0;
  cbar = 0.0;
  HTailArea = VTailArea = 0.0;
  HTailArm  = VTailArm  = 0.0;
  lbarh = lbarv = 0.0;
  vbarh = vbarv = 0.0;
  WingIncidence = 0.0;

  bind();

  Debug(0);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FGAircraft::~FGAircraft()
{
  Debug(1);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool FGAircraft::InitModel(void)
{
  if (!FGModel::InitModel()) return false;

  vForces.InitMatrix();
  vMoments.InitMatrix();

  /* int form = 1;
  double pos = 10.0;
  FDMExec->GetFCS()->SetDePos(form, pos); */

  return true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool FGAircraft::Run(bool Holding)
{
  if (FGModel::Run(Holding)) return true;
  if (Holding) return false;

  RunPreFunctions();

  vForces = in.AeroForce;
  vForces += in.PropForce;
  vForces += in.GroundForce;
  vForces += in.ExternalForce;
  vForces += in.BuoyantForce;

  /* FGColumnVector3 addedLift = FDMExec->GetAuxiliary()->resultLift(); //fontionne avec l'erreur
  vForces += addedLift; */

  vMoments = in.AeroMoment;
  vMoments += in.PropMoment;
  vMoments += in.GroundMoment;
  vMoments += in.ExternalMoment;
  vMoments += in.BuoyantMoment;

  FGColumnVector3 addedMoment = FDMExec->GetAuxiliary()->resultMoment(); //fontionne avec l'erreur
  vMoments += addedMoment;

  //std::cout << "Moments = " << addedMoment << std::endl;

  /* double roll_error = 0.0;
  double integral_roll = 0.0;
  double derivative_roll = 0.0;
  double previous_roll_error = 0.0;

  // Gains PID 
  double roll_p_gain = 0.6;
  double roll_i_gain = 0.4;
  double roll_d_gain = 0.3;

  double roll_angle = FDMExec->GetPropagate()->GetEuler(1);
  double desired_roll_angle = 0.0;

  roll_error = desired_roll_angle - roll_angle;

  // Terme proportionnel
  double p_term = roll_p_gain * roll_error;

  // Terme intégral
  integral_roll += roll_error;

  // Terme dérivatif
  derivative_roll = roll_error - previous_roll_error;

  double aileron_command = p_term + (roll_i_gain * integral_roll) + (roll_d_gain * derivative_roll);

  previous_roll_error = roll_error;

  aileron_command = std::max(0.0, std::min(1.0, aileron_command)); */

  //FDMExec->GetFCS()->SetDaCmd(aileron_command);

  RunPostFunctions();

  return false;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool FGAircraft::Load(Element* el)
{
  string element_name;
  Element* element;

  if (!FGModel::Upload(el, true)) return false;

  if (el->FindElement("wingarea"))
    WingArea = el->FindElementValueAsNumberConvertTo("wingarea", "FT2");
  if (el->FindElement("wingspan"))
    WingSpan = el->FindElementValueAsNumberConvertTo("wingspan", "FT");
  if (el->FindElement("chord"))
    cbar = el->FindElementValueAsNumberConvertTo("chord", "FT");
  if (el->FindElement("wing_incidence"))
    WingIncidence = el->FindElementValueAsNumberConvertTo("wing_incidence", "RAD");
  if (el->FindElement("htailarea"))
    HTailArea = el->FindElementValueAsNumberConvertTo("htailarea", "FT2");
  if (el->FindElement("htailarm"))
    HTailArm = el->FindElementValueAsNumberConvertTo("htailarm", "FT");
  if (el->FindElement("vtailarea"))
    VTailArea = el->FindElementValueAsNumberConvertTo("vtailarea", "FT2");
  if (el->FindElement("vtailarm"))
    VTailArm = el->FindElementValueAsNumberConvertTo("vtailarm", "FT");

  // Find all LOCATION elements that descend from this METRICS branch of the
  // config file. This would be CG location, eyepoint, etc.

  element = el->FindElement("location");
  while (element) {
    element_name = element->GetAttributeValue("name");

    if (element_name == "AERORP") vXYZrp = element->FindElementTripletConvertTo("IN");
    else if (element_name == "EYEPOINT") vXYZep = element->FindElementTripletConvertTo("IN");
    else if (element_name == "VRP") vXYZvrp = element->FindElementTripletConvertTo("IN");

    element = el->FindNextElement("location");
  }

  // calculate some derived parameters
  if (cbar != 0.0) {
    lbarh = HTailArm/cbar;
    lbarv = VTailArm/cbar;
    if (WingArea != 0.0) {
      vbarh = HTailArm*HTailArea / (cbar*WingArea);
      vbarv = VTailArm*VTailArea / (WingSpan*WingArea);
    }
  }

  PostLoad(el, FDMExec);

  Debug(2);

  return true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void FGAircraft::bind(void)
{
  typedef double (FGAircraft::*PMF)(int) const;
  PropertyManager->Tie("metrics/Sw-sqft", this, &FGAircraft::GetWingArea, &FGAircraft::SetWingArea);
  PropertyManager->Tie("metrics/bw-ft", this, &FGAircraft::GetWingSpan);
  PropertyManager->Tie("metrics/cbarw-ft", this, &FGAircraft::Getcbar);
  PropertyManager->Tie("metrics/iw-rad", this, &FGAircraft::GetWingIncidence);
  PropertyManager->Tie("metrics/iw-deg", this, &FGAircraft::GetWingIncidenceDeg);
  PropertyManager->Tie("metrics/Sh-sqft", this, &FGAircraft::GetHTailArea);
  PropertyManager->Tie("metrics/lh-ft", this, &FGAircraft::GetHTailArm);
  PropertyManager->Tie("metrics/Sv-sqft", this, &FGAircraft::GetVTailArea);
  PropertyManager->Tie("metrics/lv-ft", this, &FGAircraft::GetVTailArm);
  PropertyManager->Tie("metrics/lh-norm", this, &FGAircraft::Getlbarh);
  PropertyManager->Tie("metrics/lv-norm", this, &FGAircraft::Getlbarv);
  PropertyManager->Tie("metrics/vbarh-norm", this, &FGAircraft::Getvbarh);
  PropertyManager->Tie("metrics/vbarv-norm", this, &FGAircraft::Getvbarv);
  PropertyManager->Tie("metrics/aero-rp-x-in", this, eX, (PMF)&FGAircraft::GetXYZrp, &FGAircraft::SetXYZrp);
  PropertyManager->Tie("metrics/aero-rp-y-in", this, eY, (PMF)&FGAircraft::GetXYZrp, &FGAircraft::SetXYZrp);
  PropertyManager->Tie("metrics/aero-rp-z-in", this, eZ, (PMF)&FGAircraft::GetXYZrp, &FGAircraft::SetXYZrp);
  PropertyManager->Tie("metrics/eyepoint-x-in", this, eX, (PMF)&FGAircraft::GetXYZep);
  PropertyManager->Tie("metrics/eyepoint-y-in", this, eY,(PMF)&FGAircraft::GetXYZep);
  PropertyManager->Tie("metrics/eyepoint-z-in", this, eZ, (PMF)&FGAircraft::GetXYZep);
  PropertyManager->Tie("metrics/visualrefpoint-x-in", this, eX, (PMF)&FGAircraft::GetXYZvrp);
  PropertyManager->Tie("metrics/visualrefpoint-y-in", this, eY, (PMF)&FGAircraft::GetXYZvrp);
  PropertyManager->Tie("metrics/visualrefpoint-z-in", this, eZ, (PMF)&FGAircraft::GetXYZvrp);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    The bitmasked value choices are as follows:
//    unset: In this case (the default) JSBSim would only print
//       out the normally expected messages, essentially echoing
//       the config files as they are read. If the environment
//       variable is not set, debug_lvl is set to 1 internally
//    0: This requests JSBSim not to output any messages
//       whatsoever.
//    1: This value explicity requests the normal JSBSim
//       startup messages
//    2: This value asks for a message to be printed out when
//       a class is instantiated
//    4: When this value is set, a message is displayed when a
//       FGModel object executes its Run() method
//    8: When this value is set, various runtime state variables
//       are printed out periodically
//    16: When set various parameters are sanity checked and
//       a message is printed out when they go out of bounds

void FGAircraft::Debug(int from)
{
  if (debug_lvl <= 0) return;

  if (debug_lvl & 1) { // Standard console startup message output
    if (from == 2) { // Loading
      cout << endl << "  Aircraft Metrics:" << endl;
      cout << "    WingArea: " << WingArea  << endl;
      cout << "    WingSpan: " << WingSpan  << endl;
      cout << "    Incidence: " << WingIncidence << endl;
      cout << "    Chord: " << cbar << endl;
      cout << "    H. Tail Area: " << HTailArea << endl;
      cout << "    H. Tail Arm: " << HTailArm << endl;
      cout << "    V. Tail Area: " << VTailArea << endl;
      cout << "    V. Tail Arm: " << VTailArm << endl;
      cout << "    Eyepoint (x, y, z): " << vXYZep << endl;
      cout << "    Ref Pt (x, y, z): " << vXYZrp << endl;
      cout << "    Visual Ref Pt (x, y, z): " << vXYZvrp << endl;
    }
  }
  if (debug_lvl & 2 ) { // Instantiation/Destruction notification
    if (from == 0) cout << "Instantiated: FGAircraft" << endl;
    if (from == 1) cout << "Destroyed:    FGAircraft" << endl;
  }
  if (debug_lvl & 4 ) { // Run() method entry print for FGModel-derived objects
  }
  if (debug_lvl & 8 ) { // Runtime state variables
  }
  if (debug_lvl & 16) { // Sanity checking
  }
  if (debug_lvl & 64) {
    if (from == 0) { // Constructor
    }
  }
}

void FGAircraft::virage(double time1, double time2, double angle)
{
  // Variables de PID
  double roll_error = 0.0;
  double integral_roll = 0.0;
  double derivative_roll = 0.0;
  double previous_roll_error = 0.0;
  double dt = 0.08333;
 
  // Gains PID
  double roll_p_gain = 20.0;
  double roll_i_gain = 0.0;
  double roll_d_gain = 0.0;
 
  double roll_angle = FDMExec->GetPropagate()->GetEuler(1);//* 180/PI;
  //std::cout << "angle_roll  = " << roll_angle<<std::endl;
  double desired_roll_angle = (angle * PI)  / 180;
  //std::cout << "angle_voulue  = " << desired_roll_angle<<std::endl;
  roll_error = desired_roll_angle - roll_angle;
  //std::cout << "erreur  = " << roll_error<<std::endl;
 
  // Terme proportionnel
  double p_term = roll_p_gain * roll_error;
  //std::cout << "terme_p  = " << p_term<<std::endl;
 
  // Terme intégral
  error.push_back(roll_error);
  for ( double elem : error) {integral_roll += elem;}
  integral_roll *= dt;
  //std::cout << "terme int  = " << integral_roll<<std::endl;
 
 
 
  // Terme dérivatif
  derivative_roll = (roll_error - previous_roll_error) / dt;
 
  double aileron_command = p_term + (roll_i_gain * integral_roll) + (roll_d_gain * derivative_roll);
 
 
  if (aileron_command < -1 ){aileron_command = -1;}
 
  if (aileron_command > 1 ){aileron_command = 1;}
 
  //std::cout << "aileron_command  = " << aileron_command<<std::endl;
 
  previous_roll_error = roll_error;
 
  double T  = FDMExec->GetSimTime();
 
  /* if (time1 <= T && T< time2)
  {
    FDMExec->GetFCS()->SetDaCmd(aileron_command);
    //FDMExec->GetFCS()->SetRollTrimCmd(trim_command);
  } */

  FDMExec->GetFCS()->SetDaCmd(aileron_command);
  FDMExec->GetFCS()->SetDeCmd(-0.25);
}

} // namespace JSBSim
