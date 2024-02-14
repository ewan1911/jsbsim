/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Module:       FGAuxiliary.cpp
 Author:       Tony Peden, Jon Berndt
 Date started: 01/26/99
 Purpose:      Calculates additional parameters needed by the visual system, etc.
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
This class calculates various auxiliary parameters.

REFERENCES
  Anderson, John D. "Introduction to Flight", 3rd Edition, McGraw-Hill, 1989
                    pgs. 112-126
HISTORY
--------------------------------------------------------------------------------
01/26/99   JSB   Created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INCLUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <iostream>
#include <fstream>

#include "FGAuxiliary.h"
#include "initialization/FGInitialCondition.h"
#include "FGFDMExec.h"
#include "input_output/FGPropertyManager.h"
#include "FGInertial.h"
#include "FGAtmosphere.h"

#define PI = 3.141593

using namespace std;

namespace JSBSim {

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


int points = 4; //nombre de points de part et d'autre du centre. Ici arbitraire.
double width = in.Wingspan*0.3048;;
double dw = width/(2*points);

FGAuxiliary::FGAuxiliary(FGFDMExec* fdmex) : FGModel(fdmex)
{

  // AJOUTS POINTEUR POUR 
  std::shared_ptr<JSBSim::FGPropagate> propagatePtr = fdmex->GetPropagate(); // VOIR CHATGPT JSBSIM CODE EXPL
  Propagate = propagatePtr.get();

  
  Name = "FGAuxiliary";
  pt = FGAtmosphere::StdDaySLpressure;     // ISA SL pressure
  tat = FGAtmosphere::StdDaySLtemperature; // ISA SL temperature
  tatc = RankineToCelsius(tat);

  vcas = veas = 0.0;
  qbar = qbarUW = qbarUV = 0.0;
  Mach = MachU = 0.0;
  alpha = beta = 0.0;
  adot = bdot = 0.0;
  gamma = Vt = Vground = 0.0;
  psigt = 0.0;
  hoverbmac = hoverbcg = 0.0;
  Re = 0.0;
  Nx = Ny = Nz = 0.0;

  vPilotAccel.InitMatrix();
  vPilotAccelN.InitMatrix();
  vAeroUVW.InitMatrix();
  vAeroPQR.InitMatrix();
  vMachUVW.InitMatrix();
  vEulerRates.InitMatrix();

  bind();

  Debug(0);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool FGAuxiliary::InitModel(void)
{
  if (!FGModel::InitModel()) return false;

  pt = in.Pressure;
  tat = in.Temperature;
  tatc = RankineToCelsius(tat);

  vcas = veas = 0.0;
  qbar = qbarUW = qbarUV = 0.0;
  Mach = MachU = 0.0;
  alpha = beta = 0.0;
  adot = bdot = 0.0;
  gamma = Vt = Vground = 0.0;
  psigt = 0.0;
  hoverbmac = hoverbcg = 0.0;
  Re = 0.0;
  Nz = Ny = 0.0;

  vPilotAccel.InitMatrix();
  vPilotAccelN.InitMatrix();
  vAeroUVW.InitMatrix();
  vAeroPQR.InitMatrix();
  vMachUVW.InitMatrix();
  vEulerRates.InitMatrix();

    // ON LOAD DANS INITMODEL CAR LA FONCTION EST APPELEE QU'UNE SEULE FOIS PENDANT LA COMPILATION 
  //
  //
  //  Grid = [  [Hauteurs de la grid ( 128 valeurs)],  [Longueur de la grid (257 valeurs)],  [Largeu de la grid (257 valeurs )]    ]
  //  
  //  u, v et w sont des tableaux à trois dimension  
  // 

  
  loaduwind();
  loadvwind();
  loadwwind();
  loadgrid();

  return true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FGAuxiliary::~FGAuxiliary()
{
  Debug(1);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool FGAuxiliary::Run(bool Holding)
{
  if (FGModel::Run(Holding)) return true; // return true if error returned from base class
  if (Holding) return false;

  // Rotation

  vEulerRates(eTht) = in.vPQR(eQ)*in.CosPhi - in.vPQR(eR)*in.SinPhi;
  if (in.CosTht != 0.0) {
    vEulerRates(ePsi) = (in.vPQR(eQ)*in.SinPhi + in.vPQR(eR)*in.CosPhi)/in.CosTht;
    vEulerRates(ePhi) = in.vPQR(eP) + vEulerRates(ePsi)*in.SinTht;
  }

  // Combine the wind speed with aircraft speed to obtain wind relative speed
  vAeroPQR = in.vPQR - in.TurbPQR;
  vAeroUVW = in.vUVW - in.Tl2b * in.TotalWindNED;

  alpha = beta = adot = bdot = 0;
  double AeroU2 = vAeroUVW(eU)*vAeroUVW(eU);
  double AeroV2 = vAeroUVW(eV)*vAeroUVW(eV);
  double AeroW2 = vAeroUVW(eW)*vAeroUVW(eW);
  double mUW = AeroU2 + AeroW2;

  double Vt2 = mUW + AeroV2;
  Vt = sqrt(Vt2);

  if ( Vt > 0.001 ) {
    beta = atan2(vAeroUVW(eV), sqrt(mUW));

    if ( mUW >= 1E-6 ) {
      alpha = atan2(vAeroUVW(eW), vAeroUVW(eU));
      double Vtdot = (vAeroUVW(eU)*in.vUVWdot(eU) + vAeroUVW(eV)*in.vUVWdot(eV) + vAeroUVW(eW)*in.vUVWdot(eW))/Vt;
      adot = (vAeroUVW(eU)*in.vUVWdot(eW) - vAeroUVW(eW)*in.vUVWdot(eU))/mUW;
      bdot = (in.vUVWdot(eV)*Vt - vAeroUVW(eV)*Vtdot)/(Vt*sqrt(mUW));
    }
  }

  UpdateWindMatrices();

  //printf("Je suis dans auxiliary \n");

  Re = Vt * in.Wingchord / in.KinematicViscosity;

  double densityD2 = 0.5*in.Density;

  qbar = densityD2 * Vt2;
  qbarUW = densityD2 * (mUW);
  qbarUV = densityD2 * (AeroU2 + AeroV2);
  Mach = Vt / in.SoundSpeed;
  MachU = vMachUVW(eU) = vAeroUVW(eU) / in.SoundSpeed;
  vMachUVW(eV) = vAeroUVW(eV) / in.SoundSpeed;
  vMachUVW(eW) = vAeroUVW(eW) / in.SoundSpeed;

  // Position

  Vground = sqrt( in.vVel(eNorth)*in.vVel(eNorth) + in.vVel(eEast)*in.vVel(eEast) );

  psigt = atan2(in.vVel(eEast), in.vVel(eNorth));
  if (psigt < 0.0) psigt += 2*M_PI;
  gamma = atan2(-in.vVel(eDown), Vground);

  tat = in.Temperature*(1 + 0.2*Mach*Mach); // Total Temperature, isentropic flow
  tatc = RankineToCelsius(tat);

  pt = PitotTotalPressure(Mach, in.Pressure);

  if (abs(Mach) > 0.0) {
    vcas = VcalibratedFromMach(Mach, in.Pressure);
    veas = sqrt(2 * qbar / FGAtmosphere::StdDaySLdensity);
  }
  else
    vcas = veas = 0.0;

  vPilotAccel.InitMatrix();
  vNcg = in.vBodyAccel/in.StandardGravity;
  // Nz is Acceleration in "g's", along normal axis (-Z body axis)
  Nz = -vNcg(eZ);
  Ny =  vNcg(eY);
  Nx =  vNcg(eX);
  vPilotAccel = in.vBodyAccel + in.vPQRidot * in.ToEyePt;
  vPilotAccel += in.vPQRi * (in.vPQRi * in.ToEyePt);

  vNwcg = mTb2w * vNcg;
  vNwcg(eZ) = 1.0 - vNwcg(eZ);

  vPilotAccelN = vPilotAccel / in.StandardGravity;

  // VRP computation
  vLocationVRP = in.vLocation.LocalToLocation( in.Tb2l * in.VRPBody );

  // Recompute some derived values now that we know the dependent parameters values ...
  hoverbcg = in.DistanceAGL / in.Wingspan;

  FGColumnVector3 vMac = in.Tb2l * in.RPBody;
  hoverbmac = (in.DistanceAGL - vMac(3)) / in.Wingspan;


  // FONCTION RAJOUTEES

  double dist_long = GetLongitudeRelativePosition() * 0.3048;

  double dist_lat = GetLatitudeRelativePosition() * 0.3048;

  double dist_rel = GetDistanceRelativePosition() * 0.3048;

  double alt = Propagate->GetAltitudeASL()*0.3048;

  double lon_deg = Propagate->GetLongitudeDeg();

  double lat_deg = Propagate->GetGeodLatitudeDeg();

  double gride = grid[0][0];

  rechercheNoeuds(alt, dist_lat, dist_long, 4000.0, lon_deg);
  //std::cout << "-----------------------------------------------------------------------------------" <<std::endl;
  //std::cout << "Long = " << lon_deg << " [°]" <<std::endl;
  //std::cout << "Lat = " << lat_deg << " [°]" <<std::endl;
  //std::cout << "--------" <<std::endl;
  //std::cout << "avancement depuis la longitude (East) depuis la position initiale =  " << dist_long << " [m]" << std::endl;
  //std::cout << "avancement depuis la latitude (North)  depuis la position initiale = " << dist_lat << " [m]" <<std::endl;
  //std::cout << "Distance parcourue depuis la position initiale = " << dist_rel << " [m]" <<std::endl;
  //std::cout << "Altitude = " << alt << " [m]" <<std::endl;






// Hello I'm Simon and I'm from Belgium 






  return false;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::PitotTotalPressure(double mach, double pressure) const
{
  constexpr double SHRatio = FGAtmosphere::SHRatio;
  constexpr double a = (SHRatio-1.0) / 2.0;
  constexpr double b = SHRatio / (SHRatio-1.0);
  constexpr double c = 2.0*b;
  constexpr double d = 1.0 / (SHRatio-1.0);
  const double coeff = pow(0.5*(SHRatio+1.0), b)
                     * pow((SHRatio+1.0)/(SHRatio-1.0), d);

  if (mach < 0) return pressure;
  if (mach < 1)    //calculate total pressure assuming isentropic flow
    return pressure*pow((1.0 + a*mach*mach), b);
  else {
    // Shock in front of pitot tube, we'll assume its normal and use the
    // Rayleigh Pitot Tube Formula, i.e. the ratio of total pressure behind the
    // shock to the static pressure in front of the normal shock assumption
    // should not be a bad one -- most supersonic aircraft place the pitot probe
    // out front so that it is the forward most point on the aircraft.
    // The real shock would, of course, take on something like the shape of a
    // rounded-off cone but, here again, the assumption should be good since the
    // opening of the pitot probe is very small and, therefore, the effects of
    // the shock curvature should be small as well. AFAIK, this approach is
    // fairly well accepted within the aerospace community

    // The denominator below is zero for Mach ~ 0.38, for which
    // we'll never be here, so we're safe

    return pressure*coeff*pow(mach, c)/pow(c*mach*mach-1.0, d);
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Based on the formulas in the US Air Force Aircraft Performance Flight Testing
// Manual (AFFTC-TIH-99-01). In particular sections 4.6 to 4.8.

double FGAuxiliary::MachFromImpactPressure(double qc, double pressure) const
{
  constexpr double SHRatio = FGAtmosphere::SHRatio;
  constexpr double a = 2.0/(SHRatio-1.0);
  constexpr double b = (SHRatio-1.0)/SHRatio;
  constexpr double c = 2.0/b;
  constexpr double d = 0.5*a;
  const double coeff = pow(0.5*(SHRatio+1.0), -0.25*c)
                     * pow(0.5*(SHRatio+1.0)/SHRatio, -0.5*d);

  double A = qc / pressure + 1;
  double M = sqrt(a*(pow(A, b) - 1.0));  // Equation (4.12)

  if (M > 1.0)
    for (unsigned int i = 0; i<10; i++)
      M = coeff*sqrt(A*pow(1 - 1.0 / (c*M*M), d));  // Equation (4.17)

  return M;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::VcalibratedFromMach(double mach, double pressure) const
{
  double qc = PitotTotalPressure(mach, pressure) - pressure;
  return in.StdDaySLsoundspeed * MachFromImpactPressure(qc, FGAtmosphere::StdDaySLpressure);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::MachFromVcalibrated(double vcas, double pressure) const
{
  constexpr double StdDaySLpressure = FGAtmosphere::StdDaySLpressure;
  double qc = PitotTotalPressure(vcas / in.StdDaySLsoundspeed, StdDaySLpressure) - StdDaySLpressure;
  return MachFromImpactPressure(qc, pressure);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// From Stevens and Lewis, "Aircraft Control and Simulation", 3rd Ed., the
// transformation from body to wind axes is defined (where "a" is alpha and "B"
// is beta):
//
//   cos(a)*cos(B)     sin(B)    sin(a)*cos(B)
//  -cos(a)*sin(B)     cos(B)   -sin(a)*sin(B)
//  -sin(a)              0       cos(a)
//
// The transform from wind to body axes is then,
//
//   cos(a)*cos(B)  -cos(a)*sin(B)  -sin(a)
//          sin(B)          cos(B)     0
//   sin(a)*cos(B)  -sin(a)*sin(B)   cos(a)

void FGAuxiliary::UpdateWindMatrices(void)
{
  double ca, cb, sa, sb;

  ca = cos(alpha);
  sa = sin(alpha);
  cb = cos(beta);
  sb = sin(beta);

  mTw2b(1,1) =  ca*cb;
  mTw2b(1,2) = -ca*sb;
  mTw2b(1,3) = -sa;
  mTw2b(2,1) =  sb;
  mTw2b(2,2) =  cb;
  mTw2b(2,3) =  0.0;
  mTw2b(3,1) =  sa*cb;
  mTw2b(3,2) = -sa*sb;
  mTw2b(3,3) =  ca;

  mTb2w = mTw2b.Transposed();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::GetNlf(void) const
{
  if (in.Mass != 0)
    return (in.vFw(3))/(in.Mass*slugtolb);
  else
    return 0.;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::GetLongitudeRelativePosition(void) const
{
  return in.vLocation.GetDistanceTo(FDMExec->GetIC()->GetLongitudeRadIC(),
                                    in.vLocation.GetGeodLatitudeRad())* fttom;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::GetLatitudeRelativePosition(void) const
{
  return in.vLocation.GetDistanceTo(in.vLocation.GetLongitude(),
                                    FDMExec->GetIC()->GetGeodLatitudeRadIC())* fttom;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::GetDistanceRelativePosition(void) const
{
  auto ic = FDMExec->GetIC();
  return in.vLocation.GetDistanceTo(ic->GetLongitudeRadIC(),
                                    ic->GetGeodLatitudeRadIC())* fttom;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void FGAuxiliary::bind(void)
{
  typedef double (FGAuxiliary::*PMF)(int) const;
  typedef double (FGAuxiliary::*PF)(void) const;
  PropertyManager->Tie("propulsion/tat-r", this, &FGAuxiliary::GetTotalTemperature);
  PropertyManager->Tie("propulsion/tat-c", this, &FGAuxiliary::GetTAT_C);
  PropertyManager->Tie("propulsion/pt-lbs_sqft", this, &FGAuxiliary::GetTotalPressure);
  PropertyManager->Tie("velocities/vc-fps", this, &FGAuxiliary::GetVcalibratedFPS);
  PropertyManager->Tie("velocities/vc-kts", this, &FGAuxiliary::GetVcalibratedKTS);
  PropertyManager->Tie("velocities/ve-fps", this, &FGAuxiliary::GetVequivalentFPS);
  PropertyManager->Tie("velocities/ve-kts", this, &FGAuxiliary::GetVequivalentKTS);
  PropertyManager->Tie("velocities/vtrue-fps", this, &FGAuxiliary::GetVtrueFPS);
  PropertyManager->Tie("velocities/vtrue-kts", this, &FGAuxiliary::GetVtrueKTS);
  PropertyManager->Tie("velocities/machU", this, &FGAuxiliary::GetMachU);
  PropertyManager->Tie("velocities/p-aero-rad_sec", this, eX, (PMF)&FGAuxiliary::GetAeroPQR);
  PropertyManager->Tie("velocities/q-aero-rad_sec", this, eY, (PMF)&FGAuxiliary::GetAeroPQR);
  PropertyManager->Tie("velocities/r-aero-rad_sec", this, eZ, (PMF)&FGAuxiliary::GetAeroPQR);
  PropertyManager->Tie("velocities/phidot-rad_sec", this, ePhi, (PMF)&FGAuxiliary::GetEulerRates);
  PropertyManager->Tie("velocities/thetadot-rad_sec", this, eTht, (PMF)&FGAuxiliary::GetEulerRates);
  PropertyManager->Tie("velocities/psidot-rad_sec", this, ePsi, (PMF)&FGAuxiliary::GetEulerRates);
  PropertyManager->Tie("velocities/u-aero-fps", this, eU, (PMF)&FGAuxiliary::GetAeroUVW);
  PropertyManager->Tie("velocities/v-aero-fps", this, eV, (PMF)&FGAuxiliary::GetAeroUVW);
  PropertyManager->Tie("velocities/w-aero-fps", this, eW, (PMF)&FGAuxiliary::GetAeroUVW);
  PropertyManager->Tie("velocities/vt-fps", this, &FGAuxiliary::GetVt);
  PropertyManager->Tie("velocities/mach", this, &FGAuxiliary::GetMach);
  PropertyManager->Tie("velocities/vg-fps", this, &FGAuxiliary::GetVground);
  PropertyManager->Tie("accelerations/a-pilot-x-ft_sec2", this, eX, (PMF)&FGAuxiliary::GetPilotAccel);
  PropertyManager->Tie("accelerations/a-pilot-y-ft_sec2", this, eY, (PMF)&FGAuxiliary::GetPilotAccel);
  PropertyManager->Tie("accelerations/a-pilot-z-ft_sec2", this, eZ, (PMF)&FGAuxiliary::GetPilotAccel);
  PropertyManager->Tie("accelerations/n-pilot-x-norm", this, eX, (PMF)&FGAuxiliary::GetNpilot);
  PropertyManager->Tie("accelerations/n-pilot-y-norm", this, eY, (PMF)&FGAuxiliary::GetNpilot);
  PropertyManager->Tie("accelerations/n-pilot-z-norm", this, eZ, (PMF)&FGAuxiliary::GetNpilot);
  PropertyManager->Tie("accelerations/Nx", this, &FGAuxiliary::GetNx);
  PropertyManager->Tie("accelerations/Ny", this, &FGAuxiliary::GetNy);
  PropertyManager->Tie("accelerations/Nz", this, &FGAuxiliary::GetNz);
  PropertyManager->Tie("forces/load-factor", this, &FGAuxiliary::GetNlf);
  PropertyManager->Tie("aero/alpha-rad", this, (PF)&FGAuxiliary::Getalpha);
  PropertyManager->Tie("aero/beta-rad", this, (PF)&FGAuxiliary::Getbeta);
  PropertyManager->Tie("aero/mag-beta-rad", this, (PF)&FGAuxiliary::GetMagBeta);
  PropertyManager->Tie("aero/alpha-deg", this, inDegrees, (PMF)&FGAuxiliary::Getalpha);
  PropertyManager->Tie("aero/beta-deg", this, inDegrees, (PMF)&FGAuxiliary::Getbeta);
  PropertyManager->Tie("aero/mag-beta-deg", this, inDegrees, (PMF)&FGAuxiliary::GetMagBeta);
  PropertyManager->Tie("aero/Re", this, &FGAuxiliary::GetReynoldsNumber);
  PropertyManager->Tie("aero/qbar-psf", this, &FGAuxiliary::Getqbar);
  PropertyManager->Tie("aero/qbarUW-psf", this, &FGAuxiliary::GetqbarUW);
  PropertyManager->Tie("aero/qbarUV-psf", this, &FGAuxiliary::GetqbarUV);
  PropertyManager->Tie("aero/alphadot-rad_sec", this, (PF)&FGAuxiliary::Getadot);
  PropertyManager->Tie("aero/betadot-rad_sec", this, (PF)&FGAuxiliary::Getbdot);
  PropertyManager->Tie("aero/alphadot-deg_sec", this, inDegrees, (PMF)&FGAuxiliary::Getadot);
  PropertyManager->Tie("aero/betadot-deg_sec", this, inDegrees, (PMF)&FGAuxiliary::Getbdot);
  PropertyManager->Tie("aero/h_b-cg-ft", this, &FGAuxiliary::GetHOverBCG);
  PropertyManager->Tie("aero/h_b-mac-ft", this, &FGAuxiliary::GetHOverBMAC);
  PropertyManager->Tie("flight-path/gamma-rad", this, &FGAuxiliary::GetGamma);
  PropertyManager->Tie("flight-path/gamma-deg", this, inDegrees, (PMF)&FGAuxiliary::GetGamma);
  PropertyManager->Tie("flight-path/psi-gt-rad", this, &FGAuxiliary::GetGroundTrack);

  PropertyManager->Tie("position/distance-from-start-lon-mt", this, &FGAuxiliary::GetLongitudeRelativePosition);
  PropertyManager->Tie("position/distance-from-start-lat-mt", this, &FGAuxiliary::GetLatitudeRelativePosition);
  PropertyManager->Tie("position/distance-from-start-mag-mt", this, &FGAuxiliary::GetDistanceRelativePosition);
  PropertyManager->Tie("position/vrp-gc-latitude_deg", &vLocationVRP, &FGLocation::GetLatitudeDeg);
  PropertyManager->Tie("position/vrp-longitude_deg", &vLocationVRP, &FGLocation::GetLongitudeDeg);
  PropertyManager->Tie("position/vrp-radius-ft", &vLocationVRP, &FGLocation::GetRadius);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double FGAuxiliary::BadUnits(void) const
{
  cerr << "Bad units" << endl; return 0.0;
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

void FGAuxiliary::Debug(int from)
{
  if (debug_lvl <= 0) return;

  if (debug_lvl & 1) { // Standard console startup message output
    if (from == 0) { // Constructor

    }
  }
  if (debug_lvl & 2 ) { // Instantiation/Destruction notification
    if (from == 0) cout << "Instantiated: FGAuxiliary" << endl;
    if (from == 1) cout << "Destroyed:    FGAuxiliary" << endl;
  }
  if (debug_lvl & 4 ) { // Run() method entry print for FGModel-derived objects
  }
  if (debug_lvl & 8 ) { // Runtime state variables
  }
  if (debug_lvl & 16) { // Sanity checking
    if (Mach > 100 || Mach < 0.00)
      cout << "FGPropagate::Mach is out of bounds: " << Mach << endl;
    if (qbar > 1e6 || qbar < 0.00)
      cout << "FGPropagate::qbar is out of bounds: " << qbar << endl;
  }
  if (debug_lvl & 64) {
    if (from == 0) { // Constructor
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   
//   On import les donnés de u, v, w et la grid
//       
//       
//       

void FGAuxiliary::loaduwind() 
{
  const int dim1 = 128;  
  const int dim2 = 257;  
  const int dim3 = 257;  

  std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/u_aplati.csv");

  if (in.is_open()) {
      std::cout << "File u opened successfully." << std::endl;
  } else {
      std::cerr << "Error uopening file zebi." << std::endl;
      return;
  }

  for (int i = 0; i < dim1; ++i) {
      for (int j = 0; j < dim2; ++j) {
          for (int k = 0; k < dim3; ++k) {
              in >> u[i][j][k];
          }
      }
  }

  in.close();
}

void FGAuxiliary::loadvwind() 
{
  const int dim1 = 128;  
  const int dim2 = 257;  
  const int dim3 = 257;  

  std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/v_aplati.csv");

  if (in.is_open()) {
      std::cout << "File v opened successfully." << std::endl;
  } else {
      std::cerr << "Error v opening file zebi." << std::endl;
      return;
  }

  for (int i = 0; i < dim1; ++i) {
      for (int j = 0; j < dim2; ++j) {
          for (int k = 0; k < dim3; ++k) {
              in >> v[i][j][k];
          }
      }
  }

  in.close();
}

void FGAuxiliary::loadwwind() 
{
  const int dim1 = 128;  
  const int dim2 = 257;  
  const int dim3 = 257;  

  std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/w_aplati.csv");

  if (in.is_open()) {
      std::cout << "File w opened successfully." << std::endl;
  } else {
      std::cerr << "Error w  opening file zebi." << std::endl;
      return;
  }

  for (int i = 0; i < dim1; ++i) {
      for (int j = 0; j < dim2; ++j) {
          for (int k = 0; k < dim3; ++k) {
              in >> w[i][j][k];
          }
      }
  }

  in.close();
}

void FGAuxiliary::loadgrid() 
{
  std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/grid.csv");

  if (in.is_open()) {
      std::cout << "File grid opened successfully." << std::endl;
  } else {
      std::cerr << "Error grid  opening file zebi." << std::endl;
      return;
  }

  for (int i = 0; i < 3; ++i) {

    if (i == 0 )
    {
      for (int j = 0; j < 128; ++j){
        in >> grid[i][j];
      }
    }

    if (i == 1 )
    {
      for (int j = 0; j < 257; ++j){
        in >> grid[i][j];
      }
    }

    if (i == 2 )
    {
      for (int j = 0; j < 257; ++j){
        in >> grid[i][j];
      }
    }
  }

  in.close();
}

// on test
// re test encore 
// Ici on rentre dans la boite bigflow à la position initiale : ( 1000 [m] (altitude), 0 [m] (Longueur), 4000 [m] (largeur) )
// Donc je suis en face de la boite au milieu 
//
//Print les infos de où on est dans la grid de vent 
//
void FGAuxiliary::rechercheNoeuds(double x, double y, double z, double refz, double longi) //x,y et z selon la règle de la main droite
{    

    double z_bis = z;
    if (longi>0) // POUR SAVOIR SI ON TOURNE A DROITE OU A GAUCHE 
    {
      z *= -1;
    
    }
    // Parcourir le tableau hauteur x
    int indice1x = -1, indice2x = -1;
    for (int i = 0; i < 127 - 1; ++i) {
        if (grid[0][i] <= x && x <= grid[0][i + 1]) {
            indice1x = i;
            indice2x = i + 1;
            break;
        }
    }

    // Parcourir le tableau longueur y
    int indice1y = -1, indice2y = -1;
    for (int i = 0; i < 256 - 1; ++i) {
        if (grid[1][i] <= y && y <= grid[1][i + 1]) {
            indice1y = i;
            indice2y = i + 1;
            break;
        }
    }

    // Parcourir le tableau largeur  z
    int indice1z = -1, indice2z = -1;
    for (int i = 0; i < 256 - 1; ++i) {
        if (grid[2][i] <= z + refz  && z + refz <= grid[2][i + 1]) {
            indice1z = i;
            indice2z = i + 1;
            break;
        }
    }



    // Afficher les indices et les distances
    std::cout << "-----------------------------------------------------------------------------" << std::endl;
    std::cout << "POSITION INITIALE : Hauteur = 1000 [m], Longueur = 0 [m], Largeur = 4000 [m]" << std::endl;

    if (indice1x != -1 && indice2x != -1) {
        std::cout << "------------------HAUTEUR------------------" << std::endl;
        std::cout << "La hauteur " << x << " [m] se trouve entre les noeuds " << indice1x << " ("<<grid[0][indice1x] << "[m]) et " << indice2x<< " ("<<grid[0][indice2x] << "[m])cd " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice1x << " = " << x - grid[0][indice1x] << "[m] " << std::endl;
        std::cout << "Distance par rapport au noeud" << indice2x << " = " << grid[0][indice2x] - x << "[m] "<< std::endl;
    } else {
        std::cout << "La valeur " << x << " n'est pas présente dans le tableau." << std::endl;
    }

    
    if (indice1y != -1 && indice2y != -1) {
        std::cout << "----------------AVANT/ARRIERE----------------" << std::endl;
        std::cout << "J'ai avance de " << y << " [m] en avant depuis l'entree dans BigFlow  "  << std::endl;
        std::cout << "Je me trouve  entre les noeuds " << indice1y << " ("<<grid[1][indice1y] << "[m]) et " << indice2y << " ("<<grid[1][indice2y] << "[m]) " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice1y << " = " << y - grid[1][indice1y] << "[m] " << std::endl;
        std::cout << "Distance par rapport au noeud" << indice2y << " = " << grid[1][indice2y] - y << "[m] "<< std::endl;
    } else {
        std::cout << "La valeur " << y << " n'est pas présente dans le tableau." << std::endl;
    }

    
    if (indice1z != -1 && indice2z != -1) {
        std::cout << "-----------------DROITE/GAUCHE---------------" << std::endl;
        if (longi > 0 )
        {
          std::cout << "J'avance vers la droite" << std::endl;
        }
        else {
          std::cout << "J'avance vers la gauche" << std::endl;
        }
        std::cout << "J'ai avance de " << z_bis << " [m] depuis l'entree dans BigFlow " << std::endl;
        std::cout << "Je me trouve  entre les noeuds " << indice1z << " ("<<grid[2][indice1z] << "[m]) et " << indice2z << " ("<<grid[2][indice2z] << "[m]) " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice1z << " = " <<  (z + refz) - grid[2][indice1z] << "[m] " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice2z << " = " << grid[2][indice2z] - (z + refz) << "[m] "<< std::endl;
    } else {
        std::cout << "La valeur " << z + refz << " n'est pas présente dans le tableau." << std::endl;
    }
}

void FGAuxiliary::discretisation(double x, double y, double z, int n){
  double D;
  double positions[2*n+1][3];
  
  double yaw = Propagate->GetEuler(ePhi);
  double pitch = Propagate->GetEuler(eTht);
  double roll = Propagate->GetEuler(ePsi);

  double theta = 3.14159/2 + yaw;
  double phi = 3.14159/2 - roll;

  for (int i = 0; i < 2*n+1; i++)
  {
    if (i<=n)
    {
      D = (n-i)*dw;
      positions[i][0] = y - D*sin(phi)*cos(theta);
      positions[i][1] = z - D*sin(phi)*sin(theta);
      positions[i][2] = x - D*cos(phi);
    } else {
      D = (i-n)*dw;
      positions[i][0] = y + D*sin(phi)*cos(theta);
      positions[i][1] = z + D*sin(phi)*sin(theta);
      positions[i][2] = x + D*cos(phi);
    }
    
  }
  std::cout << "yaw: " << yaw << " pitch: " << pitch << " roll: " << roll << std::endl;
  std::cout << "Position tip gauche : (" << positions[0][0] << "," << positions[0][1] << "," << positions[0][2] << ")" << std::endl;
  std::cout << "(" << positions[1][0] << "," << positions[1][1] << "," << positions[1][2] << ")" << std::endl;
  std::cout << "(" << positions[2][0] << "," << positions[2][1] << "," << positions[2][2] << ")" << std::endl;
  std::cout << "(" << positions[3][0] << "," << positions[3][1] << "," << positions[3][2] << ")" << std::endl;
  std::cout << "Position CG : (" << positions[n][0] << "," << positions[n][1] << "," << positions[n][2] << ")" << std::endl;
  std::cout << "(" << positions[5][0] << "," << positions[5][1] << "," << positions[5][2] << ")" << std::endl;
  std::cout << "(" << positions[6][0] << "," << positions[6][1] << "," << positions[6][2] << ")" << std::endl;
  std::cout << "(" << positions[7][0] << "," << positions[7][1] << "," << positions[7][2] << ")" << std::endl;
  std::cout << "Position tip droit : (" << positions[2*n][0] << "," << positions[2*n][1] << "," << positions[2*n][2] << ")" << std::endl;
  std::cout << "-----------------------------------------------------------------------------" << std::endl;

}

FGColumnVector3 velCG;
FGMatrix33 TranfoNED2B;


void FGAuxiliary::dynamics(double **vBoite, int n) { //n le nombre d'éléments de par et d'autre CG
  double velFlowWing[2*n+1][3];
  double Clift[2*n+1];
  double lift[2*n+1];
  double alpha_e[2*n+1];

  double rho = (FDMExec->GetAtmosphere()->GetDensity())*515.378818;
  double b = in.Wingspan*0.3048;
	double S   = (FDMExec->GetAircraft()->GetWingArea())*0.092903;
  double AR = b*b/S;
  double c = in.Wingchord*0.3048; //constant pour le moment
  double U_inf;

  velCG = in.uUVW*0.3048; //in BODY frame (?)
  TransfoNED2B = in.Tl2b;

  for (int i = 0; i < 2*n+1; i++)
  {
    velFlowWing[i] = velCG - TranfoNED2B*vBoite[i];
    Clift[i] = 2*PI*(AR/(AR+2));
    alpha_e[i] = atan2(velFlowWing[i][2], velFlowWing[i][0]);
    U_inf = sqrt(velFlowWing[i][2]*velFlowWing[i][2] + velFlowWing[i][0]*velFlowWing[i][0]);
    lift[i] = 0.5*rho*U_inf*U_inf*Clift[i]*c;
  }

  double rollMoment = 0.0;
  for (int i = 0; i < 2*n+1; i++)
  {
    if (i<n)
    {
      rollMoment += lift[i]*(n-i)*dw;
    } else {
      rollMoment -= lift[i]*(i-n)*dw; //lift partie droite de l'aile contribue négativement au rolling moment
    }
  }
}

} // namespace JSBSim
