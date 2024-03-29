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
#include "FGAircraft.h"
#include "initialization/FGInitialCondition.h"
#include "FGFDMExec.h"
#include "input_output/FGPropertyManager.h"
#include "FGInertial.h"
#include "FGAtmosphere.h"

#include "FGFCS.h"

#define PI 3.141593

using namespace std;

namespace JSBSim {

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


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
  std::cout << "-----------------------------------------TEST--------------------------------------" <<std::endl;


  return true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FGAuxiliary::~FGAuxiliary()
{
  Debug(1);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int points;

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
  ajouterDonnees("Zzz_Alt",alt);

  double lon_deg = Propagate->GetLongitudeDeg();

  double lat_deg = Propagate->GetGeodLatitudeDeg();

  double gride = grid[0][0];

  points = 40; //nombre de points de part et d'autre du centre. Ici arbitraire.

  int vBoite[5][3] = {{0, 0, 3}, {0, 0, 2}, {0, 0, 1}, {0, 0, 0}, {0, 0, 0}};
  

  /* rechercheNoeuds(alt, dist_lat, dist_long, 4000.0,100.0, lon_deg);
  //std::cout << "-----------------------------------------------------------------------------------" <<std::endl;
  //std::cout << "Long = " << lon_deg << " [°]" <<std::endl;
  //std::cout << "Lat = " << lat_deg << " [°]" <<std::endl;
  //std::cout << "--------" <<std::endl;
  //std::cout << "avancement depuis la longitude (East) depuis la position initiale =  " << dist_long << " [m]" << std::endl;
  //std::cout << "avancement depuis la latitude (North)  depuis la position initiale = " << dist_lat << " [m]" <<std::endl;
  //std::cout << "Distance parcourue depuis la position initiale = " << dist_rel << " [m]" <<std::endl;
  //std::cout << "Altitude = " << alt << " [m]" <<std::endl;

  discretisation(alt, dist_lat, dist_long, points);
  dynamics(vBoite, points); */

  //getRollMoment(alt, dist_lat, dist_long, lon_deg, points, 4000.0, 100.0);


  double East_init = 5000.0; //position de départ dans la boite
  double North_init = 3000.0;

  double East_target = 7000.0; //objectif de position a atteindre
  double North_target = 6000.0;

  int box = 1; //Si box = 1, on est dans la boite. HEREEEEEEEEEEEEEEEEEEE

  double East_pos;
  double North_pos;

  if (lon_deg < 0)
  {
    dist_long *= -1;
  }

  if (lat_deg < 0)
  {
    dist_lat *= -1;
  }

  if (box == 1) //Dans boite
  {
    East_pos = dist_long + East_init; // on applique l'offset de la boite
    North_pos = dist_lat + North_init;
    getRollMoment(alt, dist_lat, dist_long, lon_deg, lat_deg, points, East_init, North_init);
  } 
  else
  {
    East_pos = dist_long; 
    North_pos = dist_lat;
  }

  ajouterDonnees("Zzz_North", North_pos);
  ajouterDonnees("Zzz_East", East_pos);
  
  //goTo(East_target, North_target, East_pos, North_pos);
  autopilot(East_target, North_target, East_pos, North_pos);

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

  //std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/u_aplati.csv");
  std::ifstream in("/Users/Simon/Documents/Aaa_Thesis/jsbsim-code/src/models/atmosphere/u_aplati_1.csv");

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

  //std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/v_aplati.csv");
  std::ifstream in("/Users/Simon/Documents/Aaa_Thesis/jsbsim-code/src/models/atmosphere/v_aplati_1.csv");

  std::cout << "TEST" << std::endl;

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

  //std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/w_aplati.csv");
  std::ifstream in("/Users/Simon/Documents/Aaa_Thesis/jsbsim-code/src/models/atmosphere/w_aplati_1.csv");

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
  //std::ifstream in("C:/Users/test/Documents/Master_2/Master_Thesis/jsbsim-master/src/models/atmosphere/grid.csv");
  std::ifstream in("/Users/Simon/Documents/Aaa_Thesis/jsbsim-code/src/models/atmosphere/grid.csv");

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
double* FGAuxiliary::rechercheNoeuds(double hauteur, double longueur, double largeur, double refz,double ref_long, double longi, double lat) //x,y et z selon la règle de la main droite
{    
    // Parcourir le tableau hauteur 
    int indice1h = -1, indice2h = -1;
    for (int i = 0; i < 127 - 1; ++i) {
        if (grid[0][i] <= hauteur && hauteur <= grid[0][i + 1]) {
            indice1h = i;
            indice2h = i + 1;
            break;
        }
    }

    //double z_bis = largeur;
    /* if (lat<0) // SI LONG < 0 ALORS ON TOURNE VERS LA GAUCHE OU ATTENDRE 1 SEC ?
    {
      longueur *= -1;
    } */

    // Parcourir le tableau longueur 
    int indice1y = -1, indice2y = -1;
    for (int i = 0; i < 256 - 1; ++i) {
        if (grid[1][i] <= (longueur + ref_long) && (longueur + ref_long) <= grid[1][i + 1]) {
            indice1y = i;
            indice2y = i + 1;
            break;
        }
    }


    double z_bis = largeur;
    /* if (longi<0) // SI LONG < 0 ALORS ON TOURNE VERS LA GAUCHE OU ATTENDRE 1 SEC ?
    {
      largeur *= -1;
    } */
    // LE METTRE DEJA DANS LA BOITE SINON PROBLEME AVEC L INTERPOLATION
    // Parcourir le tableau largeur  
    int indice1z = -1, indice2z = -1;
    for (int i = 0; i < 256 - 1; ++i) {
        if (grid[2][i] <= largeur + refz  && largeur + refz <= grid[2][i + 1]) {
            indice1z = i;
            indice2z = i + 1;
            break;
        }
    }



    // Afficher les indices et les distances
    /* std::cout << "-----------------------------------------------------------------------------" << std::endl;
    std::cout << "POSITION INITIALE : Hauteur = 1000 [m], Longueur = 100 [m], Largeur = 4000 [m]" << std::endl;
    
    std::cout << "------------------HAUTEUR------------------" << std::endl;
    if (indice1h != -1 && indice2h != -1) {
        std::cout << "La hauteur " << hauteur << " [m] se trouve entre les noeuds " << indice1h << " ("<<grid[0][indice1h] << "[m]) et " << indice2h<< " ("<<grid[0][indice2h] << "[m])cd " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice1h << " = " << hauteur - grid[0][indice1h] << "[m] " << std::endl;
        std::cout << "Distance par rapport au noeud" << indice2h << " = " << grid[0][indice2h] - hauteur << "[m] "<< std::endl;
    } else {
        std::cout << "La valeur " << hauteur << " n'est pas présente dans le tableau." << std::endl;
    }

    std::cout << "----------------AVANT/ARRIERE----------------" << std::endl;
    if (indice1y != -1 && indice2y != -1) {
        std::cout << "J'ai avance de " << longueur  << " [m] en avant depuis l'entree dans BigFlow  "  << std::endl;
        std::cout << "Je me trouve  entre les noeuds " << indice1y << " ("<<grid[1][indice1y] << "[m]) et " << indice2y << " ("<<grid[1][indice2y] << "[m]) " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice1y << " = " << (longueur + ref_long) - grid[1][indice1y] << "[m] " << std::endl;
        std::cout << "Distance par rapport au noeud" << indice2y << " = " << grid[1][indice2y] - (longueur+ref_long) << "[m] "<< std::endl;
    } else {
        std::cout << "La valeur " << longueur + ref_long << " n'est pas présente dans le tableau." << std::endl;
    }

    std::cout << "-----------------DROITE/GAUCHE---------------" << std::endl;
    if (indice1z != -1 && indice2z != -1) {
        
        if (longi > 0 )
        {
          std::cout << "J'avance vers la droite" << std::endl;
        }
        else {
          std::cout << "J'avance vers la gauche" << std::endl;
        }
        std::cout << "J'ai avance de " << z_bis << " [m] depuis l'entree dans BigFlow " << std::endl;
        std::cout << "Je me trouve  entre les noeuds " << indice1z << " ("<<grid[2][indice1z] << "[m]) et " << indice2z << " ("<<grid[2][indice2z] << "[m]) " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice1z << " = " <<  (largeur + refz) - grid[2][indice1z] << "[m] " << std::endl;
        std::cout << "Distance par rapport au noeud " << indice2z << " = " << grid[2][indice2z] - (largeur+ refz) << "[m] "<< std::endl;
    } else {
        std::cout << "La valeur " << largeur + refz << " n'est pas présente dans le tableau." << std::endl;
    } */

 
     // RATIO 
    double delta_h = hauteur - grid[0][indice1h];
    double delta_h_1_0 = grid[0][indice2h] - grid[0][indice1h];
    double ratio_h = delta_h / delta_h_1_0;

    double delta_long = (longueur + ref_long) - grid[1][indice1y];
    double delta_long_1_0 = grid[1][indice2y] - grid[1][indice1y];
    double ratio_long = delta_long / delta_long_1_0;


    double delta_larg =  (largeur + refz) - grid[2][indice1z];
    double delta_larg_1_0 = grid[2][indice2z] - grid[2][indice1z];
    double ratio_larg = delta_larg / delta_larg_1_0;

    //########################################################### u ########################################################################
    double uc000 = u[indice1h][indice1y][indice1z];
    double uc100 = u[indice1h][indice1y][indice2z];
    double uc101 = u[indice2h][indice1y][indice2z];
    double uc001 = u[indice2h][indice1y][indice1z];
    double uc011 = u[indice2h][indice2y][indice1z];
    double uc111 = u[indice2h][indice2y][indice2z];
    double uc110 = u[indice1h][indice2y][indice2z];
    double uc010 = u[indice1h][indice2y][indice1z]; 

    // ALONG Xu 
    double uc00 = uc000 * (1-ratio_larg) + uc100 * ratio_larg;
    double uc01 = uc001 * (1-ratio_larg) + uc101 * ratio_larg;
    double uc10 = uc010 * (1-ratio_larg) + uc110 * ratio_larg;
    double uc11 = uc011 * (1-ratio_larg) + uc111 * ratio_larg;

    //ALONG Yu
    double uc0 = uc00 * (1-ratio_long) + uc10 * ratio_long;
    double uc1 = uc01 * (1-ratio_long) + uc11 * ratio_long;

    // ALONG Zu
    double uc = uc0 * (1-ratio_h) + uc1 * ratio_h;

    //########################################################### v ########################################################################
    double vc000 = v[indice1h][indice1y][indice1z];
    double vc100 = v[indice1h][indice1y][indice2z];
    double vc101 = v[indice2h][indice1y][indice2z];
    double vc001 = v[indice2h][indice1y][indice1z];
    double vc011 = v[indice2h][indice2y][indice1z];
    double vc111 = v[indice2h][indice2y][indice2z];
    double vc110 = v[indice1h][indice2y][indice2z];
    double vc010 = v[indice1h][indice2y][indice1z]; 

    // ALONG Xv 
    double vc00 = vc000 * (1-ratio_larg) + vc100 * ratio_larg;
    double vc01 = vc001 * (1-ratio_larg) + vc101 * ratio_larg;
    double vc10 = vc010 * (1-ratio_larg) + vc110 * ratio_larg;
    double vc11 = vc011 * (1-ratio_larg) + vc111 * ratio_larg;

    //ALONG Yv
    double vc0 = vc00 * (1-ratio_long) + vc10 * ratio_long;
    double vc1 = vc01 * (1-ratio_long) + vc11 * ratio_long;

    // ALONG Zv
    double vc = vc0 * (1-ratio_h) + vc1 * ratio_h;

    //########################################################### w ########################################################################
    double wc000 = w[indice1h][indice1y][indice1z];
    double wc100 = w[indice1h][indice1y][indice2z];
    double wc101 = w[indice2h][indice1y][indice2z];
    double wc001 = w[indice2h][indice1y][indice1z];
    double wc011 = w[indice2h][indice2y][indice1z];
    double wc111 = w[indice2h][indice2y][indice2z];
    double wc110 = w[indice1h][indice2y][indice2z];
    double wc010 = w[indice1h][indice2y][indice1z]; 

    // ALONG Xw 
    double wc00 = wc000 * (1-ratio_larg) + wc100 * ratio_larg;
    double wc01 = wc001 * (1-ratio_larg) + wc101 * ratio_larg;
    double wc10 = wc010 * (1-ratio_larg) + wc110 * ratio_larg;
    double wc11 = wc011 * (1-ratio_larg) + wc111 * ratio_larg;

    //ALONG Yw
    double wc0 = wc00 * (1-ratio_long) + wc10 * ratio_long;
    double wc1 = wc01 * (1-ratio_long) + wc11 * ratio_long;

    // ALONG Zw
    double wc = wc0 * (1-ratio_h) + wc1 * ratio_h;

    /* std::cout << "-----------------vitesses (u,v,w) [m/s]--------------------" << std::endl;
    std::cout << "(" << uc << ", " << vc << ", " << wc << ")" <<std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl; */

    double velocities[3] = {uc, vc, wc};

    return velocities;
}

/* void FGAuxiliary::discretisation(double x, double y, double z, int n){
  double width = in.Wingspan*0.3048;
  double dw = width/(2*n);
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

} */

FGColumnVector3 velCG;
FGColumnVector3 velCGBox;
FGColumnVector3 vPointNED;
FGColumnVector3 vPointBody;
FGColumnVector3 vFlow;
FGMatrix33 TransfoNED2B;


/* void FGAuxiliary::dynamics(int vBoite[5][3], int n) { //n le nombre d'éléments de par et d'autre CG
  double lift[2*n+1];

  double rho = (FDMExec->GetAtmosphere()->GetDensity())*515.378818;
  double b = in.Wingspan*0.3048;
  double c = in.Wingchord*0.3048; //constant pour le moment
	double S   = b*c;
  double AR = b*b/S;

  double width = in.Wingspan*0.3048;
  double dw = width/(2*n);

  double a_e;
  double C_l;
  double U_inf;

  velCG = in.vUVW*0.3048; //in BODY frame (?)
  TransfoNED2B = in.Tl2b;

  for (int i = 0; i < 2*n+1; i++)
  {
    vPoint(1) = vBoite[i][1];
    vPoint(2) = -vBoite[i][2];
    vPoint(3) = -vBoite[i][0];
    vFlow = velCG - TransfoNED2B*vPoint;
    
    a_e = atan2(vFlow(3), vFlow(1));
    C_l = 2*3.141593*(AR/(AR+2))*a_e;
    U_inf = sqrt(vFlow(3)*vFlow(3) + vFlow(1)*vFlow(1));

    lift[i] = 0.5*rho*U_inf*U_inf*C_l*c;
  }

  double rollMoment = 0.0;
  for (int i = 0; i < 2*n+1; i++)
  {
    if (i<n)
    {
      rollMoment -= lift[i]*(n-i)*dw;
    } else {
      rollMoment += lift[i]*(i-n)*dw; //lift partie droite de l'aile contribue négativement au rolling moment mais il y a un moins dans la formule.
    }
  }

  std::cout << "-----------------------------------------------------------------------------" << std::endl;
  std::cout << "Rolling moment = " << rollMoment << " Nm" << std::endl;
  std::cout << "-----------------------------------------------------------------------------" << std::endl;

} */

FGColumnVector3 boxMoment;
FGColumnVector3 liftForce;

void FGAuxiliary::getRollMoment(double hauteur, double longueur, double largeur, double longi, double lat, int n, double largeur_0, double longueur_0){
  double width = in.Wingspan*0.3048;
  double dw = width/(2*n);
  double D;
  double positions[2*n+1][3];
  double vBoite[2*n+1][3];
  
  double yaw = Propagate->GetEuler(3); //ordre yaw, pitch, roll
  double pitch = Propagate->GetEuler(2);
  double roll = Propagate->GetEuler(1);

  double theta = 3.14159/2 + yaw;
  double phi = 3.14159/2 - roll;

  double lift[2*n+1];

  double rho = (FDMExec->GetAtmosphere()->GetDensity())*515.378818;
  double b = in.Wingspan*0.3048;
  double c = in.Wingchord*0.3048; //constant pour le moment
	double S = b*c;
  double AR = b*b/S;

  double a_e;
  double C_l;
  double U_inf;

  double largeur_init = 1890.0;
  double longueur_init = 100.0;

  velCG = in.vUVW*0.3048; //in BODY frame (?)
  TransfoNED2B = in.Tl2b;

  for (int i = 0; i < 2*n+1; i++)
  {
    if (i<=n)
    {
      D = (n-i)*dw;
      positions[i][0] = longueur - D*sin(phi)*cos(theta);
      positions[i][1] = largeur - D*sin(phi)*sin(theta);
      positions[i][2] = hauteur - D*cos(phi);
    } else {
      D = (i-n)*dw;
      positions[i][0] = longueur + D*sin(phi)*cos(theta);
      positions[i][1] = largeur + D*sin(phi)*sin(theta);
      positions[i][2] = hauteur + D*cos(phi);
    }
  }

  for (int i = 0; i < 2*n+1; i++)
  {
    double* vel = rechercheNoeuds(positions[i][2], positions[i][0], positions[i][1], largeur_0, longueur_0, longi, lat);
    vBoite[i][0] = vel[0];
    vBoite[i][1] = vel[1];
    vBoite[i][2] = vel[2];
  }

  for (int i = 0; i < 2*n+1; i++)
  {
    vPointNED(1) =  vBoite[i][1];
    vPointNED(2) = -vBoite[i][2];
    vPointNED(3) = -vBoite[i][0];

    vPointBody = TransfoNED2B*vPointNED;

    if (i == n)
    {
      velCGBox = vPointBody;
    }
    

    vFlow(1) = velCG(1) + vPointBody(1); 
    vFlow(2) = velCG(2) + vPointBody(2);
    vFlow(3) = velCG(3) + vPointBody(3); //Frame pointe vers le bas, d'où le - plus bas pour a_e 

    /* std::cout << "yaw = " << yaw << " roll = " << roll << " pitch = " << pitch << std::endl;
    std::cout << "vPointNED = " << vPointNED << " vPointBody = " << vPointBody << std::endl;
    std::cout << "velCG = " << velCG << " vFlow = " << vFlow << std::endl;
    //std::cout << TransfoNED2B << std::endl; */

    a_e = atan2(-vFlow(3), vFlow(1));
    C_l = 2*3.141593*(AR/(AR+2))*a_e;
    U_inf = sqrt(vFlow(3)*vFlow(3) + vFlow(1)*vFlow(1));

    lift[i] = 0.5*rho*U_inf*U_inf*C_l*c*dw;
  }

  double rollMoment = 0.0;
  double totalLift = 0.0;

  for (int i = 0; i < 2*n+1; i++)
  {
    totalLift += lift[i]; //lift local de l'élément fois la surface de cet élément.
    if (i<n)
    {
      rollMoment += lift[i]*(n-i)*dw;
    } else {
      rollMoment -= lift[i]*(i-n)*dw; //lift partie droite de l'aile contribue négativement au rolling moment mais il y a un moins dans la formule.
    }
  }

  boxMoment(1) = rollMoment;
  boxMoment(2) = 0.0;
  boxMoment(3) = 0.0;

  liftForce(1) = 0.0;
  liftForce(2) = 0.0;
  liftForce(3) = totalLift;

  /* std::cout << "-----------------------------------------------------------------------------" << std::endl;
  //std::cout << lift[0] << " " << lift[1] << " " << lift[2] << " " << lift[3] << " " << lift[4] << std::endl;
  std::cout << "yaw = " << yaw << " roll = " << roll << " pitch = " << pitch << std::endl;
  std::cout << "velCG = " << velCG << " vFlow = " << vFlow << std::endl; */
  //std::cout << "Rolling moment = " << rollMoment << " Nm" << std::endl;
}

FGColumnVector3 FGAuxiliary::resultMoment() {
  return boxMoment;
}

FGColumnVector3 FGAuxiliary::getCGWinds(){
  return velCGBox;
}

void FGAuxiliary::goTo(double x2, double y2, double x1, double y1) {
  double psi1 = FDMExec->GetPropagate()->GetEuler(3); //Yaw de l'avion
  double psi2; //Angle entre axe nord et droite vers point visé

  double errorPsi;

  if (x1 <= x2)
  {
    if (y1 <= y2)
    {
      psi2 = atan2(x2-x1, y2-y1);
      if (psi1 >= PI+psi2 && psi1 < 2*PI)
      {
        errorPsi = -((2*PI - psi1) + psi2); //négatif car on veut qu'il aille à droite
      } else 
      {
        errorPsi = psi1 - psi2;
      }
      
    } 
    else 
    {
      psi2 = PI - atan2(x2-x1, y1-y2);
      if (psi1 >= PI+psi2 && psi1 < 2*PI)
      {
        errorPsi = -((2*PI - psi1) + psi2);
      } 
      else 
      {
        errorPsi = psi1 - psi2;
      }
    }
  } 
  else 
  {
    if (y1 <= y2)
    {
      psi2 = 2*PI - atan2(x1-x2, y2-y1);
      if (psi1 >= 0 && psi1 < psi2-PI)
      {
        errorPsi = psi1 + (2*PI - psi2);
      }
      else
      {
        errorPsi = psi1 - psi2;
      }
    } else 
    {
      psi2 = PI + atan2(x1-x2, y1-y2);
      if (psi1 >= 0 && psi1 < psi2-PI)
      {
        errorPsi = psi1 + (2*PI - psi2);
      }
      else
      {
        errorPsi = psi1 - psi2;
      }
    }
  }

  errorInt += errorPsi;
  
  double gainP = 0.4;
  double gainI = 0.0;
  double gainD = 30.0;

  double P = -gainP*errorPsi;
  double I = -gainI*errorInt;
  double D;
  if (prevError == 0.0)
  {
    D = 0.0;
  } else 
  {
    if (psi1 >= psi2)
    {
      D = gainD*(errorPsi - prevError); //ramene vers la droite (<0)
    } else {
      D = -gainD*(errorPsi - prevError); //ramene vers la gauche (>0)
    }
  }

  prevError = errorPsi;

  double GoTo = P + I + D;

  //////////////// ROLL LIMITER //////////////////

  double phi = FDMExec->GetPropagate()->GetEuler(1); //Roll de l'avion
  double maxPhi = PI/4;

  double gainRollP = -0.025;
  double gainRollD = -5.0;

  double P_Roll = gainRollP*phi;
  double D_Roll;

  if (prevError_Roll == 0)
  {
    D_Roll = 0;
  } else {
    D_Roll = gainRollD*phi;
  }

  prevError_Roll = phi;

  double rollLimiter = P_Roll + D_Roll; 
  
  double ailerons = GoTo + rollLimiter; //Final command for the ailerons

  FDMExec->GetFCS()->SetDaCmd(ailerons);

  double distance = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

  if (distance <= 99.9)
  {
    FDMExec->GetFCS()->SetDeCmd(-0.15);
  }
  

  double a_R = FDMExec->GetFCS()->GetDaRPos();
  double T = FDMExec->GetSimTime();
  double rolleee = FDMExec->GetPropagate()->GetEuler(1);
  ajouterDonnees("Zzz_aR", a_R);
  ajouterDonnees("Zzz_Time", T);
  ajouterDonnees("Zzz_Roll", rolleee);

  
  //std::cout << "-----------------------------------------------------------------------------" << std::endl;
  //std::cout << x1 << " " << psi1 << " " << psi2 << " " << " delta Yaw = " << errorPsi << " " << ailerons << std::endl;
  //std::cout << T << " " << errorInt << " prev: " << prevError << std::endl;
  //std::cout << "P: " << P << " I: " << I << " D: " << D << std::endl;
}

void FGAuxiliary::initialiserFichier(const std::string& nomFichier) {
    // Ouvrir un fichier en mode création
    std::ofstream fichier(nomFichier);

    // Vérifier si le fichier est ouvert avec succès
    if (fichier.is_open()) {
        std::cout << "Fichier cree avec succes." << std::endl;

        // Fermer le fichier
        fichier.close();
    } else {
        std::cerr << "Impossible de creer le fichier." << std::endl;
    }
}

void FGAuxiliary::ajouterDonnees(const std::string& nomFichier,double valeur) {
    // Ouvrir le fichier en mode ajout
    const std::string chemin = "/Users/Simon/Documents/Aaa_Thesis/git_jsbsim/jsbsim/";

    std::ofstream fichier(chemin + nomFichier + ".txt", std::ios::app);

    // Vérifier si le fichier est ouvert avec succès
    if (fichier.is_open()) {
        // Ajouter des données au fichier
        fichier  << valeur << std::endl;

        // Fermer le fichier
        fichier.close();
        //std::cout << "Donnees ajoutees avec succes." << std::endl;
    } else {
        std::cerr << "Impossible d'ouvrir le fichier pour ajout." << std::endl;
    }
}

void FGAuxiliary::autopilot(double x_2, double y_2, double x_1, double y_1){
  double time = FDMExec->GetSimTime();
  double dist = sqrt((x_2-x_1)*(x_2-x_1) + (y_2-y_1)*(y_2-y_1));
  double rollInst = boxMoment(1);
  double altInst;
  double timeInst;

  double rollLimit = 120.0;
  double updraftLimit = 1.5;
  double deltaTime = 4.0;
  double turn_angle = 6.0; //degres de roll

  FGColumnVector3 boxWind = getCGWinds();
  double updraft = boxWind(3) * -1; // - car en NED le updraft est négatif. Je le remets positif pour plus de clareté et ne pas se tromper.

  std::cout << time << std::endl;
  
  if (turn == 1) //si le virage est en cours
  {
    if (direction == 1)
    {
      FDMExec->GetAircraft()->virage(50.0, 150.0, turn_angle);
      direction = 1;
      turn = 1;
    } else if (direction == -1) {
      FDMExec->GetAircraft()->virage(50.0, 150.0, -turn_angle);
      direction = -1;
      turn = 1;
    }
    if ((time-timeInit) >= 300)
    {
      goTo(x_2, y_2, x_1, y_1); //changement de target
      turn = 0; //on indique que le virage est désamorcé
      direction = 0;
      waitTime = time;
      std::cout << "OUT OF TURN" << std::endl;
    }
  }
  else { // turn == 0
    if (time - waitTime < 200.0 && time >= 200.0)
    {
      //std::cout << time - waitTime << " " << x_1 << " " << y_1 << std::endl;
      turn = 0;
      goTo(x_2, y_2, x_1, y_1);
    } 
    else //on a attendu
    {
      //std::cout << "I HAVE WAITED" << std::endl;
      if ((abs(rollInst) > rollLimit || updraft > 0.5) && trigger == 0 && ok == 0)
      {
        trigger = 1;
        triggerTime = time;
        turn = 0;
        goTo(x_2, y_2, x_1, y_1); //on continue vers la target
      } 
      else if ((abs(rollInst) > rollLimit || updraft > 0.5) && trigger == 1 && (time - triggerTime) < deltaTime && ok == 0){
        trigger = 1;
        std::cout << time - triggerTime << std::endl;
        rollTest += rollInst; //permet de confirmer que le roll va dans un sens ou l'autre et pas juste en un point qui peut fausser le jugement.
        goTo(x_2, y_2, x_1, y_1);
      } 
      else if (((abs(rollInst) > rollLimit || updraft > 0.5) && trigger == 1 && (time - triggerTime) >= deltaTime) || ok == 1) { // on va se rapprocher du centre de la plume puis tourner
        ok = 1;
        if (rollTest < 0.0) //on dit qu'on va à droite
        {
          //std::cout << newTarget << std::endl;
          if (updraft < updraftLimit && newTarget == 0)
          {
            direction = 1;
            double yawInst = FDMExec->GetPropagate()->GetEuler(3);
            std::cout << x_1 << " " << y_1 << std::endl;
            goToCenter(x_1, y_1, yawInst, direction, 0);
            goTo(newEastTarget, newNorthTarget, x_1, y_1);
            newTarget = 1;
            //std::cout << newEastTarget << " " << newNorthTarget << std::endl;
          } 
          else if (updraft < updraftLimit && newTarget == 1 && sqrt(pow(newEastTarget - x_1,2) + pow(newNorthTarget - y_1,2)) > 200.0) 
          {
            direction = 1;
            newTarget = 1;
            goTo(newEastTarget, newNorthTarget, x_1, y_1);
            //std::cout << newEastTarget << " " << newNorthTarget << std::endl;
          } 
          else //on entame le virage
          { 
            FDMExec->GetAircraft()->virage(50.0, 150.0, turn_angle);
            direction = 1;
            turn = 1;
            trigger = 0;
            rollTest = 0.0;
            altInit = Propagate->GetAltitudeASL()*0.3048; //on set la hauteur initale du virage
            timeInit = time; //on set le time inital du virage. On jouera avec un des deux dans futures conditions.
            newTarget = 0;
            std::cout << "IN RIGHT TURN AT " << x_1 << " " << y_1 << std::endl;
            goToCenter(x_1, y_1, 0.0, direction, 1);
            ok = 0;
          }
        } 
        else 
        { // à gauche
          if (updraft < updraftLimit && newTarget == 0)
          {
            direction = -1;
            double yawInst = FDMExec->GetPropagate()->GetEuler(3);
            std::cout << x_1 << " " << y_1 << std::endl;
            goToCenter(x_1, y_1, yawInst, direction, 0);
            goTo(newEastTarget, newNorthTarget, x_1, y_1);
            newTarget = 1;
            //std::cout << newEastTarget << " " << newNorthTarget << std::endl;
          } 
          else if (updraft < updraftLimit && newTarget == 1 && sqrt(pow(newEastTarget - x_1,2) + pow(newNorthTarget - y_1,2)) > 200.0) 
          {
            direction = -1;
            newTarget = 1;
            goTo(newEastTarget, newNorthTarget, x_1, y_1);
            //std::cout << newEastTarget << " " << newNorthTarget << std::endl;
          } 
          else //on entame le virage
          { 
            FDMExec->GetAircraft()->virage(50.0, 150.0, turn_angle);
            direction = -1;
            turn = 1;
            trigger = 0;
            rollTest = 0.0;
            altInit = Propagate->GetAltitudeASL()*0.3048; //on set la hauteur initale du virage
            timeInit = time; //on set le time inital du virage. On jouera avec un des deux dans futures conditions.
            newTarget = 0;
            std::cout << "IN LEFT TURN AT " << x_1 << " " << y_1 << std::endl;
            goToCenter(x_1, y_1, 0.0, direction, 1);
            ok = 0;
          }
        }
      } else {
        goTo(x_2, y_2, x_1, y_1); //on continue vers la target
        turn = 0;
        direction = 0;
        trigger = 0;
      }
    }
  }
  /* if (x_1 >= 5300 && x_1 <= 5700)
  {
    std::cout << "East: " << x_1 << " North : " << y_1 << std::endl;
    std::cout << "time : " << time << " updraft : " << updraft << " RM = " << rollInst << std::endl;
    std::cout << "TEST, turn = " << turn << " direction = " << direction << " trigger : " << trigger << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;
  } */
}

void FGAuxiliary::goToCenter(double x, double y, double yaw, int direction, int reset){
  if (yaw >= 0.0 && yaw < PI/2 && reset == 0)
  {
    if (direction == 1)
    {
      newEastTarget  = x + 500.0 * cos(yaw);
      newNorthTarget = y - 500.0 * sin(yaw);
    } else {
      newEastTarget  = x - 500.0 * cos(yaw);
      newNorthTarget = y + 500.0 * sin(yaw);
    }
  } else if (yaw >= PI/2 && yaw < PI && reset == 0) 
  {
    if (direction == 1)
    {
      newEastTarget  = x - 500.0 * cos(PI - yaw);
      newNorthTarget = y - 500.0 * sin(PI - yaw);
    } else {
      newEastTarget  = x + 500.0 * cos(PI - yaw);
      newNorthTarget = y + 500.0 * sin(PI - yaw);
    }
  } else if (yaw >= PI && yaw < 3*PI/2 && reset == 0) 
  {
    if (direction == 1)
    {
      newEastTarget  = x - 500.0 * cos(yaw - PI);
      newNorthTarget = y + 500.0 * sin(yaw - PI);
    } else {
      newEastTarget  = x + 500.0 * cos(yaw - PI);
      newNorthTarget = y - 500.0 * sin(yaw - PI);
    }
  } else if (yaw >= 3*PI/2 && yaw < 2*PI && reset == 0) 
  {
    if (direction == 1)
    {
      newEastTarget  = x + 500.0 * cos(2*PI - yaw);
      newNorthTarget = y + 500.0 * sin(2*PI - yaw);
    } else {
      newEastTarget  = x - 500.0 * cos(2*PI - yaw);
      newNorthTarget = y - 500.0 * sin(2*PI - yaw);
    }
  } else if (reset == 1) {
      newEastTarget  = 0.0;
      newNorthTarget = 0.0;
  }
  std::cout << "YES" << std::endl;
  std::cout << newEastTarget << " " << newNorthTarget << std::endl;
  //newCoord = {newEastTarget, newNorthTarget};
  //return newTarget;
}

} // namespace JSBSim
