/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Header:       FGAuxiliary.h
 Author:       Jon Berndt
 Date started: 01/26/99

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

HISTORY
--------------------------------------------------------------------------------
11/22/98   JSB   Created
  1/1/00   TP    Added calcs and getters for VTAS, VCAS, VEAS, Vground, in knots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SENTRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FGAUXILIARY_H
#define FGAUXILIARY_H

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INCLUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "FGModel.h"
#include "math/FGLocation.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FORWARD DECLARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

namespace JSBSim {

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS DOCUMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/** Encapsulates various uncategorized scheduled functions.
    Pilot sensed accelerations are calculated here. This is used
    for the coordinated turn ball instrument. Motion base platforms sometimes
    use the derivative of pilot sensed accelerations as the driving parameter,
    rather than straight accelerations.

    The theory behind pilot-sensed calculations is presented:

    For purposes of discussion and calculation, assume for a minute that the
    pilot is in space and motionless in inertial space. She will feel
    no accelerations. If the aircraft begins to accelerate along any axis or
    axes (without rotating), the pilot will sense those accelerations. If
    any rotational moment is applied, the pilot will sense an acceleration
    due to that motion in the amount:

    [wdot X R]  +  [w X (w X R)]
    Term I          Term II

    where:

    wdot = omegadot, the rotational acceleration rate vector
    w    = omega, the rotational rate vector
    R    = the vector from the aircraft CG to the pilot eyepoint

    The sum total of these two terms plus the acceleration of the aircraft
    body axis gives the acceleration the pilot senses in inertial space.
    In the presence of a large body such as a planet, a gravity field also
    provides an accelerating attraction. This acceleration can be transformed
    from the reference frame of the planet so as to be expressed in the frame
    of reference of the aircraft. This gravity field accelerating attraction
    is felt by the pilot as a force on her tushie as she sits in her aircraft
    on the runway awaiting takeoff clearance.

    In JSBSim the acceleration of the body frame in inertial space is given
    by the F = ma relation. If the vForces vector is divided by the aircraft
    mass, the acceleration vector is calculated. The term wdot is equivalent
    to the JSBSim vPQRdot vector, and the w parameter is equivalent to vPQR.

    @author Tony Peden, Jon Berndt
*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLASS DECLARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

// FONCTION RAJOUTEES
class FGPropagate;

class JSBSIM_API FGAuxiliary : public FGModel {
public:
  /** Constructor
      @param Executive a pointer to the parent executive object */
  explicit FGAuxiliary(FGFDMExec* Executive);

  /// Destructor
  ~FGAuxiliary();

  bool InitModel(void) override;

  /** Runs the Auxiliary routines; called by the Executive
      Can pass in a value indicating if the executive is directing the
      simulation to Hold.
      @param Holding if true, the executive has been directed to hold the sim
                     from advancing time. Some models may ignore this flag, such
                     as the Input model, which may need to be active to listen
                     on a socket for the "Resume" command to be given.  @return
                     false if no error */
  bool Run(bool Holding) override;

// GET functions

  /** Compute the total pressure in front of the Pitot tube. It uses the
  *   Rayleigh formula for supersonic speeds (See "Introduction to Aerodynamics
  *   of a Compressible Fluid - H.W. Liepmann, A.E. Puckett - Wiley & sons
  *   (1947)" §5.4 pp 75-80)
  *   @param mach      The Mach number
  *   @param pressure  Pressure in psf
  *   @return The total pressure in front of the Pitot tube in psf */
  double PitotTotalPressure(double mach, double pressure) const;

  /** Compute the Mach number from the differential pressure (qc) and the
  *   static pressure. Based on the formulas in the US Air Force Aircraft
  *   Performance Flight Testing Manual (AFFTC-TIH-99-01).
  *   @param qc        The differential/impact pressure
  *   @param pressure  Pressure in psf
  *   @return The Mach number */
  double MachFromImpactPressure(double qc, double p) const;

  /** Calculate the calibrated airspeed from the Mach number. Based on the
  *   formulas in the US Air Force Aircraft Performance Flight Testing
  *   Manual (AFFTC-TIH-99-01).
  *   @param mach      The Mach number
  *   @param pressure  Pressure in psf
  *   @return The calibrated airspeed (CAS) in ft/s
  * */
  double VcalibratedFromMach(double mach, double pressure) const;

  /** Calculate the Mach number from the calibrated airspeed.Based on the
  *   formulas in the US Air Force Aircraft Performance Flight Testing
  *   Manual (AFFTC-TIH-99-01).
  *   @param vcas      The calibrated airspeed (CAS) in ft/s
  *   @param pressure  Pressure in psf
  *   @return The Mach number
  * */
  double MachFromVcalibrated(double vcas, double pressure) const;

  // Atmospheric parameters GET functions
  /** Returns Calibrated airspeed in feet/second.*/
  double GetVcalibratedFPS(void) const { return vcas; }
  /** Returns Calibrated airspeed in knots.*/
  double GetVcalibratedKTS(void) const { return vcas*fpstokts; }
  /** Returns equivalent airspeed in feet/second. */
  double GetVequivalentFPS(void) const { return veas; }
  /** Returns equivalent airspeed in knots. */
  double GetVequivalentKTS(void) const { return veas*fpstokts; }
  /** Returns the true airspeed in feet per second. */
  double GetVtrueFPS() const { return Vt; }
  /** Returns the true airspeed in knots. */
  double GetVtrueKTS() const { return Vt * fpstokts; }

  /** Returns the total pressure.
      Total pressure is freestream total pressure for
      subsonic only. For supersonic it is the 1D total pressure
      behind a normal shock. */
  double GetTotalPressure(void) const { return pt; }

  /** Returns the total temperature.
    The total temperature ("tat", isentropic flow) is calculated:
    @code
    tat = in.Temperature*(1 + 0.2*Mach*Mach)
    @endcode
    (where "in.Temperature" is standard temperature calculated by the atmosphere
    model) */

  double GetTotalTemperature(void) const { return tat; }
  double GetTAT_C(void) const { return tatc; }

  double GetPilotAccel(int idx)  const { return vPilotAccel(idx);  }
  double GetNpilot(int idx)      const { return vPilotAccelN(idx); }
  double GetAeroPQR(int axis)    const { return vAeroPQR(axis);    }
  double GetEulerRates(int axis) const { return vEulerRates(axis); }

  const FGColumnVector3& GetPilotAccel (void) const { return vPilotAccel;  }
  const FGColumnVector3& GetNpilot     (void) const { return vPilotAccelN; }
  const FGColumnVector3& GetNcg        (void) const { return vNcg;         }
  double GetNcg                     (int idx) const { return vNcg(idx);    }
  double GetNlf                        (void) const;
  const FGColumnVector3& GetAeroPQR    (void) const { return vAeroPQR;     }
  const FGColumnVector3& GetEulerRates (void) const { return vEulerRates;  }
  const FGColumnVector3& GetAeroUVW    (void) const { return vAeroUVW;     }
  const FGLocation&      GetLocationVRP(void) const { return vLocationVRP; }

  double GetAeroUVW (int idx) const { return vAeroUVW(idx); }
  double Getalpha   (void) const { return alpha;      }
  double Getbeta    (void) const { return beta;       }
  double Getadot    (void) const { return adot;       }
  double Getbdot    (void) const { return bdot;       }
  double GetMagBeta (void) const { return fabs(beta); }

  double Getalpha   (int unit) const { if (unit == inDegrees) return alpha*radtodeg;
                                       else return BadUnits(); }
  double Getbeta    (int unit) const { if (unit == inDegrees) return beta*radtodeg;
                                       else return BadUnits(); }
  double Getadot    (int unit) const { if (unit == inDegrees) return adot*radtodeg;
                                       else return BadUnits(); }
  double Getbdot    (int unit) const { if (unit == inDegrees) return bdot*radtodeg;
                                       else return BadUnits(); }
  double GetMagBeta (int unit) const { if (unit == inDegrees) return fabs(beta)*radtodeg;
                                       else return BadUnits(); }

  /** Calculates and returns the wind-to-body axis transformation matrix.
      @return a reference to the wind-to-body transformation matrix.
      */
  const FGMatrix33& GetTw2b(void) const { return mTw2b; }

  /** Calculates and returns the body-to-wind axis transformation matrix.
      @return a reference to the wind-to-body transformation matrix.
      */
  const FGMatrix33& GetTb2w(void) const { return mTb2w; }

  double Getqbar          (void) const { return qbar;       }
  double GetqbarUW        (void) const { return qbarUW;     }
  double GetqbarUV        (void) const { return qbarUV;     }
  double GetReynoldsNumber(void) const { return Re;         }

  /** Gets the magnitude of total vehicle velocity including wind effects in
      feet per second. */
  double GetVt            (void) const { return Vt;         }

  /** Gets the ground speed in feet per second.
      The magnitude is the square root of the sum of the squares (RSS) of the 
      vehicle north and east velocity components.
      @return The magnitude of the vehicle velocity in the horizontal plane. */
  double GetVground       (void) const { return Vground;    }

  /** Gets the Mach number. */
  double GetMach          (void) const { return Mach;       }

  /** The mach number calculated using the vehicle X axis velocity. */
  double GetMachU         (void) const { return MachU;      }

  /** The longitudinal acceleration in g's of the aircraft center of gravity. */
  double GetNx            (void) const { return Nx;         }

  /** The lateral acceleration in g's of the aircraft center of gravity. */
  double GetNy            (void) const { return Ny;         }

  /** The vertical acceleration in g's of the aircraft center of gravity. */
  double GetNz            (void) const { return Nz;         }

  const FGColumnVector3& GetNwcg(void) const { return vNwcg; }

  double GetHOverBCG(void) const { return hoverbcg; }
  double GetHOverBMAC(void) const { return hoverbmac; }

  double GetGamma(void)              const { return gamma;         }
  double GetGroundTrack(void)        const { return psigt;         }

  double GetGamma(int unit) const {
    if (unit == inDegrees) return gamma*radtodeg;
    else return BadUnits();
  }

  double GetLongitudeRelativePosition (void) const;
  double GetLatitudeRelativePosition  (void) const;
  double GetDistanceRelativePosition  (void) const;

  void SetAeroPQR(const FGColumnVector3& tt) { vAeroPQR = tt; }

  FGColumnVector3 resultMoment();
  FGColumnVector3 getCGWinds();

  double newEastTarget;
  double newNorthTarget;

  struct Inputs {
    double Pressure;
    double Density;
    double Temperature;
    double StdDaySLsoundspeed;
    double SoundSpeed;
    double KinematicViscosity;
    double DistanceAGL;
    double Wingspan;
    double Wingchord;
    double StandardGravity;
    double Mass;
    FGMatrix33 Tl2b;
    FGMatrix33 Tb2l;
    FGColumnVector3 vPQR;
    FGColumnVector3 vPQRi;
    FGColumnVector3 vPQRidot;
    FGColumnVector3 vUVW;
    FGColumnVector3 vUVWdot;
    FGColumnVector3 vVel;
    FGColumnVector3 vBodyAccel;
    FGColumnVector3 ToEyePt;
    FGColumnVector3 RPBody;
    FGColumnVector3 VRPBody;
    FGColumnVector3 vFw;
    FGLocation vLocation;
    double CosTht;
    double SinTht;
    double CosPhi;
    double SinPhi;
    FGColumnVector3 TotalWindNED;
    FGColumnVector3 TurbPQR;
  } in;

private:
  double vcas, veas;
  double pt, tat, tatc; // Don't add a getter for pt!

  FGMatrix33 mTw2b;
  FGMatrix33 mTb2w;

  FGColumnVector3 vPilotAccel;
  FGColumnVector3 vPilotAccelN;
  FGColumnVector3 vNcg;
  FGColumnVector3 vNwcg;
  FGColumnVector3 vAeroPQR;
  FGColumnVector3 vAeroUVW;
  FGColumnVector3 vEulerRates;
  FGColumnVector3 vMachUVW;
  FGLocation vLocationVRP;

  double Vt, Vground;
  double Mach, MachU;
  double qbar, qbarUW, qbarUV;
  double Re; // Reynolds Number = V*c/mu
  double alpha, beta;
  double adot,bdot;
  double psigt, gamma;
  double Nx, Ny, Nz;

  double hoverbcg, hoverbmac;

  void UpdateWindMatrices(void);

  void CalculateRelativePosition(void);

  void bind(void);
  double BadUnits(void) const;
  void Debug(int from) override;


   // FONCTION RAJOUTEES
  FGPropagate* Propagate;


  void loadgrid();
  double grid[3][257];
  
  void loaduwind();
  double u[128][257][257];

  void loadvwind();
  double v[128][257][257];

  void loadwwind();
  double w[128][257][257];

  const int n = 10;

  double* rechercheNoeuds(double hauteur, double longueur, double largeur, double refz,double ref_long, double longi, double lat);
  void discretisation(double x, double y, double z, int n);
  void dynamics(int vBoite[5][3], int n);

  void getRollMoment(double hauteur, double longueur, double largeur, double longi, double lat, int n, double largeur_0, double longueur_0);

  void initialiserFichier(const std::string& nomFichier);
  void ajouterDonnees(const std::string& nomFichier,double valeur);
  void goTo(double x2, double y2, double x1, double y1);
  void autopilot(double x_1, double y_1, double x_2, double y_2);
  void goToCenter(double x, double y, double yaw, int direction, int reset);

  void sideslipController(double x, double y);


  //AUTOPILOT 2
  void autopilot2(double updraft, double time, double x, double y, double x_t, double y_t); //x_t and y_t target position
  void inThermal(double updraft, double time, double x, double y, double x_t, double y_t);
  void thermalCentering(double updraft, double time, double x, double y, double x_t, double y_t);
  void turning(double angle);
  void exitStrategy(double updraft, double time, double x, double y, double x_t, double y_t);
  double getPsiError(double x1, double y1, double x2, double y2);

  double prevUpdraft;

  double timeUp;
  double timeDown;
  int triggerDown;
  int triggerUp;

  double oldRollTurn = 0.0;
  double integralRoll = 0.0;

  int exit;
  double exitTime;

  double alt;

  int triggerRoll;
  double timeTriggerRoll;
  int closer;

  double altMax = 0.0;

  double East_target;
  double North_target;
  int n_target = 1;
  ///////////////

  double errorInt;
  double prevError;
  double prevError_Roll;
  int turn = 0;
  int direction = 0; //1->right -1->left
  double altInit = 0.0;
  double timeInit = 0.0;
  int trigger = 0; // double car on va soustraire avec time
  double triggerTime = 0.0;
  double triggerDirection = 0; //1->right -1->left
  double rollTest = 0.0;

  int newTarget = 0;
  //double newCoord[2] = {0.0, 0.0};

  double waitTime = 0.0;

  int ok = 0;

  //sideslip controller
  double oldEast = 0.0;
  double oldNorth = 0.0;
  double oldSideslip = 0.0;
};

} // namespace JSBSim

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif
