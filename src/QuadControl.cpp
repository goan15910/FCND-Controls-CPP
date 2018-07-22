#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath> 

#define G_CONST 9.8f

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  float arm = L / 1.414;

  //float c_bar = (mass * 9.81f) / kappa;
  float c_bar = collThrustCmd;
  float p_bar = momentCmd.x / arm;
  float q_bar = momentCmd.y / arm;
  float r_bar = -momentCmd.z / kappa; // negative sign for reactive moment

  // the four angular speed square
  float Thrust0 = (c_bar + p_bar + q_bar + r_bar) / 4;
  float Thrust1 = (c_bar - p_bar + q_bar - r_bar) / 4;
  float Thrust2 = (c_bar + p_bar - q_bar - r_bar) / 4;
  float Thrust3 = (c_bar - p_bar - q_bar + r_bar) / 4;

  /**/
  cmd.desiredThrustsN[0] = CONSTRAIN(Thrust0, minMotorThrust, maxMotorThrust); // front left
  cmd.desiredThrustsN[1] = CONSTRAIN(Thrust1, minMotorThrust, maxMotorThrust); // front right
  cmd.desiredThrustsN[2] = CONSTRAIN(Thrust2, minMotorThrust, maxMotorThrust); // rear left
  cmd.desiredThrustsN[3] = CONSTRAIN(Thrust3, minMotorThrust, maxMotorThrust); // rear right
  /**/

  /*
  cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
  cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
  cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
  cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right
  */

  //printf("collThrust: %f\n", collThrustCmd);
  //printf("momentCmd: (%f, %f, %f)\n", momentCmd.x, momentCmd.y, momentCmd.z);
  //printf("(c_bar, p_bar, q_bar, r_bar): (%f, %f, %f, %f)\n", c_bar, p_bar, q_bar, r_bar);
  //printf("motor: (%f, %f, %f, %f)\n", Thrust0, Thrust1, Thrust2, Thrust3);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  V3F rateError = pqrCmd - pqr;

  V3F rawCmd;
  rawCmd.x = Ixx * kpPQR.x * rateError.x;
  rawCmd.y = Iyy * kpPQR.y * rateError.y;
  rawCmd.z = Izz * kpPQR.z * rateError.z;

  momentCmd = rawCmd;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  //float c_d = collThrustCmd / (mass + 1e-7);

  if(collThrustCmd > 0.0) {
    V3F targetR = accelCmd * mass / (collThrustCmd + 1e-7);
    float targetR13 = targetR[0]; // cos(pitch) * cos(yaw)
    float targetR23 = targetR[1]; // -cos(pitch) * sin(yaw)
    float minCosSqr = pow(cosf(maxTiltAngle), 2.);
    float maxSin = sinf(maxTiltAngle);
    targetR13 = CONSTRAIN(targetR13, minCosSqr, 1.);
    targetR23 = CONSTRAIN(targetR23, -maxSin, maxSin);

    float R13Error = targetR13 - R(0, 2);
    float R23Error = targetR23 - R(1, 2);

    float b_c_x_dot = kpBank * R13Error;
    float b_c_y_dot = kpBank * R23Error;

    pqrCmd.x = (R(1, 0) * b_c_x_dot - R(0, 0) * b_c_y_dot) / (R(2, 2) + 1e-7);
    pqrCmd.y = (R(1, 1) * b_c_x_dot - R(0, 1) * b_c_y_dot) / (R(2, 2) + 1e-7);

    printf("accelCmd: (%f, %f)\n", \
            accelCmd.x, \
            accelCmd.y);

    printf("Target (R13, R23): (%f, %f)\n", \
            targetR13, \
            targetR23);

    //printf("Error (R13, R23): (%f, %f)\n", \
           R13Error, \
           R23Error);
  }
  else {
      printf("negative thrust command");
  }

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}


float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, \
                                   Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  /*
  float zDotCmd = kpPosZ * (posZCmd - posZ) + velZCmd;
  zDotCmd = CONSTRAIN(zDotCmd, -maxDescentRate, maxAscentRate);
  accelZCmd += kpVelZ * (zDotCmd - velZ);
  */

  accelZCmd += kpVelZ * (velZCmd - velZ) + kpPosZ * (posZCmd - posZ);
  thrust = mass * (accelZCmd - G_CONST) / (R(2, 2) + 1e-7);
  thrust = abs(thrust);
  
  //printf("accelZ: %f\n", accelZCmd);
  printf("thrust: %f\n", thrust);

  /////////////////////////////// END STUDENT CODE ////////////////////////////
  
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  float velMag = velCmd.mag();
  if(velMag > maxSpeedXY) {
    velCmd = velCmd.norm() * maxSpeedXY;
  }

  V3F posError = posCmd - pos;
  V3F velError = velCmd - vel;

  accelCmd += kpPosXY * posError + kpVelXY * velError;
  accelCmd.z = 0; // make sure nothing done to accelCmd.z

  printf("PosError: %f\n", posError.x);
  //printf("Pos: %f\n", pos.x);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b]. 
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  yaw = ValidateYaw(yaw);
  float yawError = yawCmd - yaw;
  yawError = ValidateYaw(yawError);

  // printf("yawError: %d\n", int(yawError * M_1_PI * 180));

  yawRateCmd = kpYaw * yawError;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, \
                                        estPos.z, estVel.z, estAtt, \
                                        curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, \
                                      estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  printf("\n============================\n");
  //printf("curTrajPoint (position): (%f, %f, %f)\n", curTrajPoint.position.x, curTrajPoint.position.y, curTrajPoint.position.z);
  //printf("curTrajPoint (velocity): (%f, %f, %f)\n", curTrajPoint.velocity.x, curTrajPoint.velocity.y, curTrajPoint.velocity.z);
  //printf("curTrajPoint (acceleration): (%f, %f, %f)\n", curTrajPoint.accel.x, curTrajPoint.accel.y, curTrajPoint.accel.z);
  //printf("estPos: (%f, %f, %f)\n", estPos.x, estPos.y, estPos.z);
  //printf("estVel: (%f, %f, %f)\n", estVel.x, estVel.y, estVel.z);
  //printf("estAtt: (%f, %f, %f)\n", estAtt[0], estAtt[1], estAtt[2]);
  //printf("diffPos: (%f, %f, %f)\n", curTrajPoint.position.x - estPos.x, curTrajPoint.position.y - estPos.y, curTrajPoint.position.z - estPos.z);
  //printf("diffVel: (%f, %f, %f)\n", curTrajPoint.velocity.x - estVel.x, curTrajPoint.velocity.y - estVel.y, curTrajPoint.velocity.z - estVel.z);
  //printf("collThrustCmd: %f\n", collThrustCmd);
  //printf("desAcc: (%f, %f, %f)\n", desAcc.x, desAcc.y, desAcc.z);
  //printf("desOmega: (%f, %f, %f)\n", desOmega.x, desOmega.y, desOmega.z);
  //printf("desMoment: (%f, %f, %f)\n", desMoment.x, desMoment.y, desMoment.z);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}

V3F QuadControl::SaturateCommand(V3F cmd, V3F maxValue)
{
  V3F validCmd;
  for(int i = 0; i < 3; i++) {
    validCmd[i] = (cmd[i] > maxValue[i]) ? cmd[i] : cmd.norm()[i] * maxValue[i];
  }

  return validCmd;
}

float QuadControl::ValidateYaw(float yaw)
{
  // let yaw be in [-pi, pi]
  yaw = fmodf(yaw, 2 * M_PI);
  if (yaw > M_PI) {
    yaw -= 2 * M_PI;
  }
  else if (yaw < -M_PI) {
    yaw += 2 * M_PI;
  }

  return yaw;
}
