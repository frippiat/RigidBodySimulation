// ----------------------------------------------------------------------------
// RigidSolver.hpp
//
//  Created on: 18 Dec 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: Simple Rigid Body Solver (DO NOT DISTRIBUTE!)
//
// Copyright 2020-2024 Kiwon Um
// ----------------------------------------------------------------------------

#ifndef _RIGIDSOLVER_HPP_
#define _RIGIDSOLVER_HPP_

#include <glm/ext/matrix_transform.hpp>
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include <cmath>
#include <vector>
#include <iostream>

struct Quaternion {
  tReal w, x, y, z;

  Quaternion(tReal w = 1, tReal x = 0, tReal y = 0, tReal z = 0)
    : w(w), x(x), y(y), z(z) {}

  Quaternion operator*(const Quaternion& q) const {
    return Quaternion(
      w * q.w - x * q.x - y * q.y - z * q.z,
      w * q.x + x * q.w + y * q.z - z * q.y,
      w * q.y - x * q.z + y * q.w + z * q.x,
      w * q.z + x * q.y - y * q.x + z * q.w
    );
  }

  Quaternion normalized() const {
    tReal mag = std::sqrt(w * w + x * x + y * y + z * z);
    return Quaternion(w / mag, x / mag, y / mag, z / mag);
  }

  Mat3f toMatrix() const {
    tReal xx = x * x, yy = y * y, zz = z * z;
    tReal xy = x * y, xz = x * z, yz = y * z;
    tReal wx = w * x, wy = w * y, wz = w * z;

    return Mat3f(
      1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy),
      2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx),
      2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)
    );
  }
};

struct BodyAttributes {
  BodyAttributes() :
    X(0, 0, 0), R(Mat3f::I()), P(0, 0, 0), L(0, 0, 0),
    V(0, 0, 0), omega(0, 0, 0), F(0, 0, 0), tau(0, 0, 0), orientation(1, 0, 0, 0) {}

  glm::mat4 worldMat() const {
    return glm::mat4(           // column-major
      R(0,0), R(1,0), R(2,0), 0,
      R(0,1), R(1,1), R(2,1), 0,
      R(0,2), R(1,2), R(2,2), 0,
      X[0],   X[1],   X[2],   1);
  }

  tReal M;                      // mass
  Mat3f I0, I0inv;              // inertia tensor and its inverse in body space
  Mat3f Iinv;                   // inverse of inertia tensor

  // rigid body state
  Vec3f X;                      // position
  Mat3f R;                      // rotation
  Quaternion orientation;       // quaternion orientation
  Vec3f P;                      // linear momentum
  Vec3f L;                      // angular momentum

  // auxiliary quantities
  Vec3f V;                      // linear velocity
  Vec3f omega;                  // angular velocity

  // force and torque
  Vec3f F;                      // force
  Vec3f tau;                    // torque

  // mesh's vertices in body space
  std::vector<Vec3f> vdata0;
};

class Box : public BodyAttributes {
public:
  explicit Box(
    const tReal w=1.0, const tReal h=1.0, const tReal d=1.0, const tReal dens=10.0,
    const Vec3f v0=Vec3f(0, 0, 0), const Vec3f omega0=Vec3f(0, 0, 0)) :
    width(w), height(h), depth(d)
  {
    V = v0;                     // initial velocity
    omega = omega0;             // initial angular velocity

    M = dens * w * h * d;       // mass calculation
    I0 = Mat3f(                 // inertia tensor for box
      (1.0 / 12.0) * M * (h * h + d * d), 0, 0,
      0, (1.0 / 12.0) * M * (w * w + d * d), 0,
      0, 0, (1.0 / 12.0) * M * (w * w + h * h)
    );
    I0inv = I0.invert();
    Iinv = I0inv;

    // vertices data (8 vertices)
    vdata0.push_back(Vec3f(-0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h, -0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h, -0.5*d));

    vdata0.push_back(Vec3f(-0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h,  0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h,  0.5*d));
  }

  tReal width, height, depth;
};

class RigidSolver {
public:
  explicit RigidSolver(
    BodyAttributes *body0=nullptr, const Vec3f g=Vec3f(0, -0.98, 0)) :
    body(body0), _g(g), _step(0), _sim_t(0) {}

  void init(BodyAttributes *body0) {
    body = body0;
    _step = 0;
    _sim_t = 0;
  }

  void step(const tReal dt) {
    std::cout << "t=" << _sim_t << " (dt=" << dt << ")" << std::endl;

    computeForceAndTorque();

    // Linear momentum update
    body->P += body->F * dt;
    body->V = body->P / body->M;
    body->X += body->V * dt;

    // Angular momentum update
    body->L += body->tau * dt;
    body->omega = body->Iinv * body->L;

    // Update orientation using quaternion
    Quaternion dq(0, body->omega[0], body->omega[1], body->omega[2]);
    dq = dq * body->orientation * 0.5;
    body->orientation = (Quaternion(
  body->orientation.w + dq.w * dt,
  body->orientation.x + dq.x * dt,
  body->orientation.y + dq.y * dt,
  body->orientation.z + dq.z * dt
)).normalized();

    body->R = body->orientation.toMatrix();

    ++_step;
    _sim_t += dt;
  }

  BodyAttributes *body;

private:
  void computeForceAndTorque() {
    body->F = body->M * _g; // Gravity

    if (_step == 1) {
      Vec3f instantForce(0.15, 0.25, 0.03);
      body->F += instantForce;
      Vec3f r = body->vdata0[0] - body->X; // Relative position of 0th vertex
      body->tau = Vec3f(
  r[1] * instantForce[2] - r[2] * instantForce[1],
  r[2] * instantForce[0] - r[0] * instantForce[2],
  r[0] * instantForce[1] - r[1] * instantForce[0]
);

    }
  }

  Vec3f _g;                     // gravity
  tIndex _step;                 // step count
  tReal _sim_t;                 // simulation time
};

#endif  /* _RIGIDSOLVER_HPP_ */
