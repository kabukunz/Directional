// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2020 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLYVECTOR_TO_RAW_H
#define DIRECTIONAL_POLYVECTOR_TO_RAW_H

#include <iostream>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/PI.h>
#include <igl/local_basis.h>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>



namespace directional
  {
  
  IGL_INLINE void multiply_polynomials(Eigen::RowVectorXcd& p1,
                                       const Eigen::RowVectorXcd& p2)
  {
    Eigen::RowVectorXcd newp = Eigen::RowVectorXcd::Zero(p1.size()+p2.size());
    for (int i=0;i<p1.size();i++)
      for (int j=0;j<p2.size();j++)
        newp(i+j)+=p1(i)*p2(j);
    
    p1 =newp;
  }
  
  IGL_INLINE void multiply_polynomials(Eigen::MatrixXcd& p1,
                                       const Eigen::MatrixXcd& p2)
  {
    assert(p1.rows()==p2.rows());
    Eigen::MatrixXcd newp = Eigen::MatrixXcd::Zero(p1.rows(),p1.cols()+p2.cols());
    for (int i=0;i<p1.cols();i++)
      for (int j=0;j<p2.cols();j++)
        newp.col(i+j).array()+=p1.col(i).array()*p2.col(j).array();
    
    p1 =newp;
  }
  
  IGL_INLINE void polynomial_eval_from_roots(const Eigen::RowVectorXcd& roots,
                                             const std::complex<double>& evalPoint,
                                             std::complex<double>& polyValue){
    polyValue=std::complex<double>(1.0,0.0);
    for (int i=0;i<roots.size();i++)
      polyValue*=(evalPoint-roots(i));
  }
  
  IGL_INLINE void polynomial_eval_from_roots(const Eigen::MatrixXcd& roots,
                                             const Eigen::VectorXcd& evalPoints,
                                             Eigen::VectorXcd& polyValues){
    assert(evalPoints.rows()==roots.rows());
    polyValues=Eigen::VectorXcd::Constant(evalPoints.rows(),std::complex<double>(1.0,0.0));
    for (int i=0;i<roots.cols();i++)
      polyValues.array()*=(evalPoints-roots.col(i)).array();  //rowwise-product faster?;
    
  }
  
  IGL_INLINE void polynomial_eval(const Eigen::RowVectorXcd& coeffs,
                                  const std::complex<double>& evalPoint,
                                  std::complex<double>& polyValue){
    polyValue=std::complex<double>(0.0,0.0);
    std::complex<double> z=1;
    for (int i=0;i<coeffs.size();i++){
      polyValue+=coeffs(i)*z;
      z*=evalPoint;
    }
    std::cout<<"z: "<<z<<std::endl;
    polyValue+=z;  //the biggest and unit power
  }
  
  IGL_INLINE void polynomial_eval(const Eigen::MatrixXcd& coeffs,
                                  const Eigen::VectorXcd& evalPoints,
                                  Eigen::VectorXcd& polyValues){
    polyValues=Eigen::VectorXcd::Zero(coeffs.rows());
    Eigen::VectorXcd z=Eigen::VectorXcd::Constant(coeffs.rows(),1.0);
    for (int i=0;i<coeffs.cols();i++){
      polyValues.array()+=coeffs.col(i).array()*z.array();
      z.array()*=evalPoints.array();
    }
    //std::cout<<"z: "<<z<<std::endl;
    polyValues.array()+=z.array();  //the biggest and unit power
  }
  
  // Converts a field in PolyVector representation to raw represenation.
  // Inputs:
  //  B1, B2:           #F by 3 matrices representing the local base of each face.
  //  polyVectorField:  #F by N complex PolyVectors
  //  N:                The degree of the field.
  // Outputs:
  //  raw:              #F by 3*N matrix with all N explicit vectors of each directional in raw format xyzxyz
  //    returns true if succeeded
  IGL_INLINE bool polyvector_to_raw(const Eigen::MatrixXd& B1,
                                    const Eigen::MatrixXd& B2,
                                    const Eigen::MatrixXcd& pvField,
                                    const int N,
                                    Eigen::MatrixXd& rawField,
                                    bool signSymmetry=true,
                                    const double rootTolerance=1e-8)
  {
    using namespace std;
    using namespace Eigen;
    rawField.resize(B1.rows(), 3 * N);
    if (N%2!=0) signSymmetry=false;  //by definition
    MatrixXcd actualPVField;
    int actualN;
    if (signSymmetry){
      actualPVField.resize(pvField.rows(),N/2);
      for (int i=0;i<N;i+=2)
        actualPVField.col(i/2)=pvField.col(i);
      actualN=N/2;
    } else{
      actualPVField=pvField;
      actualN=N;
    }
    
    MatrixXcd roots(actualPVField.rows(),actualN);
    roots.col(0).array()=(-actualPVField.col(0).array()).pow(1.0/(double)actualN);
    for (int i=1;i<actualN;i++)
      roots.col(i).array()=roots.col(i-1).array()*std::exp(std::complex<double>(0,2.0*igl::PI/(double)actualN));
    
    MatrixXd rootError=MatrixXd::Constant(actualPVField.rows(),actualN, 1000.0);
    double maxError=1000.0;
    int currRoot=0;
    MatrixXcd mostRoots(actualPVField.rows(),actualN-1);
    int maxIterations=1000;
    int currIteration=0;
    do{
      mostRoots.block(0,0, mostRoots.rows(),currRoot)=roots.block(0,0, roots.rows(),currRoot);
      mostRoots.block(0,currRoot, mostRoots.rows(),actualN-currRoot-1)=roots.block(0,currRoot+1, roots.rows(),actualN-currRoot-1);
      VectorXcd numerator; polynomial_eval(actualPVField, roots.col(currRoot), numerator);
      VectorXcd denominator; polynomial_eval_from_roots(mostRoots, roots.col(currRoot), denominator);
      roots.col(currRoot).array()-=numerator.array()/denominator.array();
      rootError.col(currRoot)=numerator.cwiseAbs();
      maxError=rootError.lpNorm<Infinity>();  //optimize by only the relevant column?
      currRoot=(currRoot+1)%actualN;
      currIteration++;
    }while ((maxError>rootTolerance)&&(currIteration<maxIterations));
    
    if (currIteration>=maxIterations)
      return false;
    
    if (signSymmetry)
      roots=roots.cwiseSqrt();
  
    for (int f=0;f<roots.rows();f++){
      RowVectorXcd rowRoots=roots.row(f);
      std::sort(rowRoots.data(), rowRoots.data() + rowRoots.size(), [](std::complex<double> a, std::complex<double> b){return arg(a) < arg(b);});
      roots.row(f)=rowRoots;
    }
    
    if (signSymmetry){
      MatrixXcd actualRoots(roots.rows(), 2*roots.cols());
      actualRoots<<roots, -roots;
      roots=actualRoots;
    }
   
    for (int i = 0; i < N; i++)
      for (int f=0;f<pvField.rows();f++){
        std::complex<double> root = roots(f,i);
        rawField.block<1, 3>(f, 3 * i) = B1.row(f) * root.real() + B2.row(f) * root.imag();
      }
      
    return true;
  }
  
  }

#endif
