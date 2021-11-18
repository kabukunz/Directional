// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2017 Daniele Panozzo <daniele.panozzo@gmail.com>, Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLYVECTOR_FIELD_H
#define DIRECTIONAL_POLYVECTOR_FIELD_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Polynomials>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <igl/barycenter.h>
#include <igl/speye.h>
#include <igl/eigs.h>
#include <iostream>
#include <directional/circumcircle.h>

namespace directional
{

#define UNIFORM_WEIGHTS 0
#define BARYCENTRIC_WEIGHTS 1
#define INV_COT_WEIGHTS 2

  struct PolyVectorData{
  public:
    
    //User parameters
    Eigen::VectorXi constFaces;   //list of faces where there are (partial) constraints. The faces can repeat to constrain more vectors
    Eigen::MatrixXd constVectors; //corresponding to constFaces.
    
    int N;                        //Degree of field
    bool signSymmetry;            //whever field enforces a ssign symmetry (only when N is even, otherwise by default set to false)
    double wSmooth;               //Weight of smoothness
    double wRoSy;                 //Weight of rotational-symmetry
    Eigen::VectorXd wAlignment;   //Weight of alignment per each of the constfaces. "-1" means a fixed face
    
    int lapType;                  //Choice of weights (from UNIFORM_WEIGHTS,BARYCENTRIC_WEIGHTS or INV_COT_WEIGHTS)
    

    Eigen::SparseMatrix<std::complex<double>> smoothMat;
    Eigen::SparseMatrix<std::complex<double>> roSyMat;
    Eigen::SparseMatrix<std::complex<double>> alignmentMat;     //(soft) alignment energy.
    Eigen::SparseMatrix<std::complex<double>> reducMat;         //reducing the fixed dofs (for instance with sign symmetry or fixed partial constraints)
    Eigen::VectorXd reducRhs;                                   //The uncompressed PV coeffs are redrucMat*true_dofs+reducRhs
    Eigen::VectorXd alignRhs;                                   //encoding the soft constraints
    
    //Mass and stiffness matrices
    Eigen::SparseMatrix<std::complex<double>> WSmooth, WAlign, WRoSy, M;
    
    PolyVectorData():signSymmetry(true), lapType(BARYCENTRIC_WEIGHTS){}
    ~PolyVectorData(){}
  };
  
  
  // Precalculate the operators according to the user-prescribed parameters. Must be called whenever any of them changes
  // Inputs:
  //  V:      #V by 3 vertex coordinates.
  //  F:      #F by 3 face vertex indices.
  //  EV:     #E by 2 matrix of edges (vertex indices)
  //  EF:     #E by 2 matrix of oriented adjacent faces
  //  B1, B2: #F by 3 matrices representing the local base of each face.
  //  PolyVectorData (must fill constFaces, N, wSmooth, wOrth, wAlignment; lapType is default BARYCENTRIC and signSymmetry is "true"):     User parameters for the requested variation of the algorithm.
  // Outputs:
  //  PolyVectorData:       Updated structure with all operators
  IGL_INLINE void polyvector_precompute(const Eigen::MatrixXd& V,
                                        const Eigen::MatrixXi& F,
                                        const Eigen::MatrixXi& EV,
                                        const Eigen::MatrixXi& EF,
                                        const Eigen::MatrixXd& B1,
                                        const Eigen::MatrixXd& B2,
                                        PolyVectorData& pvData)
  {
    
    using namespace std;
    using namespace Eigen;
    
    
    //Building the smoothness matrices, with an energy term for each inner edge and degree
    int rowCounter=0;
    std::vector< Triplet<complex<double> > > dTriplets, WTriplets;
    if (pvData.N%2!=0) pvData.signSymmetry=false;  //it has to be for odd N
    
    
    /************Smoothness matrices****************/
    
    //Stiffness weights
    VectorXd stiffnessWeights=VectorXd::Zero(EF.rows());
    
    if (pvData.lapType==UNIFORM){
      stiffnessWeights.setOnes();
    } else {
      VectorXd faceCenter;
      double radius;  //stub
      if (pvData.lapType==BARYCENTRIC_WEIGHTS)
        igl::barycenter(V,F,faceCenter);
      if (pvData.lapType==INV_COT_WEIGHTS)
        directional::circumcircle(V,F,faceCenter,radius);
      for (int i=0;i<EF.rows;i++){
        if ((EF(i,0)==-1)||(EF(i,1)==-1))
          continue;  //boundary edge
        
        RowVector3d midEdge = (V.row(EV(i,0))+V.row(EV(i,1)))/2.0;
        dualLength = (faceCenter.row(EF(i,0))-midEdge).norm()+(faceCenter.row(EF(i,1))-midEdge).norm();   //TODO: that's only an approximation of barycentric height...
        primalLength = (V.row(EV(i,0))-V.row(EV(i,1))).norm();
        if (dualLength>10e-9)  //smaller than 10e-9might happen for inv cot weights
          stiffnessWeights(i)=primalLength/dualLength;
      }
    }
    
    for (int n = 0; n < pvData.N; n++)
    {
      for (int i=0;i<EF.rows();i++){
        
        if ((EF(i,0)==-1)||(EF(i,1)==-1))
          continue;  //boundary edge
        
        // Compute the complex representation of the common edge
        RowVector3d e = V.row(EV(i,1)) - V.row(EV(i,0));
        RowVector2d vef = Vector2d(e.dot(B1.row(EF(i,0))), e.dot(B2.row(EF(i,0)))).normalized();
        complex<double> ef(vef(0), vef(1));
        Vector2d veg = Vector2d(e.dot(B1.row(EF(i,1))), e.dot(B2.row(EF(i,1)))).normalized();
        complex<double> eg(veg(0), veg(1));
        
        // Add the term conj(f)^n*ui - conj(g)^n*uj to the differential matrix
        dTriplets.push_back(Triplet<complex<double> >(rowCounter, n*F.rows()+EF(i,0), pow(conj(ef), pvData.N-n)));
        dTriplets.push_back(Triplet<complex<double> >(rowCounter, n*F.rows()+EF(i,1), -1.*pow(conj(eg), pvData.N-n)));
        
        //stiffness weights
        WSmoothTriplets.push_back(Triplet<complex<double> >(rowCounter, rowCounter, stiffnessWeights(i)));
        rowCounter++;
      }
      
      for (int i=0;i<F.rows();i++)
        MTriplets.push_back(Triplet<complex<double>>(n*F.rows()+i, n*F,rows()+i, doubleAreas(i)/2.0));
    }
    
    pvData.smoothMat.resize(rowCounter, pvData.N*F.rows());
    pvData.smoothMat.setFromTriplets(dTriplets.begin(), dTriplets.end());
    
    pvData.WSmooth.resize(rowCounter, rowCounter);
    pvData.WSmooth.setFromTriplets(WSmoothTriplets.begin(), WSmoothTriplets.end());
    
    pvData.M.resize(pvData.N*F.rows(), pvData.N*F.rows());
    pvData.M.setFromTriplets(MTriplets.begin(), MTriplets.end());
    
    //creating reduction transformation
    VectorXi numFaceConstraints = EIgen::VectorXi::Zero(F.rows());
    int realN = (pvData.signSymmetry ? N/2 : N);
    MatrixXd faceConstraints(F.rows(),realN);
    std::vector<MatrixXcd> localFaceReducMats; localFaceReducMats.resize(F.rows());
    std::vector<VectorXcd> localFaceReducRhs;  localFaceReducRhs.resize(F.rows());
    
    for (int i=0;i<F.rows();i++){
      localFaceReducMats[i]=MatrixXcd::Identity(realN);
      localFaceReducRhs[i]=VectorXcd::Zero(realN);
    }
    
    /*************Reduction matrices******************/
    for (int i=0;i<pvData.constFaces.size();i++){
      if (wAlignment>=0.0)
        continue;  //here we only handle the reduction caused by a fixed dof
      if (numFaceConstraints(pvData.constFaces(i))==realN)
        continue; //overconstrained; we ignore any further constraints on that face
      
      complex<double> constVectorComplex=std::complex<double>(constVectors(i).dot(B1.row(pvData.constFaces(i))), constVectors(i).dot(B2.row(pvData.constFaces(i))));
      constVectorComplex = (pvData.signSymmetry ? constVectorComplex*constVectorComplex : constVectorComplex)
      faceconstraints(pvData.constFaces(i), numFaceConstraints(pvData.constFaces(i))++) = constVectorComplex;
      
      MatrixXcd singleReducMat=MatrixXcd::Zero(realN,realN-1);
      VectorXcd singleReducRhs=VectorXcd::Zero(realN);
      for (int j=0;j<realN-1;j++){
        singleReducMat(j,j)=-constVectorComplex;
        singleReducMat(j+1,j)=complex<double>(1.0,0.0);
      }
      singleReducRhs(realN-1) = -constVectorComplex;
      
      localFaceReducMats[pvData.constFaces(i)] = localFaceReducMats(pvData.constFaces(i))*singleReducMat;
      localFaceReducRhs[pvData.constFaces(i)] = localFaceReducMats(pvData.constFaces(i))*singleReducRhs + localFaceReducRhs(pvData.constFaces(i));
    }
    
   
    colCounter=0;
    pvData.reducRhs.resize(pvData.N*F.rows());
    vector<Triplets<complex<double>>> reducMatTriplets;
    int jump = (pvData.signSymmetry ? 2 : 1);
    for (int i=0;i<F.rows();i++){
      for (int j=0;j<pvData.N;j+=jump)
        for (int k=0;k<localFaceReducMats[i].cols();k++)
          reducMatTriplets.push_back(complex<double>(j*F.rows()+i, colCounter+k, localFaceReducMats[i](j,k)));
      pvData.reducRhs.segment(colCounter, localFaceReducMats[i].cols()) = localFaceReducRhs[i];
      colCounter+=localFaceReducMats[i].cols();
    }
    
    pvData.reducMat.resize(pvData.N*F.rows(), colCounter);
    pvData.reducMat.setFromTriplets(reducMatTriplets.begin(), reducMatTriplets.end());
    
    
    /****************rotational-symmetry matrices********************/
    vector<Triplet<complex<double>>> roSyTriplets, WRoSyTriplets;
    for (int i=F.rows();i<pvData.N*F.rows;i++){
      roSyTriplets.push_back(complex<double>(i,i,1.0));
      WRosyTriplets.push_back(complex<double>(i,i,doubleAreas(i%F.rows())));
    }
    
    pvData.roSyMat.resize(N*F.rows(), N*F.rows());
    pvData.roSyMat.setFromTriplets(roSyTriplets.begin(), roSyTriplets.end());
    
    pvData.WRoSy.resize(N*F.rows(), N*F.rows());
    pvData.WRoSy.setFromTriplets(WRoSyTriplets.begin(), WRoSyTriplets.end());
    
    
    /*****************Soft alignment matrices*******************/
    rowCounter=0;
    vector<Triplet<complex<double>>> alignTriplets;
    vector<VectorXd> alignRhsList;
    vector<Triplet<complex<double>>> WAlignTriplets;
    for (int i=0;i<pvData.constFaces.size();i++){
      if (wAlignment(i)>0.0)
        continue;  //here we only handle soft alignments
      if (numFaceConstraints(pvData.constFaces(i))==realN)
        continue; //overconstrained; we ignore any further constraints on that face
      
      complex<double> constVectorComplex=std::complex<double>(constVectors(i).dot(B1.row(pvData.constFaces(i))), constVectors(i).dot(B2.row(pvData.constFaces(i))));
      constVectorComplex = (pvData.signSymmetry ? constVectorComplex*constVectorComplex : constVectorComplex)
      faceconstraints(pvData.constFaces(i), numFaceConstraints(pvData.constFaces(i))++) = constVectorComplex;
      
      MatrixXcd singleReducMat=MatrixXcd::Zero(realN,realN-1);
      VectorXcd singleReducRhs=VectorXcd::Zero(realN);
      for (int j=0;j<realN-1;j++){
        singleReducMat(j,j)=-constVectorComplex;
        singleReducMat(j+1,j)=complex<double>(1.0,0.0);
      }
      singleReducRhs(realN-1) = -constVectorComplex;
      
      MatrixXcd IAiA = MatrixXcd::Identity(realN) - A*(A.completeOrthogonalDecomposition().pseudoInverse());
      singleReducRhs = IAiA*singleReducRhs;
      for (int j=0;j<pvData.N;j++)
        for (int k=0;k<pvData.N;k+=jump)
          alignTriplets.push_back(Triplet<complex<double>>(rowCounter+j, k*F.rows()+i, IAiA(j,k)));
      
      alignRhsList.push_back(singleReducRhs);
      for (int j=0;j<singleReducRhs.size();j++)
        WAlignTriplets.push_back(complex<double>(rowCounter+j, rowCounter+j, wAlignment(i)*doubleareas(i)/2.0));
      rowCounter+=realN;
    }
    
    pvData.alignRhs.resize(rowCounter);
    for (int i=0;i<alignRhsList.size();i++)
      pvData.alignRhs.segment(i*realN,realN)=alignRhsList[i];
    
    pvData.alignMat.resize(rowCounter, N*F.rows());
    pvData.alignMat.setFromTriplets(alignTriplets.begin(), alignTriplets.end());
    
    pvData.WAlign.resize(rowCounter,rowCounter);
    pvData.WAlign.setFromTriplets(WAlignTriplets.begin(), WAlignTriplets.end());
    
  }
  
  
  // Computes a polyvector on the entire mesh
  // Inputs:
  //  PolyVectorData: The data structure which should have been initialized with polyvector_precompute()
  // Outputs:
  //  polyVectorField: #F by N The output interpolated field, in polyvector (complex polynomial) format.
  IGL_INLINE void polyvector_field(const PolyVectorData& pvData,
                                   Eigen::MatrixXcd& polyVectorField)
  {
    using namespace std;
    using namespace Eigen;
    
    //forming total energy matrix;
    SparseMatrix<complex<double>> totalUnreducedLhs = pvData.wSmooth * pvData.smoothMat.adjoint()*pvData.WSmooth*pvData.smoothMat + pvData.wRosy*pvData.roSyMat.adjoint()*pvData.WRosy*pvData.rosyMat + pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignMat;
    VectorXcd totalUnreducedRhs = -pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignRhs;
    
    //TODO: make sparse matrix have a zero row in case no soft alignments
    
    SparseMatrix<complex<double>> totalLhs = totalUnreducedLhs*pvData.reducMat;
    VectorXcd totalRhs = -totalUnreducedLhs*pvData.reducRhs + totalUnreducedRhs;
    
    if (pvData.constFaces == 0)  //alignmat should be empty and the reduction matrix should be only sign symmetry, if applicable
    {
      
    }
    
   
    /*if (bc.size() == 0)
    {
      //extracting first eigenvector into the field
      //Have to use reals because libigl does not currently support complex eigs.
      SparseMatrix<double> M; igl::speye(2*B1.rows(), 2*B1.rows(), M);
      //creating a matrix of only the N-rosy interpolation
      SparseMatrix<std::complex<double> > AfullNRosy(Afull.rows()/N,Afull.cols()/N);
      std::vector<Triplet<std::complex<double> > > AfullNRosyTriplets;
      for (int k=0; k<Afull.outerSize(); ++k)
        for (SparseMatrix<std::complex<double> >::InnerIterator it(Afull,k); it; ++it)
        {
          if ((it.row()<Afull.rows()/N)&&(it.col()<Afull.cols()/N))
            AfullNRosyTriplets.push_back(Triplet<std::complex<double> > (it.row(), it.col(), it.value()));
        }
      
      AfullNRosy.setFromTriplets(AfullNRosyTriplets.begin(), AfullNRosyTriplets.end());
      
      SparseMatrix<std::complex<double>> LComplex =AfullNRosy.adjoint()*AfullNRosy;
      SparseMatrix<double> L(2*B1.rows(),2*B1.rows());
      std::vector<Triplet<double> > LTriplets;
      for (int k=0; k<LComplex.outerSize(); ++k)
        for (SparseMatrix<std::complex<double>>::InnerIterator it(LComplex,k); it; ++it)
        {
          LTriplets.push_back(Triplet<double>(it.row(), it.col(), it.value().real()));
          LTriplets.push_back(Triplet<double>(it.row(), LComplex.cols()+it.col(), -it.value().imag()));
          LTriplets.push_back(Triplet<double>(LComplex.rows()+it.row(), it.col(), it.value().imag()));
          LTriplets.push_back(Triplet<double>(LComplex.rows()+it.row(), LComplex.cols()+it.col(), it.value().real()));
        }
      L.setFromTriplets(LTriplets.begin(), LTriplets.end());
      Eigen::MatrixXd U;
      Eigen::VectorXd S;
      igl::eigs(L,M,10,igl::EIGS_TYPE_SM,U,S);
      int smallestIndex; S.minCoeff(&smallestIndex);
      
      polyVectorField=MatrixXcd::Constant(B1.rows(), N, complex<double>());
      
      polyVectorField.col(0) = U.block(0,smallestIndex,U.rows()/2,1).cast<std::complex<double> >().array()*std::complex<double>(1,0)+
      U.block(U.rows()/2,smallestIndex,U.rows()/2,1).cast<std::complex<double> >().array()*std::complex<double>(0,1);
      //cout<<"polyVectorField.col(0): "<<polyVectorField.col(0)<<endl;
      //cout<<"B1: "<<B1<<endl;
      return;
    }
    
    MatrixXcd constValuesMat(b.rows(),N);
    assert((b.cols()==3*N)||(b.cols()==3));
    if (b.cols()==3)  //N-RoSy constraint
    {
      constValuesMat.setZero();
      for (int i=0;i<b.rows();i++){
        complex<double> bComplex=complex<double>(b.row(i).dot(B1.row(bc(i))), b.row(i).dot(B2.row(bc(i))));
        constValuesMat(i,0)=-pow(bComplex, N);
      }
    } else {
      for (int i=0;i<b.rows();i++){
        RowVectorXcd poly,roots(N);
        for (int n=0;n<N;n++){
          RowVector3d vec=b.block(i,3*n,1,3);
          roots(n)=complex<double>(vec.dot(B1.row(bc(i))), vec.dot(B2.row(bc(i))));
        }
        roots_to_monicPolynomial(roots, poly);
        constValuesMat.row(i)<<poly.head(N);
      }
    }
    
    VectorXi constIndices(N*bc.size());
    VectorXcd constValues(N*b.size());
    
    for (int n=0;n<N;n++){
      constIndices.segment(bc.rows()*n,bc.rows())=bc.array()+n*B1.rows();
      constValues.segment(b.rows()*n,b.rows())=constValuesMat.col(n);
    }
    
    VectorXcd torhs(N*B1.rows(),1);
    torhs.setZero();
    for (int i=0;i<constIndices.size();i++)
      torhs(constIndices(i))=constValues(i);
    
    VectorXcd rhs=-AVar.adjoint()*Afull*torhs;
    VectorXcd varFieldVector=solver.solve(rhs);
    assert(solver.info() == Success);
    
    VectorXcd polyVectorFieldVector(N*B1.rows());
    VectorXi varMask=VectorXi::Constant(N*B1.rows(),1);
    for (int i=0;i<constIndices.size();i++)
      varMask(constIndices(i))=0;
    
    VectorXi full2var=VectorXi::Constant(N*B1.rows(),-1);
    int varCounter=0;
    for (int i=0;i<N*B1.rows();i++)
      if (varMask(i))
        full2var(i)=varCounter++;
    
    assert(varCounter==N*(B1.rows()-bc.size()));
    
    for (int i=0;i<constIndices.size();i++)
      polyVectorFieldVector(constIndices(i))=constValues(i);
    
    for (int i=0;i<N*B1.rows();i++)
      if (full2var(i)!=-1)
        polyVectorFieldVector(i)=varFieldVector(full2var(i));
    
    //converting to matrix form
    polyVectorField.conservativeResize(B1.rows(),N);
    for (int n=0;n<N;n++)
      polyVectorField.col(n)=polyVectorFieldVector.segment(n*B1.rows(),B1.rows());*/
    
  }
  
  
  // minimal version without auxiliary data
  IGL_INLINE void polyvector_field(const Eigen::MatrixXd& V,
                                   const Eigen::MatrixXi& F,
                                   const Eigen::VectorXi& bc,
                                   const Eigen::MatrixXd& b,
                                   const int N,
                                   Eigen::MatrixXcd& polyVectorField)
  {
    Eigen::MatrixXi EV, xi, EF;
    igl::edge_topology(V, F, EV, xi, EF);
    Eigen::MatrixXd B1, B2, xd;
    igl::local_basis(V, F, B1, B2, xd);
    Eigen::SparseMatrix<std::complex<double>> Afull, AVar;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    polyvector_precompute(V,F,EV,EF,B1,B2,bc,N, solver,Afull,AVar);
    polyvector_field(B1, B2, bc, b, solver, Afull, AVar, N, polyVectorField);
  }
}
#endif
