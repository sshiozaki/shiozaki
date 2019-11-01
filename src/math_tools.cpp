/**
 * @file math_tools.cpp
 * @brief mathTool class
 * @author T. Otani
 */

#include "math_tools.h"
#include <omp.h>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

// #################################################################
/**
 * @brief calc norm of vector
 */
double mathTool::vectorNorm(const int &nump,const DOUBLEARRAY1 &x)
{
  double norm=0e0;

  #pragma omp parallel for reduction(+:norm)
  for(int i=0;i<nump;i++){
    norm += fabs(x[i]);
  }
  return norm;
}

// #################################################################
/**
 * @brief calc inner product
 */
double mathTool::innerProduct(const int &nump,const DOUBLEARRAY1 &x,const DOUBLEARRAY1 &y)
{
  double dot_p=0e0;

  #pragma omp parallel for reduction(+:dot_p)
  for(int i=0;i<nump;i++){
    dot_p += x[i] * y[i];
  }
  return dot_p;
}

// #################################################################
/**
 * @brief calc cross product
 */
void mathTool::crossProduct(const double (&a)[3],const double (&b)[3],double (&c)[3],double &dc)
{
  dc = 0e0;

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];

  for(int j=0;j<3;j++) dc += (c[j]*c[j]);

  dc = sqrt(dc);
}
// #################################################################
/**
 * @brief calc determinant
 */
double mathTool::calcDeterminant_3x3(const DOUBLEARRAY2 &a)
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}
// #################################################################
/**
 * @brief calc determinant
 */
double mathTool::calcDeterminant_3x3(const double (&a)[3][3])
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}
// #################################################################
/**
 * @brief calc inverse matrix
 */
void mathTool::calcInverseMatrix_3x3(DOUBLEARRAY2 &inv_a,const DOUBLEARRAY2 &a)
{
  double det = calcDeterminant_3x3(a);

  inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
  inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];

  inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
  inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];

  inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
  inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) inv_a[i][j] = inv_a[i][j] / det;
  }
}
// #################################################################
/**
 * @brief calc inverse matrix
 */
void mathTool::calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3])
{
  double det;

  det = calcDeterminant_3x3(a);

  inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
  inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
  inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
  inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) inv_a[i][j] = inv_a[i][j] / det;
  }
}

// #################################################################
/**
 * @brief calc matrix x matrrix4
 */
void mathTool::calcMatrix_x_matrix4(double (&ans)[4][4],const double (&a)[4][4],const double (&b)[4][4])
{
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      ans[i][j] = 0e0;
      for(int k=0;k<4;k++) ans[i][j] += a[i][k] * b[k][j];
    }
  }

}

/*************************************************************/
/*  実対称行列の固有値・固有ベクトル（ヤコビ法）              */
/*      n : 次数                                             */
/*      ct : 最大繰り返し回数                                */
/*      eps : 収束判定条件                                   */
/*      A : 対象とする行列                                   */
/*      A1, A2 : 作業域（nxnの行列），A1の対角要素が固有値   */
/*      X1, X2 : 作業域（nxnの行列），X1の各列が固有ベクトル */
/*      return : =0 : 正常                                   */
/*               =1 : 収束せず                               */
/*      coded by Y.Suganuma                                  */
/*************************************************************/
int mathTool::Jacobi3x3(const int &ct, const double &eps, double (&A)[3][3], double (&A1)[3][3],double (&X1)[3][3])
{
  int n=3;
  double A2[3][3],X2[3][3];
  double max, s, t, v, sn, cs;
  int i1, i2, k = 0, ind = 1, p = 0, q = 0;
          // 初期設定
  for (i1 = 0; i1 < n; i1++) {
    for (i2 = 0; i2 < n; i2++) {
      A1[i1][i2] = A[i1][i2];
      X1[i1][i2] = 0.0;
    }
    X1[i1][i1] = 1.0;
  }
          // 計算
  while (ind > 0 && k < ct) {
            // 最大要素の探索
    max = 0.0;
    for (i1 = 0; i1 < n; i1++) {
      for (i2 = 0; i2 < n; i2++) {
        if (i2 != i1) {
          if (fabs(A1[i1][i2]) > max) {
            max = fabs(A1[i1][i2]);
            p   = i1;
            q   = i2;
          }
        }
      }
    }
            // 収束判定
              // 収束した
    if (max < eps)
      ind = 0;
              // 収束しない
    else {
                // 準備
      s  = -A1[p][q];
      t  = 0.5 * (A1[p][p] - A1[q][q]);
      v  = fabs(t) / sqrt(s * s + t * t);
      sn = sqrt(0.5 * (1.0 - v));
      if (s*t < 0.0)
        sn = -sn;
      cs = sqrt(1.0 - sn * sn);
                // Akの計算
      for (i1 = 0; i1 < n; i1++) {
        if (i1 == p) {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == p)
              A2[p][p] = A1[p][p] * cs * cs + A1[q][q] * sn * sn -
                                       2.0 * A1[p][q] * sn * cs;
            else if (i2 == q)
              A2[p][q] = 0.0;
            else
              A2[p][i2] = A1[p][i2] * cs - A1[q][i2] * sn;
          }
        }
        else if (i1 == q) {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == q)
              A2[q][q] = A1[p][p] * sn * sn + A1[q][q] * cs * cs +
                                       2.0 * A1[p][q] * sn * cs;
            else if (i2 == p)
              A2[q][p] = 0.0;
            else
              A2[q][i2] = A1[q][i2] * cs + A1[p][i2] * sn;
          }
        }
        else {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == p)
              A2[i1][p] = A1[i1][p] * cs - A1[i1][q] * sn;
            else if (i2 == q)
              A2[i1][q] = A1[i1][q] * cs + A1[i1][p] * sn;
            else
              A2[i1][i2] = A1[i1][i2];
          }
        }
      }
                // Xkの計算
      for (i1 = 0; i1 < n; i1++) {
        for (i2 = 0; i2 < n; i2++) {
          if (i2 == p)
            X2[i1][p] = X1[i1][p] * cs - X1[i1][q] * sn;
          else if (i2 == q)
            X2[i1][q] = X1[i1][q] * cs + X1[i1][p] * sn;
          else
            X2[i1][i2] = X1[i1][i2];
        }
      }
                // 次のステップへ
      k++;
      for (i1 = 0; i1 < n; i1++) {
        for (i2 = 0; i2 < n; i2++) {
          A1[i1][i2] = A2[i1][i2];
          X1[i1][i2] = X2[i1][i2];
        }
      }
    }
  }

  return ind;
}

