/*

  MIT License

  Copyright (c) 2025 Yafiyogi

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "boost/numeric/ublas/operation.hpp"

#include "yy_matrix_util.hpp"

#include "yy_ekf.hpp"

namespace yafiyogi::yy_maths {

ekf::ekf(vector & p_diag,
         matrix & p_h) noexcept:
  m_n(p_h.size2()),
  m_m(p_h.size1()),
  m_x(zero_vector{m_n}),
  m_P(zero_matrix{m_n, m_n}),
  m_Q(matrix{identity_matrix{m_n, m_n}} * EPS),
  m_R(matrix{identity_matrix{m_m, m_m}} * EPS),
  m_F(identity_matrix{m_n, m_n}),
  m_H(p_h)
{
  const size_type size = p_diag.size();

  for(size_type i = 0; i < size; ++i)
  {
    m_P(i, i) = p_diag(i);
  }
}

ekf::ekf(ekf && other) noexcept:
  m_n(other.m_n),
  m_m(other.m_m),
  m_x(),
  m_P(),
  m_Q(),
  m_R(),
  m_F(),
  m_H()
{
  other.m_n = 0;
  other.m_m = 0;

  m_x.swap(other.m_x);
  m_P.swap(other.m_P);
  m_Q.swap(other.m_Q);
  m_R.swap(other.m_R);
  m_F.swap(other.m_F);
  m_H.swap(other.m_H);
}

ekf & ekf::operator=(ekf && other) noexcept
{
  if(this != &other)
  {
    m_n = other.m_n;
    other.m_n = 0;
    m_m = other.m_m;
    other.m_m = 0;

    m_x = vector{};
    m_x.swap(other.m_x);
    m_P = matrix{};
    m_P.swap(other.m_P);
    m_Q = matrix{};
    m_Q.swap(other.m_Q);
    m_R = matrix{};
    m_R.swap(other.m_R);
    m_F = matrix{};
    m_F.swap(other.m_F);
    m_H = matrix{};
    m_H.swap(other.m_H);
  }
  return *this;
}

void ekf::predict() noexcept
{
  namespace bnu = boost::numeric::ublas;

  matrix FP{m_n, m_n};
  bnu::axpy_prod(m_F, m_P, FP, true);

  matrix Ft{bnu::trans(m_F)};

  matrix FPFt{m_n, m_n};
  bnu::axpy_prod(FP, Ft, FPFt, true);
  FPFt += m_Q;

  m_P.swap(FPFt);
}

bool ekf::update(const vector & p_z, // observations
                 const vector & p_hx) noexcept
{
  namespace bnu = boost::numeric::ublas;

  // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
  matrix HP{m_m, m_n};
  bnu::axpy_prod(m_H, m_P, HP, true);

  matrix Ht{bnu::trans(m_H)};
  matrix HpHtR{m_m, m_m};
  bnu::axpy_prod(HP, Ht, HpHtR, true);
  HpHtR += m_R;

  matrix HPHtRinv{m_m, m_m};
  if(!invert(HpHtR, HPHtRinv))
  {
    return false;
  }

  matrix PHt{m_n, m_m};
  bnu::axpy_prod(m_P, Ht, PHt, true);

  matrix G{m_n, m_m};
  bnu::axpy_prod(PHt, HPHtRinv, G, true);

  // \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k))
  vector z_hx{p_z};
  z_hx -= p_hx;

  bnu::axpy_prod(G, z_hx, m_x, false);

  // P_k = (I - G_k H_k) P_k
  matrix GH{m_n, m_n};
  bnu::axpy_prod(G, m_H, GH, true);
  GH *= -1; // negate
  GH += identity_matrix{GH.size1()};

  matrix GHP{m_n, m_n};
  bnu::axpy_prod(GH, m_P, GHP, true);
  m_P.swap(GHP);

  return true;
}

} // namespace yafiyogi::yy_maths
