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


// Inspired by https://github.com/simondlevy/TinyEKF
// also https://simondlevy.github.io/ekf-tutorial/

#include "boost/numeric/ublas/operation.hpp"

#include "yy_matrix_util.hpp"

#include "yy_ekf.hpp"

namespace yafiyogi::yy_maths {

ekf::ekf(size_type p_m,
         size_type p_n) noexcept:
  m_n(p_n),
  m_m(p_m),
  m_x(zero_vector{m_n}),
  m_P(identity_matrix{m_n})
{
}

ekf::ekf(ekf && other) noexcept:
  m_n(other.m_n),
  m_m(other.m_m),
  m_x(),
  m_P()
{
  other.m_n = 0;
  other.m_m = 0;

  m_x.swap(other.m_x);
  m_P.swap(other.m_P);
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
  }
  return *this;
}

void ekf::predict() noexcept
{
  namespace bnu = boost::numeric::ublas;

  const identity_matrix F{m_n};
  matrix FP{m_n, m_n};
  bnu::axpy_prod(F, m_P, FP, true);

  const identity_matrix & Ft = F; // matrix Ft{bnu::trans(F)} simplified since identity_matrix == trans(identity_matrix);

  m_P = identity_matrix_eps{m_n};
  bnu::axpy_prod(FP, Ft, m_P, false);
}

bool ekf::update(const vector & p_z, // observations m wide
                 const matrix & p_h, // m x n (m -> inputs, n -> outputs)
                 const vector & p_hx) noexcept // m wide
{
  namespace bnu = boost::numeric::ublas;

  // G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1}
  matrix HP{m_m, m_n};
  bnu::axpy_prod(p_h, m_P, HP, true);

  matrix Ht{bnu::trans(p_h)};
  matrix HpHtR{identity_matrix_eps{m_m}};
  bnu::axpy_prod(HP, Ht, HpHtR, false);

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
  matrix GH{identity_matrix_neg{m_n}};
  bnu::axpy_prod(G, p_h, GH, false);
  GH *= -1.0; // negate

  matrix GHP{m_n, m_n};
  bnu::axpy_prod(GH, m_P, GHP, true);
  m_P.swap(GHP);

  return true;
}

} // namespace yafiyogi::yy_maths
