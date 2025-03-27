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

#pragma once

#include "yy_matrix.hpp"
#include "yy_diagonal_matrix.hpp"

namespace yafiyogi::yy_maths {

class ekf final
{
  public:
    using value_type = double;
    static constexpr value_type EPS = 1e-4;

    using matrix = yy_maths::matrix<value_type>;
    using identity_matrix = yy_maths::identity_matrix<value_type>;
    using diagonal_matrix_eps = diagonal_matrix_fixed<value_type, EPS>;
    using diagonal_matrix_neg = diagonal_matrix_fixed<value_type, -1.0>;
    using diagonal_matrix_type = diagonal_matrix<value_type>;
    using zero_matrix = yy_maths::zero_matrix<value_type>;
    using vector = yy_maths::vector<value_type>;
    using zero_vector = yy_maths::zero_vector<value_type>;
    using size_type = matrix::size_type;

    ekf(size_type p_m,
        size_type p_n) noexcept;
    ekf(size_type p_m,
        size_type p_n,
        vector & p_r) noexcept;

    constexpr ekf() noexcept = default;
    ekf(const ekf & other) noexcept = default;
    ekf(ekf && other) noexcept;

    ekf & operator=(const ekf & other) noexcept = default;
    ekf & operator=(ekf && other) noexcept;

    void predict() noexcept;
    bool update(const vector & p_z, // observations m wide
                const matrix & p_h, // m x n (m -> inputs, n -> outputs)
                const vector & p_hx) noexcept; // m wide

    const vector & X() const noexcept
    {
      return m_x;
    }

    const value_type & X(size_type idx) const noexcept
    {
      return m_x(idx);
    }

    constexpr size_type N() const noexcept
    {
      return m_n;
    }

    constexpr size_type M() const noexcept
    {
      return m_m;
    }

  private:
    size_type m_n = 0; // Number of outputs
    size_type m_m = 0; // Number of inputs
    vector m_x{}; // State vector.
    matrix m_P{}; // Prediction error covariance
    diagonal_matrix_type m_R{}; // Measurement noise.
};

} // namespace yafiyogi::yy_maths
