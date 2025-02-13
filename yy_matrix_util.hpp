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

#pragma once

#include "yy_matrix.hpp"

namespace yafiyogi::yy_maths {

namespace matrix_util_detail {
// From https://web.archive.org/web/20231002021242/http://jean-pierre.moreau.pagesperso-orange.fr:80/Cplus/choles_cpp.txt
// and https://github.com/simondlevy/TinyEKF/blob/master/src/tinyekf.h

template<typename T>
constexpr bool choldc1(matrix<T> & a,
                       vector<T> & p) noexcept
{
  using matrix_type = matrix<T>;
  using value_type = matrix_type::value_type;
  using size_type = matrix_type::size_type;

  const size_type size = a.size1();

  for(size_type i = 0; i < size; ++i)
  {
    for(size_type j = i; j < size; ++j)
    {
      value_type sum = a(i, j);

      for(int k = static_cast<int>(i) - 1; k >= 0; --k)
      {
        sum -= a(i, static_cast<size_type>(k)) * a(j, static_cast<size_type>(k));
      }

      if(i == j)
      {
        if(sum <= 0)
        {
          return false; /* error */
        }
        p(i) = sqrt(sum);
      }
      else
      {
        a(j, i) = sum / p(i);
      }
    }
  }

  return true; // success
}

template<typename T>
constexpr bool choldcsl(const matrix<T> & A,
                        matrix<T> & a,
                        vector<T> & p) noexcept
{
  using matrix_type = matrix<T>;
  using value_type = matrix_type::value_type;
  using size_type = matrix_type::size_type;

  a = A;

  if(!choldc1(a, p))
  {
    return false;
  }

  const size_type size = A.size1();
  for(size_type i = 0; i < size; ++i)
  {
    a(i, i) = 1 / p[i];
    for(size_type j = i + 1; j < size; ++j)
    {
      value_type sum = 0.0;
      for(size_type k = i; k < j; ++k)
      {
        sum -= a(j, k) * a(k, i);
      }
      a(j, i) = sum / p[j];
    }
  }

  return true;
}

template<typename T>
constexpr bool cholsl(const matrix<T> & A,
                      matrix<T> & a,
                      vector<T> & p) noexcept
{
  using matrix_type = matrix<T>;
  using value_type = matrix_type::value_type;
  using size_type = matrix_type::size_type;

  if(!choldcsl(A,a,p))
  {
    return false;
  }

  const size_type size = A.size1();
  for(size_type i = 0; i < size; ++i)
  {
    for(size_type j = i + 1; j < size; ++j)
    {
      a(i, j) = 0.0;
    }
  }

  for(size_type i = 0; i < size; ++i)
  {
    auto & a_ii = a(i, i);

    a_ii *= a_ii;
    for(size_type k = i + 1; k < size; ++k)
    {
      value_type a_ki = a(k, i);

      a_ii += a_ki * a_ki;
    }

    for(size_type j = i + 1; j < size; ++j)
    {
      for(size_type k = j; k < size; ++k)
      {
        a(i, j) += a(k, i) * a(k, j);
      }
    }
  }

  for(size_type i = 0; i < size; ++i)
  {
    for(size_type j = 0; j < i; ++j)
    {
      a(i, j) = a(j, i);
    }
  }

  return true; // success
}

} // namespace matrix_util_detail

template<typename T>
constexpr bool invert(const matrix<T> & A,
                      matrix<T> & a) noexcept
{
  if((A.size1() != A.size2())
     || (A.size1() != a.size2())
     || (a.size1() != a.size2()))
  {
    return false;
  }

  vector<T> tmp{A.size1()};

  return matrix_util_detail::cholsl(A, a, tmp);
}

} // namespace yafiyogi::yy_maths
