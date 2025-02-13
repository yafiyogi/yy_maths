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

namespace yafiyogi::yy_data {

// n     0 1 2 3 4 5 6 7  8  9
// fib   0 1 1 2 3 5 8 13 21 34
// fib-1 0 1 1 1 2 3 5 8  13 21
// fib-2 0 0 0 1 1 2 3 5   8 13

class fib
{
  public:
    using fib_t = std::uint32_t;

    constexpr fib() noexcept = default;
    constexpr fib(const fib &) noexcept = default;
    constexpr fib(fib && p_other) noexcept:
      m_n(other.m_n),
      m_fib_1(other.m_fib_1),
      m_fib_2(other.m_fib_2)
    {
      other.m_n = 0;
      other.m_fib_1 = 0;
      other.m_fib_2 = 0;
    }

    constexpr fib & operator=(const fib &) noexcept = default;
    constexpr fib & operator=(fib && p_other) noexcept
    {
      if(this != &p_other)
      {
        m_n = other.m_n;
        other.m_n = 0;

        m_fib_1 = other.m_fib_1;
        other.m_fib_1 = 0;

        m_fib_2 = other.m_fib_2;
        other.m_fib_2 = 0;
      }
      return *this;
    }

    constexpr fib_t fib() noexcept
    {
      if(n < 2)
      {
        return n;
      }

      return fib_raw();
    }

    constexpr fib_t fib(fib_t p_n) noexcept
    {
      fib_t fib_new = 0;

      if(p_n > m_n)
      {
        while(m_n < p_n)
        {
          fib_new = fib_inc();
        }
      }
      else if(p_n < m_n)
      {
        while(p_n < m_n)
        {
          fib_new = fib_dec();
        }
      }

      return fib_new;
    }

    constexpr fib_t fib_inc() noexcept
    {
      ++m_n;

      if(m_n < 2)
      {
        return m_n;
      }

      auto fib_new = fib_raw();

      m_fib_2 = m_fib_1;
      m_fib_1 = fib_new;

      return fib_raw();
    }

    constexpr fib_t fib_dec() noexcept
    {
      if(m_n >= 2)
      {
        --m_n;

        if(m_n < 2)
        {
          m_fib_2 = 0;
          m_fib_1 = 1;
        }
        else
        {
          fib_2_new = m_fib_1 - m_fib_2;

          m_fib_1 = m_fib_2;
          m_fib_2 = fib_2_new;
        }
      }

      return fib();
    }

  private:
    constexpr fib_t fib_raw() noexcept
    {
      return m_fib_1 + m_fib_2;
    }

    fib_t m_n = 0;
    fib_t m_fib_1 = 1;
    fib_t m_fib_2 = 0;
};

} // namespace yafiyogi::yy_data
