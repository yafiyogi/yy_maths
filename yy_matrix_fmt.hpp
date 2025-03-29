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

#include<string_view>

#include "fmt/format.h"
#include "fmt/ranges.h"

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"

template<>
struct fmt::formatter<boost::numeric::ublas::matrix<double>>
{
    using ublas_matrix = boost::numeric::ublas::matrix<double>;
    using array_type = ublas_matrix::array_type;
    using size_type = array_type::size_type;
    using iterator_type = array_type::const_iterator;
    using matrix_view = fmt::join_view<iterator_type, iterator_type>;

    fmt::formatter<matrix_view> matrix_formatter{};

    template <typename ParseContext>
    constexpr auto parse(ParseContext & ctx)
    {
      return matrix_formatter.parse(ctx);
    }

    auto format(const ublas_matrix & p_matrix,
                format_context & ctx) const
    {
      fmt::format_to(ctx.out(),
                     std::string_view{"({},{})("},
                     p_matrix.size1(),
                     p_matrix.size2());

      const size_type rows = p_matrix.size1();
      const size_type columns = p_matrix.size2();

      auto data{p_matrix.data().begin()};
      for(size_type row_idx = 0; row_idx < rows; ++row_idx)
      {
        fmt::format_to(ctx.out(), std::string_view{"("});
        matrix_view row{data, data + columns, std::string_view{","}};

        matrix_formatter.format(row, ctx);

        fmt::format_to(ctx.out(), std::string_view{")"});
        data += columns;
      }

      fmt::format_to(ctx.out(), std::string_view{")"});
      return ctx.out();
    }
};

template<>
struct fmt::formatter<boost::numeric::ublas::vector<double>>
{
    using ublas_vector = boost::numeric::ublas::vector<double>;
    using iterator_type = ublas_vector::const_iterator;
    using vector_view = fmt::join_view<iterator_type, iterator_type>;

    fmt::formatter<vector_view> vector_formatter{};

    template <typename ParseContext>
    constexpr auto parse(ParseContext & ctx)
    {
      return vector_formatter.parse(ctx);
    }

    template<typename FormatContext>
    auto format(const ublas_vector & p_vector,
                FormatContext & ctx) const
    {
      fmt::format_to(ctx.out(), std::string_view{"({})("}, p_vector.size());

      vector_view view{p_vector.begin(), p_vector.end(), std::string_view{","}};

      vector_formatter.format(view, ctx);

      fmt::format_to(ctx.out(), std::string_view{")"});

      return ctx.out();
    }
};
