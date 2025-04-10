/*

  Copyright (c) 2000-2010
  Joerg Walter, Mathias Koch, Gunter Winkler, David Bellot
  Copyright (c) 2014, Athanasios Iliopoulos

  Distributed under the Boost Software License, Version 1.0. (See
  accompanying file LICENSE_1_0.txt or copy at
  http://www.boost.org/LICENSE_1_0.txt)

  The authors gratefully acknowledge the support of
  GeNeSys mbH & Co. KG in producing this work.

  Boost Software License - Version 1.0 - August 17th, 2003

  Permission is hereby granted, free of charge, to any person or organization
  obtaining a copy of the software and accompanying documentation covered by
  this license (the "Software") to use, reproduce, display, distribute,
  execute, and transmit the Software, and to prepare derivative works of the
  Software, and to permit third-parties to whom the Software is furnished to
  do so, all subject to the following:

  The copyright notices in the Software and this entire statement, including
  the above license grant, this restriction and the following disclaimer,
  must be included in all copies of the Software, in whole or in part, and
  all derivative works of the Software, unless such copies or derivative
  works are solely in the form of machine-executable object code generated by
  a source language processor.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
  SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
  FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.
*/

/*

  MIT License

  Parts copyright (c) 2025 Yafiyogi

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

#include "boost/numeric/ublas/matrix.hpp"

// diagonal_matrix_fixed is copied from identity_matrix in boost/numeric/ublas/matrix.hpp
// Amended by me.

namespace yafiyogi::yy_maths {
namespace bnu = boost::numeric::ublas;

template<typename T,
         T Value,
         typename ALLOC = std::allocator<T>>
class diagonal_matrix_fixed:
      public bnu::matrix_container<diagonal_matrix_fixed<T, Value, ALLOC> > {

    typedef const T *const_pointer;
    typedef diagonal_matrix_fixed<T, Value, ALLOC> self_type;
  public:
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
    using bnu::matrix_container<self_type>::operator ();
#endif
    typedef typename boost::allocator_size_type<ALLOC>::type size_type;
    typedef typename boost::allocator_difference_type<ALLOC>::type difference_type;
    typedef T value_type;
    typedef const T &const_reference;
    typedef T &reference;
    typedef const bnu::matrix_reference<const self_type> const_closure_type;
    typedef bnu::matrix_reference<self_type> closure_type;
    typedef bnu::sparse_tag storage_category;
    typedef bnu::unknown_orientation_tag orientation_category;

    // Construction and destruction
    BOOST_UBLAS_INLINE
    diagonal_matrix_fixed ():
      bnu::matrix_container<self_type> (),
      m_size (0), size_common_(0) {}
    BOOST_UBLAS_INLINE
    diagonal_matrix_fixed (size_type size):
      bnu::matrix_container<self_type> (),
      m_size (size), size_common_(size) {}
    BOOST_UBLAS_INLINE
    diagonal_matrix_fixed (size_type size1, size_type size2):
      bnu::matrix_container<self_type> (),
      m_size (std::min(size1, size2)), size_common_(std::min(size1, size2)) {}
    BOOST_UBLAS_INLINE
    diagonal_matrix_fixed (const diagonal_matrix_fixed &m):
      bnu::matrix_container<self_type> (),
      m_size (m.m_size), size_common_(m.size_common_) {}

    // Accessors
    BOOST_UBLAS_INLINE
    size_type size1 () const {
      return m_size;
    }
    BOOST_UBLAS_INLINE
    size_type size2 () const {
      return m_size;
    }

    // Resizing
    BOOST_UBLAS_INLINE
    void resize (size_type size, bool /*preserve*/ = true) {
      m_size = size;
    }
    BOOST_UBLAS_INLINE
    void resize (size_type size1, size_type size2, bool /*preserve*/ = true) {
      m_size = std::min(size1, size2);
    }

    // Element access
    BOOST_UBLAS_INLINE
    const_reference operator () (size_type i, size_type j) const {
      if (i == j)
        return m_value;
      else
        return m_zero;
    }

    // Assignment
    BOOST_UBLAS_INLINE
    diagonal_matrix_fixed &operator = (const diagonal_matrix_fixed &m) {
      m_size = m.m_size;
      size_common_ = m.size_common_;
      return *this;
    }
    BOOST_UBLAS_INLINE
    diagonal_matrix_fixed &assign_temporary (diagonal_matrix_fixed &m) {
      swap (m);
      return *this;
    }

    // Swapping
    BOOST_UBLAS_INLINE
    constexpr void swap (diagonal_matrix_fixed &m) noexcept {
      if (this != &m) {
        std::swap (m_size, m.m_size);
        std::swap (size_common_, m.size_common_);
      }
    }
    BOOST_UBLAS_INLINE
    friend constexpr void swap (diagonal_matrix_fixed &m1, diagonal_matrix_fixed &m2) noexcept {
      m1.swap (m2);
      m1.swap (m2);
    }

    // Iterator types
  private:
    // Use an index
    typedef size_type const_subiterator_type;

  public:
    class const_iterator1;
    class const_iterator2;
    typedef bnu::reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
    typedef bnu::reverse_iterator_base2<const_iterator2> const_reverse_iterator2;

    // Element lookup
    BOOST_UBLAS_INLINE
    const_iterator1 find1 (int rank, size_type i, size_type j) const {
      if (rank == 1) {
        i = (std::max) (i, j);
        i = (std::min) (i, j + 1);
      }
      return const_iterator1 (*this, i);
    }
    BOOST_UBLAS_INLINE
    const_iterator2 find2 (int rank, size_type i, size_type j) const {
      if (rank == 1) {
        j = (std::max) (j, i);
        j = (std::min) (j, i + 1);
      }
      return const_iterator2 (*this, j);
    }

    class const_iterator1:
      public bnu::container_const_reference<diagonal_matrix_fixed>,
      public bnu::bidirectional_iterator_base<bnu::sparse_bidirectional_iterator_tag,
                                              const_iterator1, value_type> {
      public:
        typedef typename diagonal_matrix_fixed::value_type value_type;
        typedef typename diagonal_matrix_fixed::difference_type difference_type;
        typedef typename diagonal_matrix_fixed::const_reference reference;
        typedef typename diagonal_matrix_fixed::const_pointer pointer;

        typedef const_iterator2 dual_iterator_type;
        typedef const_reverse_iterator2 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator1 ():
          bnu::container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator1 (const self_type &m, const const_subiterator_type &it):
          bnu::container_const_reference<self_type> (m), it_ (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator1 &operator ++ () {
          BOOST_UBLAS_CHECK (it_ < (*this) ().size1 (), bnu::bad_index ());
          ++it_;
          return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator1 &operator -- () {
          BOOST_UBLAS_CHECK (it_ > 0, bnu::bad_index ());
          --it_;
          return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const {
          return m_value;
        }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 begin () const {
          return const_iterator2 ((*this) (), it_);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 cbegin () const {
          return begin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 end () const {
          return const_iterator2 ((*this) (), it_ + 1);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 cend () const {
          return end ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 rbegin () const {
          return const_reverse_iterator2 (end ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 crbegin () const {
          return rbegin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 rend () const {
          return const_reverse_iterator2 (begin ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 crend () const {
          return rend ();
        }
#endif

        // Indices
        BOOST_UBLAS_INLINE
        size_type index1 () const {
          return it_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
          return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator1 &operator = (const const_iterator1 &it) {
          bnu::container_const_reference<self_type>::assign (&it ());
          it_ = it.it_;
          return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const const_iterator1 &it) const {
          BOOST_UBLAS_CHECK (&(*this) () == &it (), bnu::external_logic ());
          return it_ == it.it_;
        }

      private:
        const_subiterator_type it_;
    };

    typedef const_iterator1 iterator1;

    BOOST_UBLAS_INLINE
    const_iterator1 begin1 () const {
      return const_iterator1 (*this, 0);
    }
    BOOST_UBLAS_INLINE
    const_iterator1 cbegin1 () const {
      return begin1 ();
    }
    BOOST_UBLAS_INLINE
    const_iterator1 end1 () const {
      return const_iterator1 (*this, size_common_);
    }
    BOOST_UBLAS_INLINE
    const_iterator1 cend1 () const {
      return end1 ();
    }

    class const_iterator2:
      public bnu::container_const_reference<diagonal_matrix_fixed>,
      public bnu::bidirectional_iterator_base<bnu::sparse_bidirectional_iterator_tag,
                                              const_iterator2, value_type> {
      public:
        typedef typename diagonal_matrix_fixed::value_type value_type;
        typedef typename diagonal_matrix_fixed::difference_type difference_type;
        typedef typename diagonal_matrix_fixed::const_reference reference;
        typedef typename diagonal_matrix_fixed::const_pointer pointer;

        typedef const_iterator1 dual_iterator_type;
        typedef const_reverse_iterator1 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator2 ():
          bnu::container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator2 (const self_type &m, const const_subiterator_type &it):
          bnu::container_const_reference<self_type> (m), it_ (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator2 &operator ++ () {
          BOOST_UBLAS_CHECK (it_ < (*this) ().size_common_, bnu::bad_index ());
          ++it_;
          return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator2 &operator -- () {
          BOOST_UBLAS_CHECK (it_ > 0, bnu::bad_index ());
          --it_;
          return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const {
          return m_value;
        }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 begin () const {
          return const_iterator1 ((*this) (), it_);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 cbegin () const {
          return begin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 end () const {
          return const_iterator1 ((*this) (), it_ + 1);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 cend () const {
          return end ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 rbegin () const {
          return const_reverse_iterator1 (end ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 crbegin () const {
          return rbegin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 rend () const {
          return const_reverse_iterator1 (begin ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 crend () const {
          return rend ();
        }
#endif

        // Indices
        BOOST_UBLAS_INLINE
        size_type index1 () const {
          return it_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
          return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator2 &operator = (const const_iterator2 &it) {
          bnu::container_const_reference<self_type>::assign (&it ());
          it_ = it.it_;
          return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const const_iterator2 &it) const {
          BOOST_UBLAS_CHECK (&(*this) () == &it (), bnu::external_logic ());
          return it_ == it.it_;
        }

      private:
        const_subiterator_type it_;
    };

    typedef const_iterator2 iterator2;

    BOOST_UBLAS_INLINE
    const_iterator2 begin2 () const {
      return const_iterator2 (*this, 0);
    }
    BOOST_UBLAS_INLINE
    const_iterator2 cbegin2 () const {
      return begin2 ();
    }
    BOOST_UBLAS_INLINE
    const_iterator2 end2 () const {
      return const_iterator2 (*this, size_common_);
    }
    BOOST_UBLAS_INLINE
    const_iterator2 cend2 () const {
      return end2 ();
    }

    // Reverse iterators

    BOOST_UBLAS_INLINE
    const_reverse_iterator1 rbegin1 () const {
      return const_reverse_iterator1 (end1 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator1 crbegin1 () const {
      return rbegin1 ();
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator1 rend1 () const {
      return const_reverse_iterator1 (begin1 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator1 crend1 () const {
      return rend1 ();
    }

    BOOST_UBLAS_INLINE
    const_reverse_iterator2 rbegin2 () const {
      return const_reverse_iterator2 (end2 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator2 crbegin2 () const {
      return rbegin2 ();
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator2 rend2 () const {
      return const_reverse_iterator2 (begin2 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator2 crend2 () const {
      return rend2 ();
    }

    // Serialization
    template<typename Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */){

      // we need to copy to a collection_size_type to get a portable
      // and efficient serialization
      boost::serialization::collection_size_type s1 (m_size);
      boost::serialization::collection_size_type s2 (m_size);

      // serialize the sizes
      ar & boost::serialization::make_nvp("size1",s1)
        & boost::serialization::make_nvp("size2",s2);

      // copy the values back if loading
      if (Archive::is_loading::value) {
        m_size = std::min(s1, s2);
      }
    }

  private:
    size_type m_size;
    size_type size_common_;
    static const value_type m_zero;
    static const value_type m_value;
};

template<typename T, T Value, typename ALLOC>
const typename diagonal_matrix_fixed<T, Value, ALLOC>::value_type diagonal_matrix_fixed<T, Value, ALLOC>::m_zero{};
template<typename T, T Value, typename ALLOC>
const typename diagonal_matrix_fixed<T, Value, ALLOC>::value_type diagonal_matrix_fixed<T, Value, ALLOC>::m_value{Value};

// diagonal_matrix is copied from identity_matrix in boost/numeric/ublas/matrix.hpp
// Amended by me.

template<typename T,
         typename ALLOC = std::allocator<T>>
class diagonal_matrix:
      public bnu::matrix_container<diagonal_matrix<T, ALLOC> > {

    typedef const T *const_pointer;
    typedef diagonal_matrix<T, ALLOC> self_type;
  public:
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
    using bnu::matrix_container<self_type>::operator ();
#endif
    typedef typename boost::allocator_size_type<ALLOC>::type size_type;
    typedef typename boost::allocator_difference_type<ALLOC>::type difference_type;
    typedef T value_type;
    typedef const T &const_reference;
    typedef T &reference;
    typedef const bnu::matrix_reference<const self_type> const_closure_type;
    typedef bnu::matrix_reference<self_type> closure_type;
    typedef bnu::sparse_tag storage_category;
    typedef bnu::unknown_orientation_tag orientation_category;
    typedef bnu::vector<value_type> vector_type;

    // Construction and destruction
    BOOST_UBLAS_INLINE
    diagonal_matrix ():
      bnu::matrix_container<self_type> (),
      m_size (0) {}
    BOOST_UBLAS_INLINE
    diagonal_matrix (size_type size):
      bnu::matrix_container<self_type> (),
      m_size (size) {}
    BOOST_UBLAS_INLINE
    diagonal_matrix (size_type size1, size_type size2):
      bnu::matrix_container<self_type> (),
      m_size (std::min(size1, size2)) {}
    BOOST_UBLAS_INLINE
    diagonal_matrix (const diagonal_matrix &m):
      bnu::matrix_container<self_type> (),
      m_size (m.m_size),
      m_values (m.m_values){}
    BOOST_UBLAS_INLINE
    diagonal_matrix (const vector_type & p_values):
      bnu::matrix_container<self_type> (),
      m_size (p_values.size()),
      m_values(p_values) {}

    // Accessors
    BOOST_UBLAS_INLINE
    size_type size1 () const {
      return m_size;
    }
    BOOST_UBLAS_INLINE
    size_type size2 () const {
      return m_size;
    }

    // Resizing
    BOOST_UBLAS_INLINE
    void resize (size_type size, bool /*preserve*/ = true) {
      m_size = size;
    }
    BOOST_UBLAS_INLINE
    void resize (size_type size1, size_type size2, bool /*preserve*/ = true) {
      m_size = std::min(size1, size2);
    }

    // Element access
    BOOST_UBLAS_INLINE
    const_reference operator () (size_type i, size_type j) const {
      if((m_size != 0) && (i == j))
      {
        return m_values(i);
      }
      return m_zero;
    }

    // Assignment
#ifdef BOOST_UBLAS_MOVE_SEMANTICS

    /*! @note "pass by value" the key idea to enable move semantics */
    BOOST_UBLAS_INLINE
    diagonal_matrix &operator = (diagonal_matrix m) {
      assign_temporary(m);
      return *this;
    }
#else
    BOOST_UBLAS_INLINE
    diagonal_matrix &operator = (const diagonal_matrix &m) {
      m_size = m.m_size;
      m_values = m.m_values;
      return *this;
    }
#endif
    BOOST_UBLAS_INLINE
    diagonal_matrix &assign_temporary (diagonal_matrix &m) {
      swap (m);
      return *this;
    }

    // Swapping
    BOOST_UBLAS_INLINE
    constexpr void swap (diagonal_matrix &m) noexcept {
      if (this != &m) {
        std::swap (m_size, m.m_size);
        std::swap (m_values, m.m_values);
      }
    }
    BOOST_UBLAS_INLINE
    constexpr friend void swap (diagonal_matrix &m1, diagonal_matrix &m2) noexcept {
      m1.swap (m2);
    }

    // Iterator types
  private:
    // Use an index
    typedef size_type const_subiterator_type;

  public:
    class const_iterator1;
    class const_iterator2;
    typedef bnu::reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
    typedef bnu::reverse_iterator_base2<const_iterator2> const_reverse_iterator2;

    // Element lookup
    BOOST_UBLAS_INLINE
    const_iterator1 find1 (int rank, size_type i, size_type j) const {
      if (rank == 1) {
        i = (std::max) (i, j);
        i = (std::min) (i, j + 1);
      }
      return const_iterator1 (*this, i);
    }
    BOOST_UBLAS_INLINE
    const_iterator2 find2 (int rank, size_type i, size_type j) const {
      if (rank == 1) {
        j = (std::max) (j, i);
        j = (std::min) (j, i + 1);
      }
      return const_iterator2 (*this, j);
    }

    class const_iterator1:
      public bnu::container_const_reference<diagonal_matrix>,
      public bnu::bidirectional_iterator_base<bnu::sparse_bidirectional_iterator_tag,
                                              const_iterator1, value_type> {
      public:
        typedef typename diagonal_matrix::value_type value_type;
        typedef typename diagonal_matrix::difference_type difference_type;
        typedef typename diagonal_matrix::const_reference reference;
        typedef typename diagonal_matrix::const_pointer pointer;

        typedef const_iterator2 dual_iterator_type;
        typedef const_reverse_iterator2 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator1 ():
          bnu::container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator1 (const self_type &m, const const_subiterator_type &it):
          bnu::container_const_reference<self_type> (m), it_ (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator1 &operator ++ () {
          BOOST_UBLAS_CHECK (it_ < (*this) ().size1 (), bnu::bad_index ());
          ++it_;
          return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator1 &operator -- () {
          BOOST_UBLAS_CHECK (it_ > 0, bnu::bad_index ());
          --it_;
          return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const {
          return (*this)().m_values(it_);
        }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 begin () const {
          return const_iterator2 ((*this) (), it_);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 cbegin () const {
          return begin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 end () const {
          return const_iterator2 ((*this) (), it_ + 1);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 cend () const {
          return end ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 rbegin () const {
          return const_reverse_iterator2 (end ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 crbegin () const {
          return rbegin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 rend () const {
          return const_reverse_iterator2 (begin ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 crend () const {
          return rend ();
        }
#endif

        // Indices
        BOOST_UBLAS_INLINE
        size_type index1 () const {
          return it_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
          return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator1 &operator = (const const_iterator1 &it) {
          bnu::container_const_reference<self_type>::assign (&it ());
          it_ = it.it_;
          return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const const_iterator1 &it) const {
          BOOST_UBLAS_CHECK (&(*this) () == &it (), bnu::external_logic ());
          return it_ == it.it_;
        }

      private:
        const_subiterator_type it_;
    };

    typedef const_iterator1 iterator1;

    BOOST_UBLAS_INLINE
    const_iterator1 begin1 () const {
      return const_iterator1 (*this, 0);
    }
    BOOST_UBLAS_INLINE
    const_iterator1 cbegin1 () const {
      return begin1 ();
    }
    BOOST_UBLAS_INLINE
    const_iterator1 end1 () const {
      return const_iterator1 (*this, m_size);
    }
    BOOST_UBLAS_INLINE
    const_iterator1 cend1 () const {
      return end1 ();
    }

    class const_iterator2:
      public bnu::container_const_reference<diagonal_matrix>,
      public bnu::bidirectional_iterator_base<bnu::sparse_bidirectional_iterator_tag,
                                              const_iterator2, value_type> {
      public:
        typedef typename diagonal_matrix::value_type value_type;
        typedef typename diagonal_matrix::difference_type difference_type;
        typedef typename diagonal_matrix::const_reference reference;
        typedef typename diagonal_matrix::const_pointer pointer;

        typedef const_iterator1 dual_iterator_type;
        typedef const_reverse_iterator1 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator2 ():
          bnu::container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator2 (const self_type &m, const const_subiterator_type &it):
          bnu::container_const_reference<self_type> (m), it_ (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator2 &operator ++ () {
          BOOST_UBLAS_CHECK (it_ < (*this) ().m_size, bnu::bad_index ());
          ++it_;
          return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator2 &operator -- () {
          BOOST_UBLAS_CHECK (it_ > 0, bnu::bad_index ());
          --it_;
          return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const {
          return (*this)().m_values(it_);
        }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 begin () const {
          return const_iterator1 ((*this) (), it_);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 cbegin () const {
          return begin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 end () const {
          return const_iterator1 ((*this) (), it_ + 1);
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 cend () const {
          return end ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 rbegin () const {
          return const_reverse_iterator1 (end ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 crbegin () const {
          return rbegin ();
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 rend () const {
          return const_reverse_iterator1 (begin ());
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 crend () const {
          return rend ();
        }
#endif

        // Indices
        BOOST_UBLAS_INLINE
        size_type index1 () const {
          return it_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const {
          return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator2 &operator = (const const_iterator2 &it) {
          bnu::container_const_reference<self_type>::assign (&it ());
          it_ = it.it_;
          return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const const_iterator2 &it) const {
          BOOST_UBLAS_CHECK (&(*this) () == &it (), bnu::external_logic ());
          return it_ == it.it_;
        }

      private:
        const_subiterator_type it_;
    };

    typedef const_iterator2 iterator2;

    BOOST_UBLAS_INLINE
    const_iterator2 begin2 () const {
      return const_iterator2 (*this, 0);
    }
    BOOST_UBLAS_INLINE
    const_iterator2 cbegin2 () const {
      return begin2 ();
    }
    BOOST_UBLAS_INLINE
    const_iterator2 end2 () const {
      return const_iterator2 (*this, m_size);
    }
    BOOST_UBLAS_INLINE
    const_iterator2 cend2 () const {
      return end2 ();
    }

    // Reverse iterators

    BOOST_UBLAS_INLINE
    const_reverse_iterator1 rbegin1 () const {
      return const_reverse_iterator1 (end1 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator1 crbegin1 () const {
      return rbegin1 ();
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator1 rend1 () const {
      return const_reverse_iterator1 (begin1 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator1 crend1 () const {
      return rend1 ();
    }

    BOOST_UBLAS_INLINE
    const_reverse_iterator2 rbegin2 () const {
      return const_reverse_iterator2 (end2 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator2 crbegin2 () const {
      return rbegin2 ();
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator2 rend2 () const {
      return const_reverse_iterator2 (begin2 ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator2 crend2 () const {
      return rend2 ();
    }

    // Serialization
    template<typename Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */){

      // we need to copy to a collection_size_type to get a portable
      // and efficient serialization
      boost::serialization::collection_size_type s1 (m_size);
      boost::serialization::collection_size_type s2 (m_size);

      // serialize the sizes
      ar & boost::serialization::make_nvp("size1",s1)
        & boost::serialization::make_nvp("size2",s2);

      // copy the values back if loading
      if (Archive::is_loading::value) {
        m_size = std::min(s1, s2);
      }
    }

  private:
    size_type m_size{};
    static const value_type m_zero;
    vector_type m_values{};
};

template<typename T, typename ALLOC>
const typename diagonal_matrix<T, ALLOC>::value_type diagonal_matrix<T, ALLOC>::m_zero{};

} // namespace yafiyogi::yy_maths
