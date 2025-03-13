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

#include "boost/numeric/ublas/matrix.hpp"

namespace boost::numeric::ublas {

template<class T, T Value = 1, class ALLOC = std::allocator<T>>
class identity_matrix_v:
      public matrix_container<identity_matrix_v<T, Value, ALLOC> > {

    typedef const T *const_pointer;
    typedef identity_matrix_v<T, Value, ALLOC> self_type;
  public:
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
    using matrix_container<self_type>::operator ();
#endif
    typedef typename boost::allocator_size_type<ALLOC>::type size_type;
    typedef typename boost::allocator_difference_type<ALLOC>::type difference_type;
    typedef T value_type;
    typedef const T &const_reference;
    typedef T &reference;
    typedef const matrix_reference<const self_type> const_closure_type;
    typedef matrix_reference<self_type> closure_type;
    typedef sparse_tag storage_category;
    typedef unknown_orientation_tag orientation_category;

    // Construction and destruction
    BOOST_UBLAS_INLINE
    identity_matrix_v ():
      matrix_container<self_type> (),
      size1_ (0), size2_ (0), size_common_ (0) {}
    BOOST_UBLAS_INLINE
    identity_matrix_v (size_type size):
      matrix_container<self_type> (),
      size1_ (size), size2_ (size), size_common_ ((std::min) (size1_, size2_)) {}
    BOOST_UBLAS_INLINE
    identity_matrix_v (size_type size1, size_type size2):
      matrix_container<self_type> (),
      size1_ (size1), size2_ (size2), size_common_ ((std::min) (size1_, size2_)) {}
    BOOST_UBLAS_INLINE
    identity_matrix_v (const identity_matrix_v &m):
      matrix_container<self_type> (),
      size1_ (m.size1_), size2_ (m.size2_), size_common_ ((std::min) (size1_, size2_)) {}

    // Accessors
    BOOST_UBLAS_INLINE
    size_type size1 () const {
      return size1_;
    }
    BOOST_UBLAS_INLINE
    size_type size2 () const {
      return size2_;
    }

    // Resizing
    BOOST_UBLAS_INLINE
    void resize (size_type size, bool /*preserve*/ = true) {
      size1_ = size;
      size2_ = size;
      size_common_ = ((std::min)(size1_, size2_));
    }
    BOOST_UBLAS_INLINE
    void resize (size_type size1, size_type size2, bool /*preserve*/ = true) {
      size1_ = size1;
      size2_ = size2;
      size_common_ = ((std::min)(size1_, size2_));
    }

    // Element access
    BOOST_UBLAS_INLINE
    const_reference operator () (size_type i, size_type j) const {
      if (i == j)
        return value_;
      else
        return zero_;
    }

    // Assignment
    BOOST_UBLAS_INLINE
    identity_matrix_v &operator = (const identity_matrix_v &m) {
      size1_ = m.size1_;
      size2_ = m.size2_;
      size_common_ = m.size_common_;
      return *this;
    }
    BOOST_UBLAS_INLINE
    identity_matrix_v &assign_temporary (identity_matrix_v &m) {
      swap (m);
      return *this;
    }

    // Swapping
    BOOST_UBLAS_INLINE
    void swap (identity_matrix_v &m) {
      if (this != &m) {
        std::swap (size1_, m.size1_);
        std::swap (size2_, m.size2_);
        std::swap (size_common_, m.size_common_);
      }
    }
    BOOST_UBLAS_INLINE
    friend void swap (identity_matrix_v &m1, identity_matrix_v &m2) {
      m1.swap (m2);
    }

    // Iterator types
  private:
    // Use an index
    typedef size_type const_subiterator_type;

  public:
    class const_iterator1;
    class const_iterator2;
    typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
    typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;

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
      public container_const_reference<identity_matrix_v>,
      public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                         const_iterator1, value_type> {
      public:
        typedef typename identity_matrix_v::value_type value_type;
        typedef typename identity_matrix_v::difference_type difference_type;
        typedef typename identity_matrix_v::const_reference reference;
        typedef typename identity_matrix_v::const_pointer pointer;

        typedef const_iterator2 dual_iterator_type;
        typedef const_reverse_iterator2 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator1 ():
          container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator1 (const self_type &m, const const_subiterator_type &it):
          container_const_reference<self_type> (m), it_ (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator1 &operator ++ () {
          BOOST_UBLAS_CHECK (it_ < (*this) ().size1 (), bad_index ());
          ++it_;
          return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator1 &operator -- () {
          BOOST_UBLAS_CHECK (it_ > 0, bad_index ());
          --it_;
          return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const {
          return value_;
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
          container_const_reference<self_type>::assign (&it ());
          it_ = it.it_;
          return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const const_iterator1 &it) const {
          BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
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
      public container_const_reference<identity_matrix_v>,
      public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
                                         const_iterator2, value_type> {
      public:
        typedef typename identity_matrix_v::value_type value_type;
        typedef typename identity_matrix_v::difference_type difference_type;
        typedef typename identity_matrix_v::const_reference reference;
        typedef typename identity_matrix_v::const_pointer pointer;

        typedef const_iterator1 dual_iterator_type;
        typedef const_reverse_iterator1 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator2 ():
          container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator2 (const self_type &m, const const_subiterator_type &it):
          container_const_reference<self_type> (m), it_ (it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator2 &operator ++ () {
          BOOST_UBLAS_CHECK (it_ < (*this) ().size_common_, bad_index ());
          ++it_;
          return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator2 &operator -- () {
          BOOST_UBLAS_CHECK (it_ > 0, bad_index ());
          --it_;
          return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const {
          return value_;
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
          container_const_reference<self_type>::assign (&it ());
          it_ = it.it_;
          return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const const_iterator2 &it) const {
          BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
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
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */){

      // we need to copy to a collection_size_type to get a portable
      // and efficient serialization
      serialization::collection_size_type s1 (size1_);
      serialization::collection_size_type s2 (size2_);

      // serialize the sizes
      ar & serialization::make_nvp("size1",s1)
        & serialization::make_nvp("size2",s2);

      // copy the values back if loading
      if (Archive::is_loading::value) {
        size1_ = s1;
        size2_ = s2;
        size_common_ = ((std::min)(size1_, size2_));
      }
    }

  private:
    size_type size1_;
    size_type size2_;
    size_type size_common_;
    static const value_type zero_;
    static const value_type value_;
};

template<class T, T Value, class ALLOC>
const typename identity_matrix_v<T, Value, ALLOC>::value_type identity_matrix_v<T, Value, ALLOC>::zero_ = T(/*zero*/);
template<class T, T Value, class ALLOC>
const typename identity_matrix_v<T, Value, ALLOC>::value_type identity_matrix_v<T, Value, ALLOC>::value_ (Value); // ISSUE: need 'one'-traits here

} // namespace boost::numeric::ublas
