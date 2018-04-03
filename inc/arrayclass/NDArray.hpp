#ifndef ARRAYITERATOR_HPP
#define ARRAYITERATOR_HPP

#include "Array.h"

namespace Array {

template <typename T>
class NDarray1 : public array1<T>
{
public:
    /**
     * @brief The const_iterator class
     */
    class abstract_iterator
    {
    public:
        /**
         * @brief value_type std::iterator_traits::value_type
         */
        typedef T value_type;

        /**
         * @brief reference std::iterator_traits::reference
         */
        typedef T& reference;

        /**
         * @brief pointer std::iterator_traits::pointer
         */
        typedef T* pointer;

        /**
         * @brief iterator_category std::iterator_traits::iterator_category
         */
        typedef std::random_access_iterator_tag iterator_category;

        /**
         * @brief difference_type std::iterator_traits::difference_type
         */
        typedef std::ptrdiff_t difference_type;

        /**
         * default and initializing constructors
         */
        abstract_iterator(T* ptr=nullptr) : _ptr(ptr) {}

        friend void swap(abstract_iterator& lhs, abstract_iterator& rhs)
        { std::swap(lhs._ptr,rhs._ptr); }

        abstract_iterator& operator=(abstract_iterator& other)
        { swap(*this, other); return *this; }

    #if __cplusplus > 199711L // for C++11, explicitly define use of default
        abstract_iterator(const abstract_iterator&) = default;

        abstract_iterator(abstract_iterator&&) = default;

        abstract_iterator& operator=(abstract_iterator&&) = default;
    #endif // C++11

    public: // arithmetic operators
        inline abstract_iterator& operator+=(const std::ptrdiff_t& n)
        { _ptr += n; return this; }

        inline const abstract_iterator operator+(const std::ptrdiff_t& n) const
        { return abstract_iterator(*this) += n; }

        inline abstract_iterator operator++()
        { _ptr++; return *this; }

        inline abstract_iterator operator++(int)
        { abstract_iterator i = *this; _ptr++; return i; }

        friend inline const abstract_iterator
        operator+( const std::ptrdiff_t& lhs, const abstract_iterator& rhs)
        { return rhs+lhs; }

        inline abstract_iterator& operator-=(const std::ptrdiff_t& n)
        { _ptr -= n; return this; }

        inline abstract_iterator operator--()
        { _ptr--; return *this; }

        inline abstract_iterator operator--(int)
        { abstract_iterator i = *this; _ptr--; return i; }

        inline const abstract_iterator operator-(const std::ptrdiff_t& n) const
        { return abstract_iterator(*this) -= n; }

        friend inline const abstract_iterator
        operator-( const std::ptrdiff_t& lhs, const abstract_iterator& rhs)
        { return rhs-lhs; }

    public: // comparission operators
        inline const bool operator==(const abstract_iterator& other)
        { return _ptr == other._ptr; }

        inline const bool operator!=(const abstract_iterator& other)
        { return !(*this == other); }

        inline const bool operator>(const abstract_iterator& other)
        { return _ptr > other._ptr; }

        inline const bool operator<(const abstract_iterator& other)
        { return _ptr < other._ptr; }

        inline const bool operator>=(const abstract_iterator& other)
        { return !operator<(other); }

        inline const bool operator<=(const abstract_iterator& other)
        { return !operator>(other); }

    protected:
        T* _ptr;
    };

    class const_iterator : public abstract_iterator
    {
    public:
        // default and initializing constructors
        const_iterator(T* ptr=nullptr) : abstract_iterator(ptr) {}

    public: // access operators
        inline const T& operator[](std::ptrdiff_t& n)
        { return *(abstract_iterator::_ptr + n); }

        inline const T& operator*()
        { return *abstract_iterator::_ptr; }

        inline const T* operator->()
        { return abstract_iterator::_ptr; }
    };

    class iterator : public abstract_iterator
    {
    public:
        // default and initializing constructors
        iterator(T* ptr=nullptr) : abstract_iterator(ptr) {}

    public: // access operators
        inline T& operator[](std::ptrdiff_t& n)
        { return *(abstract_iterator::_ptr + n); }

        inline T& operator*()
        { return *abstract_iterator::_ptr; }

        inline T* operator->()
        { return abstract_iterator::_ptr; }
    };

    NDarray1(size_t size)
    : array1<T>(size)
    {}

    inline size_t size() const
    { return array1<T>::size; }

    iterator begin()
    { return iterator(array1<T>::v); }

    iterator end()
    { return iterator(array1<T>::v + array1<T>::size); }

    abstract_iterator begin() const
    { return abstract_iterator(array1<T>::v); }

    abstract_iterator end() const
    { return abstract_iterator(array1<T>::v + array1<T>::size); }
};

} // namespace Array


#endif // ARRAYITERATOR_HPP
