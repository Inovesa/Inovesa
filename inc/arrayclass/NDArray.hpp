#ifndef ARRAYITERATOR_HPP
#define ARRAYITERATOR_HPP

#include "Array.h"

#include <initializer_list>

namespace Array {

template <typename T>
class NDarray1 : public array1<T>
{
private:
    /**
     * @brief The NDarray1::abstract_iterator class serves
     *  as base class for NDarray1::const_iterator and NDarray1::iterator.
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

        friend void swap(abstract_iterator& lhs, abstract_iterator& rhs)
        { std::swap(lhs._ptr,rhs._ptr); }

        /**
         * @brief operator =
         * @param other
         * @return
         */
        abstract_iterator& operator=(abstract_iterator& other)
        { swap(*this, other); return *this; }

        /* constructors are marked protected,
         * so that no abstract_iterator is constructed */
    protected:
        /**
         * @brief abstract_iterator default and initializing constructors
         * @param ptr
         */
        abstract_iterator(T* ptr=nullptr) : _ptr(ptr) {}

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
        /* access operators (const and non-const iterators use them) */
        inline const T& operator[](std::ptrdiff_t& n) const
        { return *(abstract_iterator::_ptr + n); }

        inline const T& operator*() const
        { return *abstract_iterator::_ptr; }

        inline const T* operator->() const
        { return abstract_iterator::_ptr; }

    protected:
        T* _ptr;
    };

public:
    /**
     * @brief The NDarray1::const_iterator class
     */
    class const_iterator : public abstract_iterator
    {
    public:
        // default and initializing constructors
        const_iterator(T* ptr=nullptr) : abstract_iterator(ptr) {}

    public: // access operators
        inline const T& operator[](std::ptrdiff_t& n)
        { return (abstract_iterator::operator[](n)); }

        inline const T& operator*()
        { return (abstract_iterator::operator*()); }

        inline const T* operator->()
        { return (abstract_iterator::operator->()); }
    };

    /**
     * @brief The NDarray1::iterator class
     */
    class iterator : public abstract_iterator
    {
    public:
        // default and initializing constructors
        iterator(T* ptr=nullptr) : abstract_iterator(ptr) {}

    public: // access operators
        inline T& operator[](std::ptrdiff_t& n)
        { return const_cast<T&>(abstract_iterator::operator[](n)); }

        inline T& operator*()
        { return const_cast<T&>(abstract_iterator::operator*()); }

        inline T* operator->()
        { return const_cast<T&>(abstract_iterator::operator->()); }
    };

    /**
     * @brief NDarray1
     */
    NDarray1() = delete;

    /**
     * @brief NDarray1
     * @param size
     * @param align
     */
    NDarray1(size_t size, size_t align=0)
    : array1<T>(size, align)
    {}

    /**
     * @brief NDarray1
     * @param il initializer list
     *
     * It migth be worth to change this to allow move construction
     * from non-NDarray data.
     */
    NDarray1(std::initializer_list<T> il)
    : array1<T>(il.size())
    { std::copy(il.begin(),il.end(),array1<T>::v); }

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
