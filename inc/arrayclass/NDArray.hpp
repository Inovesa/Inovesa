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
        virtual inline abstract_iterator
        operator+=(const std::ptrdiff_t& n)
        { _ptr += n; return *this; }

        virtual inline const abstract_iterator
        operator+(const std::ptrdiff_t& n) const
        { return abstract_iterator(*this) += n; }

        virtual inline abstract_iterator operator++()
        { _ptr++; return *this; }

        virtual inline abstract_iterator operator++(int)
        { abstract_iterator i = *this; _ptr++; return i; }

        virtual inline abstract_iterator operator-=(const std::ptrdiff_t& n)
        { _ptr -= n; return *this; }

        virtual inline abstract_iterator operator--()
        { _ptr--; return *this; }

        virtual inline abstract_iterator operator--(int)
        { abstract_iterator i = *this; _ptr--; return i; }

        virtual inline const abstract_iterator
        operator-(const std::ptrdiff_t& n) const
        { return abstract_iterator(*this) -= n; }

        friend inline const abstract_iterator&
        operator+( const std::ptrdiff_t& lhs, const abstract_iterator& rhs)
        { return rhs+lhs; }

        friend inline const abstract_iterator
        operator-( const std::ptrdiff_t& lhs, const abstract_iterator& rhs)
        { return rhs-lhs; }

    public: // comparission operators
        inline const bool operator==(const abstract_iterator& other)
        { return _ptr == other._ptr; }

        inline const bool operator!=(const abstract_iterator& other)
        { return !(*this == other); }

        virtual inline const bool operator>(const abstract_iterator& other)
        { return _ptr > other._ptr; }

        virtual inline const bool operator<(const abstract_iterator& other)
        { return _ptr < other._ptr; }

        virtual inline const bool operator>=(const abstract_iterator& other)
        { return !operator<(other); }

        virtual inline const bool operator<=(const abstract_iterator& other)
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

    class abstract_reverse_iterator: public abstract_iterator
    {
    protected:
        // default and initializing constructors
        abstract_reverse_iterator(T* ptr=nullptr)
            : abstract_iterator((ptr==nullptr)?nullptr:ptr-1) {}

    public:
        virtual inline abstract_iterator
        operator+=(const std::ptrdiff_t& n)
        { return abstract_iterator::operator -=(n); }

        virtual inline const abstract_iterator
        operator+(const std::ptrdiff_t& n) const final
        { return abstract_iterator::operator-(n); }

        virtual inline abstract_iterator operator++() final
        { return abstract_iterator::operator--(); }

        virtual inline abstract_iterator operator++(int) final
        { return abstract_iterator::operator--(0); }

        virtual inline abstract_iterator
        operator-=(const std::ptrdiff_t& n) final
        { return abstract_iterator::operator +=(n); }

        virtual inline abstract_iterator operator--() final
        { return abstract_iterator::operator ++(); }

        virtual inline abstract_iterator operator--(int) final
        { return abstract_iterator::operator ++(0); }

        virtual inline const abstract_iterator
        operator-(const std::ptrdiff_t& n) const final
        { return abstract_iterator::operator+(n); }

    public: // comparission operators
        inline const bool operator>(const abstract_iterator& other) final
        { return abstract_iterator::operator<(other); }

        inline const bool operator<(const abstract_iterator& other) final
        { return abstract_iterator::operator>(other); }

        inline const bool operator>=(const abstract_iterator& other) final
        { return abstract_iterator::operator<=(other); }

        inline const bool operator<=(const abstract_iterator& other) final
        { return abstract_iterator::operator>=(other); }
    };



public:
    /**
     * @brief The NDarray1::iterator class
     */
    class iterator : public abstract_iterator
    {
    public:
        // default and initializing constructors
        iterator(T* ptr=nullptr)
        : abstract_iterator(ptr) {}

    public: // access operators
        inline T& operator[](std::ptrdiff_t n)
        { return const_cast<T&>(abstract_iterator::operator[](n)); }

        inline T& operator*()
        { return const_cast<T&>(abstract_iterator::operator*()); }

        inline T* operator->()
        { return const_cast<T&>(abstract_iterator::operator->()); }
    };

    /**
     * @brief The NDarray1::reverse_iterator class
     */
    class reverse_iterator : public abstract_reverse_iterator
    {
    public:
        // default and initializing constructors
        reverse_iterator(T* ptr=nullptr)
        : abstract_reverse_iterator(ptr) {}

    public: // access operators
        inline T& operator[](std::ptrdiff_t n)
        { return const_cast<T&>(abstract_iterator::operator[](-n)); }

        inline T& operator*()
        { return const_cast<T&>(abstract_iterator::operator*()); }

        inline T* operator->()
        { return const_cast<T&>(abstract_iterator::operator->()); }
    };

    /**
     * @brief The NDarray1::const_iterator class
     */
    class const_iterator : public abstract_iterator
    {
    public:
        // default and initializing constructors
        const_iterator(T* ptr=nullptr)
        : abstract_iterator(ptr) {}

    public: // access operators
        inline const T& operator[](std::ptrdiff_t n)
        { return (abstract_iterator::operator[](n)); }

        inline const T& operator*()
        { return (abstract_iterator::operator*()); }

        inline const T* operator->()
        { return (abstract_iterator::operator->()); }
    };

    /**
     * @brief The NDarray1::const_reverse_iterator class
     */
    class const_reverse_iterator: public abstract_reverse_iterator
    {
    public:
        // default and initializing constructors
        const_reverse_iterator(T* ptr=nullptr)
        : abstract_reverse_iterator(ptr) {}

    public: // access operators
        inline const T& operator[](std::ptrdiff_t n)
        { return (abstract_iterator::operator[](-n)); }

        inline const T& operator*()
        { return (abstract_iterator::operator*()); }

        inline const T* operator->()
        { return (abstract_iterator::operator->()); }
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

    iterator begin() noexcept
    { return iterator(array1<T>::v); }

    iterator end() noexcept
    { return iterator(array1<T>::v + array1<T>::size); }

    reverse_iterator rbegin() noexcept
    { return reverse_iterator(array1<T>::v + array1<T>::size); }

    reverse_iterator rend() noexcept
    { return reverse_iterator(array1<T>::v); }

    const_iterator begin() const noexcept
    { return const_iterator(array1<T>::v); }

    const_iterator end() const noexcept
    { return const_iterator(array1<T>::v + array1<T>::size); }

    const_iterator cbegin() const noexcept
    { return const_iterator(array1<T>::v); }

    const_iterator cend() const noexcept
    { return const_iterator(array1<T>::v + array1<T>::size); }

    const_reverse_iterator rbegin() const noexcept
    { return const_reverse_iterator(array1<T>::v + array1<T>::size); }

    const_reverse_iterator rend() const noexcept
    { return const_reverse_iterator(array1<T>::v); }

    const_reverse_iterator crbegin() const noexcept
    { return const_reverse_iterator(array1<T>::v + array1<T>::size); }

    const_reverse_iterator crend() const noexcept
    { return const_reverse_iterator(array1<T>::v); }
};

} // namespace Array


#endif // ARRAYITERATOR_HPP
