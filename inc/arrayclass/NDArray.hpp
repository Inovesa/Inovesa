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
     *
     * @todo: allow "int + const_iterator"
     */
    class const_iterator
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
        const_iterator(T* ptr=nullptr) : _ptr(ptr) {}

    #if __cplusplus > 199711L
    public: // explicitly define use of default (for C++11)
        const_iterator(const const_iterator&) = default;

        const_iterator(const_iterator&&) = default;

        const_iterator& operator=(const const_iterator&) = default;

        const_iterator& operator=(const_iterator&&) = default;
    #endif // C++11

    public: // arithmetic operators
        inline const_iterator& operator+=(const std::ptrdiff_t& n)
        { _ptr += n; return this; }

        inline const const_iterator operator+(const std::ptrdiff_t& n) const
        { return const_iterator(*this) += n; }

        inline const_iterator operator++()
        { _ptr++; return *this; }

        inline const_iterator operator++(int)
        { const_iterator i = *this; _ptr++; return i; }

        inline const_iterator& operator-=(const std::ptrdiff_t& n)
        { _ptr -= n; return this; }

        inline const const_iterator operator-(const std::ptrdiff_t& n) const
        { return const_iterator(*this) -= n; }

        inline const_iterator operator--()
        { _ptr--; return *this; }

        inline const_iterator operator--(int)
        { const_iterator i = *this; _ptr--; return i; }

    public: // access operators
        inline T& operator[](std::ptrdiff_t& n)
        { return *(_ptr + n); }

        inline const T& operator*()
        { return *_ptr; }

        inline T* operator->()
        { return _ptr; }

    public: // comparission operators
        inline bool operator==(const const_iterator& other)
        { return _ptr == other._ptr; }

        inline bool operator!=(const const_iterator& other)
        { return !(*this == other); }

        inline bool operator>(const const_iterator& other)
        { return _ptr > other._ptr; }

        inline bool operator<(const const_iterator& other)
        { return _ptr < other._ptr; }

        inline bool operator>=(const const_iterator& other)
        { return !operator<(other); }

        inline bool operator<=(const const_iterator& other)
        { return !operator>(other); }

    protected:
        T* _ptr;
    };

    class iterator : public const_iterator
    {
    public:
        // default and initializing constructors
        iterator(T* ptr=nullptr) : const_iterator(ptr) {}

        inline T& operator*()
        { return *const_iterator::_ptr; }
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

    const_iterator begin() const
    { return const_iterator(array1<T>::v); }

    const_iterator end() const
    { return const_iterator(array1<T>::v + array1<T>::size); }
};

} // namespace Array


#endif // ARRAYITERATOR_HPP
