#ifndef ARRAYITERATOR_HPP
#define ARRAYITERATOR_HPP

#include "Array.h"

namespace Array {

template <typename T>
class NDarray1 : public array1<T>
{
public:
    class iterator
    {
    public:
        typedef T value_type;
        typedef T& reference;
        typedef T* pointer;
        typedef std::forward_iterator_tag iterator_category;
        typedef std::ptrdiff_t difference_type;
        iterator(T* ptr) : __ptr(ptr) { }
        iterator operator++() { __ptr++; return *this; }
        iterator operator++(int) { iterator i = *this; __ptr++; return i; }
        T& operator*() { return *__ptr; }
        T* operator->() { return __ptr; }
        bool operator==(const iterator& rhs) { return __ptr == rhs.__ptr; }
        bool operator!=(const iterator& rhs) { return __ptr != rhs.__ptr; }
    private:
        T* __ptr;
   };

    class const_iterator
    {
    public:
        typedef T value_type;
        typedef T& reference;
        typedef T* pointer;
        typedef std::forward_iterator_tag iterator_category;
        typedef int difference_type;
        const_iterator(T* ptr) : __ptr(ptr) { }
        const_iterator operator++() { __ptr++; return *this; }
        const_iterator operator++(int) { const_iterator i = *this; __ptr++; return i; }
        T& operator*() { return *__ptr; }
        T* operator->() { return __ptr; }
        bool operator==(const iterator& rhs) { return __ptr == rhs.__ptr; }
        bool operator!=(const iterator& rhs) { return __ptr != rhs.__ptr; }
    private:
        T* __ptr;
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
