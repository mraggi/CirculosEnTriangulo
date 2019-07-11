#pragma once

#include <iterator>
#include <type_traits>

template <class Iter>
class View
{
public:
    using iterator = Iter;
    using const_iterator = Iter;
    using value_type = typename Iter::value_type;
    using size_type = long long;

    View(Iter first, Iter last) : begin_(first), end_(last) {}

    Iter begin() const { return begin_; }
    Iter end() const { return end_; }

    value_type operator[](std::size_t i) const { return *(begin() + i); }
    value_type& operator[](std::size_t i) { return *(begin() + i); }

    size_type size() const { return end_ - begin_; }

private:
    Iter begin_;
    Iter end_;
};

inline void replace_by_bigger(double& a, double b)
{
    if (a < b)
        a = b;
}

inline void replace_by_smaller(double& a, double b)
{
    if (a > b)
        a = b;
}
