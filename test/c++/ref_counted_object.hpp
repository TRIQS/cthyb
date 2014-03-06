#pragma once
#include <triqs/utility/exceptions.hpp>
#include <iostream>

namespace triqs {
namespace utility {

 // derive the object to be counted from this little object
 struct ref_counted_object {
  int n_ref;
  ref_counted_object() { n_ref = 1; }
 };

 // A light shared_ptr
 template <typename Object> class ref_counted_object_ptr {
  public:
  Object* p;
  void decref() {
   //if (p && ((--p->n_ref) == 0)) {
    //std::cerr << "deallocating " << p << std::endl;
   // delete p;
   //}
  }
  void incref() {
   //if (p) p->n_ref++;
  }

  public:
  ref_counted_object_ptr(Object* ptr) : p(ptr) {}
  ref_counted_object_ptr(std::nullptr_t) : p(nullptr) {}
  ref_counted_object_ptr() : ref_counted_object_ptr(nullptr) {}

  template <typename... T> explicit ref_counted_object_ptr(T&&... x) {
   p = new Object{x...};
   //std::cerr << " allocating " << p << std::endl;
  }

  ~ref_counted_object_ptr() { decref(); }

  ref_counted_object_ptr(ref_counted_object_ptr&& x) {
   p = x.p;
   x.p = nullptr;
  }

  ref_counted_object_ptr(ref_counted_object_ptr const& x) {
   p = x.p;
   incref();
  }

  ref_counted_object_ptr& operator=(ref_counted_object_ptr const& x) {
   decref();
   p = x.p;
   incref();
   return *this;
  }

  ref_counted_object_ptr& operator=(ref_counted_object_ptr&& x) {
   decref();
   p = x.p;
   x.p = nullptr;
   return *this;
  }

  Object* operator->() { return p; }
  Object const * operator->() const { return p; }
  operator bool() const { return p != nullptr; }
  bool operator==(ref_counted_object_ptr x) const { return p == x.p; }
  bool operator!=(ref_counted_object_ptr x) const { return !(operator==(x)); }
 };
}
}

