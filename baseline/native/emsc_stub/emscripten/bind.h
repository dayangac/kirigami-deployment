// Minimal native stub for <emscripten/bind.h>.
//
// The authors' sources include this header from utils/Hmesh.h,
// utils/conversions.h and geometry/unit_pattern.h.  Only `emscripten::val` is
// used, and only inside the JS glue (from_js_mesh / to_js_mesh /
// hmesh_to_js_object / js_object_to_hmesh) that the native CLI never calls.
// This stub lets those declarations compile; every operation throws if it is
// ever actually reached, so a silent wrong answer is impossible.
#pragma once
#include <stdexcept>
#include <string>

namespace emscripten {

[[noreturn]] inline void val_not_available() {
  throw std::logic_error(
      "emscripten::val is not available in the native baseline build");
}

class val {
public:
  val() = default;

  static val array() { val_not_available(); }
  static val object() { val_not_available(); }
  static val undefined() { val_not_available(); }
  static val null() { val_not_available(); }

  template <typename K> val operator[](const K &) const {
    val_not_available();
  }
  template <typename T> T as() const { val_not_available(); }
  template <typename K, typename V> void set(const K &, const V &) {
    val_not_available();
  }
  template <typename R, typename... Args>
  R call(const char *, Args &&...) const {
    val_not_available();
  }
  bool hasOwnProperty(const char *) const { val_not_available(); }
};

} // namespace emscripten
